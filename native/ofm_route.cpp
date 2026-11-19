/* Parallel finite-volume OFM routing. MATLAB reference: src/explicit.m and
 * src/diffusive_face.m. Worker threads use only C++ and caller-owned numeric
 * buffers. All MATLAB/MEX API calls happen on the calling MATLAB thread.
 */
#include "mex.h"
#include <algorithm>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <string>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

class Team {
    std::vector<std::thread> workers;
    std::mutex mutex;
    std::condition_variable ready, done;
    std::function<void(size_t,size_t)> job;
    size_t generation=0, count=0, pending=0;
    bool stopping=false;
public:
    explicit Team(size_t n) {
        try {
        for(size_t i=0;i<n;++i) workers.emplace_back([this,i,n]() {
            size_t seen=0;
            for(;;) {
                std::unique_lock<std::mutex> lock(mutex);
                ready.wait(lock,[&](){return stopping || generation!=seen;});
                if(stopping) return;
                seen=generation;
                const size_t begin=count*i/n, end=count*(i+1)/n;
                lock.unlock(); job(begin,end); lock.lock();
                if(--pending==0) done.notify_one();
            }
        });
        } catch(...) {
            {std::lock_guard<std::mutex> lock(mutex); stopping=true;}
            ready.notify_all();
            for(auto &worker:workers) worker.join();
            throw;
        }
    }
    ~Team() {
        {std::lock_guard<std::mutex> lock(mutex); stopping=true;}
        ready.notify_all();
        for(auto &worker:workers) worker.join();
    }
    size_t size() const {return workers.size();}
    void run(size_t n,const std::function<void(size_t,size_t)> &f) {
        if(workers.size()==1) {f(0,n); return;}
        std::unique_lock<std::mutex> lock(mutex);
        job=f; count=n; pending=workers.size(); ++generation;
        ready.notify_all(); done.wait(lock,[&](){return pending==0;});
        job=nullptr;
    }
};
static Team *team=nullptr;
static bool registered=false;
static void release() {delete team; team=nullptr;}
static unsigned cores() {return std::max(1u,std::thread::hardware_concurrency());}
struct RoutingError {const char *id; std::string message;};
static void require(bool good,const char *message) {
    if(!good) throw RoutingError{"OFM:NativeRoutingInput",message};
}
static double scalar(const mxArray *a) {
    require(mxIsDouble(a) && !mxIsComplex(a) && !mxIsSparse(a) && mxGetNumberOfElements(a)==1,"Expected a real double scalar.");
    return mxGetScalar(a);
}
static const double *grid(const mxArray *a,size_t n) {
    require(mxIsDouble(a) && !mxIsComplex(a) && !mxIsSparse(a) && mxGetNumberOfElements(a)==n,"Grid inputs must be full real double arrays of the same size.");
    return mxGetDoubles(a);
}
static double positive(double x) {return std::max(0.0,x);}

static void route(int nout,mxArray *output[],int nin,const mxArray *input[]) {
    if(nin==1 && mxIsChar(input[0])) {
        char command[16]; mxGetString(input[0],command,sizeof(command));
        require(std::string(command)=="info" && nout<=1,"Unknown native routing command.");
        if(nout) output[0]=mxCreateDoubleScalar(cores());
        return;
    }
    require(nin==10 && nout>=1 && nout<=2,"Use [rate,flux]=ofm_route(shape,dt,depression,dx,Manning,bed,depth,source,fixedDepth,threads).");
    const size_t rows=mxGetM(input[0]), cols=mxGetN(input[0]), n=rows*cols;
    require(rows>0 && cols>0,"Routing grids must not be empty.");
    const double *shape=grid(input[0],n), *dep=grid(input[2],n), *manning=grid(input[4],n);
    const double *bed=grid(input[5],n), *old=grid(input[6],n), *source=grid(input[7],n), *fixedDepth=grid(input[8],n);
    const double interval=scalar(input[1]), dx=scalar(input[3]), requested=scalar(input[9]);
    require(std::isfinite(interval) && interval>0 && std::isfinite(dx) && dx>0,"Timestep and cell size must be finite and positive.");
    require(std::isfinite(requested) && requested>=1 && requested==std::floor(requested),"Thread count must be a positive integer.");
    const unsigned threads=static_cast<unsigned>(std::min(requested,static_cast<double>(cores())));
    if(!registered) {mexAtExit(release); registered=true;}
    if(!team || team->size()!=threads) {
        release();
        try {team=new Team(threads);} catch(const std::exception &e) {
            throw RoutingError{"OFM:NativeRoutingThreads",std::string("Cannot create routing threads: ")+e.what()};
        }
    }
    std::vector<double> h(n), csum(n), outgoing(n);
    std::vector<uint8_t> active(n), fixed(n), roof(n);
    for(size_t i=0;i<n;++i) {
        active[i]=shape[i]!=-9999; fixed[i]=active[i] && (shape[i]==4 || shape[i]==5);
        roof[i]=active[i] && (shape[i]==999 || bed[i]==999);
        if(!active[i]) {h[i]=0; continue;}
        require(std::isfinite(bed[i]) && std::isfinite(dep[i]) && dep[i]>=0 && std::isfinite(manning[i]) && manning[i]>0,"Active cells need finite bed, nonnegative depression and positive Manning coefficients.");
        require(std::isfinite(old[i]) && std::isfinite(source[i]) && (!fixed[i] || (std::isfinite(fixedDepth[i]) && fixedDepth[i]>=0)),"Depths, sources and prescribed depths must be finite.");
        h[i]=positive(old[i]+interval*source[i]);
        if(fixed[i]) h[i]=fixedDepth[i];
    }
    const size_t ne=rows*(cols-1), ns=(rows-1)*cols;
    std::vector<double> qe(ne), ce(ne), qs(ns), cs(ns);
    output[0]=mxCreateDoubleMatrix(rows,cols,mxREAL);
    double *result=mxGetDoubles(output[0]), *flux=nullptr;
    if(nout==2) {
        mwSize dims[3]={rows,cols,2}; output[1]=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
        flux=mxGetDoubles(output[1]);
    }
    const double roofCoefficient=0.6*0.5*std::sqrt(2*9.81)/(dx*dx);
    auto face=[&](size_t a,size_t b,double &q,double &c) {
        q=0; c=0;
        if(!active[a] || !active[b]) return;
        if(roof[a]!=roof[b]) {
            const size_t donor=roof[a]?a:b;
            const double depth=positive(h[donor]-dep[donor]);
            q=(roof[a]?1:-1)*roofCoefficient*std::pow(depth,1.5);
            c=1.5*roofCoefficient*std::sqrt(depth);
            return;
        }
        const double head=(bed[a]+h[a])-(bed[b]+h[b]), difference=std::abs(head);
        if(difference<=1e-12) return;
        const size_t donor=head>=0?a:b;
        const double depth=positive(std::min(h[donor]-dep[donor],std::max(bed[a]+h[a],bed[b]+h[b])-std::max(bed[a],bed[b])));
        if(depth<=0) return;
        q=(head>=0?1:-1)*std::pow(depth,5.0/3.0)*std::sqrt(difference/dx)/(manning[donor]*dx);
        c=std::abs(q)/difference+(5.0/3.0)*std::abs(q)/depth;
    };
    double remaining=interval;
    size_t substeps=0;
    while(remaining>0) {
        if(++substeps>100000) throw RoutingError{"OFM:Substeps","Routing exceeded 100000 substeps. Reduce DELTA_T or inspect the input."};
        team->run(ne+ns,[&](size_t begin,size_t end) {
            for(size_t k=begin;k<end;++k) {
                if(k<ne) face(k,k+rows,qe[k],ce[k]);
                else {
                    const size_t j=k-ne, row=j%(rows-1), col=j/(rows-1), a=row+col*rows;
                    face(a,a+1,qs[j],cs[j]);
                }
            }
        });
        team->run(n,[&](size_t begin,size_t end) {
            for(size_t i=begin;i<end;++i) {
                const size_t row=i%rows, col=i/rows;
                double c=0, out=0;
                // Same east/west/south/north accumulation order as explicit.m.
                if(col+1<cols) {c+=ce[i]; out+=positive(qe[i]);}
                if(col>0) {c+=ce[i-rows]; out+=positive(-qe[i-rows]);}
                if(row+1<rows) {const size_t s=row+col*(rows-1); c+=cs[s]; out+=positive(qs[s]);}
                if(row>0) {const size_t s=row-1+col*(rows-1); c+=cs[s]; out+=positive(-qs[s]);}
                csum[i]=c; outgoing[i]=out;
            }
        });
        double maxC=0, availableTime=remaining;
        for(size_t i=0;i<n;++i) {
            maxC=std::max(maxC,csum[i]);
            if(active[i] && !fixed[i] && outgoing[i]>0)
                availableTime=std::min(availableTime,0.9*positive(h[i]-dep[i])/outgoing[i]);
        }
        double dt=std::min(remaining,availableTime);
        if(maxC>0) dt=std::min(dt,0.4/maxC);
        if(!std::isfinite(dt) || dt<=0 || remaining-dt==remaining)
            throw RoutingError{"OFM:SubstepUnderflow","Native routing cannot advance. Inspect depths, roughness and grid size."};
        team->run(n,[&](size_t begin,size_t end) {
            for(size_t i=begin;i<end;++i) {
                if(!active[i] || fixed[i]) continue;
                const size_t row=i%rows, col=i/rows;
                double change=0;
                if(col+1<cols) change-=qe[i];
                if(col>0) change+=qe[i-rows];
                if(row+1<rows) change-=qs[row+col*(rows-1)];
                if(row>0) change+=qs[row-1+col*(rows-1)];
                h[i]+=dt*change;
            }
        });
        if(flux) {
            team->run(ne+ns,[&](size_t begin,size_t end) {
                for(size_t k=begin;k<end;++k) {
                    if(k<ne) flux[k]+=dt*qe[k]/interval;
                    else {const size_t j=k-ne; flux[n+j%(rows-1)+(j/(rows-1))*rows]+=dt*qs[j]/interval;}
                }
            });
        }
        remaining=positive(remaining-dt);
    }
    for(size_t i=0;i<n;++i) result[i]=active[i]?(h[i]-old[i])/interval:0;
}

void mexFunction(int nout,mxArray *output[],int nin,const mxArray *input[]) {
    try {route(nout,output,nin,input);}
    catch(const RoutingError &e) {mexErrMsgIdAndTxt(e.id,"%s",e.message.c_str());}
    catch(const std::exception &e) {mexErrMsgIdAndTxt("OFM:NativeRouting","%s",e.what());}
}

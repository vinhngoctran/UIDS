/* Apple Metal routing backend, called on the MATLAB client. The GPU stores
 * water in compensated float pairs; CPU double checks verify the budget.
 * MATLAB reference: explicit.m, diffusive_face.m; CPU backend: ofm_route.cpp.
 */
#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#include "mex.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <string>
#include "ofm_metal_kernels.hpp"

struct Pair {float hi,lo;};
struct Params {uint32_t rows,cols,n,ne,ns,writeFlux; float dx,roofCoefficient; Pair dt,weight;};
static Pair split(double x) {float hi=static_cast<float>(x); return {hi,static_cast<float>(x-static_cast<double>(hi))};}
static double value(Pair x) {return static_cast<double>(x.hi)+x.lo;}
struct Error {const char *id; std::string message;};
static void need(bool condition,const char *message) {
    if(!condition) throw Error{"OFM:MetalInput",message};
}
static void checkBuffer(id<MTLBuffer> buffer) {
    if(!buffer) throw Error{"OFM:MetalMemory","Metal could not allocate routing buffers. Reduce the grid size or use CPU routing."};
}
static id<MTLDevice> device;
static id<MTLCommandQueue> queue;
static id<MTLComputePipelineState> facePipeline,limitPipeline,updatePipeline;
static bool registered=false;
struct Sum {
    double sum=0,correction=0;
    void add(double value) {double y=value-correction,t=sum+y;correction=(t-sum)-y;sum=t;}
};
static void release() {facePipeline=nil; limitPipeline=nil; updatePipeline=nil; queue=nil; device=nil;}
static std::string message(NSError *e) {return e?std::string([[e localizedDescription] UTF8String]):"Unknown Metal error.";}
static void initialize() {
    if(facePipeline && limitPipeline && updatePipeline) return;
    if(!registered) {mexAtExit(release); registered=true;}
    device=MTLCreateSystemDefaultDevice();
    if(!device || ![device hasUnifiedMemory]) throw Error{"OFM:MetalDevice","This backend requires an Apple silicon Metal GPU with unified memory."};
    queue=[device newCommandQueue];
    if(!queue) throw Error{"OFM:MetalDevice","Cannot create a Metal command queue."};
    MTLCompileOptions *options=[MTLCompileOptions new]; options.fastMathEnabled=NO;
    NSError *error=nil;
    NSString *source=[NSString stringWithUTF8String:ofmMetalSource];
    id<MTLLibrary> library=[device newLibraryWithSource:source options:options error:&error];
    if(!library) throw Error{"OFM:MetalCompile",message(error)};
    facePipeline=[device newComputePipelineStateWithFunction:[library newFunctionWithName:@"faces"] error:&error];
    if(!facePipeline) throw Error{"OFM:MetalCompile",message(error)};
    limitPipeline=[device newComputePipelineStateWithFunction:[library newFunctionWithName:@"limits"] error:&error];
    if(!limitPipeline) throw Error{"OFM:MetalCompile",message(error)};
    updatePipeline=[device newComputePipelineStateWithFunction:[library newFunctionWithName:@"update"] error:&error];
    if(!updatePipeline) throw Error{"OFM:MetalCompile",message(error)};
}
static id<MTLBuffer> buffer(size_t bytes) {
    id<MTLBuffer> b=[device newBufferWithLength:std::max(size_t(8),bytes) options:MTLResourceStorageModeShared];
    checkBuffer(b); return b;
}
static const double *grid(const mxArray *a,size_t n) {
    need(mxIsDouble(a) && !mxIsComplex(a) && !mxIsSparse(a) && mxGetNumberOfElements(a)==n,"Use full real double grids of the same size.");
    return mxGetDoubles(a);
}
static double scalar(const mxArray *a) {need(mxGetNumberOfElements(a)==1,"Expected scalar input."); return *grid(a,1);}
static void dispatch(id<MTLComputeCommandEncoder> encoder,id<MTLComputePipelineState> pipeline,size_t n,size_t width) {
    [encoder setComputePipelineState:pipeline];
    [encoder dispatchThreadgroups:MTLSizeMake((n+width-1)/width,1,1) threadsPerThreadgroup:MTLSizeMake(width,1,1)];
}
static id<MTLCommandBuffer> newCommand() {
    id<MTLCommandBuffer> command=[queue commandBuffer];
    if(!command) throw Error{"OFM:MetalMemory","Cannot allocate a Metal command buffer."};
    return command;
}
static id<MTLComputeCommandEncoder> newEncoder(id<MTLCommandBuffer> command) {
    id<MTLComputeCommandEncoder> encoder=[command computeCommandEncoder];
    if(!encoder) throw Error{"OFM:MetalMemory","Cannot allocate a Metal compute encoder."};
    return encoder;
}
static void complete(id<MTLCommandBuffer> command) {
    [command commit]; [command waitUntilCompleted];
    if(command.status!=MTLCommandBufferStatusCompleted) throw Error{"OFM:MetalExecution",message(command.error)};
}
static void route(int nout,mxArray *out[],int nin,const mxArray *in[]) {
    if(nin==1 && mxIsChar(in[0])) {
        char command[16]; mxGetString(in[0],command,sizeof(command));
        need(std::string(command)=="info" && nout<=1,"Unknown Metal command.");
        initialize();
        if(nout) {
            const char *fields[]={"Name","Available","Precision"};
            out[0]=mxCreateStructMatrix(1,1,3,fields);
            mxSetField(out[0],0,"Name",mxCreateString([[device name] UTF8String]));
            mxSetField(out[0],0,"Available",mxCreateLogicalScalar(true));
            mxSetField(out[0],0,"Precision",mxCreateString("Compensated FP32 storage; FP32 hydraulic powers"));
        }
        return;
    }
    need(nin==9 && nout>=1 && nout<=3,"Use [rate,flux,diagnostics]=ofm_metal(shape,dt,depression,dx,Manning,bed,depth,source,fixedDepth).");
    const size_t rows=mxGetM(in[0]),cols=mxGetN(in[0]),n=rows*cols,ne=rows*(cols-1),ns=(rows-1)*cols;
    need(rows>0 && cols>0 && n<=UINT32_MAX && ne+ns<=UINT32_MAX,"GPU grid is empty or exceeds 32-bit indexing.");
    const double *shape=grid(in[0],n),*dep=grid(in[2],n),*rough=grid(in[4],n),*bed=grid(in[5],n);
    const double *old=grid(in[6],n),*source=grid(in[7],n),*fixed=grid(in[8],n);
    const double interval=scalar(in[1]),dx=scalar(in[3]);
    need(std::isfinite(interval) && interval>0 && std::isfinite(dx) && dx>0,"Timestep and cell size must be positive and finite.");
    initialize();
    id<MTLBuffer> hBuf=buffer(n*sizeof(Pair)),depBuf=buffer(n*sizeof(Pair));
    id<MTLBuffer> dzBuf=buffer((ne+ns)*sizeof(Pair)),roughBuf=buffer(n*sizeof(float));
    id<MTLBuffer> flagsBuf=buffer(n*sizeof(uint32_t)),qcBuf=buffer((ne+ns)*2*sizeof(float));
    id<MTLBuffer> boundaryBuf=buffer(n*sizeof(Pair)),fluxBuf=buffer((nout>=2?2*n:1)*sizeof(Pair));
    const size_t maxWidth=std::min({size_t(256),static_cast<size_t>(limitPipeline.maxTotalThreadsPerThreadgroup),
        static_cast<size_t>(facePipeline.maxTotalThreadsPerThreadgroup),static_cast<size_t>(updatePipeline.maxTotalThreadsPerThreadgroup)});
    size_t width=1; while(2*width<=maxWidth) width*=2;
    need(width>=32,"Unsupported GPU reduction width.");
    const size_t groups=(n+width-1)/width;
    id<MTLBuffer> partialBuf=buffer(groups*2*sizeof(float));
    Pair *h=static_cast<Pair *>(hBuf.contents),*d=static_cast<Pair *>(depBuf.contents),*dz=static_cast<Pair *>(dzBuf.contents);
    uint32_t *flags=static_cast<uint32_t *>(flagsBuf.contents);
    float *r=static_cast<float *>(roughBuf.contents);
    std::memset(boundaryBuf.contents,0,n*sizeof(Pair));
    if(nout>=2) std::memset(fluxBuf.contents,0,2*n*sizeof(Pair));
    Sum applied;
    for(size_t i=0;i<n;++i) {
        flags[i]=shape[i]!=-9999?1:0;
        if(!flags[i]) {h[i]={0,0}; d[i]={0,0}; r[i]=1; continue;}
        need(std::isfinite(bed[i]) && std::isfinite(dep[i]) && dep[i]>=0 && std::isfinite(rough[i]) && rough[i]>0 &&
             std::isfinite(old[i]) && std::isfinite(source[i]),"Invalid active-cell bed, storage, roughness or forcing.");
        if(shape[i]==4 || shape[i]==5) flags[i]|=2;
        if(shape[i]==999 || bed[i]==999) flags[i]|=4;
        double water=std::max(0.0,old[i]+interval*source[i]);
        if(flags[i]&2) {need(std::isfinite(fixed[i]) && fixed[i]>=0,"Invalid prescribed boundary depth."); water=fixed[i];}
        h[i]=split(water); d[i]=split(dep[i]); r[i]=static_cast<float>(rough[i]);
        need(std::isfinite(h[i].hi) && std::isfinite(d[i].hi) && std::isfinite(r[i]) && r[i]>0,"Inputs exceed Metal floating-point range.");
        applied.add(value(h[i]));
    }
    for(size_t k=0;k<ne+ns;++k) {
        size_t a,b;
        if(k<ne) {a=k;b=a+rows;}
        else {size_t j=k-ne; a=j%(rows-1)+(j/(rows-1))*rows;b=a+1;}
        dz[k]=(flags[a]&1) && (flags[b]&1)?split(bed[a]-bed[b]):Pair{0,0};
        need(std::isfinite(dz[k].hi),"Terrain differences exceed Metal floating-point range.");
    }
    Params p={static_cast<uint32_t>(rows),static_cast<uint32_t>(cols),static_cast<uint32_t>(n),
        static_cast<uint32_t>(ne),static_cast<uint32_t>(ns),static_cast<uint32_t>(nout>=2),
        static_cast<float>(dx),static_cast<float>(0.6*0.5*std::sqrt(2*9.81)/(dx*dx)),{0,0},{0,0}};
    need(std::isfinite(p.dx) && p.dx>0 && std::isfinite(p.roofCoefficient),"Cell size exceeds Metal floating-point range; select CPU routing.");
    auto encodeLimits=[&](id<MTLCommandBuffer> command) {
        if(ne+ns) {
            id<MTLComputeCommandEncoder> e=newEncoder(command);
            [e setBuffer:hBuf offset:0 atIndex:0];[e setBuffer:depBuf offset:0 atIndex:1];
            [e setBuffer:dzBuf offset:0 atIndex:2];[e setBuffer:roughBuf offset:0 atIndex:3];
            [e setBuffer:flagsBuf offset:0 atIndex:4];[e setBuffer:qcBuf offset:0 atIndex:5];
            [e setBytes:&p length:sizeof(p) atIndex:6]; dispatch(e,facePipeline,ne+ns,width);[e endEncoding];
        }
        id<MTLComputeCommandEncoder> e=newEncoder(command);
        [e setBuffer:hBuf offset:0 atIndex:0];[e setBuffer:depBuf offset:0 atIndex:1];
        [e setBuffer:flagsBuf offset:0 atIndex:2];[e setBuffer:qcBuf offset:0 atIndex:3];
        [e setBuffer:partialBuf offset:0 atIndex:4];[e setBytes:&p length:sizeof(p) atIndex:5];
        dispatch(e,limitPipeline,n,width);[e endEncoding];
    };
    id<MTLCommandBuffer> first=newCommand(); encodeLimits(first); complete(first);
    double remaining=interval; size_t steps=0;
    while(remaining>0) {
        if(++steps>100000) throw Error{"OFM:Substeps","GPU routing exceeded 100000 substeps. Reduce DELTA_T or check the input."};
        const float *limits=static_cast<float *>(partialBuf.contents);
        double cmax=0,dt=remaining;
        for(size_t j=0;j<groups;++j) {
            if(!std::isfinite(limits[2*j]) || std::isnan(limits[2*j+1])) throw Error{"OFM:MetalState","GPU returned nonfinite hydraulic limits."};
            cmax=std::max(cmax,static_cast<double>(limits[2*j])); dt=std::min(dt,static_cast<double>(limits[2*j+1]));
        }
        if(cmax>0) dt=std::min(dt,.4/cmax);
        if(!std::isfinite(dt) || dt<=0 || remaining-dt==remaining) throw Error{"OFM:SubstepUnderflow","GPU routing cannot advance; check input or select CPU routing."};
        p.dt=split(dt); p.weight=split(dt/interval);
        if(p.dt.hi<=0) throw Error{"OFM:SubstepUnderflow","GPU timestep is below its numerical range; select CPU routing."};
        id<MTLCommandBuffer> command=newCommand();
        id<MTLComputeCommandEncoder> e=newEncoder(command);
        [e setBuffer:hBuf offset:0 atIndex:0];[e setBuffer:flagsBuf offset:0 atIndex:1];
        [e setBuffer:qcBuf offset:0 atIndex:2];[e setBuffer:boundaryBuf offset:0 atIndex:3];
        [e setBuffer:fluxBuf offset:0 atIndex:4];[e setBytes:&p length:sizeof(p) atIndex:5];
        dispatch(e,updatePipeline,n,width);[e endEncoding];
        remaining=std::max(0.0,remaining-dt);
        if(remaining>0) encodeLimits(command);
        complete(command);
    }
    out[0]=mxCreateDoubleMatrix(rows,cols,mxREAL); double *result=mxGetDoubles(out[0]);
    const Pair *boundary=static_cast<Pair *>(boundaryBuf.contents);
    Sum final,boundaryTotal; double minDepth=INFINITY;
    for(size_t i=0;i<n;++i) {
        if(!(flags[i]&1)) {result[i]=0;continue;}
        double water=value(h[i]);
        if(!std::isfinite(water) || water< -1e-10) throw Error{"OFM:MetalState","GPU produced an invalid surface depth."};
        result[i]=(water-old[i])/interval;
        final.add(water); boundaryTotal.add(value(boundary[i])); minDepth=std::min(minDepth,water);
    }
    if(!std::isfinite(minDepth)) minDepth=0;
    const double residual=final.sum-applied.sum-boundaryTotal.sum;
    if(std::abs(residual)>1e-10*std::max(1.0,std::abs(applied.sum)))
        throw Error{"OFM:MetalBalance","GPU surface water budget failed. Use CPU routing for this case."};
    if(nout>=2) {
        mwSize dims[3]={rows,cols,2};out[1]=mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
        const Pair *f=static_cast<Pair *>(fluxBuf.contents);double *destination=mxGetDoubles(out[1]);
        for(size_t i=0;i<2*n;++i) destination[i]=value(f[i]);
    }
    if(nout>=3) {
        const char *fields[]={"Substeps","BalanceResidual","BoundaryDepth","MinimumDepth","Precision","Device"};
        out[2]=mxCreateStructMatrix(1,1,6,fields);
        mxSetField(out[2],0,"Substeps",mxCreateDoubleScalar(steps));
        mxSetField(out[2],0,"BalanceResidual",mxCreateDoubleScalar(residual));
        mxSetField(out[2],0,"BoundaryDepth",mxCreateDoubleScalar(boundaryTotal.sum));
        mxSetField(out[2],0,"MinimumDepth",mxCreateDoubleScalar(minDepth));
        mxSetField(out[2],0,"Precision",mxCreateString("Compensated FP32 storage; FP32 hydraulic powers"));
        mxSetField(out[2],0,"Device",mxCreateString([[device name] UTF8String]));
    }
}
void mexFunction(int nout,mxArray *out[],int nin,const mxArray *in[]) {
    @autoreleasepool {
        try {route(nout,out,nin,in);}
        catch(const Error &e) {mexErrMsgIdAndTxt(e.id,"%s",e.message.c_str());}
        catch(const std::exception &e) {mexErrMsgIdAndTxt("OFM:Metal","%s",e.what());}
    }
}

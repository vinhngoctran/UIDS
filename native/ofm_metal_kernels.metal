#include <metal_stdlib>
using namespace metal;
// Compensated two-float storage, about 48 significant bits. Hydraulic
// powers use FP32; additions/transfers retain the low storage component.
typedef float2 ds;
inline ds norm(float a,float b) {float s=a+b; return ds(s,b-(s-a));}
inline ds add(ds a,ds b) {
    float s=a.x+b.x, v=s-a.x;
    float e=(a.x-(s-v))+(b.x-v)+a.y+b.y;
    return norm(s,e);
}
inline ds neg(ds a) {return -a;}
inline ds sub(ds a,ds b) {return add(a,neg(b));}
inline ds mul(ds a,ds b) {
    float p=a.x*b.x;
    float e=fma(a.x,b.x,-p)+a.x*b.y+a.y*b.x+a.y*b.y;
    return norm(p,e);
}
inline bool less(ds a,ds b) {return a.x<b.x || (a.x==b.x && a.y<b.y);}
inline ds minimum(ds a,ds b) {return less(a,b)?a:b;}
inline ds maximum(ds a,ds b) {return less(a,b)?b:a;}
inline ds pos(ds a) {return maximum(ds(0),a);}
inline float val(ds a) {return a.x+a.y;}
struct Params {uint rows,cols,n,ne,ns,writeFlux; float dx,roofCoefficient; ds dt,weight;};
// Per-cell data: depression. Bed differences for faces are precomputed in
// double on the host, preventing cancellation at large terrain elevations.
kernel void faces(device const ds *h [[buffer(0)]],device const ds *dep [[buffer(1)]],
    device const ds *dz [[buffer(2)]],device const float *rough [[buffer(3)]],
    device const uint *flags [[buffer(4)]],device float2 *qc [[buffer(5)]],
    constant Params &p [[buffer(6)]],uint k [[thread_position_in_grid]]) {
    if(k>=p.ne+p.ns) return;
    uint a,b;
    if(k<p.ne) {a=k; b=a+p.rows;}
    else {uint j=k-p.ne; a=j%(p.rows-1)+(j/(p.rows-1))*p.rows; b=a+1;}
    if(!(flags[a]&1) || !(flags[b]&1)) {qc[k]=float2(0); return;}
    bool roofA=flags[a]&4, roofB=flags[b]&4;
    if(roofA!=roofB) {
        uint donor=roofA?a:b; float depth=val(pos(sub(h[donor],dep[donor])));
        qc[k]=float2((roofA?1.f:-1.f)*p.roofCoefficient*pow(depth,1.5f),
                     1.5f*p.roofCoefficient*sqrt(depth));
        return;
    }
    ds head=add(dz[k],sub(h[a],h[b]));
    float gradient=abs(val(head));
    if(gradient<=1e-12f) {qc[k]=float2(0); return;}
    bool forward=!less(head,ds(0)); uint donor=forward?a:b;
    ds aboveBed=maximum(add(h[a],minimum(ds(0),dz[k])),sub(h[b],maximum(ds(0),dz[k])));
    float depth=val(pos(minimum(sub(h[donor],dep[donor]),aboveBed)));
    if(depth<=0) {qc[k]=float2(0); return;}
    float q=(forward?1.f:-1.f)*pow(depth,5.f/3.f)*sqrt(gradient/p.dx)/(rough[donor]*p.dx);
    qc[k]=float2(q,abs(q)/gradient+(5.f/3.f)*abs(q)/depth);
}
kernel void limits(device const ds *h [[buffer(0)]],device const ds *dep [[buffer(1)]],
    device const uint *flags [[buffer(2)]],device const float2 *qc [[buffer(3)]],
    device float2 *partial [[buffer(4)]],constant Params &p [[buffer(5)]],
    uint i [[thread_position_in_grid]],uint lane [[thread_index_in_threadgroup]],
    uint group [[threadgroup_position_in_grid]],uint width [[threads_per_threadgroup]]) {
    threadgroup float2 shared[256];
    float c=0,out=0;
    if(i<p.n && (flags[i]&1)) {
        uint r=i%p.rows, col=i/p.rows;
        if(col+1<p.cols) {c+=qc[i].y; out+=max(0.f,qc[i].x);}
        if(col>0) {c+=qc[i-p.rows].y; out+=max(0.f,-qc[i-p.rows].x);}
        if(r+1<p.rows) {uint j=p.ne+r+col*(p.rows-1); c+=qc[j].y; out+=max(0.f,qc[j].x);}
        if(r>0) {uint j=p.ne+r-1+col*(p.rows-1); c+=qc[j].y; out+=max(0.f,-qc[j].x);}
    }
    float available=INFINITY;
    if(i<p.n && (flags[i]&1) && !(flags[i]&2) && out>0)
        available=.9f*val(pos(sub(h[i],dep[i])))/out;
    shared[lane]=float2(c,available);
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for(uint stride=width/2;stride>0;stride/=2) {
        if(lane<stride) shared[lane]=float2(max(shared[lane].x,shared[lane+stride].x),
                                           min(shared[lane].y,shared[lane+stride].y));
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    if(lane==0) partial[group]=shared[0];
}
kernel void update(device ds *h [[buffer(0)]],device const uint *flags [[buffer(1)]],
    device const float2 *qc [[buffer(2)]],device ds *boundary [[buffer(3)]],
    device ds *flux [[buffer(4)]],constant Params &p [[buffer(5)]],
    uint i [[thread_position_in_grid]]) {
    if(i>=p.n) return;
    uint r=i%p.rows, col=i/p.rows;
    if(flags[i]&1) {
        ds change=ds(0);
        if(col+1<p.cols) change=sub(change,ds(qc[i].x,0));
        if(col>0) change=add(change,ds(qc[i-p.rows].x,0));
        if(r+1<p.rows) change=sub(change,ds(qc[p.ne+r+col*(p.rows-1)].x,0));
        if(r>0) change=add(change,ds(qc[p.ne+r-1+col*(p.rows-1)].x,0));
        ds transfer=mul(p.dt,change);
        if(flags[i]&2) boundary[i]=sub(boundary[i],transfer);
        else h[i]=add(h[i],transfer);
    }
    if(p.writeFlux) {
        if(col+1<p.cols) flux[i]=add(flux[i],mul(p.weight,ds(qc[i].x,0)));
        if(r+1<p.rows) flux[p.n+i]=add(flux[p.n+i],mul(p.weight,ds(qc[p.ne+r+col*(p.rows-1)].x,0)));
    }
}

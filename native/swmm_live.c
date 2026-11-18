/* Direct MATLAB/EPA SWMM coupling. EPA's engine is kept open for the run.
 * No simulation results or inflow series are saved/read by advance().
 * Including the pinned engine main module exposes its private deadline so
 * each adaptive hydraulic step stops exactly at an OFM synchronization time.
 */
#include "mex.h"
#include "swmm5.c"

static int opened = 0, started = 0, registered = 0;
static int *nodes = NULL, *links = NULL, *outfalls = NULL;
static mwSize nn = 0, nl = 0, no = 0;

static void release(void)
{
    if (started) { swmm_end(); started = 0; }
    if (opened) { swmm_close(); opened = 0; if (mexIsLocked()) mexUnlock(); }
    if (nodes) mxFree(nodes);
    if (links) mxFree(links);
    if (outfalls) mxFree(outfalls);
    nodes = links = outfalls = NULL;
    nn = nl = no = 0;
}

static void fail(const char *id, const char *message)
{
    release();
    mexErrMsgIdAndTxt(id,"%s",message);
}

static void check(int code, const char *operation)
{
    if (code) {
        char detail[1024], message[1200];
        swmm_getError(detail,sizeof(detail));
        snprintf(message,sizeof(message),"SWMM %s failed (%d): %s",operation,code,detail);
        fail("OFM:SWMM",message);
    }
}

static int *indices(const mxArray *ids, int type)
{
    mwSize k, count = mxGetNumberOfElements(ids);
    int *result = (int *)mxCalloc(count ? count : 1,sizeof(int));
    mexMakeMemoryPersistent(result);
    for (k=0; k<count; ++k) {
        char *id = mxArrayToString(mxGetCell(ids,k));
        mwSize j;
        result[k] = swmm_getIndex(type,id);
        if (result[k] < 0) {
            char message[512];
            snprintf(message,sizeof(message),"SWMM object '%s' was not found.",id);
            mxFree(id); mxFree(result);
            fail("OFM:SWMMObject",message);
        }
        mxFree(id);
        for (j=0; j<k; ++j) if (result[j] == result[k]) {
            mxFree(result);
            fail("OFM:SWMMObject","Duplicate SWMM identifiers in an exchange vector.");
        }
    }
    return result;
}

static mxArray *vector(mwSize count)
{
    return mxCreateDoubleMatrix(1,count,mxREAL);
}

static mxArray *state(const double *flood, const double *outVolume)
{
    const char *names[] = {"TimeSeconds","NodeDepth","NodeHead","NodeInflow",
        "NodeOverflow","NodeFloodVolume","LinkFlow","OutFlow","OutDepth",
        "OutVolume","SewerVolume","EngineVersion"};
    mxArray *s = mxCreateStructMatrix(1,1,12,names);
    mxArray *depth=vector(nn), *head=vector(nn), *inflow=vector(nn), *overflow=vector(nn);
    mxArray *fv=vector(nn), *flow=vector(nl), *out=vector(no), *outDepth=vector(no), *ov=vector(no);
    mwSize k;
    double volume = 0;
    for (k=0; k<nn; ++k) {
        mxGetDoubles(depth)[k]=swmm_getValue(swmm_NODE_DEPTH,nodes[k]);
        mxGetDoubles(head)[k]=swmm_getValue(swmm_NODE_HEAD,nodes[k]);
        mxGetDoubles(inflow)[k]=swmm_getValue(swmm_NODE_INFLOW,nodes[k]);
        mxGetDoubles(overflow)[k]=swmm_getValue(swmm_NODE_OVERFLOW,nodes[k]);
        mxGetDoubles(fv)[k]=flood ? flood[k] : 0;
    }
    for (k=0; k<nl; ++k) mxGetDoubles(flow)[k]=swmm_getValue(swmm_LINK_FLOW,links[k]);
    for (k=0; k<no; ++k) {
        mxGetDoubles(out)[k]=swmm_getValue(swmm_NODE_INFLOW,outfalls[k]);
        mxGetDoubles(outDepth)[k]=swmm_getValue(swmm_NODE_DEPTH,outfalls[k]);
        mxGetDoubles(ov)[k]=outVolume ? outVolume[k] : 0;
    }
    for (k=0; k<(mwSize)Nobjects[NODE]; ++k) volume += Node[k].newVolume;
    for (k=0; k<(mwSize)Nobjects[LINK]; ++k) volume += Link[k].newVolume;
    mxSetField(s,0,"TimeSeconds",mxCreateDoubleScalar(NewRoutingTime/1000.0));
    mxSetField(s,0,"NodeDepth",depth); mxSetField(s,0,"NodeHead",head);
    mxSetField(s,0,"NodeInflow",inflow); mxSetField(s,0,"NodeOverflow",overflow);
    mxSetField(s,0,"NodeFloodVolume",fv); mxSetField(s,0,"LinkFlow",flow);
    mxSetField(s,0,"OutFlow",out); mxSetField(s,0,"OutDepth",outDepth);
    mxSetField(s,0,"OutVolume",ov);
    mxSetField(s,0,"SewerVolume",mxCreateDoubleScalar(volume*UCF(VOLUME)));
    mxSetField(s,0,"EngineVersion",mxCreateDoubleScalar(swmm_getVersion()));
    return s;
}

static void validate_ids(const mxArray *ids)
{
    mwSize k;
    if (!mxIsCell(ids)) mexErrMsgIdAndTxt("OFM:SWMMArgument","Identifiers must be cell arrays of character vectors.");
    for (k=0; k<mxGetNumberOfElements(ids); ++k)
        if (!mxGetCell(ids,k) || !mxIsChar(mxGetCell(ids,k)))
            mexErrMsgIdAndTxt("OFM:SWMMArgument","Identifiers must be character vectors.");
}

void mexFunction(int nout, mxArray *output[], int nin, const mxArray *input[])
{
    char command[32];
    mwSize k;
    if (nin < 1 || !mxIsChar(input[0]) || nout > 1)
        mexErrMsgIdAndTxt("OFM:SWMMArgument","Expected a command and at most one output.");
    mxGetString(input[0],command,sizeof(command));
    if (!registered) { mexAtExit(release); registered=1; }
    if (strcmp(command,"close")==0) {
        release();
        return;
    }
    if (strcmp(command,"open")==0) {
        char *path, *report;
        if (opened) mexErrMsgIdAndTxt("OFM:SWMMSession","A SWMM session is already open.");
        if (nin != 6 || !mxIsChar(input[1]) || !mxIsChar(input[2]))
            mexErrMsgIdAndTxt("OFM:SWMMArgument","open expects input, report, node, link and outfall identifiers.");
        for (k=3; k<6; ++k) validate_ids(input[k]);
        path=mxArrayToString(input[1]); report=mxArrayToString(input[2]);
        opened=1; mexLock();
        { int code=swmm_open(path,report,""); mxFree(path); mxFree(report); check(code,"open"); }
        if (FlowUnits != CMS) fail("OFM:SWMMUnits","SWMM must use FLOW_UNITS CMS.");
        nn=mxGetNumberOfElements(input[3]); nl=mxGetNumberOfElements(input[4]); no=mxGetNumberOfElements(input[5]);
        nodes=indices(input[3],swmm_NODE);
        links=indices(input[4],swmm_LINK);
        outfalls=indices(input[5],swmm_NODE);
        for (k=0; k<nn; ++k) {
            if (Node[nodes[k]].type == OUTFALL)
                fail("OFM:SWMMObject","An exchange node cannot be an outfall.");
            /* OFM stores surface ponding. SWMM must release flooding once. */
            Node[nodes[k]].pondedArea=0;
        }
        for (k=0; k<no; ++k) if (Node[outfalls[k]].type != OUTFALL)
            fail("OFM:SWMMObject","An outfall identifier refers to another node type.");
        swmm_setValue(swmm_NOREPORT,0,1);
        started=1; check(swmm_start(0),"start");
        if (nout) output[0]=state(NULL,NULL);
        return;
    }
    if (!opened || !started) mexErrMsgIdAndTxt("OFM:SWMMSession","SWMM is not open.");
    if (strcmp(command,"advance")==0) {
        double target, *q, *flood, *outVolume, deadline, previous;
        const double tolerance=1e-6; /* milliseconds */
        if (nin != 3 || !mxIsDouble(input[1]) || mxIsSparse(input[1]) || mxIsComplex(input[1]) || mxGetNumberOfElements(input[1]) != 1 ||
            !mxIsDouble(input[2]) || mxIsSparse(input[2]) || mxIsComplex(input[2]) || mxGetNumberOfElements(input[2]) != nn)
            mexErrMsgIdAndTxt("OFM:SWMMArgument","advance expects target seconds and one double inflow per node.");
        target=mxGetScalar(input[1])*1000.0;
        if (!mxIsFinite(target) || target <= NewRoutingTime || target > TotalDuration+tolerance)
            mexErrMsgIdAndTxt("OFM:SWMMTime","Target time must advance within the simulation duration.");
        q=mxGetDoubles(input[2]);
        for (k=0; k<nn; ++k) if (!mxIsFinite(q[k]) || q[k]<0)
            mexErrMsgIdAndTxt("OFM:SWMMInflow","Surface inflows must be finite and nonnegative.");
        for (k=0; k<nn; ++k) swmm_setValue(swmm_NODE_LATFLOW,nodes[k],q[k]);
        flood=(double *)mxCalloc(nn ? nn : 1,sizeof(double));
        outVolume=(double *)mxCalloc(no ? no : 1,sizeof(double));
        deadline=RoutingDuration;
        RoutingDuration=MIN(target,TotalDuration);
        while (NewRoutingTime < RoutingDuration-tolerance) {
            double elapsed, seconds;
            previous=NewRoutingTime;
            check(swmm_step(&elapsed),"step");
            seconds=(NewRoutingTime-previous)/1000.0;
            if (seconds <= 0) fail("OFM:SWMMTime","SWMM did not advance to the synchronization time.");
            for (k=0; k<nn; ++k) flood[k] += MAX(0,Node[nodes[k]].overflow)*UCF(FLOW)*seconds;
            for (k=0; k<no; ++k) outVolume[k] += Node[outfalls[k]].inflow*UCF(FLOW)*seconds;
        }
        RoutingDuration=deadline;
        if (fabs(NewRoutingTime-target)>tolerance)
            fail("OFM:SWMMTime","SWMM and OFM clocks differ at synchronization.");
        if (nout) output[0]=state(flood,outVolume);
        mxFree(flood); mxFree(outVolume);
        return;
    }
    mexErrMsgIdAndTxt("OFM:SWMMArgument","Unknown SWMM command.");
}

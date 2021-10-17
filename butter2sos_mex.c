/*=============================*
 * BUTTER2SOS MEX API
 *=============================*/

#include "mex.h"
#include "matrix.h"
#include "butter2sos_design.h"

// To build: mex -I C:\...\dir butter2sos_mex.c butter2sos_design.c

void mexFunction(int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[])
{
    int order, type, idx, nbiquads, stgIdx, coefIdx;
    float fc, fs, *mat;
    double *buff;

    // indicate help if no arguments are provided for the function.
    if (!nrhs) {
        mexPrintf("***********************************************\n");
        mexPrintf("BUTTER2SOS_MEX is the BUTTER2SOS MATLAB API.\n");
        mexPrintf("\t sos = butter2sos_mex(N, Fc, Fs, type);\n");
        mexPrintf("\t type (scalar) indicates 0=LPF, 1=HPF, 2=APF.\n");
        mexPrintf("\t N (double) is the order.\n\t Fc (double) is the corner frequency.\n");
        mexPrintf("\t Fs (double) is the sampling rate.\n\t sos (double) is the SOS matrix.\n");
        mexPrintf("***********************************************\n");
        return;
    } else if ((nrhs < 4) || !mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3])) {
        mexErrMsgTxt("BUTTER2SOS_MEX requires 4 input arguments indicating order, Fc, Fs, and type.");
        return;
    }
    order = (int)(*mxGetPr(prhs[0]));
    fc = (regular_t)(*mxGetPr(prhs[1]));
    fs = (regular_t)(*mxGetPr(prhs[2]));
    type = (int)(*mxGetPr(prhs[3]));

    nbiquads = (int)ceil(order/2);
    if (order % 2) {
        nbiquads += 1;
    }

    // generate the SOS matrix.
    mat = (regular_t*)butter(order, fc, fs, type);
    if (mat[0] == NULL) {
        mexErrMsgTxt("BUTTER2SOS_MEX butter function failure. Unable to construct SOS matrix.\n");
        return;
    }


    buff = mxMalloc(N_SOSCOEFFS * nbiquads * sizeof(double)); // output buffer for MATLAB.
    
    // re-format the output buffer since MATLAB expects column major data.
    idx = 0;
    for (coefIdx = 0; coefIdx < N_SOSCOEFFS; coefIdx++) {
        for (stgIdx = 0; stgIdx < nbiquads; stgIdx++) {
            buff[idx++] = (double)mat[coefIdx+(N_SOSCOEFFS*stgIdx)];
        }
    }

    // free the buffer allocated by butter2sos.
    free(mat);

    // create the output for MATLAB [0 by 0, going to swap pointer with the allocated buffer]
    plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    if (plhs[0] == NULL) {
        mexErrMsgTxt("BUTTER2SOS_MEX could not allocate memory.\n");
        return;
    }

    // set pointer and indicate dimensions.
    mxSetPr(plhs[0], buff);
    mxSetM(plhs[0], biquads);
    mxSetN(plhs[0], N_SOSCOEFFS);

    // don't free buff since that points to the output.
    return;
}
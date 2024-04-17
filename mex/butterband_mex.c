/*=============================*
 * BUTTERBAND MEX API
 *=============================*/

#include "mex.h"
#include "matrix.h"
#include "butter2sos_design.h"

// To build: mex -I C:\...\dir butterband_mex.c butter2sos_design.c

void mexFunction(int nlhs, mxArray *plhs[], const int nrhs, const mxArray *prhs[])
{
    int order, idx, numbiquads, stgIdx, coefIdx, type;
    double fc1, fc2, fs, *mat;
    double *buff;

    // indicate help if no arguments are provided for the function.
    if (!nrhs) {
        mexPrintf("***********************************************\n");
        mexPrintf("BUTTERBAND_MEX is the BUTTERBAND MATLAB API.\n");
        mexPrintf("\t sos = butterband_mex(N, Flo, Fhi, Fs, type);\n");
        mexPrintf("\t type (scalar) indicates 0=bandpass, 1=bandstop.\n");
        mexPrintf("\t N (double) is the order.\n\t Fc (double) is the corner frequency.\n");
        mexPrintf("\t Flo (scalar) is the lower corner frequency.\n");
        mexPrintf("\t Fhi (scalar) is the upper corner frequency.\n");
        mexPrintf("\t Fs (double) is the sampling rate.\n\t sos (double) is the SOS matrix.\n");
        mexPrintf("***********************************************\n");
        return;
    } else if ((nrhs < 5) || !mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3]) || !mxIsDouble(prhs[4])) {
        mexErrMsgTxt("BUTTERBAND_MEX requires 5 input arguments indicating order, frequencies, Fs, and type.");
        return;
    }
    order = (int)(*mxGetPr(prhs[0]));
    fc1 = (real64_t)(*mxGetPr(prhs[1]));
    fc2 = (real64_t)(*mxGetPr(prhs[2]));
    fs = (real64_t)(*mxGetPr(prhs[3]));
    type = (int)(*mxGetPr(prhs[4]));
    
    if (order <= 0) {
        mexErrMsgTxt("BUTTERBAND_MEX requires a positive finite scalar order.\n");
        return;
    } else if ((fc1 > fs/2) || (fc2 > fs/2) || (fc1 <= 0) || (fc2 <= 0)) {
        mexErrMsgTxt("BUTTERBAND_MEX requires corner frequencies in (0, Nq] range.\n");
        return;
    }
    numbiquads = order; // same number of stages as the requested order.

    // generate the SOS matrix.
    mat = (real64_t*)butterband(order, fc1, fc2, fs, type);
    if (mat[0] == NULL) {
        mexErrMsgTxt("BUTTERBAND_MEX butter function failure. Unable to construct SOS matrix.\n");
        return;
    }

    buff = mxMalloc(N_SOSCOEFFS * numbiquads * sizeof(double)); // output buffer for MATLAB.
    
    // re-format the output buffer since MATLAB expects column major data.
    idx = 0;
    for (coefIdx = 0; coefIdx < N_SOSCOEFFS; coefIdx++) {
        for (stgIdx = 0; stgIdx < numbiquads; stgIdx++) {
            buff[idx++] = (double)mat[coefIdx+(N_SOSCOEFFS*stgIdx)];
        }
    }

    // free the buffer allocated by BUTTERBAND.
    free(mat);

    // create the output for MATLAB [0 by 0, going to swap pointer with the allocated buffer]
    plhs[0] = mxCreateNumericMatrix(0, 0, mxDOUBLE_CLASS, mxREAL);
    if (plhs[0] == NULL) {
        mexErrMsgTxt("BUTTERBAND_MEX could not allocate memory.\n");
        return;
    }

    // set pointer and indicate dimensions.
    mxSetPr(plhs[0], buff);
    mxSetM(plhs[0], numbiquads);
    mxSetN(plhs[0], N_SOSCOEFFS);

    // don't free buff since that points to the output.
    return;
}
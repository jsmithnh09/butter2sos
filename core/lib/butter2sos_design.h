/*
 * =======================================================================
 *  buttersos_design.h
 * =======================================================================
 */
 
#ifndef BUTTERSOS_DESIGN_H
#define BUTTERSOS_DESIGN_H

#include <stdio.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <complex.h>


typedef double real64_t;
typedef complex double complex64_t;

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

#define HALF_PI   (real64_t) M_PI / 2.0
#define TWO_PI    (real64_t) M_PI * 2.0
#define REG_PI    (real64_t) M_PI
#define N_SOSCOEFFS 6
#define BUTTER2SOS_DEBUG_PRINT printf("Line %d\n", __LINE__);

enum butter_filttype_e {FILT_LPF, FILT_HPF, FILT_APF};

// print the sos matrix to stdout
void printsosmatrix(const real64_t* matrix, const int nstages);
 
// principal design function (LPF/HPF/APF)
real64_t* butter(const int order, const real64_t fc, real64_t fs, const int type);

// band-based design function
real64_t* butterband(const int order, real64_t flo, real64_t fhi, real64_t fs, const int type);

// julia API for LPF/HPF/APF
void jl_butter(const int order, const double fc, const double fs, const int type, double* sos);

// julia API for bandpass/bandstop
void jl_butterband(const int order, const double flo, const double fhi, const double fs, const int type, double *mat);

#endif

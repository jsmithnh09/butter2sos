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


// linear space generation.
real64_t* linspace(const int start, const int stop, const int step, int *outsize);

// pole seed around the unit circle.
complex64_t* seedpoles(const int order, int* numpoles);

// print the sos matrix to stdout
void printsosmatrix(const real64_t* matrix, const int nstages);
 
// create the sos matrix
real64_t* mksosmatrix(const int order, const int type);

// print the complex array
void printcarray(const complex64_t* x, const int len);

// complex division
complex64_t compdiv(complex64_t x, complex64_t y);

// euler identity
complex64_t euler(real64_t x);

// complex subtraction
complex64_t compsub(complex64_t z1, complex64_t z2);

// complex addition
complex64_t compadd(complex64_t z1, complex64_t z2);

// complex square root
complex64_t compsqrt(complex64_t x);

// complex multiplication
complex64_t compmult(complex64_t x, complex64_t y);

// check if a singularity is real
int isreal(complex64_t x);

// pole sorting based on proximity to the unit circle.
void polesort(complex64_t* poles, int numpoles, int order);

// warp lowpass to bandpass prototype
void bpfwarp(complex64_t* poles, int* numpoles, complex64_t* zeros, int* numzeros, real64_t* gain, const real64_t* bwidth, const real64_t* Wn);

// warp lowpass to bandstop prototype
void bsfwarp(complex64_t* poles, int* numpoles, complex64_t* zeros, real64_t* gain, const real64_t* bwidth, const real64_t* Wn);

// warp lowpass to highpass prototype
void hpfwarp(complex64_t* poles, const int numpoles, real64_t* zeros, int* nzeros, real64_t* gain, const real64_t omega);

// warp lowpass to lowpass prototype
void lpfwarp(complex64_t* poles, int numpoles, real64_t* zeros, int* nzeros, real64_t* gain, real64_t omega);

// bilinear transform band-based singularities, (BPF/BSF.)
void bilinear_band_s2z(complex64_t* poles, const int* numpoles, complex64_t* zeros, int* numzeros, real64_t* gain, const real64_t* fs);

// bilinear transform in-place singularities.
void bilinear_s2z(complex64_t* poles, const int numpoles, const int numzeros, real64_t* gain, const real64_t fs);

// principal design function (LPF/HPF/APF)
real64_t* butter(const int order, const real64_t fc, real64_t fs, const int type);

// band-based design function
real64_t* butterband(const int order, real64_t flo, real64_t fhi, real64_t fs, const int type);

// julia API function
void jl_butter(const int order, const double fc, const double fs, const int type, double* sos);

// julia API for bandpass/bandstop
void jl_butterband(const int order, const double flo, const double fhi, const double fs, const int type, double *mat);

#endif

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


typedef float regular_t;
typedef complex float complex_t;

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

#define BUTTER2SOS_PRECISION 32
#define HALF_PI   (regular_t) M_PI / 2.0
#define TWO_PI    (regular_t) M_PI * 2.0
#define REG_PI    (regular_t) M_PI
#define N_SOSCOEFFS 6
#define BUTTER2SOS_DEBUG_PRINT printf("Line %d\n", __LINE__);


// linear space generation.
regular_t* linspace(const int start, const int stop, const int step, int *outsize);

// pole seed around the unit circle.
complex_t* seedpoles(const int order, int* numpoles);

// print the sos matrix to stdout
void printsosmatrix(const regular_t* matrix, const int nstages);
 
// create the sos matrix
regular_t* mksosmatrix(const int order, const int type);

// print the complex array
void printcarray(const complex_t* x, const int len);

// complex division
complex_t compdiv(complex_t x, complex_t y);

// euler identity
complex_t euler(regular_t x);

// complex subtraction
complex_t compsub(complex_t z1, complex_t z2);

// complex addition
complex_t compadd(complex_t z1, complex_t z2);

// complex square root
complex_t compsqrt(complex_t x);

// complex multiplication
complex_t compmult(complex_t x, complex_t y);

// check if a singularity is real
int isreal(complex_t x);

// pole sorting based on proximity to the unit circle.
void polesort(complex_t* poles, int numpoles, int order);

// warp lowpass to bandpass prototype
void bpfwarp(complex_t* poles, int* numpoles, complex_t* zeros, int* numzeros, regular_t* gain, const regular_t* bwidth, const regular_t* Wn);

// warp lowpass to bandstop prototype
void bsfwarp(complex_t* poles, int* numpoles, complex_t* zeros, regular_t* gain, const regular_t* bwidth, const regular_t* Wn);

// warp lowpass to highpass prototype
void hpfwarp(complex_t* poles, const int numpoles, regular_t* zeros, int* nzeros, regular_t* gain, const regular_t omega);

// warp lowpass to lowpass prototype
void lpfwarp(complex_t* poles, int numpoles, regular_t* zeros, int* nzeros, regular_t* gain, regular_t omega);

// bilinear transform band-based singularities, (BPF/BSF.)
void bilinear_band_s2z(complex_t* poles, const int* numpoles, complex_t* zeros, int* numzeros, regular_t* gain, const regular_t* fs);

// bilinear transform in-place singularities.
void bilinear_s2z(complex_t* poles, const int numpoles, const int numzeros, regular_t* gain, const regular_t fs);

// principal design function (LPF/HPF/APF)
regular_t* butter(const int order, const regular_t fc, regular_t fs, const int type);

// band-based design function
regular_t* butterband(const int order, regular_t flo, regular_t fhi, regular_t fs, const int type);

// julia API function
void jl_butter(const int order, const float fc, const float fs, const int type, float* sos);

// julia API for bandpass/bandstop
void jl_butterband(const int order, const float flo, const float fhi, const float fs, const int type, float *mat);

#endif

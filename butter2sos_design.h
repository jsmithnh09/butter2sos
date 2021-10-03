/*
 * =======================================================================
 *  buttersos_design.h
 * =======================================================================
 */
 
#ifndef BUTTERSOS_DESIGN_H
#define BUTTERSOS_DESIGN_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>


typedef float regular_t;
typedef complex float complex_t;

#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

#define PRECISION 32
#define HALF_PI   (regular_t) M_PI / 2.0
#define TWO_PI    (regular_t) M_PI * 2.0
#define REG_PI    (regular_t) M_PI
#define N_SOSCOEFFS 6


// linear space generation.
static regular_t* linspace(const int start, const int stop, const int step, int *outsize);

// pole seed around the unit circle.
static complex_t* seedpoles(const int order, int* numpoles);

// print the sos matrix to stdout
void printsosmatrix(const regular_t* matrix, const int nstages);
 
// create the sos matrix
static regular_t* mksosmatrix(const int order, const int type);

// print the complex array
void printcarray(const complex_t* x, const int len);

// euler identity
static complex_t euler(regular_t x);
 
// complex division
static complex_t compdiv(complex_t x, complex_t y);

// complex multiplication
static complex_t compmult(complex_t x, complex_t y);

// warp lowpass to bandpass prototype
static void bpfwarp(complex_t* poles, const int npoles, regular_t* zeros, int* nzeros, regular_t* gain, const regular_t* bwidth, const regular_t* Wn);

// warp lowpass to highpass prototype
static void hpfwarp(complex_t* poles, const int numpoles, regular_t* zeros, int* nzeros, regular_t* gain, const regular_t omega);

// warp lowpass to lowpass prototype
static void lpfwarp(complex_t* poles, int numpoles, regular_t* zeros, int* nzeros, regular_t* gain, regular_t omega);

// bilinear transform in-place singularities.
static void bilinear_s2z(complex_t* poles, const int numpoles, const int numzeros, regular_t* gain, const regular_t fs);

// principal design function
regular_t* butter(const int order, const regular_t fc, regular_t fs, const int type);


#endif

/*
 * ==============================================================================
 * buttersos_design.c
 *
 * Designs lowpass, highpass, and allpass filters at specified sampling
 * rates. Even/odd order filters supported.
 *
 * References:
 *
 * [1] J.G. Proakis and D.G. Manolakis, Digital Signal Processing, Prentice
 *    Hall, 2007, chapter 10, section 3.
 * [2] A.V. Oppenheim and R.W. Schafer, Digital Signal Processing, Prentice
 *    Hall, 1975, chapter 5, sections 1 through 3.
 * [3] J.O. Smith III, Digital Filters with Audio Applications, BookSurge
 *    Publishing, 2007, appendix I.
 *
 * ==============================================================================
 */

/*
 * Jordan R. Smith
 *
 */

#include "butter2sos_design.h"

/******************************************************************
 * Function linspace                                              *
 *    Create a linearly spaced element vector. This is for        *
 *    generating poles around the unit circle.                    *
 *  Inputs:                                                       *
 *    start (int) is the starting element                         *
 *    end (int) is the ending element                             *
 *    step (int) is the step between elements.                    *
 *    outsize (int) is a pointer for indicating the array length. *
 *  Output is the pointer to the linear array, or NULL on failure.*
 ******************************************************************/

regular_t* linspace(const int start, const int stop, const int step, int *outsize)
{
    int m = round((int)(stop-start)/step);
    regular_t* y;
    y = (regular_t*)malloc(sizeof(regular_t)*(m+1));
    *outsize = (int) m+1;
    for (int d = 0; d < m+1; d++)
    {
        y[d] = (regular_t)start + (regular_t)d*step;
    }
    return &(y[0]);
}


/**********************************************************
 *  Function compadd                                      *
 *    Generates a complex number using the equation       *
 *      z1 + z2 = (re(z1) + re(z2)) + (im(z1) + im(z2))i  *
 **********************************************************/

inline complex_t compadd(complex_t z1, complex_t z2)
{
  return (complex_t)((creal(z1) + creal(z2)) + (cimag(z1) + cimag(z2))*I);
}

/**********************************************************
 *  Function compsub                                      *
 *    Generates a complex number using the equation       *
 *      z1 - z2 = (re(z1) - re(z2)) + (im(z1) - im(z2))i  *
 **********************************************************/

inline complex_t compsub(complex_t z1, complex_t z2)
{
  return (complex_t)((creal(z1) - creal(z2)) + (cimag(z1) - cimag(z2))*I);
}

/*************************************************************
 *  Function compsqrt                                        *
 *    complex square root using the equation                 *
 *    sqrt(z) = sqrt(r)*(cos(phi/2) + sin(phi/2)i), where    *
 *      r = abs(z), phi = angle(z).                          *
 *************************************************************/

inline complex_t compsqrt(complex_t x)
{
  #if BUTTER2SOS_PRECISION == 32
    return (complex_t) csqrtf(x);
  #else
    return (complex_t) csqrt(x);
  #endif
}

/**********************************************************************
 *  Function euler                                                    *
 *    Generates a complex number using the identity                   *
 *    e^jx = cos(x) + j*sin(x).                                       *
 *  Inputs:                                                           *
 *    x (regular_t) is the real number "x" to convert to a complex.   *
 *  Outputs a complex number using euler's equation.                  *
 **********************************************************************/

inline complex_t euler(regular_t x)
{
  #if BUTTER2SOS_PRECISION == 32
    return (complex_t) (cosf(x) + sinf(x)*I);
  #else
    return (complex_t) (cos(x) + sin(x)*I);
  #endif
}

/*************************************************
 *  Function isreal                              *
 *    determine if a singularity is purely real. *
 *  Inputs:                                      *
 *    x (complex_t*) is a complex value.         *
 *  Outputs:                                     *
 *    isreal (int) is 0 if its not purely real.  *
 * Notes:                                        *
 * Comparing against a threshold to eliminate    *
 * any imaginary residue.                        *
 *************************************************/
inline int isreal(complex_t x)
{
 #if BUTTER2SOS_PRECISION == 32
  return (int)(cimag(x) < 100*FLT_EPSILON);
 #else
  return (int)(cimag(x) < 100*DBL_EPSILON);
 #endif
}


/**********************************************************************
 * Function seedpoles                                                 *
 *    Generates S-plane poles on the unit circle. If the              *
 *    order is odd, the last pole is located at (-1, 0) pre-warping.  *
 *  Inputs:                                                           *
 *    order (int) indicates the filter order.                         *
 *  Output is a complex_t array of poles around the unit circle.      *
 **********************************************************************/

complex_t* seedpoles(const int order, int* numpoles)
{
  int* len = (int*)malloc(sizeof(int));
  
  // positions around the unit circle [1]
  regular_t* points = linspace((int)1, order-1, (int)2, len);
  for (int pInd = 0; pInd < *len; pInd++)
  {
    points[pInd] = REG_PI * (points[pInd]/(2*order)) + HALF_PI;
  }
  
  // allocate the complex positions around the circle. 
  int Npoles = (order%2) ? (int)2*(*len)+1 : (int)2*(*len);
  complex_t* poles = (complex_t*)malloc(sizeof(complex_t) * Npoles);

  
  // iterate over each pair and calculate the complex values.
  int posInd = 0;
  for (int pInd = 0; pInd < *len; pInd++)
  {
    poles[2*pInd] = euler(points[pInd]); // generate the pole position.
    poles[2*pInd+1] = conj(poles[2*pInd]); // add the conjugate.
    posInd+=2;
  }
  
  // append the pole at (-1, 0) if odd order, [1].
  if (order % 2)
  {
    poles[posInd] = (regular_t)-1.0 + 0.0*I;
  }

  // indicate the length.
  *numpoles = (int)Npoles;
  
  // free the intermediary values.
  free(points);
  free(len);
  return &(poles[0]);
  
}

/******************************************************
 *  Function mksosmatrix                              *
 *    Allocates a flat SOS matrix given the order N.  *
 *  Inputs:                                           *
 *    order (int) is the order of the filter.         *
 *    type (int) indicates 0 for LPF/APF and 1 for HPF*
 *  Outputs:                                          *
 *    matrix (regular_t*) is the matrix, [(N/2)x6].   *
 ******************************************************/

regular_t* mksosmatrix(const int order, const int type)
{
  // determine the number of SOS stages.
  #if BUTTER2SOS_PRECISION == 32
    int N = (int)ceilf((regular_t)order/2);
  #else
    int N = (int)ceil((regular_t)order/2);
  #endif
  int b0, a0, sign;

  regular_t* matrix = (regular_t*)malloc(sizeof(regular_t)* N * N_SOSCOEFFS);
  for (int grpInd = 0; grpInd < N; grpInd++)
  {
    b0 = grpInd*N_SOSCOEFFS;
    a0 = grpInd*N_SOSCOEFFS + 3;

    matrix[b0] = (regular_t)1.0;
    matrix[a0++] = (regular_t)1.0;
    matrix[a0++] = (regular_t)0.0;
    matrix[a0] = (regular_t)0.0;

    
    if (!type)
    {
      // LPF/APF type -> zeros at -1, or [1, +2, 1].
      sign = (regular_t)1.0;
    }
    else
    {
      // HPF type -> zeros at +1, or [1, -2, 1].
      sign = (regular_t)-1.0;
    }

    matrix[b0++] = 1.0;
    if ((grpInd == 0) && (order % 2)) 
    {
      // first stage is going to be [1, +/-1, 0] if odd order for dangling pole.
      matrix[b0++] = sign * 1.0;
      matrix[b0++] = 0.0;
    }
    else
    {
      // if not first stage or odd order, we have [1, +/-2, 1] for zero roots.
      matrix[b0++] = sign * 2.0;  
      matrix[b0++] = 1.0;
    }
  } 
  return &(matrix[0]);
}

/*******************************************
 *  Function printsosmatrix                *
 *    prints out the SOS matrix to stdout. *
 *  Inputs:                                *
 *    matrix (regular_t*) is a pointer to  *
 *        the SOS matrix.                  *
 *    nstages (int) indicates the number   *
 *        of biquads to print.             *
 *******************************************/

void printsosmatrix(const regular_t* matrix, int nstages)
{
  // a0 term ignored, implicitly 1.
  printf("#      b0              b1              b2              a0              a1              a2\n");
  for (int grpInd = 0; grpInd < nstages; grpInd++) {
    printf("%- 16.7e%- 16.7e%- 16.7e%- 16.7e%- 16.7e%- 16.7e\n", matrix[grpInd*N_SOSCOEFFS], \
            matrix[grpInd*N_SOSCOEFFS+1], matrix[grpInd*N_SOSCOEFFS+2], \
            matrix[grpInd*N_SOSCOEFFS+3], \
            matrix[grpInd*N_SOSCOEFFS+4], matrix[grpInd*N_SOSCOEFFS+5]);
  }
}

/*******************************************
 *  Function printcarray                   *
 *    prints out a complex array.          *
 *******************************************/

void printcarray(const complex_t* x, const int len)
{
  regular_t re, im;
  printf("[1x%d] complex\n", len);
  for (int cInd = 0; cInd < len; cInd++)
  {
    re = creal(x[cInd]);
    im = cimag(x[cInd]);
    
    printf("%p : %2.5f + %2.5fi\n", (void*)(&x + cInd), re, im);
  }
}
      
      

/************************************************************************
 *  Function lpfwarp                                                    *
 *    Warps the pole positions given a normalized rotational frequency. *
 *  Inputs:                                                             *
 *    poles (complex_t*) is the array of poles to transform.            *
 *    np (int) are the number of poles in the array.                    *
 *    gain (regular_t) is the gain change associated with moving poles. *
 ************************************************************************/

void lpfwarp(complex_t* poles, int numpoles, regular_t* zeros, int* nzeros, regular_t* gain, regular_t omega)
{
  for (int pInd = 0; pInd < numpoles; pInd++)
  {
    poles[pInd] = (complex_t) (creal(poles[pInd]) * omega) + (cimag(poles[pInd] * omega))*I; // w .* p
  }
  // k = k * w^(order), where order = Np - Nz (no zeros in butterworth design.)
  #if BUTTER2SOS_PRECISION == 32
    *gain *= powf(omega, numpoles);
  #else
    *gain *= pow(omega, numpoles); 
  #endif
  
  // no zeros generated from an LPF warp.
  *nzeros = 0;
  zeros = (regular_t*)NULL;
}


/***********************************************************
 *  Function compdiv (complex division)                    *
 *    Performs complex division (X/Y).                     *
 *  Inputs:                                                *
 *    x (complex_t) is the numerator                       *
 *    y (complex_t) is the denominator                     *
 *  Outputs:                                               *
 *    z (complex_t) is the result of the division.         *
 ***********************************************************/

complex_t compdiv(complex_t x, complex_t y)
{
  /* for X = a+bi and Y = c+di, Z = X/Y is equal to:
   *  ((ac + bd) / (c^2 + d^2)) + ((bc - ad)/(c^2 + d^2)i
   */
  regular_t num1 = creal(x)*creal(y) + cimag(x)*cimag(y);
  regular_t num2 = cimag(x)*creal(y) - creal(x)*cimag(y);
  #if BUTTER2SOS_PRECISION == 32
    regular_t den = powf(creal(y), 2) + powf(cimag(y), 2);
  #else
    regular_t den = pow(creal(y), 2) + pow(cimag(y), 2);
  #endif
  return (complex_t)(num1/den) + (num2/den)*I;
}


/***********************************************************
 *  Function compmult (complex division)                   *
 *    Performs complex multiplication (X/Y).               *
 *  Inputs:                                                *
 *    x (complex_t) is the first complex number.           *
 *    y (complex_t) is the second complex number.          *
 *  Outputs:                                               *
 *    z (complex_t) is the result of the multiplication.   *
 ***********************************************************/

complex_t compmult(complex_t x, complex_t y)
{
  /* for X = a+bi and Y = c+di, Z = X*Y is equal to:
   *  (ac - bd) + (ad + bc)i
   */
  regular_t re = creal(x)*creal(y) - cimag(x)*cimag(y);
  regular_t im = creal(x)*cimag(y) + cimag(x)*creal(y);
  return (complex_t) re + im*I;
}


/**********************************************************************
 *  Function hpfwarp                                                  *
 *    Warps the pole positions given a normalized rotation frequency. *
 *    This function also returns zeros at the origin if any at Inf.   *
 *  Inputs/Outputs:                                                   *
 *    poles (complex_t*) is the complex pole array (modified).        *
 *    npoles (int) is the length of the pole array.                   *
 *    gain (regular_t*) is a pointer to the gain of the filter.       *
 **********************************************************************/

void hpfwarp(complex_t* poles, const int numpoles, regular_t* zeros, int* nzeros, regular_t* gain, const regular_t omega)
{
  complex_t wc = omega + (regular_t)0.0*I;              // omega + 0i, since frequency is real.
  regular_t gshift = (regular_t)0.0;                    // gain change from HPF warping.
  complex_t accum = (regular_t)1.0 + 0*I;               // accumulated gain change in poles.
  complex_t negone = (regular_t)-1.0 + 0*I;             // negative one.
  complex_t posone = (regular_t)1.0 + 0.0*I;            // positive one.
  complex_t interm = (regular_t)0.0 + 0.0*I;            // intermediate value from -1*p[i].
  
  for (int pInd = 0; pInd < numpoles; pInd++)
  {
    interm = compmult(negone, poles[pInd]);
    accum = compmult(accum, interm); // prod(-p).
    poles[pInd] = compdiv(wc, poles[pInd]); // w ./ p
  }
  
  /* since no zeros pre-transform, the gain approx. for real(prod(-z)/prod(-p)) is
   * equal to real(1/prod(-p)).
   */
  gshift = (regular_t)creal(compdiv(posone, accum));
  *gain *= gshift; // k = k * real(prod(-z)/prod(-p)).
  
  // allocating the zeros at the origin.
  *nzeros = numpoles;
  zeros = (regular_t*)malloc(sizeof(regular_t)*numpoles);
  for (int zInd = 0; zInd < numpoles; zInd++)
  {
    zeros[zInd] = (regular_t)0.0; // zero-padding.
  }
}

/************************************************************************
 *  Function bpfwarp                                                    *
 *    Warps the pole positions based on the corner frequencies.         *
 *  Inputs:                                                             *
 *    poles (complex_t*) is the complex pole array.                     *
 *    numpoles (int) are the number of poles. This will double in size. *
 *    zeros (regular_t*) is the zero array.                             *
 *    numzeros (int*) is the number of zeros. This grows to zeros(np,1) *
 *    gain (regular_t*) is the gain approximation.                      *
 *    bwidth (regular_t*) is the bandwidth between corner frequencies.  *
 *    Wn (regular_t*) is the center normalized frequency.               *
 ************************************************************************/

void bpfwarp(complex_t* poles, int* numpoles, complex_t* zeros, int* numzeros, regular_t* gain, const regular_t* bwidth, const regular_t* Wn)
{
  int order = *numpoles;                // need original order for the gain approximation.
  complex_t bw2 = (*bwidth)/2 + 0.0*I;  // half bandwidth, bw/2.
  complex_t bsq, bsum, pold, cWn;       // intermediate complex values for warp calculations.
  cWn = (regular_t)*Wn + 0.0*I;
  complex_t Wn2 = compmult(cWn, cWn); // complex Wn^2.
  
  // p * bw/2 for the original data.
  for (int pInd = 0; pInd < (*numpoles); pInd++)
  {
    poles[pInd] = (complex_t)compmult(poles[pInd], bw2);  
  }
  
  // [1,2,3,4] -> [1, 0, 2, 0, 3, 0, 4, 0]. expand the array.
  int oldInd = (*numpoles)-1;
  for (int pInd = 2*(*numpoles)-2; pInd > 0; pInd-=2)
  {
    poles[pInd] = poles[oldInd];
    poles[pInd+1] = (complex_t)0.0 + 0.0*I;
    poles[oldInd--] = (complex_t)0.0 + 0.0*I;
  }
  
  // compute the p +/- sqrt(p^2 - Wn^2) for each group.
  for (int pInd = 0; pInd < (*numpoles); pInd++)
  {
    // extract pole intermediate values for BPF warping.
    pold = poles[2*pInd];
    bsum = (complex_t)compsub(compmult(pold, pold), Wn2); // p^2 - Wn^2
    bsq = (complex_t)csqrt(bsum); // sqrt(p^2 - Wn^2)
    
    // p(i)   = p(i) + sqrt(p(i)^2 - Wn^2).
    // p(i+1) = p(i) - sqrt(p(i)^2 - Wn^2).
    
    poles[2*pInd] = (complex_t)compadd(pold, bsq);
    poles[2*pInd+1] = (complex_t)compsub(pold,  bsq);
  }
  
  // indicate the new number of poles from this process.
  *numpoles = (int)2*(*numpoles);

  // k = k * bw^(order)
  *gain *= (regular_t)pow(*bwidth, order);
}


/************************************************************
 *  Function polesort                                       *
 *    Sorts poles based on proximity to the unit circle.    *
 *  Inputs:                                                 *
 *    poles (complex_t*) is a pointer to the pole array.    *
 *    numpoles (int*) are the number of poles in the array. *
 *    order (int*) indicates the order of the filter. If    *
 *      odd, the last two poles are real from               *
 *      BPF warping.                                        *
 ************************************************************/

void polesort(complex_t* poles, int numpoles, int order)
{
  int len = (order % 2) ? numpoles-2 : numpoles;
  
  // bubble sorting the poles into descending distance from unit circle.
  complex_t temp;
  for (int outerInd = 0; outerInd < len-1; outerInd++)
  {
    for (int compInd = 0; compInd < len-outerInd-1; compInd++)
    {
      if (cabs(poles[compInd]) < cabs(poles[compInd+1]))
      {
        temp = poles[compInd];
        poles[compInd] = poles[compInd+1];
        poles[compInd+1] = temp;
      }
    }
  }
  
  // flip the real poles in the case of an odd order bandpass.
  if ((order % 2) && (creal(poles[len]) < creal(poles[len+1])))
  {
    temp = poles[len];
    poles[len] = poles[len+1];
    poles[len+1] = temp;
  }
}

/**************************************************************
 *  Function bsfwarp                                          *
 *    Warps the poles and zeros to a bandstop filter type.    *
 *  Inputs:                                                   *
 *    poles (complex_t*) are the seeded poles.                *
 *    numpoles (int*) is the number of poles in the array.    *
 *    zeros (complex_t*) is a pointer to the zeros array.     *
 *    numzeros (int*) is the number of zeros in the array.    *
 *    gain (regular_t*) is the gain of the filter.            *
 *    bwidth (regular_t*) is the bandwidth of the stop region.*
 *    Wn (regular_t*) is the normalized center frequency.     *
 * Outputs:                                                   *
 *    numpoles will be modified to double in size.            *
 *    numzeros will be the same length as numpoles with       *
 *      repeating complex conjugates.                         *
 **************************************************************/

void bsfwarp(complex_t* poles, int* numpoles, complex_t* zeros, regular_t* gain, const regular_t* bwidth, const regular_t* Wn)
{
  complex_t bw2 = (*bwidth / 2.0) + 0.0*I;
  complex_t pold, bsum, bsq;
  complex_t Wn2 = compmult((complex_t)*Wn+0.0*I, (complex_t)*Wn+0.0*I); // complex Wn^2.
  
  // gain shift is just 1, since k * real(prod(-z)/prod(-p)) = 1
  // when the poles are around the unit circle and no zeros.
  *gain = (regular_t)1.0;
  
  // warp the poles akin to the HPF transformation, by inverting.
  for (int pInd = 0; pInd < *numpoles; pInd++)
  {
    poles[pInd] = compdiv(bw2, poles[pInd]); // (bw/2) ./ p
  }
  
  // expanding [1, 2, 3, 4] -> [1, 0, 2, 0, 3, 0, 4, 0].
  int oldInd = (*numpoles)-1;
  for (int pInd = 2*(*numpoles)-2; pInd > 0; pInd-=2)
  {
    poles[pInd] = poles[oldInd];
    poles[pInd+1] = (complex_t)0.0 + 0.0*I;
    poles[oldInd--] = (complex_t)0.0 + 0.0*I;
  }
  
  // zeros already same length as poles, (0+/-j*Wn).
  for (int zInd = 0; zInd < *numpoles; zInd++)
  {
    zeros[2*zInd] = (complex_t)0.0 + (*Wn)*I;
    zeros[2*zInd+1] = (complex_t)0.0 - (*Wn)*I;
  }
  
  // calculating p +/- sqrt(p^2 - Wn^2)
  for (int pInd = 0; pInd < *numpoles; pInd++)
  {
    pold = poles[2*pInd];
    bsum = compsub(compmult(pold, pold), Wn2); // p^2 - Wn^2
    bsq = csqrt(bsum);
    
    // p(i) = p(old) + sqrt(p(old)^2 - Wn^2)
    // p(i+1) = p(old) - sqrt(p(old)^2 - Wn^2)
    poles[2*pInd] = compadd(pold, bsq);
    poles[2*pInd+1] = compsub(pold, bsq);
  }
}

/************************************************************************
 *  Function bilinear_band_s2z                                          *
 *    Warps the pole and zero positions for bandpass/bandstop filters.  *
 *  Inputs:                                                             *
 *    poles (complex_t*) is the complex pole array to warp.             *
 *    npoles (int*) is the number of poles.                             *
 *    zeros (complex_t*) is the complex zero array to warp.             *
 *    nzeros (int*) are the number of zeros.                            *
 *    gain (regular_t*) is the pointer to the gain of the filter.       *
 *                                                                      *
 ************************************************************************/

void bilinear_band_s2z(complex_t* poles, const int* numpoles, complex_t* zeros, int* numzeros, regular_t* gain, const regular_t* fs)
{
  complex_t fs2 = (regular_t)2*(*fs) + 0.0*I;  // 2/T equivalent.
  complex_t pos1 = (regular_t)1.0 + 0.0*I;  // bilinear operation for (z+1)/(z-1)
  complex_t pgain = (regular_t)1.0 + 0.0*I;
  complex_t zgain = (regular_t)1.0 + 0.0*I; // zero singularity gain also starts at 1.
  complex_t num, den, warp, sub;            // intermediate values.

  // prod(fs2 - p)
  for (int pInd = 0; pInd < *numpoles; pInd++)
  {
    sub = compsub(fs2, poles[pInd]);
    pgain = compmult(pgain, sub);
  }

  // prod(fs2 - z)
  for (int zInd = 0; zInd < *numzeros; zInd++)
  {
    sub = compsub(fs2, zeros[zInd]);
    zgain = compmult(zgain, sub);
  }

  // k = k * real(prod(fs2-p)/prod(fs2-z))
  *gain *= (regular_t)creal(compdiv(zgain, pgain)); // calculates gain change.
  
  for (int pInd = 0; pInd < *numpoles; pInd++)
  {
    warp = compdiv(poles[pInd], fs2);
    num = compadd(pos1, warp); // 1 + z
    den = compsub(pos1, warp); // 1 - z
    poles[pInd] = compdiv(num, den); // (1+z)/(1-z)
  }

  // bilinear xform zeros.
  for (int zInd = 0; zInd < *numzeros; zInd++)
  {
    warp = compdiv(zeros[zInd], fs2);
    num = compadd(pos1, warp);
    den = compsub(pos1, warp);
    zeros[zInd] = compdiv(num, den);
  }
}


/***************************************************************
 *  Function butterband                                        *
 *    This is the principal entry point for designing bandpass *
 *    or bandstop SOS filter matrices.                         *
 *                                                             *
 * Inputs:                                                     *
 *  order (int) is the butterworth filter order.               *
 *  Flo (regular_t) is the lower corner frequency.             *
 *  Fhi (regular_t) is the upper corner frequency.             *
 *  Fs (regular_t) is the sampling rate to design at.          *
 *  type (int) indicates 0=BPF, 1=BSF.                         *
 * Outputs:                                                    *
 *  matrix (regular_t*) is a pointer to the SOS matrix.        *
 ***************************************************************/

regular_t* butterband(const int order, regular_t flo, regular_t fhi, regular_t fs, const int type)
{
  
  // check if the corner frequencies need to swap, or if they're equal.
  regular_t swap;
  if (flo > fhi)
  {
    swap = flo;
    flo = fhi;
    fhi = swap;
  }
  else if (flo == fhi) 
  {
    printf("ERROR: corner frequencies are equivalent.\n");
    return (regular_t*)NULL;
  }

  // check if the upper frequency exceeds nyquist.
  if (fhi >= fs/2.0)
  {
    printf("ERROR: upper corner frequency exceeds Nyquist.\n");
    return (regular_t*)NULL;
  }

  int npoles, nzeros = 0;
  regular_t* mat;
  complex_t* poles = (complex_t*)NULL;
  complex_t* zeros = (complex_t*)NULL;
  regular_t gain, w1, w2, bwidth, Wn;

  gain = 1.0;
  w1 = 2.0 * (flo/fs);
  w2 = 2.0 * (fhi/fs);

  #if BUTTER2SOS_PRECISION == 32
    w1 = (regular_t)tanf(HALF_PI * w1) * 4.0;
    w2 = (regular_t)tanf(HALF_PI * w2) * 4.0;
  #else
    N = (int)ceil((regular_t)order/2);
    w1 = (regular_t)tan(HALF_PI * w1) * 4.0;
    w2 = (regular_t)tan(HALF_PI * w2) * 4.0;
  #endif

  bwidth = w2 - w1;     // normalized bandwidth
  Wn = sqrt(w1 * w2);   // center frequency
  fs = 2.0;             // normalized Wq = 1, resetting Fs.

  // seed the poles around the unit circle.
  poles = seedpoles(order, &npoles);

  // make the SOS matrix.
  mat = (regular_t*)malloc(sizeof(regular_t)*N_SOSCOEFFS*order);

  // double the length for the BPF/BSF to manage.
  poles = (complex_t*)realloc(poles, sizeof(complex_t)*2*npoles);
  if (!type)
  {

    // bandpass warp.
    bpfwarp(poles, &npoles, zeros, &nzeros, &gain, &bwidth, &Wn);

    // since no zeros were seeded, create S-plane zeros, (making same length as poles.)
    zeros = (complex_t*)malloc(sizeof(complex_t)*npoles);
    nzeros = npoles/2.0;

    // bilinear warps poles to [+1,+1,....-1,-1].
    for (int zInd = 0; zInd < npoles; zInd++)
    {
      zeros[zInd] = zInd < nzeros ? (complex_t)0.0 + 0.0*I : (complex_t)-1.0 + 0.0;
    }
  }
  else
  {
    // bandstop filter warp. BSF will place complex conjugates in the zeros array.
    nzeros = npoles*2;
    zeros = (complex_t*)malloc(sizeof(complex_t)*nzeros);
    bsfwarp(poles, &npoles, zeros, &gain, &bwidth, &Wn);
    npoles *= 2; // only change the number after having used the loop in BSF warp.

  }

  printf("S-plane poles:\n");
  printcarray(poles, npoles);
  printf("S-plane zeros:\n");
  printcarray(zeros, nzeros);


  // bilinear transform.
  bilinear_band_s2z(poles, &npoles, zeros, &nzeros, &gain, &fs);

  // sort the poles after warping to Z-plane.
  polesort(poles, npoles, order);

  /*
   * sos matrix quantization. We either have real poles or
   * repeating complex conjugates based on the filter type.
   */
  int pInd, zInd, posInd = 0;
  for (int stgInd = 0; stgInd < order; stgInd++)
  {
    pInd = 2*stgInd;
    zInd = 2*stgInd;               // normal order for zeros.
    posInd = (stgInd*N_SOSCOEFFS); // position in the stage.

    // zero check
    if (isreal(zeros[zInd]) && isreal(zeros[zInd+1])) 
    {
      if ((creal(zeros[zInd]) == -1) && (creal(zeros[zInd+1]) == -1)) 
      {
        // [-1, -1] -> [1, 2, 1]
        mat[posInd++] = 1.0;
        mat[posInd++] = 2.0;
        mat[posInd++] = 1.0;
      }
      else if((creal(zeros[zInd]) == 1) && (creal(zeros[zInd+1]) == -1)) 
      {
        // [1, -1] -> [1, 0, -1]
        mat[posInd++] = 1.0;
        mat[posInd++] = 0.0;
        mat[posInd++] = -1.0;
      }
      else if((creal(zeros[zInd]) == 1) && (creal(zeros[zInd+1]) == 1)) 
      {
        // [1, 1] -> [1, -2, 1]
        mat[posInd++] = 1.0;
        mat[posInd++] = -2.0;
        mat[posInd++] = 1.0;
      }
    }
    else 
    {
      // complex conjugates.
      mat[posInd++] = 1.0;
      mat[posInd++] = -2.0 * creal(zeros[zInd]);
      mat[posInd++] = creal(compmult(zeros[zInd], conj(zeros[zInd])));
    }

    // a0 implicitly 1.
    mat[posInd++] = 1.0;


    // pole check
    if (isreal(poles[pInd]) && isreal(poles[pInd+1])) 
    {
      mat[posInd++] = -1.0*creal(poles[pInd]) + -1.0*creal(poles[pInd]); // -p1 + -p2
      mat[posInd] = (-1.0*creal(poles[pInd])) * (-1.0*creal(poles[pInd])); // -p1 * -p2
    }
    else 
    {
      mat[posInd++] = -2.0*creal(poles[pInd]);
      mat[posInd] = creal(compmult(poles[pInd], conj(poles[pInd])));
    }
  }

  // encorporate the gain into the last stage.
  mat[(N_SOSCOEFFS*(order-1))] *= gain;
  mat[(N_SOSCOEFFS*(order-1))+1] *= gain;
  mat[(N_SOSCOEFFS*(order-1))+2] *= gain;

  // release the complex array.
  free(poles);
  free(zeros);

  return &(mat[0]);
}

/**********************************************************************
 *  Function bilinear_s2z                                             *
 *    Bilinear transformation from S to Z-plane. The transform is     *
 *    recognized as:                                                  *
 *                                                                    *
 *      H(z) = H(s) where s = 2/T * (z+1)/(z-1).                      *
 *                                                                    *
 *  Inputs:                                                           *
 *    poles (complex_t*) are the poles that are being re-mapped.      *
 *    numpoles (int) are the number of poles in the pole array.       *
 *    numzeros (int) indicates if zeros were generated based on       *
 *      the filter type, (LP/HP/APF).                                 *
 *    gain (regular_t*) is a scalar gain that will be updated based   *
 *      on the transformation.                                        *
 *    fs (regular_t) is the sampling rate we're mapping to.           *
 *                                                                    *
 **********************************************************************/

void bilinear_s2z(complex_t* poles, const int numpoles, const int numzeros, regular_t* gain, const regular_t fs)
{
  complex_t fs2 = (regular_t)2*fs + 0.0*I;  // 2/T equivalent in bilinear xform.
  complex_t pgain = (regular_t)1.0 + 0.0*I; // complex +1, (starting value for pole gain shift.)
  complex_t num, den, warp, psub;           // intermediate values for bilinear xform.
  regular_t zgain;                          // gain-shift for zeros.
  
  // gain adjustment from poles is "prod(Fs-p)", pre-bilinear xform.
  for (int pInd = 0; pInd < numpoles; pInd++)
  {
    // (a+bi) - (c+di) = (a - c) + (b - d)*i. Fs-p -> (Fs-real(p)) + (0 - imag(p))i, since Fs is real.
    psub = (complex_t)((creal(fs2) - creal(poles[pInd])) + (cimag(fs2) - cimag(poles[pInd]))*I);
    pgain = compmult(pgain, psub);
  }
  
  // bilinear transform, S to Z plane.
  for (int pInd = 0; pInd < numpoles; pInd++)
  {
    // xform => (1+p/Fs)./(1-p/Fs), or (1+(p/Fs))./(1-(p/Fs)).
    // (a + c) + (b + d)i, 1 + (p/Fs) -> (1 + Re(p/Fs)) + (0 + Im(p/Fs))i
    // (a - c) + (b - d)i, 1 - (p/Fs) -> (1 - Re(p/Fs)) + (0 - Im(p/Fs))i
    warp = compdiv(poles[pInd], fs2);
    num = ((regular_t)1.0 + creal(warp)) + cimag(warp)*I;
    den = ((regular_t)1.0 - creal(warp)) + ((regular_t)0.0 - cimag(warp))*I;
    poles[pInd] = compdiv(num, den);
  }
  
  // check HPF vs. LPF/APF design cases.
  if (numzeros != 0)
  {
    /* for the HPF case, zeros on S-plane move from 0 to 1 on the Z-plane. Then,
     * for determining the gain adjustment of the zeros, there are an equal number
     * of zeros and poles. This means the "prod(Fs-z|s)" term is simply fs2^npoles.
     */
     #if BUTTER2SOS_PRECISION == 32
        zgain = powf((regular_t)2*fs, (regular_t)numpoles);
     #else
        zgain = pow((regular_t)2*fs, (regular_t)numpoles);
     #endif
  }
  else
  {
    // the lowpass/allpass case has no zeros generated; the "prod(Fs-z)" term is simply 1.
    zgain = 1.0;
  }

  // k = k * real(prod(Fs-z)/prod(Fs-p))
  *gain *= creal(compdiv((complex_t)(zgain + 0.0*I), pgain));
}

/****************************************************************
 * Function butter                                              *
 *    This is the principal entry point for designing the       *
 *    butterworth SOS filter matrix.                            *
 *                                                              *
 * Inputs:                                                      *
 *  order (int) is the butterworth filter order.                *
 *  Fc (regular_t) is the corner frequency (-3 dB point.)       *
 *  Fs (regular_t) is the discrete sampling rate to design at.  *
 *  type (int) indicates 0=LPF, 1=HPF, 2=APF.                   *
 * Outputs:                                                     *
 *  matrix (regular_t*) is a pointer to the SOS matrix.         *
 ****************************************************************/

regular_t* butter(const int order, const regular_t fc, regular_t fs, const int type)
{

  if (fc >= fs/2) {
    printf("ERROR: corner frequency must be less than Nyquist rate.\n");
    return (regular_t*)NULL;
  }
  if ((type < 0) || (type > 2)) {
    printf("ERROR: Unrecognized filter type. 0=LPF, 1=HPF, 2=APF.\n");
    return (regular_t*)NULL;
  }

  int N, npoles, nzeros;
  regular_t* mat;
  complex_t* poles;
  regular_t gain, zeros, omega, warp;

  gain = 1.0;             // starting with unity gain.
  omega = 2.0 * (fc/fs);  // normalizing frequency.

  #if BUTTER2SOS_PRECISION == 32
    N = (int)ceilf((regular_t)order/2);
    warp = (regular_t)tanf(HALF_PI * omega) * 4.0;
  #else
    N = (int)ceil((regular_t)order/2);
    warp = (regular_t)tan(HALF_PI * omega) * 4.0;
  #endif

  // since omega is now normalized about Nq = 1.set Fs=2.
  fs = 2.0;

  // APF defaults to LPF formulation.
  if (type == 2) 
  {
    mat = mksosmatrix(order, (int)0);
  }
  else 
  {
    mat = mksosmatrix(order, (int)type);
  }

  // seed the poles around the unit circle.
  poles = seedpoles(order, &npoles);

  // warp the poles based on the filter type.
  if ((type == 0) || (type == 2))
  {
    lpfwarp(poles, npoles, &zeros, &nzeros, &gain, warp);
  }
  else 
  {
    hpfwarp(poles, npoles, &zeros, &nzeros, &gain, warp);
  }

  // bilinear transform S to Z.
  bilinear_s2z(poles, npoles, nzeros, &gain, fs);

  // build the SOS polynomials and incorporate the gain.
  int posInd, poleInd = 0;
  poleInd = (order % 2) ? npoles-2 : npoles-1; // cascade "up" in reverse order.
  for (int stageInd = 0; stageInd < N; stageInd++)
  {
    posInd = (int)(N_SOSCOEFFS*stageInd) + 4; // moving to a2 position.
    // first stage would contain dangling pole if odd order.
    if ((stageInd == 0) && (order % 2))
    {
      /* since it's already known the last pole from the seedpoles
       * function is -1+0i (S-plane), don't need to increment pole index.
       * sosmatrix is already [1 0 0 1 0 0], S-plane pole gives [1, -pole, 0],
       * or just modifying a1.
       */
      mat[posInd] = -1.0 * creal(poles[npoles-1]);
      continue;
    }

    // place the position index at a1, (nstage*ncoeffs + 4)
    mat[posInd++] = -2.0 * creal(poles[poleInd]);
    mat[posInd] = creal(compmult(poles[poleInd], conj(poles[poleInd])));
    poleInd -= 2; 
  }

  free(poles); // no longer need complex pole positions.

  // since poles are now quantized, loop over if allpass was specified.
  if (type == 2)
  {
    // allpass flips the poles into the zeros for phase relationship.
    posInd = 0;
    for (int stageInd = 0; stageInd < N; stageInd++)
    {
      posInd = (int)N_SOSCOEFFS*stageInd;
      if ((stageInd == 0) && (order % 2))
      {
        // first stage is going to be [-pole 1 0, 1 -pole 0] for odd order.
        mat[posInd] = mat[posInd+4];
        mat[posInd+1] = 1.0;
        mat[posInd+2] = 0.0;
      }
      else
      {
        // now [a2 a1 1 1 a1 a2].
        mat[posInd] = mat[posInd+5];
        mat[posInd+1] = mat[posInd+4];
        mat[posInd+2] = 1.0;
      }
    }
  }
  else
  {
    // encorporate the gain into the last stage, if its not an allpass.
    mat[(N_SOSCOEFFS*(N-1))+0] *= gain;
    mat[(N_SOSCOEFFS*(N-1))+1] *= gain;
    mat[(N_SOSCOEFFS*(N-1))+2] *= gain;
  }
  return &(mat[0]);
}

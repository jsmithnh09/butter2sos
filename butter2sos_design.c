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
 *    generating poles on the S-plane axis.                       *
 *  Inputs:                                                       *
 *    start (int) is the starting element                         *
 *    end (int) is the ending element                             *
 *    step (int) is the step between elements.                    *
 *    outsize (int) is a pointer for indicating the array length. *
 *  Output is the pointer to the linear array, or NULL on failure.*
 ******************************************************************/

real64_t* linspace(const int start, const int stop, const int step, int *outsize)
{
    int m = round((int)(stop-start)/step);
    real64_t* y;
    y = (real64_t*)malloc(sizeof(real64_t)*(m+1));
    *outsize = (int) m+1;
    for (int d = 0; d < m+1; d++)
    {
        y[d] = (real64_t)start + (real64_t)d*step;
    }
    return &(y[0]);
}


/**********************************************************
 *  Function compadd                                      *
 *    Generates a complex number using the equation       *
 *      z1 + z2 = (re(z1) + re(z2)) + (im(z1) + im(z2))i  *
 **********************************************************/

inline complex64_t compadd(complex64_t z1, complex64_t z2)
{
  return (complex64_t)((creal(z1) + creal(z2)) + (cimag(z1) + cimag(z2))*I);
}

/**********************************************************
 *  Function compsub                                      *
 *    Generates a complex number using the equation       *
 *      z1 - z2 = (re(z1) - re(z2)) + (im(z1) - im(z2))i  *
 **********************************************************/

inline complex64_t compsub(complex64_t z1, complex64_t z2)
{
  return (complex64_t)((creal(z1) - creal(z2)) + (cimag(z1) - cimag(z2))*I);
}

/*************************************************************
 *  Function compsqrt                                        *
 *    complex square root using the equation                 *
 *    sqrt(z) = sqrt(r)*(cos(phi/2) + sin(phi/2)i), where    *
 *      r = abs(z), phi = angle(z).                          *
 *************************************************************/

inline complex64_t compsqrt(complex64_t x)
{
  return (complex64_t) csqrt(x);
}

/**********************************************************************
 *  Function euler                                                    *
 *    Generates a complex number using the identity                   *
 *    e^jx = cos(x) + j*sin(x).                                       *
 *  Inputs:                                                           *
 *    x (real64_t) is the real number "x" to convert to a complex.   *
 *  Outputs a complex number using euler's equation.                  *
 **********************************************************************/

inline complex64_t euler(real64_t x)
{
  return (complex64_t) (cos(x) + sin(x)*I);
}

/*************************************************
 *  Function isreal                              *
 *    determine if a singularity is purely real. *
 *  Inputs:                                      *
 *    x (complex64_t*) is a complex value.       *
 *  Outputs:                                     *
 *    isreal (int) is 0 if its not purely real.  *
 * Notes:                                        *
 * Comparing against a threshold to eliminate    *
 * any imaginary residue.                        *
 *************************************************/
inline int isreal(complex64_t x)
{
  return (int)(cimag(x) < 100*DBL_EPSILON);
}


/**********************************************************************
 * Function seedpoles                                                 *
 *    Generates S-plane poles on the re/im axis.  If the              *
 *    order is odd, the last pole is located at (-1, 0) pre-warping.  *
 *  Inputs:                                                           *
 *    order (int) indicates the filter order.                         *
 *  Output is a complex64_t array of poles on the S-plane.              *
 **********************************************************************/

complex64_t* seedpoles(const int order, int* numpoles)
{
  int* len = (int*)malloc(sizeof(int));
  
  // determine positions on the Re/Im axis [1]
  real64_t* points = linspace((int)1, order-1, (int)2, len);
  for (int pInd = 0; pInd < *len; pInd++)
  {
    points[pInd] = REG_PI * (points[pInd]/(2*order)) + HALF_PI;
  }
  
  // allocate the complex positions.
  int Npoles = (order%2) ? (int)2*(*len)+1 : (int)2*(*len);
  complex64_t* poles = (complex64_t*)malloc(sizeof(complex64_t) * Npoles);

  
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
    poles[posInd] = (real64_t)-1.0 + 0.0*I;
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
 *    matrix (real64_t*) is the matrix, [(N/2)x6].   *
 ******************************************************/

real64_t* mksosmatrix(const int order, const int type)
{
  // determine the number of SOS stages.
  int N = (int)ceil((real64_t)order/2);
  int b0, a0, sign;

  real64_t* matrix = (real64_t*)malloc(sizeof(real64_t)* N * N_SOSCOEFFS);
  for (int grpInd = 0; grpInd < N; grpInd++)
  {
    b0 = grpInd*N_SOSCOEFFS;
    a0 = grpInd*N_SOSCOEFFS + 3;

    matrix[b0] = (real64_t)1.0;
    matrix[a0++] = (real64_t)1.0;
    matrix[a0++] = (real64_t)0.0;
    matrix[a0] = (real64_t)0.0;

    
    if (!type)
    {
      // LPF/APF type -> zeros at -1, or [1, +2, 1].
      sign = (real64_t)1.0;
    }
    else
    {
      // HPF type -> zeros at +1, or [1, -2, 1].
      sign = (real64_t)-1.0;
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
 *    matrix (real64_t*) is a pointer to  *
 *        the SOS matrix.                  *
 *    nstages (int) indicates the number   *
 *        of biquads to print.             *
 *******************************************/

void printsosmatrix(const real64_t* matrix, int nstages)
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

void printcarray(const complex64_t* x, const int len)
{
  real64_t re, im;
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
 *    poles (complex64_t*) is the array of poles to transform.            *
 *    np (int) are the number of poles in the array.                    *
 *    gain (real64_t) is the gain change associated with moving poles. *
 ************************************************************************/

void lpfwarp(complex64_t* poles, int numpoles, real64_t* zeros, int* nzeros, real64_t* gain, real64_t omega)
{
  for (int pInd = 0; pInd < numpoles; pInd++)
  {
    poles[pInd] = (complex64_t) (creal(poles[pInd]) * omega) + (cimag(poles[pInd] * omega))*I; // w .* p
  }
  // k = k * w^(order), where order = Np - Nz (no zeros in butterworth design.)
  *gain *= pow(omega, numpoles); 
  
  // no zeros generated from an LPF warp.
  *nzeros = 0;
  zeros = (real64_t*)NULL;
}


/***********************************************************
 *  Function compdiv (complex division)                    *
 *    Performs complex division (X/Y).                     *
 *  Inputs:                                                *
 *    x (complex64_t) is the numerator                       *
 *    y (complex64_t) is the denominator                     *
 *  Outputs:                                               *
 *    z (complex64_t) is the result of the division.         *
 ***********************************************************/

complex64_t compdiv(complex64_t x, complex64_t y)
{
  /* for X = a+bi and Y = c+di, Z = X/Y is equal to:
   *  ((ac + bd) / (c^2 + d^2)) + ((bc - ad)/(c^2 + d^2)i
   */
  real64_t num1 = creal(x)*creal(y) + cimag(x)*cimag(y);
  real64_t num2 = cimag(x)*creal(y) - creal(x)*cimag(y);
  real64_t den = pow(creal(y), 2) + pow(cimag(y), 2);
  return (complex64_t)(num1/den) + (num2/den)*I;
}


/***********************************************************
 *  Function compmult (complex division)                   *
 *    Performs complex multiplication (X/Y).               *
 *  Inputs:                                                *
 *    x (complex64_t) is the first complex number.           *
 *    y (complex64_t) is the second complex number.          *
 *  Outputs:                                               *
 *    z (complex64_t) is the result of the multiplication.   *
 ***********************************************************/

complex64_t compmult(complex64_t x, complex64_t y)
{
  /* for X = a+bi and Y = c+di, Z = X*Y is equal to:
   *  (ac - bd) + (ad + bc)i
   */
  real64_t re = creal(x)*creal(y) - cimag(x)*cimag(y);
  real64_t im = creal(x)*cimag(y) + cimag(x)*creal(y);
  return (complex64_t) re + im*I;
}


/**********************************************************************
 *  Function hpfwarp                                                  *
 *    Warps the pole positions given a normalized rotation frequency. *
 *    This function also returns zeros at the origin if any at Inf.   *
 *  Inputs/Outputs:                                                   *
 *    poles (complex64_t*) is the complex pole array (modified).        *
 *    npoles (int) is the length of the pole array.                   *
 *    gain (real64_t*) is a pointer to the gain of the filter.       *
 **********************************************************************/

void hpfwarp(complex64_t* poles, const int numpoles, real64_t* zeros, int* nzeros, real64_t* gain, const real64_t omega)
{
  complex64_t wc = omega + (real64_t)0.0*I;              // omega + 0i, since frequency is real.
  real64_t gshift = (real64_t)0.0;                    // gain change from HPF warping.
  complex64_t accum = (real64_t)1.0 + 0*I;               // accumulated gain change in poles.
  complex64_t negone = (real64_t)-1.0 + 0*I;             // negative one.
  complex64_t posone = (real64_t)1.0 + 0.0*I;            // positive one.
  complex64_t interm = (real64_t)0.0 + 0.0*I;            // intermediate value from -1*p[i].
  
  for (int pInd = 0; pInd < numpoles; pInd++)
  {
    interm = compmult(negone, poles[pInd]);
    accum = compmult(accum, interm); // prod(-p).
    poles[pInd] = compdiv(wc, poles[pInd]); // w ./ p
  }
  
  /* since no zeros pre-transform, the gain approx. for real(prod(-z)/prod(-p)) is
   * equal to real(1/prod(-p)).
   */
  gshift = (real64_t)creal(compdiv(posone, accum));
  *gain *= gshift; // k = k * real(prod(-z)/prod(-p)).
  
  // allocating the zeros at the origin.
  *nzeros = numpoles;
  zeros = (real64_t*)malloc(sizeof(real64_t)*numpoles);
  for (int zInd = 0; zInd < numpoles; zInd++)
  {
    zeros[zInd] = (real64_t)0.0; // zero-padding.
  }
}

/************************************************************************
 *  Function bpfwarp                                                    *
 *    Warps the pole positions based on the corner frequencies.         *
 *  Inputs:                                                             *
 *    poles (complex64_t*) is the complex pole array.                     *
 *    numpoles (int) are the number of poles. This will double in size. *
 *    zeros (real64_t*) is the zero array.                             *
 *    numzeros (int*) is the number of zeros. This grows to zeros(np,1) *
 *    gain (real64_t*) is the gain approximation.                      *
 *    bwidth (real64_t*) is the bandwidth between corner frequencies.  *
 *    Wn (real64_t*) is the center normalized frequency.               *
 ************************************************************************/

void bpfwarp(complex64_t* poles, int* numpoles, complex64_t* zeros, int* numzeros, real64_t* gain, const real64_t* bwidth, const real64_t* Wn)
{
  int order = *numpoles;                // need original order for the gain approximation.
  complex64_t bw2 = (*bwidth)/2 + 0.0*I;  // half bandwidth, bw/2.
  complex64_t bsq, bsum, pold, cWn;       // intermediate complex values for warp calculations.
  cWn = (real64_t)*Wn + 0.0*I;
  complex64_t Wn2 = compmult(cWn, cWn); // complex Wn^2.
  
  // p * bw/2 for the original data.
  for (int pInd = 0; pInd < (*numpoles); pInd++)
  {
    poles[pInd] = (complex64_t)compmult(poles[pInd], bw2);  
  }
  
  // [1,2,3,4] -> [1, 0, 2, 0, 3, 0, 4, 0]. expand the array.
  int oldInd = (*numpoles)-1;
  for (int pInd = 2*(*numpoles)-2; pInd > 0; pInd-=2)
  {
    poles[pInd] = poles[oldInd];
    poles[pInd+1] = (complex64_t)0.0 + 0.0*I;
    poles[oldInd--] = (complex64_t)0.0 + 0.0*I;
  }
  
  // compute the p +/- sqrt(p^2 - Wn^2) for each group.
  for (int pInd = 0; pInd < (*numpoles); pInd++)
  {
    // extract pole intermediate values for BPF warping.
    pold = poles[2*pInd];
    bsum = (complex64_t)compsub(compmult(pold, pold), Wn2); // p^2 - Wn^2
    bsq = (complex64_t)csqrt(bsum); // sqrt(p^2 - Wn^2)
    
    // p(i)   = p(i) + sqrt(p(i)^2 - Wn^2).
    // p(i+1) = p(i) - sqrt(p(i)^2 - Wn^2).
    
    poles[2*pInd] = (complex64_t)compadd(pold, bsq);
    poles[2*pInd+1] = (complex64_t)compsub(pold,  bsq);
  }
  
  // indicate the new number of poles from this process.
  *numpoles = (int)2*(*numpoles);

  // k = k * bw^(order)
  *gain *= (real64_t)pow(*bwidth, order);
}


/************************************************************
 *  Function polesort                                       *
 *    Sorts poles based on proximity to the unit circle.    *
 *  Inputs:                                                 *
 *    poles (complex64_t*) is a pointer to the pole array.    *
 *    numpoles (int*) are the number of poles in the array. *
 *    order (int*) indicates the order of the filter. If    *
 *      odd, the last two poles are real from               *
 *      BPF warping.                                        *
 ************************************************************/

void polesort(complex64_t* poles, int numpoles, int order)
{
  int len = (order % 2) ? numpoles-2 : numpoles;
  
  // bubble sorting the poles into descending distance from unit circle.
  complex64_t temp;
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
 *    poles (complex64_t*) are the seeded poles.                *
 *    numpoles (int*) is the number of poles in the array.    *
 *    zeros (complex64_t*) is a pointer to the zeros array.     *
 *    numzeros (int*) is the number of zeros in the array.    *
 *    gain (real64_t*) is the gain of the filter.            *
 *    bwidth (real64_t*) is the bandwidth of the stop region.*
 *    Wn (real64_t*) is the normalized center frequency.     *
 * Outputs:                                                   *
 *    numpoles will be modified to double in size.            *
 *    numzeros will be the same length as numpoles with       *
 *      repeating complex conjugates.                         *
 **************************************************************/

void bsfwarp(complex64_t* poles, int* numpoles, complex64_t* zeros, real64_t* gain, const real64_t* bwidth, const real64_t* Wn)
{
  complex64_t bw2 = (*bwidth / 2.0) + 0.0*I;
  complex64_t pold, bsum, bsq;
  complex64_t Wn2 = compmult((complex64_t)*Wn+0.0*I, (complex64_t)*Wn+0.0*I); // complex Wn^2.
  
  // gain shift is just 1, since k * real(prod(-z)/prod(-p)) = 1
  // when the poles are around the unit circle and no zeros.
  *gain = (real64_t)1.0;
  
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
    poles[pInd+1] = (complex64_t)0.0 + 0.0*I;
    poles[oldInd--] = (complex64_t)0.0 + 0.0*I;
  }
  
  // zeros already same length as poles, (0+/-j*Wn).
  for (int zInd = 0; zInd < *numpoles; zInd++)
  {
    zeros[2*zInd] = (complex64_t)0.0 + (*Wn)*I;
    zeros[2*zInd+1] = (complex64_t)0.0 - (*Wn)*I;
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
 *    poles (complex64_t*) is the complex pole array to warp.             *
 *    npoles (int*) is the number of poles.                             *
 *    zeros (complex64_t*) is the complex zero array to warp.             *
 *    nzeros (int*) are the number of zeros.                            *
 *    gain (real64_t*) is the pointer to the gain of the filter.       *
 *                                                                      *
 ************************************************************************/

void bilinear_band_s2z(complex64_t* poles, const int* numpoles, complex64_t* zeros, int* numzeros, real64_t* gain, const real64_t* fs)
{
  complex64_t fs2 = (real64_t)2*(*fs) + 0.0*I;  // 2/T equivalent.
  complex64_t pos1 = (real64_t)1.0 + 0.0*I;  // bilinear operation for (z+1)/(z-1)
  complex64_t pgain = (real64_t)1.0 + 0.0*I;
  complex64_t zgain = (real64_t)1.0 + 0.0*I; // zero singularity gain also starts at 1.
  complex64_t num, den, warp, sub;            // intermediate values.

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
  *gain *= (real64_t)creal(compdiv(zgain, pgain)); // calculates gain change.
  
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
 *  Flo (real64_t) is the lower corner frequency.             *
 *  Fhi (real64_t) is the upper corner frequency.             *
 *  Fs (real64_t) is the sampling rate to design at.          *
 *  type (int) indicates 0=BPF, 1=BSF.                         *
 * Outputs:                                                    *
 *  matrix (real64_t*) is a pointer to the SOS matrix.        *
 ***************************************************************/

real64_t* butterband(const int order, real64_t flo, real64_t fhi, real64_t fs, const int type)
{
  
  // check if the corner frequencies need to swap, or if they're equal.
  real64_t swap;
  if (flo > fhi)
  {
    swap = flo;
    flo = fhi;
    fhi = swap;
  }
  else if (flo == fhi) 
  {
    printf("ERROR: corner frequencies are equivalent.\n");
    return (real64_t*)NULL;
  }

  // check if the upper frequency exceeds nyquist.
  if (fhi >= fs/2.0)
  {
    printf("ERROR: upper corner frequency exceeds Nyquist.\n");
    return (real64_t*)NULL;
  }

  int npoles, nzeros = 0;
  real64_t* mat;
  complex64_t* poles = (complex64_t*)NULL;
  complex64_t* zeros = (complex64_t*)NULL;
  complex64_t* spoles = (complex64_t*)NULL;
  real64_t gain, w1, w2, bwidth, Wn;

  gain = 1.0;
  w1 = 2.0 * (flo/fs);
  w2 = 2.0 * (fhi/fs);

  int N = (int)ceil((real64_t)order/2);
  w1 = (real64_t)tan(HALF_PI * w1) * 4.0;
  w2 = (real64_t)tan(HALF_PI * w2) * 4.0;

  bwidth = w2 - w1;     // normalized bandwidth
  Wn = sqrt(w1 * w2);   // center frequency
  fs = 2.0;             // normalized Wq = 1, resetting Fs.

  // seed the poles around the unit circle.
  spoles = seedpoles(order, &npoles);

  // make the SOS matrix.
  mat = (real64_t*)malloc(sizeof(real64_t)*N_SOSCOEFFS*order);

  // double the length for the BPF/BSF to manage.
  poles = (complex64_t*)malloc(sizeof(complex64_t)*2*npoles);
  for (int pInd = 0; pInd < npoles; pInd++)
  {
      poles[pInd] = spoles[pInd];
  }
  free(spoles);
    
  if (!type)
  {

    // bandpass warp.
    bpfwarp(poles, &npoles, zeros, &nzeros, &gain, &bwidth, &Wn);

    // since no zeros were seeded, create S-plane zeros, (making same length as poles.)
    zeros = (complex64_t*)malloc(sizeof(complex64_t)*npoles);
    nzeros = npoles/2.0;

    // bilinear warps poles to [+1,+1,....-1,-1].
    for (int zInd = 0; zInd < npoles; zInd++)
    {
      zeros[zInd] = zInd < nzeros ? (complex64_t)0.0 + 0.0*I : (complex64_t)-1.0 + 0.0;
    }
  }
  else
  {
    // bandstop filter warp. BSF will place complex conjugates in the zeros array.
    nzeros = npoles*2;
    zeros = (complex64_t*)malloc(sizeof(complex64_t)*nzeros);
    bsfwarp(poles, &npoles, zeros, &gain, &bwidth, &Wn);
    npoles *= 2; // only change the number after having used the loop in BSF warp.

  }

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
    pInd = (npoles-2) - 2*stgInd;  // reverse pole ordering.
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

  // encorporate the gain into the first stage.
  mat[0] *= gain;
  mat[1] *= gain;
  mat[2] *= gain;

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
 *    poles (complex64_t*) are the poles that are being re-mapped.      *
 *    numpoles (int) are the number of poles in the pole array.       *
 *    numzeros (int) indicates if zeros were generated based on       *
 *      the filter type, (LP/HP/APF).                                 *
 *    gain (real64_t*) is a scalar gain that will be updated based   *
 *      on the transformation.                                        *
 *    fs (real64_t) is the sampling rate we're mapping to.           *
 *                                                                    *
 **********************************************************************/

void bilinear_s2z(complex64_t* poles, const int numpoles, const int numzeros, real64_t* gain, const real64_t fs)
{
  complex64_t fs2 = (real64_t)2*fs + 0.0*I;  // 2/T equivalent in bilinear xform.
  complex64_t pgain = (real64_t)1.0 + 0.0*I; // complex +1, (starting value for pole gain shift.)
  complex64_t num, den, warp, psub;           // intermediate values for bilinear xform.
  real64_t zgain;                          // gain-shift for zeros.
  
  // gain adjustment from poles is "prod(Fs-p)", pre-bilinear xform.
  for (int pInd = 0; pInd < numpoles; pInd++)
  {
    // (a+bi) - (c+di) = (a - c) + (b - d)*i. Fs-p -> (Fs-real(p)) + (0 - imag(p))i, since Fs is real.
    psub = (complex64_t)((creal(fs2) - creal(poles[pInd])) + (cimag(fs2) - cimag(poles[pInd]))*I);
    pgain = compmult(pgain, psub);
  }
  
  // bilinear transform, S to Z plane.
  for (int pInd = 0; pInd < numpoles; pInd++)
  {
    // xform => (1+p/Fs)./(1-p/Fs), or (1+(p/Fs))./(1-(p/Fs)).
    // (a + c) + (b + d)i, 1 + (p/Fs) -> (1 + Re(p/Fs)) + (0 + Im(p/Fs))i
    // (a - c) + (b - d)i, 1 - (p/Fs) -> (1 - Re(p/Fs)) + (0 - Im(p/Fs))i
    warp = compdiv(poles[pInd], fs2);
    num = ((real64_t)1.0 + creal(warp)) + cimag(warp)*I;
    den = ((real64_t)1.0 - creal(warp)) + ((real64_t)0.0 - cimag(warp))*I;
    poles[pInd] = compdiv(num, den);
  }
  
  // check HPF vs. LPF/APF design cases.
  if (numzeros != 0)
  {
    /* for the HPF case, zeros on S-plane move from 0 to 1 on the Z-plane. Then,
     * for determining the gain adjustment of the zeros, there are an equal number
     * of zeros and poles. This means the "prod(Fs-z|s)" term is simply fs2^npoles.
     */
    zgain = pow((real64_t)2*fs, (real64_t)numpoles);
  }
  else
  {
    // the lowpass/allpass case has no zeros generated; the "prod(Fs-z)" term is simply 1.
    zgain = 1.0;
  }

  // k = k * real(prod(Fs-z)/prod(Fs-p))
  *gain *= creal(compdiv((complex64_t)(zgain + 0.0*I), pgain));
}

/****************************************************************
 * Function butter                                              *
 *    This is the principal entry point for designing the       *
 *    butterworth SOS filter matrix.                            *
 *                                                              *
 * Inputs:                                                      *
 *  order (int) is the butterworth filter order.                *
 *  Fc (real64_t) is the corner frequency (-3 dB point.)       *
 *  Fs (real64_t) is the discrete sampling rate to design at.  *
 *  type (int) indicates 0=LPF, 1=HPF, 2=APF.                   *
 * Outputs:                                                     *
 *  matrix (real64_t*) is a pointer to the SOS matrix.         *
 ****************************************************************/

real64_t* butter(const int order, const real64_t fc, real64_t fs, const int type)
{

  if (fc >= fs/2) {
    printf("ERROR: corner frequency must be less than Nyquist rate.\n");
    return (real64_t*)NULL;
  }
  if ((type < 0) || (type > 2)) {
    printf("ERROR: Unrecognized filter type. 0=LPF, 1=HPF, 2=APF.\n");
    return (real64_t*)NULL;
  }

  int N, npoles, nzeros;
  real64_t* mat;
  complex64_t* poles;
  real64_t gain, zeros, omega, warp;

  gain = 1.0;             // starting with unity gain.
  omega = 2.0 * (fc/fs);  // normalizing frequency.
  N = (int)ceil((real64_t)order/2);
  warp = (real64_t)tan(HALF_PI * omega) * 4.0;

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

/*************************************************************
 * Function jl_butter                                        *
 *  Julia API that requires address of pre-allocated memory. *
 *************************************************************/ 

void jl_butter(const int order, const double fc, const double fs, const int type, double *mat)
{
  int N, ncoeffs;
  real64_t *sos;
  sos = butter(order, fc, fs, type);
  N = (order % 2) ? (order-1)/2 + 1 : order/2;
  ncoeffs = N * N_SOSCOEFFS;
  for (int iC = 0; iC < ncoeffs; iC++) {
    mat[iC] = (double)sos[iC];
  }
  free(sos);
}

/*************************************************************
 * Function jl_butterband                                    *
 *  Julia API that requires address of pre-allocated memory. *
 *************************************************************/
void jl_butterband(const int order, const double flo, const double fhi, const double fs, const int type, double *mat)
{
  int ncoeffs;
  real64_t *sos;
  sos = butterband(order, flo, fhi, fs, type);
  ncoeffs = order * N_SOSCOEFFS; // bandpass/bandstop order is the same number of biquads.
  for (int iC = 0; iC < ncoeffs; iC++) {
    mat[iC] = (double)sos[iC];
  }
  free(sos);
}
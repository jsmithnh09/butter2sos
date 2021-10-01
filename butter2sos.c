/*
 *==========================================
 * butter2sos.c
 *
 * Primary entrypoint for designing 
 * butterworth filters.	
 *
 * Make: gcc -std=c99 buttersos_design.c buttersos.c -o buttersos
 *
 *==========================================
 */

/*
 * Contact: Jordan R. Smith
 */

#include "butter2sos_design.h"

int main(int argc, char const *argv[])
{
	int order, type, N;
	regular_t fs, fc;
	regular_t* mat;
  
  type = 0;               // 0=LPF, 1=HPF, 2=APF.
	fs = 44100.0;           // sampling rate.
	fc = 300.0;             // 300 Hz
  order = 7;              // order specification.
  N = (order % 2) ? (order+1)/2 : order/2;

  // check filter type.
  if (type == 0) {
    printf("# %-6s: %s\n", "type", "lowpass");
  }
  else if (type == 1) {
    printf("# %-6s: %s\n", "type", "highpass");
  }
  else if (type == 2) {
    printf("# %-6s: %s\n", "type", "allpass");
  }
  printf("# %-6s: %g Hz\n", "Fs", fs);
  printf("# %-6s: %g Hz\n", "Fc", fc);
  printf("# %-6s: %d\n", "order", order);
  printf("\n");

  mat = butter(order, fc, fs, type);
  if ((void *)mat == NULL)
  {
    printf("Unable to construct the SOS matrix.\n");
    return 1;
  }
  
  // print the results and clean up.
  printf("SAMPLE RATE: %g\n\n", fs);
  printf("COEFFICIENTS:\n\n");
  printsosmatrix(mat, N); 
  free(mat);
	return 0;
}

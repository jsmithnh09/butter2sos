/*
 *==========================================
 * butter2sos.c
 *
 * Primary entrypoint for designing 
 * butterworth filters.	
 *
 * ./butter2sos -n<order> -r<sampling rate> -c<corner frequency> -t<filter type (0, 1, 2)>
 *
 * Make: gcc -std=c99 buttersos_design.c buttersos.c -o buttersos
 *
 *==========================================
 */

/*
 * Contact: Jordan R. Smith
 */

#include "butter2sos_design.h"
#include "version.h"
#ifdef _WIN32
  #include <io.h>
#else
  #include <unistd.h>
#endif

int main(int argc, char **argv)
{
  int order = 0;
  int type = 0;
  int N = 0;
  real64_t fs = 1.0; 
  real64_t fc = 0.5;
  real64_t *mat;
  
  #ifdef _WIN32
    static char usage[] = "%s v%s: <order> <corner_frequency> <sample_rate> <type 0,1,2>\n";
    if (argc == 1) {
      printf(usage, argv[0], BUTTER_SOS_VERSION);
      exit(0);
    } else if (argc != 5) {
      fprintf(stderr, "Invalid # of input arguments.\n");
      exit(1);
    }
    order = atoi(argv[1]);
    fc = atof(argv[2]);
    fs = atof(argv[3]);
    type = atoi(argv[4]);
  #else
    int c;
    static char usage[] = "%s v%s: -n <order> -c <cornerfrequency> -r <samplerate> -t <type 0,1,2>\n";
    if (argc == 1) {
      printf(usage, argv[0], BUTTER_SOS_VERSION);
      exit(0);
    }
    while ((c = getopt (argc, argv, "n:c:r:t:")) != -1) {
      switch (c) {
        case 'n':
          order = atoi(optarg);
          break;
        case 'r':
          fs = atof(optarg);
          break;
        case 'c':
          fc = atof(optarg);
          break;
        case 't':
          type = atoi(optarg);
          break;
        case '?':
          fprintf(stderr, "ERROR: unknown argument %c provided\n", optopt);
          exit(1);
        default:
          printf(usage, argv[0]);
          exit(1);
        }
    }
  #endif

  // confirm the input parameters are suitable.
  if (order <= 0) {
    fprintf(stderr, "ERROR: specified order must be a positive integer number.\n");
    exit(1);
  }
  else if (fc >= fs/2) {
    fprintf(stderr, "ERROR: specified corner frequency exceeds Nyquist.\n");
    exit(1);
  }
  else if (fs <= 0) {
    fprintf(stderr, "ERROR: specified sampling rate must be a positive value.\n");
    exit(1);
  }

  // even/odd order spec. defines # biquads
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
  } else {
    fprintf(stderr, "ERROR: Invalid filter type. 0=LPF 1=HPF 2=APF.\n");
    exit(1);
  }
  printf("# %-6s: %g Hz\n", "Fs", fs);
  printf("# %-6s: %g Hz\n", "Fc", fc);
  printf("# %-6s: %d\n", "order", order);
  printf("\n");

  mat = butter(order, fc, fs, type);
  if ((void *)mat == NULL)
  {
    fprintf(stderr, "Unable to construct the SOS matrix.\n");
    exit(1);
  }
  
  // print the results and clean up.
  printf("SAMPLE RATE: %g\n\n", fs);
  printf("COEFFICIENTS:\n\n");
  printsosmatrix(mat, N); 
  free(mat);
	exit(0);
}

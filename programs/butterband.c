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
#ifdef _WIN32
  #include <io.h>
#else
  #include <unistd.h>
#endif

int main(int argc, char **argv)
{
  int order, type, N = 0;
  real64_t fs, fc1, fc2;
  real64_t *mat;
  #ifdef _WIN32
    static char usage[] = "%s: <order> <lowband> <highband> <samplerate> <type 0,1>\n";
    if (argc == 1) {
      printf(usage, argv[0]);
      exit(0);
    }
    if (argc != 6) {
      fprintf(stderr, usage);
      fprintf(stderr, "ERROR: Invalid # of input arguments.\n");
      exit(1);
    }
    order = atoi(argv[1]);
    fc1 = atof(argv[2]);
    fc2 = atof(argv[3]);
    fs = atof(argv[4]);
    type = atoi(argv[5]);
  #else
    int c;
    static char usage[] = "%s: -n <order> -l <lowband> -h <highband> -r <samplerate> -t <type 0,1>\n";
    while ((c = getopt(argc, argv, "n:l:h:r:t:")) != -1) {
      switch (c) {
        case 'n':
          order = atoi(optarg);
          break;
        case 'r':
          fs = atof(optarg);
          break;
        case 'l':
          fc1 = atof(optarg);
          break;
        case 'h':
          fc2 = atof(optarg);
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
  else if ((fc1 >= fs/2) || (fc2 >= fs/2)) {
    fprintf(stderr, "ERROR: specified corner frequenc(ies) exceeds Nyquist.\n");
    exit(1);
  }
  else if (fs <= 0) {
    fprintf(stderr, "ERROR: specified sampling rate must be a positive value.\n");
    exit(1);
  }
  if (fc1 >= fc2) {
    fprintf(stderr, "ERROR: the lowband must be less than the highband.\n");
    exit(1);
  }

  // order matches number of biquads.
  N = order;

  // check filter type.
  if (type == 0) {
    printf("# %-8s: %s\n", "type", "bandpass");
  }
  else if (type == 1) {
    printf("# %-8s: %s\n", "type", "bandstop");
  } else {
    fprintf(stderr, "ERROR: Invalid filter type. 0=BPF, 1=BSF.\n");
    exit(1);
  }
  printf("# %-8s: %g Hz\n", "Fs", fs);
  printf("# %-8s: %g Hz\n", "lowband", fc1);
  printf("# %-8s: %g Hz\n", "highband", fc2);
  printf("# %-8s: %d\n", "order", order);
  printf("\n");

  mat = butterband(order, fc1, fc2, fs, type);
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

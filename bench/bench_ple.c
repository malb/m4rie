#include <m4rie/m4rie.h>
#include <cpucycles.h>
#include "benchmarking.h"

struct ple_params {
  rci_t k; 
  rci_t m;
  rci_t n;
  rci_t c;
  rci_t r;
  char const *algorithm;  
};

int run(void *_p, unsigned long long *data, int *data_len) {
  struct ple_params *p = (struct ple_params *)_p;
  *data_len = 2;

  gf2e *ff = gf2e_init(irreducible_polynomials[p->k][1]);
  mzed_t *A = mzed_init(ff, p->m, p->n);
  mzed_randomize(A);

  mzp_t *P = mzp_init(p->m);
  mzp_t *Q = mzp_init(p->n);

  data[1] = cpucycles();
  data[0] = walltime(0.0);

  if(strcmp(p->algorithm,"default")==0)
    p->r = _mzed_ple(A, P, Q, p->c);
  else if(strcmp(p->algorithm,"newton-john")==0)
    p->r = mzed_ple_newton_john(A, P, Q);
  else if(strcmp(p->algorithm,"naive")==0)
    p->r = mzed_ple_naive(A, P, Q);
  else
    p->r = mzed_echelonize(A, 1);
  data[1] = cpucycles() - data[1];
  data[0] = walltime(data[0]);

  mzed_free(A);
  mzp_free(P);
  mzp_free(Q);
  gf2e_free(ff);
  return 0;
}

void print_help() {
  printf("bench_ple:\n\n");
  printf("REQUIRED\n");
  printf("  e -- integer between 2 and 10\n");
  printf("  m -- integer > 0\n");
  printf("  n -- integer > 0\n");
  printf("  algorithm -- default\n");
  printf("               newton-john\n");
  printf("               naive\n");
  printf("  c -- cutoff (for 'default')\n");
  printf("\n");
  bench_print_global_options(stdout);
}

int main(int argc, char **argv) {
  global_options(&argc, &argv);

  if (argc < 4) {
    print_help();
    m4ri_die("");
  }

  struct ple_params params;

  params.k = atoi(argv[1]);
  params.m = atoi(argv[2]);
  params.n = atoi(argv[3]);
  if (argc >= 5)
    params.algorithm = argv[4];
  else
    params.algorithm = (char*)"default";
  if (argc >= 6)
    params.c = atoi(argv[5]);
  else
    params.c = 0;
  if(argc >= 7) {
    print_help();
    m4ri_die("");
  }

  srandom(17);
  unsigned long long data[2];
  run_bench(run, (void*)&params, data, 2);

  double cc_per_op = ((double)data[1])/ ( (double)params.m * (double)params.n * powl((double)params.r,__M4RIE_OMEGA-2) );

  printf("e: %2d, m: %5d, n: %5d, algorithm: %10s, cutoff: %10d, cpu cycles: %10llu, cc/(mnr^0.807): %.5lf, wall time: %lf\n", params.k, params.m, params.n, params.algorithm, params.c, data[1], cc_per_op, data[0] / 1000000.0);
}

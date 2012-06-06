#include <gf2e_cxx/finite_field_givaro.h>
#include <m4rie/m4rie.h>
#include <cpucycles.h>
#include "benchmarking.h"

using namespace M4RIE;

struct elim_params {
  rci_t k; 
  rci_t m;
  rci_t n;
  rci_t r;
  char const *algorithm;  
  char type;
};

int run_mzed(void *_p, unsigned long long *data, int *data_len) {
  struct elim_params *p = (struct elim_params *)_p;
  *data_len = 2;

  FiniteField *F = (FiniteField*)(new GFqDom<int>(2,p->k));
  gf2e *ff = gf2e_init_givgfq(F);

  mzed_t *A = mzed_init(ff,p->m,p->n);
  mzed_randomize(A);

  data[0] = walltime(0);
  data[1] = cpucycles();

  if(strcmp(p->algorithm,"newton-john")==0)
    p->r=  mzed_echelonize_newton_john(A, 1);
  else if(strcmp(p->algorithm,"naive")==0)
    p->r = mzed_echelonize_naive(A, 1);
  else if(strcmp(p->algorithm,"ple")==0)
    p->r = mzed_echelonize_ple(A, 1);
  else if(strcmp(p->algorithm,"default")==0)
    p->r = mzed_echelonize(A, 1);
  else
    m4ri_die("uknown algorithm '%s'\n.",p->algorithm);
  data[1] = cpucycles() - data[1];
  data[0] = walltime(data[0]);

  mzed_free(A);
  gf2e_free(ff);
  delete F;
  return 0;
}

int run_mzd_slice(void *_p, unsigned long long *data, int *data_len) {
  struct elim_params *p = (struct elim_params *)_p;
  *data_len = 2;

  FiniteField *F = (FiniteField*)(new GFqDom<int>(2,p->k));
  gf2e *ff = gf2e_init_givgfq(F);

  mzd_slice_t *A = mzd_slice_init(ff,p->m,p->n);
  mzd_slice_randomize(A);

  data[0] = walltime(0);
  data[1] = cpucycles();

  if(strcmp(p->algorithm,"newton-john")==0) {

    mzed_t *B = mzed_cling(NULL, A);
    p->r=  mzed_echelonize_newton_john(B, 1);
    mzed_slice(A, B);
    mzed_free(B);

  } else if(strcmp(p->algorithm,"naive")==0) {

    mzed_t *B = mzed_cling(NULL, A);
    p->r = mzed_echelonize_naive(B, 1);
    mzed_slice(A, B);
    mzed_free(B);

  } else if(strcmp(p->algorithm,"ple")==0) {

    p->r = mzd_slice_echelonize_ple(A, 1);

  } else if(strcmp(p->algorithm,"default")==0) {

    p->r = mzd_slice_echelonize(A, 1);

  } else {

    m4ri_die("uknown algorithm '%s'\n.",p->algorithm);

  }
  data[1] = cpucycles() - data[1];
  data[0] = walltime(data[0]);

  mzd_slice_free(A);
  gf2e_free(ff);
  delete F;
  return 0;
}

void print_help() {
  printf("bench_elimination:\n\n");
  printf("REQUIRED\n");
  printf("  e -- integer between 2 and 10\n");
  printf("  m -- integer > 0, number of rows\n");
  printf("  n -- integer > 0, number of columns\n");
  printf("  algorithm -- default -- let M4RIE decide\n");
  printf("               naive -- cubic Gaussian elimination\n");
  printf("               newton-john -- Newton-John tables\n");
  printf("               ple -- PLE based\n");
  printf(" type -- mzed_t or mzd_slice_t (default: mzed_t)\n");
  printf("\n");
  bench_print_global_options(stdout);
}

int main(int argc, char **argv) {
  global_options(&argc, &argv);

  if (argc < 4) {
    print_help();
    m4ri_die("");
  }

  struct elim_params params;

  params.k = atoi(argv[1]);
  params.m = atoi(argv[2]);
  params.n = atoi(argv[3]);
  if (argc >= 5)
    params.algorithm = argv[4];
  else
    params.algorithm = (char*)"default";
  if (argc >= 6) {
    if (strcmp("mzed_t",argv[5]) == 0) 
      params.type = 0;
    else if (strcmp("mzd_slice_t",argv[5]) == 0) 
      params.type = 1;
    else
      m4ri_die("unknown type '%s'\n",argv[5]);
  } else {
    params.type = 0;
  }

  srandom(17);
  unsigned long long data[2];
  if (params.type == 0)
    run_bench(run_mzed, (void*)&params, data, 2);
  else
    run_bench(run_mzd_slice, (void*)&params, data, 2);

  double cc_per_op = ((double)data[1])/ ((double)params.m * (double)params.n * powl((double)params.r,__M4RIE_OMEGA-2) );

  printf("e: %2d, m: %5d, n: %5d, type: %d, algo: %10s, cpu cycles: %10llu, cc/(mnr^0.807): %.5lf, wall time: %lf\n", params.k, params.m, params.n, params.type, params.algorithm, data[1], cc_per_op, data[0] / 1000000.0);
}

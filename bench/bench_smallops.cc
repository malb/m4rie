#include <gf2e_cxx/finite_field_givaro.h>
#include <m4rie/m4rie.h>
#include <cpucycles.h>
#include "benchmarking.h"

using namespace M4RIE;

struct smallops_params {
  rci_t k; 
  rci_t m;
  rci_t n;
  char const *algorithm;
};

int run_mzed_add(void *_p, unsigned long long *data, int *data_len) {
  struct smallops_params *p = (struct smallops_params *)_p;
  *data_len = 2;

  FiniteField *F = (FiniteField*)(new GFqDom<int>(2,p->k));
  gf2e *ff = gf2e_init_givgfq(F);

  mzed_t *A = mzed_init(ff,p->m,p->n);
  mzed_randomize(A);
  mzed_t *B = mzed_init(ff,p->m,p->n);
  mzed_randomize(B);
  mzed_t *C = mzed_init(ff,p->m,p->n);

  data[0] = walltime(0);
  data[1] = cpucycles();

  mzed_add(C, A, B);

  data[1] = cpucycles() - data[1];
  data[0] = walltime(data[0]);

  mzed_free(A);
  mzed_free(B);
  mzed_free(C);

  gf2e_free(ff);
  delete F;
  return 0;
}


int run_mzed_slice(void *_p, unsigned long long *data, int *data_len) {
  struct smallops_params *p = (struct smallops_params *)_p;
  *data_len = 2;

  FiniteField *F = (FiniteField*)(new GFqDom<int>(2,p->k));
  gf2e *ff = gf2e_init_givgfq(F);

  mzed_t *A = mzed_init(ff,p->m,p->n);
  mzed_randomize(A);
  mzd_slice_t *a = mzd_slice_init(ff,p->m,p->n);

  data[0] = walltime(0);
  data[1] = cpucycles();

  mzed_slice(a, A);

  data[1] = cpucycles() - data[1];
  data[0] = walltime(data[0]);

  mzed_free(A);
  mzd_slice_free(a);

  gf2e_free(ff);
  delete F;
  return 0;
}

int run_mzed_cling(void *_p, unsigned long long *data, int *data_len) {
  struct smallops_params *p = (struct smallops_params *)_p;
  *data_len = 2;

  FiniteField *F = (FiniteField*)(new GFqDom<int>(2,p->k));
  gf2e *ff = gf2e_init_givgfq(F);

  mzd_slice_t *a = mzd_slice_init(ff,p->m,p->n);
  mzd_slice_randomize(a);
  mzed_t *A = mzed_init(ff, p->m, p->n);

  data[0] = walltime(0);
  data[1] = cpucycles();

  mzed_cling(A, a);

  data[1] = cpucycles() - data[1];
  data[0] = walltime(data[0]);

  mzed_free(A);
  mzd_slice_free(a);

  gf2e_free(ff);
  delete F;
  return 0;
}

void print_help() {
  printf("bench_smallops:\n\n");
  printf("REQUIRED\n");
  printf("  e -- integer between 2 and 10\n");
  printf("  m -- integer > 0\n");
  printf("  n -- integer > 0\n");
  printf("  what -- mze_cling\n");
  printf("          mzed_slice\n");
  printf("          mzed_add\n");
  printf("\n");
  bench_print_global_options(stdout);
}


int main(int argc, char **argv) {
  global_options(&argc, &argv);

  if (argc < 5) {
    print_help();
    m4ri_die("");
  }

  struct smallops_params params;

  params.k = atoi(argv[1]);
  params.m = atoi(argv[2]);
  params.n = atoi(argv[3]);

  srandom(17);
  unsigned long long data[2];

  if(strcmp(argv[4],"mzed_slice") == 0) {
    run_bench(run_mzed_slice, (void*)&params, data, 2);
  } else if(strcmp(argv[4],"mzed_cling") == 0) {
    run_bench(run_mzed_cling, (void*)&params, data, 2);
  } else if(strcmp(argv[4],"mzed_add") == 0) {
    run_bench(run_mzed_add, (void*)&params, data, 2);
  }


  double cc_per_op = ((double)data[1])/ ( (double)params.m * (double)params.n );

  printf("%s: m: %5d, n: %5d, cpu cycles: %10llu, cc/(mn): %.5lf, wall time: %lf\n", argv[4], params.m, params.n, data[1], cc_per_op, data[0] / 1000000.0);
}



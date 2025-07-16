#include <m4rie/m4rie.h>
#include "cpucycles.h"
#include "benchmarking.h"

struct transpose_params {
  rci_t e;
  rci_t m;
  rci_t n;
  char const *algorithm;
};

int run_mzed (void *_p, unsigned long long *data, int *data_len) {
  struct transpose_params *p = (struct transpose_params *)_p;
  *data_len = 2;

  gf2e *ff = gf2e_init(irreducible_polynomials[p->e][1]);
  
  mzed_t *A = mzed_init(ff,p->m,p->n);
  mzed_randomize(A);
  mzed_t *B = mzed_init(ff,p->n,p->m);

  data[0] = walltime(0);
  data[1] = cpucycles();
  
  if (strcmp(p->algorithm, "optimised") == 0) {
    B = mzed_transpose(B, A);
  } else if (strcmp(p->algorithm, "naive") == 0) {
    for (rci_t row = 0; row < A->nrows; row++) {
      for (rci_t col = 0; col < A->ncols; col++) {
        mzed_write_elem(B, col, row, mzed_read_elem(A, row, col));
      }
    }
  } else {
    m4ri_die("unknown algorithm '%s'\n.",p->algorithm);
  }
  
  B = mzed_transpose(B, A);

  data[1] = cpucycles() - data[1];
  data[0] = walltime(data[0]);

  mzed_free(A);
  mzed_free(B);
  gf2e_free(ff);
  return 0;
}

void print_help() {
  printf("bench_elimination:\n\n");
  printf("REQUIRED\n");
  printf("  e -- integer between 2 and 16\n");
  printf("  m -- integer > 0, number of rows\n");
  printf("  n -- integer > 0, number of columns\n");
  printf("  algorithm -- optimised -- optimised transpose\n");
  printf("               naive -- transfer element by element\n");
  printf("\n");
  bench_print_global_options(stdout);
}

int main(int argc, char **argv) {
  global_options(&argc, &argv);
  
  if (argc < 4 || argc > 5) {
    print_help();
    m4ri_die("");
  }
  
  struct transpose_params params;

  params.e = atoi(argv[1]);
  params.m = atoi(argv[2]);
  params.n = atoi(argv[3]);
  if (argc >= 5)
    params.algorithm = argv[4];
  else
    params.algorithm = (char*)"optimised";

  srandom(17);
  unsigned long long data[2];
  run_bench(run_mzed, (void*)&params, data, 2);

  double cc_per_bit = (((double)data[1])/ ((double)params.m * (double)params.n * (double)params.e) );

  printf("e: %2d, m: %5d, n: %5d, algo: %10s, cpu cycles: %10llu, cc/bit: %.5lf, wall time: %lf\n", params.e, params.m, params.n, params.algorithm, data[1], cc_per_bit, data[0] / 1000000.0);
}
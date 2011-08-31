#include <gf2e_cxx/finite_field_givaro.h>
#include <m4rie/m4rie.h>
#include <cpucycles.h>
#include "benchmarking.h"

using namespace M4RIE;


struct mul_params {
  rci_t k; 
  rci_t m;
  rci_t n;
  char const *algorithm;  
};

int run(void *_p, unsigned long long *data, int *data_len) {
  struct mul_params *p = (struct mul_params *)_p;
  *data_len = 2;

  FiniteField *F = (FiniteField*)(new GFqDom<int>(2,p->k));
  gf2e *ff = gf2e_init_givgfq(F);

  mzed_t *A = mzed_init(ff,p->m,p->n);
  mzed_randomize(A);
  mzed_t *B = mzed_init(ff,p->n,p->m);
  mzed_randomize(B);

  mzed_t *C;

  data[0] = walltime(0);
  data[1] = cpucycles();

  if(strcmp(p->algorithm, "travolta")==0)
    C = mzed_mul_travolta(NULL, A, B);
  else if(strcmp(p->algorithm,"naive")==0)
    C = mzed_mul_naive(NULL, A, B);
  else if(strcmp(p->algorithm,"strassen")==0)
    C = mzed_mul_strassen(NULL, A, B,_mzed_strassen_cutoff(NULL, A, B));
  else if(strcmp(p->algorithm,"karatsuba")==0)
    C = mzed_mul_karatsuba(NULL, A, B);
  else
    C = mzed_mul(NULL, A, B);

  data[1] = cpucycles() - data[1];
  data[0] = walltime(data[0]);

  mzed_free(A);
  mzed_free(B);
  mzed_free(C);

  gf2e_free(ff);
  delete F;
  return 0;
}

int main(int argc, char **argv) {
  global_options(&argc, &argv);

  if (argc < 4)
    m4ri_die("syntax: e m n alg\n");

  struct mul_params params;

  params.k = atoi(argv[1]);
  params.m = atoi(argv[2]);
  params.n = atoi(argv[3]);
  if (argc == 5)
    params.algorithm = argv[4];
  else
    params.algorithm = (char*)"default";

  srandom(17);
  unsigned long long data[2];
  run_bench(run, (void*)&params, data, 2);

  double cc_per_op = ((double)data[1])/ ( (double)params.m * powl((double)params.n,1.807) );

  printf("e: %2d, m: %5d, n: %5d, algo: %10s, cpu cycles: %10llu, cc/(mn^1.807): %.5lf, wall time: %lf\n", params.k, params.m, params.n, params.algorithm, data[1], cc_per_op, data[0] / 1000000.0);
}



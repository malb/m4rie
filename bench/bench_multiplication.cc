#include <m4rie/m4rie.h>
#include <cpucycles.h>
#include "benchmarking.h"

struct mul_params {
  rci_t k; 
  rci_t m;
  rci_t n;
  char const *algorithm;  
  char type;
};

int run_mzed(void *_p, unsigned long long *data, int *data_len) {
  struct mul_params *p = (struct mul_params *)_p;
  *data_len = 2;

  gf2e *ff = gf2e_init(irreducible_polynomials[p->k][1]);

  mzed_t *A = mzed_init(ff,p->m,p->n);
  mzed_randomize(A);
  mzed_t *B = mzed_init(ff,p->n,p->m);
  mzed_randomize(B);

  mzed_t *C;

  data[0] = walltime(0);
  data[1] = cpucycles();

  if(strcmp(p->algorithm, "newton-john")==0)
    C = mzed_mul_newton_john(NULL, A, B);
  else if(strcmp(p->algorithm,"naive")==0)
    C = mzed_mul_naive(NULL, A, B);
  else if(strcmp(p->algorithm,"strassen")==0)
    C = mzed_mul_strassen(NULL, A, B,_mzed_strassen_cutoff(NULL, A, B));
  else if(strcmp(p->algorithm,"karatsuba")==0)
    C = mzed_mul_karatsuba(NULL, A, B);
  else if(strcmp(p->algorithm,"default")==0)
    C = mzed_mul(NULL, A, B);
  else
    m4ri_die("uknown algorithm '%s'\n.",p->algorithm);

  data[1] = cpucycles() - data[1];
  data[0] = walltime(data[0]);

  mzed_free(A);
  mzed_free(B);
  mzed_free(C);

  gf2e_free(ff);
  return 0;
}


int run_mzd_slice(void *_p, unsigned long long *data, int *data_len) {
  struct mul_params *p = (struct mul_params *)_p;
  *data_len = 2;

  gf2e *ff = gf2e_init(irreducible_polynomials[p->k][1]);

  mzd_slice_t *A = mzd_slice_init(ff,p->m,p->n);
  mzd_slice_randomize(A);
  mzd_slice_t *B = mzd_slice_init(ff,p->n,p->m);
  mzd_slice_randomize(B);

  mzd_slice_t *C;

  data[0] = walltime(0);
  data[1] = cpucycles();

  if(strcmp(p->algorithm,"karatsuba")==0) {
    C = mzd_slice_mul_karatsuba(NULL, A, B);
  } else if(strcmp(p->algorithm,"default")==0) {
    C = mzd_slice_mul(NULL, A, B);
  } else {
    m4ri_die("uknown algorithm '%s'\n.",p->algorithm);
  }

  data[1] = cpucycles() - data[1];
  data[0] = walltime(data[0]);

  mzd_slice_free(A);
  mzd_slice_free(B);
  mzd_slice_free(C);

  gf2e_free(ff);
  return 0;
}

void print_help() {
  printf("bench_elimination:\n\n");
  printf("REQUIRED\n");
  printf("  e -- integer between 2 and 10\n");
  printf("  m -- integer > 0, number of rows\n");
  printf("  n -- integer > 0, number of columns\n");
  printf("  algorithm -- default -- let M4RIE decide (mzed_t, mzd_slice_t)\n");
  printf("               naive -- cubic multiplication (mzed_t)\n");
  printf("               newton-john -- Newton-John tables (mzed_t) \n");
  printf("               strassen -- Strassen+Newton-John (mzed_t)\n");
  printf("               karatsuba -- Karatsuba (mzed_t)\n");
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

  struct mul_params params;

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

  double cc_per_op = ((double)data[1])/ ( (double)params.m * powl((double)params.n,1.807) );

  printf("e: %2d, m: %5d, n: %5d, type: %d, algo: %10s, cpu cycles: %10llu, cc/(mn^1.807): %.5lf, wall time: %lf\n", params.k, params.m, params.n, params.type, params.algorithm, data[1], cc_per_op, data[0] / 1000000.0);
}



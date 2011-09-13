#include <gf2e_cxx/finite_field_givaro.h>
#include <m4rie/m4rie.h>
#include <cpucycles.h>
#include "benchmarking.h"

using namespace M4RIE;

struct elim_params {
  rci_t e; 
  rci_t m;
  rci_t n;
  char const *matrix_type;  
  char const *direction;  
  char const *algorithm;  
  rci_t cutoff;
};

int run_mzed(void *_p, unsigned long long *data, int *data_len) {
  struct elim_params *p = (struct elim_params *)_p;
  *data_len = 2;

  FiniteField *F = (FiniteField*)(new GFqDom<int>(2,p->e));
  gf2e *ff = gf2e_init_givgfq(F);

  mzed_t *A = mzed_init(ff,p->m,p->m);
  mzed_randomize(A);

  const int bitmask = (1<<ff->degree)-1;
  for(rci_t i=0; i<p->m; i++) {
    while(mzed_read_elem(A, i, i) == 0) {
      mzed_write_elem(A, i, i, random()&bitmask) ;
    }
  };

  mzed_t *B = mzed_init(ff,p->m,p->n);
  mzed_randomize(B);

  data[0] = walltime(0);
  data[1] = cpucycles();
  
  if (strcmp(p->direction,"lower_left")==0) {
    if(strcmp(p->algorithm,"naive")==0)
      mzed_trsm_lower_left_naive(A, B);
    else 
      _mzed_trsm_lower_left(A, B, p->cutoff);
  }   
  data[1] = cpucycles() - data[1];
  data[0] = walltime(data[0]);
  
  mzed_free(A);
  mzed_free(B);
  gf2e_free(ff);
  delete F;
  return 0;
}

int run_mzd_slice(void *_p, unsigned long long *data, int *data_len) {
  struct elim_params *p = (struct elim_params *)_p;
  *data_len = 2;

  FiniteField *F = (FiniteField*)(new GFqDom<int>(2,p->e));
  gf2e *ff = gf2e_init_givgfq(F);

  mzd_slice_t *A = mzd_slice_init(ff,p->m,p->m);
  mzd_slice_randomize(A);

  const int bitmask = (1<<ff->degree)-1;
  for(rci_t i=0; i<p->m; i++) {
    while(mzd_slice_read_elem(A, i, i) == 0) {
      mzd_slice_write_elem(A, i, i, random()&bitmask) ;
    }
  };

  mzd_slice_t *B = mzd_slice_init(ff,p->m,p->n);
  mzd_slice_randomize(B);

  data[0] = walltime(0);
  data[1] = cpucycles();
  
  if (strcmp(p->direction,"lower_left")==0) {
    if(strcmp(p->algorithm,"naive")==0)
      mzd_slice_trsm_lower_left_naive(A, B);
    else 
      _mzd_slice_trsm_lower_left(A, B, p->cutoff);
  }   
  data[1] = cpucycles() - data[1];
  data[0] = walltime(data[0]);
  
  mzd_slice_free(A);
  mzd_slice_free(B);
  gf2e_free(ff);
  delete F;
  return 0;
}


int main(int argc, char **argv) {
  global_options(&argc, &argv);

  if (argc < 6)
    m4ri_die("syntax: e m n matrix_type direction algorithm\n");

  struct elim_params params;

  params.e = atoi(argv[1]);
  params.m = atoi(argv[2]);
  params.n = atoi(argv[3]);
  params.matrix_type = argv[4];
  params.direction = argv[5];
  if (strcmp(params.direction,"lower_left") != 0)
    m4ri_die("not implemented.");

  if (argc >= 7)
    params.algorithm = argv[6];
  else
    params.algorithm = (char*)"default";
  if (argc >= 8)
    params.cutoff = atoi(argv[7]);
  else
    params.cutoff = MZED_TRSM_CUTOFF;

  srandom(17);
  unsigned long long data[2];

  if (strcmp(params.matrix_type,"mzed_t") == 0)
    run_bench(run_mzed, (void*)&params, data, 2);
  else if(strcmp(params.matrix_type,"mzd_slice_t") == 0)
    run_bench(run_mzd_slice, (void*)&params, data, 2);
  else
    m4ri_die("unknown type '%s'",params.matrix_type);
  double cc_per_op = ((double)data[1])/ ( (double)params.m * (double)params.m * powl((double)params.n,0.807) );

  printf("e: %2d, m: %5d, n: %5d, cutoff: %4d, cpu cycles: %10llu, cc/(mmn^0.807): %.5lf, wall time: %lf\n", params.e, params.m, params.n, params.cutoff, data[1], cc_per_op, data[0] / 1000000.0);
}



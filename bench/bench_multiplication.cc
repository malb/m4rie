#include <gf2e_cxx/finite_field_givaro.h>
#include <m4rie/m4rie.h>

#include "cpucycles.h"
#include "walltime.h"

using namespace M4RIE;

int main(int argc, char **argv) {
  double wt = 0.0;
  char *algorithm;
  size_t r; 

  if (argc < 4)
    m4ri_die("syntax: e m n alg");
  int k = atoi(argv[1]);
  int m = atoi(argv[2]);
  int n = atoi(argv[3]);
  if (argc == 5)
    algorithm = argv[4];
  else
    algorithm = (char*)"default";
  FiniteField *F = (FiniteField*)(new GFqDom<int>(2,k));
  gf2e *ff = gf2e_init_givgfq(F);

  mzed_t *A = mzed_init(ff,m,n);
  mzed_randomize(A);
  mzed_t *B = mzed_init(ff,n,m);
  mzed_randomize(B);

  mzed_t *C;

  wt = walltime(&wt);
  unsigned long long t = cpucycles();
  if(strcmp(algorithm, "travolta")==0)
    C = mzed_mul_travolta(NULL, A, B);
  else if(strcmp(algorithm,"naive")==0)
    C = mzed_mul_naive(NULL, A, B);
  else
    C = mzed_mul(NULL, A, B);
  printf("m: %5d, n: %5d, cpu cycles: %llu wall time: %lf\n",m, n, r, cpucycles() - t, walltime(&wt));

  mzed_free(A);

  gf2e_free(ff);
  delete F;
}

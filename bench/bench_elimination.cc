#include <gf2e_cxx/finite_field_givaro.h>
#include <m4rie/m4rie.h>

#include "cpucycles.h"
#include "walltime.h"

using namespace M4RIE;

int main(int argc, char **argv) {
  double wt = 0.0;
  size_t r; 

  if (argc < 4)
    m4ri_die("syntax: e m n alg");
  int k = atoi(argv[1]);
  int m = atoi(argv[2]);
  int n = atoi(argv[3]);
  char *algorithm = argv[4];
  FiniteField *F = (FiniteField*)(new GFqDom<int>(2,k));
  gf2e *ff = gf2e_init_givgfq(F);

  mzed_t *A = mzed_init(ff,m,n);
  mzed_randomize(A);

  mzed_echelonize_travolta(A,1);

  wt = walltime(&wt);
  unsigned long long t = cpucycles();
  if(strcmp(algorithm,"travolta")==0)
   r=  mzed_echelonize_travolta(A, 1);
  else if(strcmp(algorithm,"naive")==0)
    r = mzed_echelonize_naive(A, 1);
  printf("m: %5d, n: %5d, r: %5d, cpu cycles: %llu wall time: %lf\n",m, n, r, cpucycles() - t, walltime(&wt));

  mzed_free(A);

  gf2e_free(ff);
  delete F;
}

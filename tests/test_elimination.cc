#include <gf2e_cxx/finite_field_givaro.h>
#include <m4rie/m4rie.h>

using namespace M4RIE;

int main(int argc, char **argv) {
  int fail = 0;
  for(size_t k=2; k<=10; k++) {
    FiniteField *F = (FiniteField*)(new GFqDom<int>(2,k));
    gf2e *ff = gf2e_init_givgfq(F);
    for(size_t i=0; i<(k*1024/(1<<k)); i++) {
      size_t m = random() & 255;
      size_t n = random() & 255;
      m = m ? (m) : 1;
      n = n ? (n) : 1;
      mzed_t *A0 = mzed_init(ff,m,n);
      mzed_randomize(A0);
      mzed_t *A1 = mzed_copy(NULL, A0);
      mzed_t *A2 = mzed_copy(NULL, A0);
      mzed_echelonize_travolta(A0,1);
      mzed_echelonize_naive(A1,1);
      mzed_echelonize(A2,1);
      printf("elim: m: %5d, n: %5d, k: %2d ... ",(int)m,(int)n,(int)k);
      if (!mzed_cmp(A0,A1) && !mzed_cmp(A1,A2)) {
        printf("pass\n");
      } else {
        printf("FAIL\n");
        fail = 1;
      }
      mzed_free(A0);
      mzed_free(A1);
      mzed_free(A2);
    }
    gf2e_free(ff);
    delete F;
  }
  return fail;
}

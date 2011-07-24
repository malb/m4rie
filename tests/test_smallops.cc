#include "testing.h"

#include <gf2e_cxx/finite_field_givaro.h>
#include <m4rie/m4rie.h>

using namespace M4RIE;

mzed_t *random_mzed_t(gf2e *ff, int m, int n) {
  mzed_t *A  = mzed_init(ff,m,n);
  mzed_randomize(A);
  return A;
}

mzd_slice_t *random_mzd_slice_t(gf2e *ff, int m, int n) {
  mzd_slice_t *A  = mzd_slice_init(ff,m,n);
  mzd_slice_randomize(A);
  return A;
}

int test_slice(gf2e *ff, int m, int n) {
  int fail_ret = 0;

  mzed_t *A = random_mzed_t(ff, m, n);
  mzd_slice_t *a = mzed_slice(NULL, A);

  mzd_slice_t *b = random_mzd_slice_t(ff, m, n);
  mzed_slice(b, A);

  m4rie_check(mzd_slice_cmp(a, b) == 0);

  mzed_t *B = mzed_cling(NULL, b);
  
  m4rie_check(mzed_cmp(A, B) == 0);

  mzed_cling(A, a);

  m4rie_check(mzed_cmp(A, B) == 0);

  mzed_free(A);
  mzed_free(B);
  mzd_slice_free(a);
  mzd_slice_free(b);
  
  return fail_ret; 
};

int test_add(gf2e *ff, int m, int n) {
  int fail_ret = 0;

  mzed_t *A = random_mzed_t(ff, m, n);
  mzed_t *B = random_mzed_t(ff, m, n);
  mzed_t *C = random_mzed_t(ff, m, n);

  mzed_add(C,A,B);

  mzed_t *D = mzed_copy(NULL, C);

  mzed_add(C,C,A);
  mzed_add(C,C,B);

  m4rie_check(mzed_is_zero(C) == 1);

  mzed_add(C,A,B);

  m4rie_check(mzed_cmp(D,C) == 0);

  mzd_slice_t *a = mzed_slice(NULL, A);
  mzd_slice_t *b = mzed_slice(NULL, B);
  mzd_slice_t *c = mzed_slice(NULL, C);
  mzd_slice_t *d = mzed_slice(NULL, D);

  mzd_slice_add(c, a, b);

  d = mzd_slice_copy(NULL, c);

  mzed_cling(D, d);

  m4rie_check( mzed_cmp(D, C) == 0 );


  mzed_free(A);
  mzed_free(B);
  mzed_free(C);
  mzed_free(D);

  mzd_slice_free(a);
  mzd_slice_free(b);
  mzd_slice_free(c);
  mzd_slice_free(d);

  return fail_ret; 
}

int test_batch(gf2e *ff, int m, int n) {
  int fail_ret = 0;
  printf("testing k: %2d, m: %4d, n: %4d ",ff->degree,m,n);

  m4rie_check( test_slice(ff, m, n) == 0); printf("."); 
  m4rie_check( test_add(ff, m, n) == 0) ; printf("."); 
  if (fail_ret == 0)
    printf(" passed\n");
  else
    printf(" FAILED\n");
  return fail_ret;
}

int main(int argc, char **argv) {

  gf2e *ff[10];

  for(int k=2; k<=10; k++) {
    FiniteField *F = (FiniteField*)(new GFqDom<int>(2,k));
    ff[k] = gf2e_init_givgfq(F);
    delete F;
  }

  int fail = 0;
  for(int k=2; k<=2; k++) {
    for(int m=1; m<=128; m++) {
      for(int n=1; n<=128; n++) {
        fail += test_batch(ff[k], m, n);
      }
    }
  }
  return fail;
}
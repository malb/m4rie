#include "testing.h"

#include <gf2e_cxx/finite_field_givaro.h>
#include <m4rie/m4rie.h>

using namespace M4RIE;

int test_slice(gf2e *ff, int m, int n) {
  int fail_ret = 0;

  mzed_t *A = random_mzed_t(ff, m, n);
  mzd_slice_t *a = mzed_slice(NULL, A);
  mzd_slice_set_canary(a);

  mzd_slice_t *b = random_mzd_slice_t(ff, m, n);
  mzed_slice(b, A);

  m4rie_check(mzd_slice_canary_is_alive(a));
  m4rie_check(mzd_slice_canary_is_alive(b));
  m4rie_check(mzed_canary_is_alive(A));

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

int test_slice_known_answers(gf2e *ff, int m, int n) {
  int fail_ret = 0;
  mzed_t *A = mzed_init(ff, m, n);
  mzed_set_canary(A);

  mzd_t *one = mzd_init(m,n);
  mzd_set_ui(one, 1);

  for(int j=0; j<ff->degree; j++) {
    mzed_set_ui(A, 1<<j);
    mzd_slice_t *a = mzed_slice(NULL, A);

    for(int i=0; i<a->depth; i++) {
      if (i!=j) {
        m4rie_check( mzd_is_zero(a->x[i]) );
      } else {
        m4rie_check( mzd_cmp(a->x[i], one) == 0 );
      }
    }
    mzed_t *AA = mzed_cling(NULL, a);
    m4rie_check( mzed_cmp(AA, A) == 0 );
    m4rie_check( mzed_canary_is_alive(A) );
    mzd_slice_free(a);
    mzed_free(AA);
  }
  mzd_free(one);
  mzed_free(A);
  return fail_ret;
}

int test_add(gf2e *ff, int m, int n) {
  int fail_ret = 0;

  mzed_t *A = random_mzed_t(ff, m, n);
  mzed_t *B = random_mzed_t(ff, m, n);
  mzed_t *C = random_mzed_t(ff, m, n);

  mzed_add(C,A,B);

  mzed_t *D = mzed_copy(NULL, C);
  mzed_set_canary(D);

  mzed_add(C,C,A);
  mzed_add(C,C,B);

  m4rie_check(mzed_is_zero(C) == 1);

  mzed_add(C,A,B);

  m4rie_check(mzed_cmp(D,C) == 0);

  mzd_slice_t *a = mzed_slice(NULL, A);
  mzd_slice_t *b = mzed_slice(NULL, B);
  mzd_slice_t *c = mzed_slice(NULL, C);
  mzd_slice_t *d = mzed_slice(NULL, D);

  mzd_slice_set_canary(a);
  mzd_slice_set_canary(b);
  mzd_slice_set_canary(c);

  mzd_slice_add(c, a, b);

  d = mzd_slice_copy(NULL, c);
  mzd_slice_set_canary(d);

  mzed_cling(D, d);

  m4rie_check( mzed_cmp(D, C) == 0 );

  m4rie_check( mzed_canary_is_alive(A));
  m4rie_check( mzed_canary_is_alive(B));
  m4rie_check( mzed_canary_is_alive(C));
  m4rie_check( mzed_canary_is_alive(D));

  m4rie_check( mzd_slice_canary_is_alive(a));
  m4rie_check( mzd_slice_canary_is_alive(b));
  m4rie_check( mzd_slice_canary_is_alive(c));
  m4rie_check( mzd_slice_canary_is_alive(d));

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

  m4rie_check( test_slice(ff, m, n) == 0);   printf("."); fflush(0);
  m4rie_check( test_add(ff, m, n) == 0) ;    printf("."); fflush(0);
  m4rie_check( test_slice_known_answers(ff, m, n) == 0); printf("."); fflush(0);

  m4rie_check( test_slice(ff, m, m) == 0);   printf("."); fflush(0);
  m4rie_check( test_add(ff, m, m) == 0) ;    printf("."); fflush(0);
  m4rie_check( test_slice_known_answers(ff, m, m) == 0); printf("."); fflush(0);

  m4rie_check( test_slice(ff, n, m) == 0);   printf("."); fflush(0);
  m4rie_check( test_add(ff, n, m) == 0) ;    printf("."); fflush(0);
  m4rie_check( test_slice_known_answers(ff, n, m) == 0); printf("."); fflush(0);

  m4rie_check( test_slice(ff, n, n) == 0);  printf("."); fflush(0);
  m4rie_check( test_add(ff, n, n) == 0) ;   printf("."); fflush(0);
  m4rie_check( test_slice_known_answers(ff, n, n) == 0); printf("."); fflush(0);

  if (fail_ret == 0)
    printf(" passed\n");
  else
    printf(" FAILED\n");
  return fail_ret;
}

int main(int argc, char **argv) {

  gf2e *ff[10];
  int fail_ret = 0;

  for(int k=2; k<=10; k++) {
    FiniteField *F = (FiniteField*)(new GFqDom<int>(2,k));
    ff[k] = gf2e_init_givgfq(F);
    delete F;
  }

  for(int k=2; k<=__M4RIE_MAX_KARATSUBA_DEGREE; k++) {
    fail_ret += test_batch(ff[k],   2, m4ri_radix/gf2e_degree_to_w(ff[k]));
    fail_ret += test_batch(ff[k],   2, 2*m4ri_radix/gf2e_degree_to_w(ff[k]));
    fail_ret += test_batch(ff[k],   2, 3*m4ri_radix/gf2e_degree_to_w(ff[k]));
    fail_ret += test_batch(ff[k],   2, 4*m4ri_radix/gf2e_degree_to_w(ff[k]));
    fail_ret += test_batch(ff[k],   4,   3);
    fail_ret += test_batch(ff[k],   1,   2);
    fail_ret += test_batch(ff[k],  10,  11);
    fail_ret += test_batch(ff[k],  20,  19);
    fail_ret += test_batch(ff[k],  32,  64);
    fail_ret += test_batch(ff[k],  63,  65);
    fail_ret += test_batch(ff[k],  64,  65);
    fail_ret += test_batch(ff[k],  64, 128);
    fail_ret += test_batch(ff[k],  65, 129);
    fail_ret += test_batch(ff[k], 201, 200);
    fail_ret += test_batch(ff[k], 217,   2);
  }

  for(int k=2; k<=10; k++) {
    gf2e_free(ff[k]);
  }

  return fail_ret;
}

#include "testing.h"
#include <gf2e_cxx/finite_field_givaro.h>

using namespace M4RIE;

int test_equality(gf2e *ff, rci_t m, rci_t n) {
  int fail_ret = 0;
  mzed_t *A0 = random_mzed_t(ff, m, n);
  mzed_t *A1 = mzed_copy(NULL, A0);
  mzed_t *A2 = mzed_copy(NULL, A0);
  mzed_t *A3 = mzed_copy(NULL, A0);

  mzed_set_canary(A1);
  mzed_set_canary(A2);
  mzed_set_canary(A3);

  const rci_t r0 = mzed_echelonize_travolta(A0,1);
  const rci_t r1 = mzed_echelonize_naive(A1,1);
  const rci_t r2 = mzed_echelonize(A2,1);
  rci_t r3 = 0;
  if (ff->degree <= __M4RIE_MAX_KARATSUBA_DEGREE) {
    r3 = mzed_echelonize_ple(A3,1);
  }
  m4rie_check( r0 == r1);
  m4rie_check( mzed_cmp(A0, A1) == 0);

  m4rie_check( r1 == r2);
  m4rie_check( mzed_cmp(A1, A2) == 0);

  if(ff->degree <=  __M4RIE_MAX_KARATSUBA_DEGREE) {
    m4rie_check( r2 == r3);
    m4rie_check( mzed_cmp(A2, A3) == 0);
    m4rie_check( r3 == r0);
    m4rie_check( mzed_cmp(A3, A0) == 0);
  } else {
    m4rie_check( r2 == r0);
    m4rie_check( mzed_cmp(A2, A0) == 0);
  }

  m4rie_check( mzed_canary_is_alive(A0) );
  m4rie_check( mzed_canary_is_alive(A1) );
  m4rie_check( mzed_canary_is_alive(A2) );
  m4rie_check( mzed_canary_is_alive(A3) );

  mzed_free(A0);
  mzed_free(A1);
  mzed_free(A2);
  mzed_free(A3);

  return fail_ret;
}

int test_batch(gf2e *ff, rci_t m, rci_t n) {
  int fail_ret = 0;
  printf("elim: k: %2d, minpoly: 0x%03x m: %5d, n: %5d ",(int)ff->degree, (unsigned int)ff->minpoly, (int)m, (int)n);

  if(m == n) {
    m4rie_check(   test_equality(ff, m, n) == 0); printf("."); fflush(0);
    printf(" ");
  } else {
    m4rie_check(   test_equality(ff, m, n) == 0); printf("."); fflush(0);
    m4rie_check(   test_equality(ff, n, m) == 0); printf("."); fflush(0);
  }

  if (fail_ret == 0)
    printf(" passed\n");
  else
    printf(" FAILED\n");

  return fail_ret;
}

int main(int argc, char **argv) {
  srandom(17);

  gf2e *ff[10];
  int fail_ret = 0;

  for(int k=2; k<=10; k++) {
    FiniteField *F = (FiniteField*)(new GFqDom<int>(2,k));
    ff[k] = gf2e_init_givgfq(F);
    delete F;
  }

  for(int k=2; k<=10; k++) {
    fail_ret += test_batch(ff[k],   2,   5);
    fail_ret += test_batch(ff[k],   5,  10);
    fail_ret += test_batch(ff[k],   1,   1);
    fail_ret += test_batch(ff[k],   1,   2);
    fail_ret += test_batch(ff[k],  11,  12);
    fail_ret += test_batch(ff[k],  21,  22);
    fail_ret += test_batch(ff[k],  13,   2);
    fail_ret += test_batch(ff[k],  32,  33);
    fail_ret += test_batch(ff[k],  63,  64);
    fail_ret += test_batch(ff[k], 127, 128);
    fail_ret += test_batch(ff[k], 200,  20);
    fail_ret += test_batch(ff[k],   1,   1);
    fail_ret += test_batch(ff[k],   1,   3);
    fail_ret += test_batch(ff[k],  11,  13);
    fail_ret += test_batch(ff[k],  21,  23);
    fail_ret += test_batch(ff[k],  13,  90);
    fail_ret += test_batch(ff[k],  32,  34);
    fail_ret += test_batch(ff[k],  63,  65);
    fail_ret += test_batch(ff[k], 127, 129);
    fail_ret += test_batch(ff[k], 200, 112);
  };

  for(int k=2; k<=10; k++) {
    gf2e_free(ff[k]);
  }

  return fail_ret;
}

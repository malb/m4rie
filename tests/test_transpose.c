#include "testing.h"

int test_equality(gf2e *ff, rci_t m, rci_t n) {
  int fail_ret = 0, bad_elem;
  mzed_t *A0 = random_mzed_t(ff, m, n);
  mzed_t *A1 = mzed_init(ff, n, m);
  
  mzed_t *A2 = mzed_copy(NULL, A0);
  mzed_t *A3 = NULL;
  
  mzed_t *A4 = mzed_copy(NULL, A0);
  mzed_t *A5 = mzed_copy(NULL, A1);
  
  mzed_t *A6 = mzed_copy(NULL, A0);
  mzed_t *A7 = random_mzed_t(ff, n, m);
  
  mzed_set_canary(A1);
  mzed_set_canary(A2);
  mzed_set_canary(A4);
  mzed_set_canary(A5);
  mzed_set_canary(A6);
  mzed_set_canary(A7);
  
  for (rci_t c = 0; c < A0->ncols; c++) {
    for (rci_t r = 0; r < A0->nrows; r++) {
      mzed_write_elem(A1, c, r, mzed_read_elem(A0, r, c));
    }
  }
  A3 = mzed_transpose(A3, A2);
  A5 = mzed_transpose(A5, A4);
  A7 = mzed_transpose(A7, A6);
  
  m4rie_check( mzed_cmp(A0, A2) == 0);
  m4rie_check( mzed_cmp(A0, A4) == 0);
  m4rie_check( mzed_cmp(A0, A6) == 0);
  m4rie_check( mzed_cmp(A1, A3) == 0);
  m4rie_check( mzed_cmp(A1, A5) == 0);
  m4rie_check( mzed_cmp(A1, A7) == 0);
  
  m4rie_check( mzed_canary_is_alive(A0) );
  m4rie_check( mzed_canary_is_alive(A1) );
  m4rie_check( mzed_canary_is_alive(A2) );
  m4rie_check( mzed_canary_is_alive(A4) );
  m4rie_check( mzed_canary_is_alive(A5) );
  m4rie_check( mzed_canary_is_alive(A6) );
  m4rie_check( mzed_canary_is_alive(A7) );
  
  mzed_free(A0);
  mzed_free(A1);
  mzed_free(A2);
  mzed_free(A3);
  mzed_free(A4);
  mzed_free(A5);
  mzed_free(A6);
  mzed_free(A7);
  
  return fail_ret;
}

int test_batch(gf2e *ff, rci_t m, rci_t n) {
  int fail_ret = 0;
  printf("elim: k: %2d, minpoly: 0x%05x m: %5d, n: %5d ",(int)ff->degree, (unsigned int)ff->minpoly, (int)m, (int)n);

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

  int runlong = parse_parameters(argc, argv);

  gf2e *ff;
  int fail_ret = 0;

  for(int k=2; k<=16; k++) {
    ff = gf2e_init(irreducible_polynomials[k][1]);

    fail_ret += test_batch(ff,   2,   5);
    fail_ret += test_batch(ff,   5,  10);
    fail_ret += test_batch(ff,   1,   1);
    fail_ret += test_batch(ff,   1,   2);
    fail_ret += test_batch(ff,  11,  12);
    fail_ret += test_batch(ff,  21,  22);
    fail_ret += test_batch(ff,  13,   2);
    fail_ret += test_batch(ff,  32,  33);
    fail_ret += test_batch(ff,  63,  64);
    if (k <= 12 || runlong) {
      fail_ret += test_batch(ff, 127, 128);
      fail_ret += test_batch(ff, 200,  20);
    }
    fail_ret += test_batch(ff,   1,   3);
    fail_ret += test_batch(ff,  11,  13);
    fail_ret += test_batch(ff,  21,  23);
    fail_ret += test_batch(ff,  13,  90);
    fail_ret += test_batch(ff,  32,  34);
    fail_ret += test_batch(ff,  63,  65);
    if (k <= 12 || runlong) {
      fail_ret += test_batch(ff, 127, 129);
      fail_ret += test_batch(ff, 200, 112);
      fail_ret += test_batch(ff,  10, 200);
    }
    gf2e_free(ff);
  }

  return fail_ret;
}

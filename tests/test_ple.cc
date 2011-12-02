#include "testing.h"
#include <gf2e_cxx/finite_field_givaro.h>

using namespace M4RIE;

int test_mzd_slice_ple(gf2e *ff, const rci_t m, const rci_t n, const rci_t r) {
  int fail_ret = 0;

  mzed_t *a = random_mzed_t_rank(ff, m, n, r);
  mzd_slice_t *A = mzed_slice(NULL, a);
  mzed_free(a);

  mzd_slice_t *LE = mzd_slice_copy(NULL, A);
  mzd_slice_t *L = mzd_slice_init(ff, m, m);
  mzd_slice_t *E = mzd_slice_init(ff, m, n);

  mzp_t *P = mzp_init(m);
  mzp_t *Q = mzp_init(n);

  rci_t rbar = mzd_slice_ple(LE, P, Q);

  m4rie_check( rbar == r);

  for(rci_t j=0; j<r; j++) {
    for(rci_t i=j; i<LE->nrows; i++) {
      mzd_slice_write_elem(L,i,j, mzd_slice_read_elem(LE,i,j));
    }
  }

  for(rci_t i=0; i<r; i++) {
    mzd_slice_write_elem(E, i, Q->values[i], 1);
    for(rci_t j=Q->values[i]+1; j< LE->ncols; j++) {
      mzd_slice_write_elem(E, i, j, mzd_slice_read_elem(LE, i, j));
    }
  }

  mzd_slice_t *B = mzd_slice_mul(NULL, L, E);

  mzd_slice_apply_p_left(A, P);

  m4rie_check( mzd_slice_cmp(A, B) == 0);

  mzd_slice_free(A);
  mzd_slice_free(B);
  mzd_slice_free(LE);
  mzd_slice_free(L);
  mzd_slice_free(E);
  mzp_free(P);
  mzp_free(Q);

  return fail_ret;
}

int test_mzed_ple(gf2e *ff, const rci_t m, const rci_t n, const rci_t r) {
  int fail_ret = 0;

  mzed_t *A  = random_mzed_t_rank(ff, m, n, r);

  /**
   * We check equality first
   */

  mzed_t *LE0 = mzed_copy(NULL, A);
  mzp_t *P0 = mzp_init(m);
  mzp_t *Q0 = mzp_init(n);

  mzed_t *LE1 = mzed_copy(NULL, A);
  mzp_t *P1 = mzp_init(m);
  mzp_t *Q1 = mzp_init(n);

  mzed_t *LE2 = mzed_copy(NULL, A);
  mzp_t *P2 = mzp_init(m);
  mzp_t *Q2 = mzp_init(n);

  mzed_set_canary(LE0);
  rci_t r0 = mzed_ple_naive(   LE0, P0, Q0);
  m4rie_check( mzed_canary_is_alive(LE0) );
  m4rie_check( r0 == r);
  
  mzed_set_canary(LE1);
  rci_t r1 = mzed_ple_travolta(LE1, P1, Q1);
  m4rie_check( mzed_canary_is_alive(LE1) );
  m4rie_check( r1 == r);

  mzed_set_canary(LE2);
  rci_t r2 = mzed_ple(         LE2, P2, Q2);
  m4rie_check( mzed_canary_is_alive(LE2) );
  m4rie_check( r2 == r);

  m4rie_check( mzed_cmp(LE0, LE1) == 0 );
  m4rie_check( mzed_cmp(LE1, LE2) == 0 );
  m4rie_check( mzed_cmp(LE2, LE0) == 0 );

  /**
   * Now we check mathematical properties. Equality has been
   * established so we only deal with LE0.
   */

  mzed_t *L = mzed_init(ff, m, m);

  mzed_set_canary(L);
  for(rci_t j=0; j<r; j++) 
    for(rci_t i=j; i<LE0->nrows; i++) 
      mzed_write_elem(L,i,j, mzed_read_elem(LE0,i,j));
  m4rie_check( mzed_canary_is_alive(L) );

  mzed_t *E = mzed_init(ff, m, n);

  mzed_set_canary(E);
  for(rci_t i=0; i<r; i++) {
    mzed_write_elem(E, i, Q0->values[i], 1);
    for(rci_t j=Q0->values[i]+1; j< LE0->ncols; j++)
      mzed_write_elem(E, i, j, mzed_read_elem(LE0, i, j));
  }
  m4rie_check( mzed_canary_is_alive(E) );

  mzed_t *B = mzed_mul(NULL, L, E);
  
  mzed_apply_p_left(A, P0);
  m4rie_check( mzed_canary_is_alive(A) );

  m4rie_check( mzed_cmp(A, B) == 0);

  mzed_free(A);
  mzed_free(B);

  mzed_free(LE0);
  mzp_free(P0);
  mzp_free(Q0);

  mzed_free(LE1);
  mzp_free(P1);
  mzp_free(Q1);

  mzed_free(LE2);
  mzp_free(P2);
  mzp_free(Q2);

  mzed_free(L);
  mzed_free(E);

  return fail_ret;
}


int test_batch(gf2e *ff, const rci_t m, const rci_t n, const rci_t r) {
  assert(r <= m);
  assert(r <= n);

  printf("ple: k: %2d, minpoly: 0x%03x m: %5d, n: %5d, r: %5d ",(int)ff->degree, (unsigned int)ff->minpoly, (int)m, (int)n, (int)r);

  int fail_ret = 0;

  if(m == n) {
    m4rie_check(   test_mzed_ple(ff, m, n, r) == 0); printf("."); fflush(0);
    printf(" ");
    if(ff->degree <= 4) {
      m4rie_check(   test_mzd_slice_ple(ff, m, n, r) == 0); printf("."); fflush(0);
      printf(" ");
    } else {
      printf("  ");
    }
    fflush(0);
  } else {
    m4rie_check(   test_mzed_ple(ff, m, n, r) == 0); printf("."); fflush(0);
    m4rie_check(   test_mzed_ple(ff, n, m, r) == 0); printf("."); fflush(0);
    if(ff->degree <= 4) {
      m4rie_check(   test_mzd_slice_ple(ff, m, n, r) == 0); printf("."); fflush(0);
      m4rie_check(   test_mzd_slice_ple(ff, n, m, r) == 0); printf("."); fflush(0);
    } else {
      printf("  ");
    }
    fflush(0);
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
    GFqDom<int> GF = GFqDom<int>(2,k);
    FiniteField *F = (FiniteField*)&GF;
    ff[k] = gf2e_init_givgfq(F);
  }

  for(int k=2; k<=10; k++) {
    fail_ret += test_batch(ff[k],   1,   1,   1);
    fail_ret += test_batch(ff[k],   1,   2,   1);
    fail_ret += test_batch(ff[k],   2,   2,   2);
    fail_ret += test_batch(ff[k],   2,   3,   2);
    fail_ret += test_batch(ff[k],  11,  12,  10);
    fail_ret += test_batch(ff[k],  21,  22,  21);
    fail_ret += test_batch(ff[k],  13,   2,   2);
    fail_ret += test_batch(ff[k],  32,  33,  31);
    fail_ret += test_batch(ff[k],  63,  64,  62);
    fail_ret += test_batch(ff[k], 127, 128, 125);
    fail_ret += test_batch(ff[k], 127, 128,  12);
    fail_ret += test_batch(ff[k], 127, 128,  37);
    fail_ret += test_batch(ff[k], 127, 128,  67);
    fail_ret += test_batch(ff[k], 200,  20,  19);
    fail_ret += test_batch(ff[k],   1,   1,   0);
    fail_ret += test_batch(ff[k],   1,   3,   1);
    fail_ret += test_batch(ff[k],  11,  13,  10);
    fail_ret += test_batch(ff[k],  21,  23,  20);
    fail_ret += test_batch(ff[k],  13,  90,  10);
    fail_ret += test_batch(ff[k],  32,  34,  31);
    fail_ret += test_batch(ff[k],  63,  65,  62);
    fail_ret += test_batch(ff[k], 127, 129, 127);
    fail_ret += test_batch(ff[k], 200, 112, 111);
  };

  for(int k=2; k<=10; k++) {
    gf2e_free(ff[k]);
  }
  if (fail_ret == 0)
    printf("success\n");

  return fail_ret;
}

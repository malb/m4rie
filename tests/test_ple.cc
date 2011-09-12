#include "testing.h"
#include <gf2e_cxx/finite_field_givaro.h>

using namespace M4RIE;


int test_ple(gf2e *ff, const rci_t m, const rci_t n, const rci_t r) {
  int fail_ret = 0;

  mzed_t *A  = random_mzed_t_rank(ff, m, n, r);
  mzed_t *LE = mzed_copy(NULL, A);
  mzed_t *L = mzed_init(ff, m, m);
  mzed_t *E = mzed_init(ff, m, n);

  mzp_t *P = mzp_init(m);
  mzp_t *Q = mzp_init(n);

  mzed_set_canary(LE);

  mzed_set_canary(L);
  mzed_set_canary(E);

  rci_t rbar = mzed_ple_naive(LE,P, Q);

  m4rie_check( rbar == r);
  m4rie_check( mzed_canary_is_alive(LE) );

  for(rci_t j=0; j<r; j++) {
    for(rci_t i=j; i<LE->nrows; i++) {
      mzed_write_elem(L,i,j, mzed_read_elem(LE,i,j));
    }
  }
  m4rie_check( mzed_canary_is_alive(L) );

  for(rci_t i=0; i<r; i++) {
    mzed_write_elem(E, i, Q->values[i], 1);
    for(rci_t j=Q->values[i]+1; j< LE->ncols; j++) {
      mzed_write_elem(E, i, j, mzed_read_elem(LE, i, j));
    }
  }
  m4rie_check( mzed_canary_is_alive(E) );

  mzed_t *B = mzed_mul(NULL, L, E);

  mzed_apply_p_left(A, P);

  m4rie_check( mzed_canary_is_alive(A) );
  m4rie_check( mzed_cmp(A, B) == 0);

  mzed_free(A);
  mzed_free(B);
  mzed_free(LE);
  mzed_free(L);
  mzed_free(E);
  mzp_free(P);
  mzp_free(Q);

  return fail_ret;
}

int test_batch(gf2e *ff, const rci_t m, const rci_t n, const rci_t r) {
  assert(r <= m);
  assert(r <= n);

  printf("ple: k: %2d, minpoly: 0x%03x m: %5d, n: %5d, r: %5d ",(int)ff->degree, (unsigned int)ff->minpoly, (int)m, (int)n, (int)r);

  int fail_ret = 0;

  if(m == n) {
    m4rie_check(   test_ple(ff, m, n, r) == 0); printf("."); fflush(0);
    printf(" ");
  } else {
    m4rie_check(   test_ple(ff, m, n, r) == 0); printf("."); fflush(0);
    m4rie_check(   test_ple(ff, n, m, r) == 0); printf("."); fflush(0);
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

  return fail_ret;
}

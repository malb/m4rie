/**
 * \file test_multiplication.cc
 * \brief Test code for multiplication routines
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

/******************************************************************************
*
*            M4RIE: Linear Algebra over GF(2^e)
*
*    Copyright (C) 2010 Martin Albrecht <martinralbrecht@googlemail.com>
*
*  Distributed under the terms of the GNU General Public License (GEL)
*  version 2 or higher.
*
*    This code is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*    General Public License for more details.
*
*  The full text of the GPL is available at:
*
*                  http://www.gnu.org/licenses/
******************************************************************************/


#include "testing.h"
#include <gf2e_cxx/finite_field_givaro.h>

using namespace M4RIE;

int test_addmul(gf2e *ff, rci_t m, rci_t n, rci_t l) {
  int fail_ret = 0;

  mzed_t *A = random_mzed_t(ff, m, l);
  mzed_t *B = random_mzed_t(ff, l, n);

  mzed_t *C0 = random_mzed_t(ff, m, n);
  mzed_t *C1 = mzed_copy(NULL, C0);
  mzed_t *C2 = mzed_copy(NULL, C0);
  mzed_t *C3 = mzed_copy(NULL, C0);
  mzed_t *C4 = mzed_copy(NULL, C0);

  mzed_set_canary(C1);
  mzed_set_canary(C2);
  mzed_set_canary(C3);
  mzed_set_canary(C4);

  mzed_addmul_travolta(C0, A, B);
  mzed_addmul_naive(C1, A, B);
  mzed_addmul_strassen(C2, A, B, 64);
  mzed_addmul(C3, A, B);

  m4rie_check( mzed_cmp(C0, C1) == 0);
  m4rie_check( mzed_cmp(C1, C2) == 0);
  m4rie_check( mzed_cmp(C2, C3) == 0);

  if (ff->degree <= __M4RIE_MAX_KARATSUBA_DEGREE) {
    mzed_addmul_karatsuba(C4, A, B);
    m4rie_check( mzed_cmp(C3, C4) == 0);
  }

  m4rie_check( mzed_canary_is_alive(A) );
  m4rie_check( mzed_canary_is_alive(B) );
  m4rie_check( mzed_canary_is_alive(C1) );
  m4rie_check( mzed_canary_is_alive(C2) );
  m4rie_check( mzed_canary_is_alive(C3) );
  m4rie_check( mzed_canary_is_alive(C4) );

  mzed_free(A);
  mzed_free(B);
  mzed_free(C0);
  mzed_free(C1);
  mzed_free(C2);
  mzed_free(C3);
  mzed_free(C4);

  return fail_ret;
}

int test_mul(gf2e *ff, rci_t m, rci_t n, rci_t l) {
  int fail_ret = 0;
  const mzed_t *A = random_mzed_t(ff, m, l);
  const mzed_t *B = random_mzed_t(ff, l, n);

  mzed_t *C0 = random_mzed_t(ff, m, n);
  mzed_t *C1 = random_mzed_t(ff, m, n);
  mzed_t *C2 = random_mzed_t(ff, m, n);
  mzed_t *C3 = random_mzed_t(ff, m, n);
  mzed_t *C4 = random_mzed_t(ff, m, n);

  mzed_mul_travolta(C0, A, B);
  mzed_mul_naive(C1, A, B);
  mzed_mul_strassen(C2, A, B, 64);
  mzed_mul(C3, A, B);

  m4rie_check( mzed_cmp(C0, C1) == 0);
  m4rie_check( mzed_cmp(C1, C2) == 0);
  m4rie_check( mzed_cmp(C2, C3) == 0);

  if (ff->degree <= __M4RIE_MAX_KARATSUBA_DEGREE) {
    mzed_mul_karatsuba(C4, A, B);
    m4rie_check( mzed_cmp(C3, C4) == 0);
  }

  m4rie_check( mzed_canary_is_alive((mzed_t*)A) );
  m4rie_check( mzed_canary_is_alive((mzed_t*)B) );
  m4rie_check( mzed_canary_is_alive(C1) );
  m4rie_check( mzed_canary_is_alive(C2) );
  m4rie_check( mzed_canary_is_alive(C3) );

  mzed_free((mzed_t*)A);
  mzed_free((mzed_t*)B);
  mzed_free(C0);
  mzed_free(C1);
  mzed_free(C2);
  mzed_free(C3);
  mzed_free(C4);

  return fail_ret;
}

int test_scalar(gf2e *ff, rci_t m, rci_t n) {
  int fail_ret = 0;

  word a = random() & ((1<<ff->degree)-1);
  while (!a)
    a = random() & ((1<<ff->degree)-1);
  mzed_t *B = random_mzed_t(ff, m, n);

  mzed_t *C0 = mzed_init(ff, m, n);
  mzed_t *C1 = random_mzed_t(ff, m, n);
  mzed_t *C2 = NULL;

  C0 = mzed_mul_scalar(C0, a, B);
  C1 = mzed_mul_scalar(C1, a, B);
  C2 = mzed_mul_scalar(C2, a, B);

  m4rie_check( mzed_cmp(C0, C1) == 0);
  m4rie_check( mzed_cmp(C1, C2) == 0);

  if (ff->degree <= __M4RIE_MAX_KARATSUBA_DEGREE) {
    mzed_t *C3 = NULL;
    mzd_slice_t *BB = mzed_slice(NULL, B);
    mzd_slice_t *CC = mzd_slice_mul_scalar(NULL, a, BB);
    C3 = mzed_cling(C3, CC);
    mzd_slice_free(BB);
    mzd_slice_free(CC);
    m4rie_check( mzed_cmp(C2, C3) == 0);
    mzed_free(C3);
  }

  const word a_inv = ff->inv[a];

  mzed_t *B0 = mzed_init(ff, m, n);
  mzed_t *B1 = random_mzed_t(ff, m, n);
  mzed_t *B2 = NULL;

  B0 = mzed_mul_scalar(B0, a_inv, C0);
  B1 = mzed_mul_scalar(B1, a_inv, C1);
  B2 = mzed_mul_scalar(B2, a_inv, C2);

  m4rie_check( mzed_cmp(B, B0) == 0);
  m4rie_check( mzed_cmp(B, B1) == 0);
  m4rie_check( mzed_cmp(B, B2) == 0);

  mzed_free(C0);
  mzed_free(C1);
  mzed_free(C2);

  mzed_free(B0);
  mzed_free(B1);
  mzed_free(B2);

  return fail_ret;
}


int test_batch(gf2e *ff, rci_t m, rci_t l, rci_t n) {
  int fail_ret = 0;
  printf("mul: k: %2d, minpoly: 0x%03x m: %5d, l: %5d, n: %5d ",(int)ff->degree, (unsigned int)ff->minpoly, (int)m, (int)l, (int)n);

  m4rie_check(test_scalar(ff, m, m) == 0); printf("."); fflush(0);
  m4rie_check(test_scalar(ff, l, l) == 0); printf("."); fflush(0);
  m4rie_check(test_scalar(ff, n, n) == 0); printf("."); fflush(0);
  m4rie_check(test_scalar(ff, m, l) == 0); printf("."); fflush(0);
  m4rie_check(test_scalar(ff, l, n) == 0); printf("."); fflush(0);
  m4rie_check(test_scalar(ff, m, n) == 0); printf("."); fflush(0);
  m4rie_check(test_scalar(ff, l, m) == 0); printf("."); fflush(0);

  if(m == l && m == n) {
    m4rie_check(   test_mul(ff, m, l, n) == 0); printf("."); fflush(0);
    m4rie_check(test_addmul(ff, m, l, n) == 0); printf("."); fflush(0);
    printf("          ");
  } else {
    m4rie_check(   test_mul(ff, m, l, n) == 0); printf("."); fflush(0);
    m4rie_check(   test_mul(ff, m, n, l) == 0); printf("."); fflush(0);
    m4rie_check(   test_mul(ff, n, m, l) == 0); printf("."); fflush(0);
    m4rie_check(   test_mul(ff, n, l, m) == 0); printf("."); fflush(0);
    m4rie_check(   test_mul(ff, l, m, n) == 0); printf("."); fflush(0);
    m4rie_check(   test_mul(ff, l, n, m) == 0); printf("."); fflush(0);
    m4rie_check(test_addmul(ff, m, l, n) == 0); printf("."); fflush(0);
    m4rie_check(test_addmul(ff, m, n, l) == 0); printf("."); fflush(0);
    m4rie_check(test_addmul(ff, n, m, l) == 0); printf("."); fflush(0);
    m4rie_check(test_addmul(ff, n, l, m) == 0); printf("."); fflush(0);
    m4rie_check(test_addmul(ff, l, m, n) == 0); printf("."); fflush(0);
    m4rie_check(test_addmul(ff, l, n, m) == 0); printf("."); fflush(0);
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
    fail_ret += test_batch(ff[k],   1,   2,   3);
    fail_ret += test_batch(ff[k],  11,  12,  13);
    fail_ret += test_batch(ff[k],  21,  22,  23);
    fail_ret += test_batch(ff[k],  13,   2,  90);
    fail_ret += test_batch(ff[k],  32,  33,  34);
    fail_ret += test_batch(ff[k],  63,  64,  65);
    fail_ret += test_batch(ff[k], 127, 128, 129);
    fail_ret += test_batch(ff[k], 200,  20, 112);
  };

  for(int k=2; k<=10; k++) {
    gf2e_free(ff[k]);
  }

  return fail_ret;
}

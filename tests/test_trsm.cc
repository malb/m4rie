/**
 * \file test_multiplication.cc
 * \brief Test code for triangular system solving with matrices (TRSM) routines
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

/******************************************************************************
*
*            M4RIE: Linear Algebra over GF(2^e)
*
*    Copyright (C) 2011 Martin Albrecht <martinralbrecht@googlemail.com>
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


mzed_t *random_mzed_t_upper_left(gf2e *ff, rci_t m) {
  const int bitmask = (1<<ff->degree)-1;
  mzed_t *U = random_mzed_t(ff, m, m);
  for(rci_t i=0; i<m; i++) {
    for(rci_t j=0; j<i; j++) {
      mzed_write_elem(U, i, j, 0);
    }
    while(mzed_read_elem(U, i, i) == 0) {
      mzed_write_elem(U, i, i, random()&bitmask) ;
    }
  }
  mzed_set_canary(U);
  return U;
}

mzed_t *random_mzed_t_lower_left(gf2e *ff, rci_t m) {
  const int bitmask = (1<<ff->degree)-1;
  mzed_t *L = random_mzed_t(ff, m, m);
  for(rci_t i=0; i<m; i++) {
    for(rci_t j=i+1; j<m; j++) {
      mzed_write_elem(L, i, j, 0);
    }
    while(mzed_read_elem(L, i, i) == 0) {
      mzed_write_elem(L, i, i, random()&bitmask) ;
    }
  }
  mzed_set_canary(L);
  return L;
}

mzd_slice_t *random_mzd_slice_t_upper_left(gf2e *ff, rci_t m) {
  const int bitmask = (1<<ff->degree)-1;
  mzd_slice_t *U = random_mzd_slice_t(ff, m, m);
  for(rci_t i=0; i<m; i++) {
    for(rci_t j=0; j<i; j++) {
      mzd_slice_write_elem(U, i, j, 0);
    }
    while(mzd_slice_read_elem(U, i, i) == 0) {
      mzd_slice_write_elem(U, i, i, random()&bitmask) ;
    }
  }
  mzd_slice_set_canary(U);
  return U;
}

mzd_slice_t *random_mzd_slice_t_lower_left(gf2e *ff, rci_t m) {
  const int bitmask = (1<<ff->degree)-1;
  mzd_slice_t *L = random_mzd_slice_t(ff, m, m);
  for(rci_t i=0; i<m; i++) {
    for(rci_t j=i+1; j<m; j++) {
      mzd_slice_write_elem(L, i, j, 0);
    }
    while(mzd_slice_read_elem(L, i, i) == 0) {
      mzd_slice_write_elem(L, i, i, random()&bitmask) ;
    }
  }
  mzd_slice_set_canary(L);
  return L;
}

int test_mzed_trsm_upper_left(gf2e *ff, rci_t m, rci_t n) {
  int fail_ret = 0;

  mzed_t *U = random_mzed_t_upper_left(ff, m);
  mzed_t *H = mzed_copy(NULL, U);
  mzed_t *B = random_mzed_t(ff, m, n);
  mzed_t *X = mzed_copy(NULL, B);

  mzed_set_canary(H);
  mzed_set_canary(X);

  mzed_trsm_upper_left(U, X);

  m4rie_check(mzed_canary_is_alive(U));
  m4rie_check(mzed_canary_is_alive(X));

  mzed_addmul(B, U, X);

  m4rie_check( mzed_canary_is_alive(B));
  m4rie_check( mzed_canary_is_alive(U));
  m4rie_check( mzed_canary_is_alive(X));

  m4rie_check(mzed_is_zero(B) == 1);
  m4rie_check(mzed_cmp(U,H) == 0);

  mzed_free(U);
  mzed_free(H);
  mzed_free(B);
  mzed_free(X);

  U = random_mzed_t(ff, m, m);
  B = random_mzed_t(ff, m, n);
  X = mzed_copy(NULL, B);

  mzed_set_canary(X);

  const int bitmask = (1<<ff->degree)-1;
  for(rci_t i=0; i<m; i++) {
    while(mzed_read_elem(U, i, i) == 0) {
      mzed_write_elem(U, i, i, random()&bitmask) ;
    }
  };
  H = mzed_copy(NULL, U);
  mzed_set_canary(H);

  mzed_trsm_upper_left(U, X);

  m4rie_check( mzed_canary_is_alive(U) );
  m4rie_check( mzed_canary_is_alive(X) );

  m4rie_check(mzed_cmp(U,H) == 0);

  mzed_set_canary(U);

  for(rci_t i=0; i<m; i++) {
    for(rci_t j=0; j<i; j++) {
      mzed_write_elem(U, i, j, 0);
    }
  }
  mzed_addmul(B, U, X);

  m4rie_check( mzed_canary_is_alive(B) );
  m4rie_check( mzed_canary_is_alive(U) );
  m4rie_check( mzed_canary_is_alive(X) );
  
  m4rie_check(mzed_is_zero(B) == 1);

  mzed_free(U);
  mzed_free(H);
  mzed_free(B);
  mzed_free(X);

  return fail_ret;
}

int test_mzed_trsm_upper_left_echelonize(gf2e *ff, rci_t m, rci_t n) {
  const int bitmask = (1<<ff->degree)-1;

  int fail_ret = 0;

  mzed_t *A = random_mzed_t(ff, m, m+n);
  mzed_t *U = mzed_init_window(A, 0, 0, m, m);
  mzed_t *B = mzed_init_window(A, 0, m, m, m+n);


  for(rci_t i=0; i<m; i++) {
    for(rci_t j=0; j<i; j++) {
      mzed_write_elem(U, i, j, 0);
    }
    while(mzed_read_elem(U, i, i) == 0) {
      mzed_write_elem(U, i, i, random()&bitmask) ;
    }
  }
  mzed_t *C = mzed_copy(NULL, A);
  mzed_echelonize(C, 1);

  mzed_trsm_upper_left(U, B);

  mzed_set_ui(U, 0);

  for(rci_t i=0; i<m; i++) {
    mzed_write_elem(U, i, i, 1);
  }

  m4rie_check(mzed_cmp(C,A) == 0);

  if (mzed_cmp(C,A) != 0) {
    mzed_print(C);
    printf("\n");
    mzed_print(A);
  }

  mzed_free(A);
  mzed_free_window(U);
  mzed_free_window(B);
  mzed_free(C);
  
  return fail_ret;
}


int test_mzed_trsm_lower_left(gf2e *ff, rci_t m, rci_t n) {
  int fail_ret = 0;

  mzed_t *L = random_mzed_t_lower_left(ff, m);
  mzed_t *H = mzed_copy(NULL, L);
  mzed_t *B = random_mzed_t(ff, m, n);
  mzed_t *X = mzed_copy(NULL, B);

  mzed_set_canary(H);
  mzed_set_canary(X);

  mzed_trsm_lower_left(L, X);

  m4rie_check(mzed_canary_is_alive(L));
  m4rie_check(mzed_canary_is_alive(X));

  mzed_addmul(B, L, X);

  m4rie_check( mzed_canary_is_alive(B));
  m4rie_check( mzed_canary_is_alive(L));
  m4rie_check( mzed_canary_is_alive(X));

  m4rie_check(mzed_is_zero(B) == 1);
  m4rie_check(mzed_cmp(L,H) == 0);

  mzed_free(L);
  mzed_free(H);
  mzed_free(B);
  mzed_free(X);

  L = random_mzed_t(ff, m, m);
  B = random_mzed_t(ff, m, n);
  X = mzed_copy(NULL, B);

  mzed_set_canary(X);

  const int bitmask = (1<<ff->degree)-1;
  for(rci_t i=0; i<m; i++) {
    while(mzed_read_elem(L, i, i) == 0) {
      mzed_write_elem(L, i, i, random()&bitmask) ;
    }
  };
  H = mzed_copy(NULL, L);
  mzed_set_canary(H);

  mzed_trsm_lower_left(L, X);

  m4rie_check( mzed_canary_is_alive(L) );
  m4rie_check( mzed_canary_is_alive(X) );

  m4rie_check(mzed_cmp(L,H) == 0);

  mzed_set_canary(L);

  for(rci_t i=0; i<m; i++) {
    for(rci_t j=i+1; j<m; j++) {
      mzed_write_elem(L, i, j, 0);
    }
  }
  mzed_addmul(B, L, X);

  m4rie_check( mzed_canary_is_alive(B) );
  m4rie_check( mzed_canary_is_alive(L) );
  m4rie_check( mzed_canary_is_alive(X) );
  
  m4rie_check(mzed_is_zero(B) == 1);

  mzed_free(L);
  mzed_free(H);
  mzed_free(B);
  mzed_free(X);

  return fail_ret;
}

int test_mzd_slice_trsm_upper_left(gf2e *ff, rci_t m, rci_t n) {
  int fail_ret = 0;

  mzd_slice_t *U = random_mzd_slice_t_upper_left(ff, m);
  mzd_slice_t *H = mzd_slice_copy(NULL, U);
  mzd_slice_t *B = random_mzd_slice_t(ff, m, n);
  mzd_slice_t *X = mzd_slice_copy(NULL, B);

  mzd_slice_set_canary(H);
  mzd_slice_set_canary(X);

  mzd_slice_trsm_upper_left(U, X);

  m4rie_check( mzd_slice_canary_is_alive(U) );
  m4rie_check( mzd_slice_canary_is_alive(X) );

  /**
   * @TODO:  mzd_slice_addmul is not 'canary safe' because mzd_addmul() isn't
   */
  mzd_slice_clear_canary(X);
  mzd_slice_addmul(B, U, X);

  m4rie_check( mzd_slice_canary_is_alive(B) );
  m4rie_check( mzd_slice_canary_is_alive(U) );

  m4rie_check(mzd_slice_is_zero(B) == 1);
  m4rie_check(mzd_slice_cmp(U,H) == 0);

  mzd_slice_free(U);
  mzd_slice_free(H);
  mzd_slice_free(B);
  mzd_slice_free(X);

  U = random_mzd_slice_t(ff, m, m);
  B = random_mzd_slice_t(ff, m, n);

  X = mzd_slice_copy(NULL, B);
  mzd_slice_set_canary(X);

  const int bitmask = (1<<ff->degree)-1;
  for(rci_t i=0; i<m; i++) {
    while(mzd_slice_read_elem(U, i, i) == 0) {
      mzd_slice_write_elem(U, i, i, random()&bitmask) ;
    }
  };
  H = mzd_slice_copy(NULL, U);
  mzd_slice_set_canary(H);

  mzd_slice_trsm_upper_left(U, X);

  m4rie_check( mzd_slice_canary_is_alive(U) );
  m4rie_check( mzd_slice_canary_is_alive(X) );
  
  m4rie_check(mzd_slice_cmp(U,H) == 0);

  for(rci_t i=0; i<m; i++) {
    for(rci_t j=0; j<i; j++) {
      mzd_slice_write_elem(U, i, j, 0);
    }
  }
  /**
   * @TODO:  mzd_slice_addmul is not 'canary safe' because mzd_addmul() isn't
   */
  mzd_slice_clear_canary(X);
  mzd_slice_addmul(B, U, X);

  m4rie_check( mzd_slice_canary_is_alive(B) );
  m4rie_check( mzd_slice_canary_is_alive(U) );

  m4rie_check(mzd_slice_is_zero(B) == 1);

  mzd_slice_free(U);
  mzd_slice_free(H);
  mzd_slice_free(B);
  mzd_slice_free(X);

  return fail_ret;
}

int test_mzd_slice_trsm_lower_left(gf2e *ff, rci_t m, rci_t n) {
  int fail_ret = 0;

  mzd_slice_t *L = random_mzd_slice_t_lower_left(ff, m);
  mzd_slice_t *H = mzd_slice_copy(NULL, L);
  mzd_slice_t *B = random_mzd_slice_t(ff, m, n);
  mzd_slice_t *X = mzd_slice_copy(NULL, B);

  mzd_slice_set_canary(H);
  mzd_slice_set_canary(X);

  mzd_slice_trsm_lower_left(L, X);

  m4rie_check( mzd_slice_canary_is_alive(L) );
  m4rie_check( mzd_slice_canary_is_alive(X) );

  /**
   * @TODO:  mzd_slice_addmul is not 'canary safe' because mzd_addmul() isn't
   */
  mzd_slice_clear_canary(X);
  mzd_slice_addmul(B, L, X);

  m4rie_check( mzd_slice_canary_is_alive(B) );
  m4rie_check( mzd_slice_canary_is_alive(L) );

  m4rie_check(mzd_slice_is_zero(B) == 1);
  m4rie_check(mzd_slice_cmp(L,H) == 0);

  mzd_slice_free(L);
  mzd_slice_free(H);
  mzd_slice_free(B);
  mzd_slice_free(X);

  L = random_mzd_slice_t(ff, m, m);
  B = random_mzd_slice_t(ff, m, n);

  X = mzd_slice_copy(NULL, B);
  mzd_slice_set_canary(X);

  const int bitmask = (1<<ff->degree)-1;
  for(rci_t i=0; i<m; i++) {
    while(mzd_slice_read_elem(L, i, i) == 0) {
      mzd_slice_write_elem(L, i, i, random()&bitmask) ;
    }
  };
  H = mzd_slice_copy(NULL, L);
  mzd_slice_set_canary(H);

  mzd_slice_trsm_lower_left(L, X);

  m4rie_check( mzd_slice_canary_is_alive(L) );
  m4rie_check( mzd_slice_canary_is_alive(X) );
  
  m4rie_check(mzd_slice_cmp(L,H) == 0);

  for(rci_t i=0; i<m; i++) {
    for(rci_t j=i+1; j<m; j++) {
      mzd_slice_write_elem(L, i, j, 0);
    }
  }
  /**
   * @TODO:  mzd_slice_addmul is not 'canary safe' because mzd_addmul() isn't
   */
  mzd_slice_clear_canary(X);
  mzd_slice_addmul(B, L, X);

  m4rie_check( mzd_slice_canary_is_alive(B) );
  m4rie_check( mzd_slice_canary_is_alive(L) );

  m4rie_check(mzd_slice_is_zero(B) == 1);

  mzd_slice_free(L);
  mzd_slice_free(H);
  mzd_slice_free(B);
  mzd_slice_free(X);

  return fail_ret;
}


int test_batch(gf2e *ff, rci_t m, rci_t n) {
  int fail_ret = 0;
  printf("trsm: k: %2d, minpoly: 0x%03x m: %5d, n: %5d ",(int)ff->degree, (unsigned int)ff->minpoly, (int)m,(int)n);

  m4rie_check(test_mzed_trsm_lower_left(ff, m, n) == 0); printf("."); fflush(0);
  m4rie_check(test_mzed_trsm_upper_left(ff, m, n) == 0); printf("."); fflush(0);
  m4rie_check(test_mzed_trsm_lower_left(ff, n, m) == 0); printf("."); fflush(0);
  m4rie_check(test_mzed_trsm_upper_left(ff, n, m) == 0); printf("."); fflush(0);

  // m4rie_check(test_mzed_trsm_upper_left_echelonize(ff, m, n) == 0); printf("."); fflush(0);
  // m4rie_check(test_mzed_trsm_upper_left_echelonize(ff, n, m) == 0); printf("."); fflush(0);

  if(ff->degree <= __M4RIE_MAX_KARATSUBA_DEGREE) {
    m4rie_check(test_mzd_slice_trsm_lower_left(ff, m, n) == 0); printf("."); fflush(0);
    m4rie_check(test_mzd_slice_trsm_lower_left(ff, n, m) == 0); printf("."); fflush(0);
    m4rie_check(test_mzd_slice_trsm_upper_left(ff, m, n) == 0); printf("."); fflush(0);
    m4rie_check(test_mzd_slice_trsm_upper_left(ff, n, m) == 0); printf("."); fflush(0);
  } else {
    printf("    "); fflush(0);
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
    fail_ret += test_batch(ff[k],   1,   1);
    fail_ret += test_batch(ff[k],   1,   2);
    fail_ret += test_batch(ff[k],  11,  12);
    fail_ret += test_batch(ff[k],  21,  22);
    fail_ret += test_batch(ff[k],  13,   2);
    fail_ret += test_batch(ff[k],  32,  33);
    fail_ret += test_batch(ff[k],  63,  64);
    fail_ret += test_batch(ff[k],  65,  1);
    fail_ret += test_batch(ff[k], 127, 128);
    fail_ret += test_batch(ff[k], 200,  20);
  };

  for(int k=2; k<=10; k++) {
    gf2e_free(ff[k]);
  }
  if (fail_ret == 0)
    printf("success\n");

  return fail_ret;
}

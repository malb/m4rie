/**
 * \file test_smallops.c
 * \brief Test code for auxilary routines
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

/******************************************************************************
*
*            M4RIE: Linear Algebra over GF(2^e)
*
*    Copyright (C) 2010-2012 Martin Albrecht <martinralbrecht@googlemail.com>
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

int test_gf2e(gf2e *ff) {
  int fail_ret = 0;
  for(word a=1; a < __M4RI_TWOPOW(ff->degree); a++) {
    word a_inv = ff->inv(ff, a);
    fail_ret += ((a == ff->inv(ff, a_inv)) ^ 1);
  }
  return fail_ret;
}

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
  
  mzd_slice_set_canary(a);
  mzd_slice_set_canary(b);
  mzd_slice_set_canary(c);

  mzd_slice_add(c, a, b);

  mzd_slice_t *d = mzd_slice_copy(NULL, c);
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

  m4rie_check( test_gf2e(ff) == 0); printf("."); fflush(0);

  if (fail_ret == 0)
    printf(" passed\n");
  else
    printf(" FAILED\n");
  return fail_ret;
}

int main(int argc, char **argv) {

  gf2e *ff;
  int fail_ret = 0;

  for(int k=2; k<=16; k++) {
    ff = gf2e_init(irreducible_polynomials[k][1]);

    fail_ret += test_batch(ff,   2,   m4ri_radix/gf2e_degree_to_w(ff));
    fail_ret += test_batch(ff,   2, 2*m4ri_radix/gf2e_degree_to_w(ff));
    fail_ret += test_batch(ff,   2, 3*m4ri_radix/gf2e_degree_to_w(ff));
    fail_ret += test_batch(ff,   2, 4*m4ri_radix/gf2e_degree_to_w(ff));
    fail_ret += test_batch(ff,   4,   3);
    fail_ret += test_batch(ff,   1,   2);
    fail_ret += test_batch(ff,  10,  11);
    fail_ret += test_batch(ff,  20,  19);
    fail_ret += test_batch(ff,  32,  64);
    fail_ret += test_batch(ff,  63,  65);
    fail_ret += test_batch(ff,  64,  65);
    fail_ret += test_batch(ff,  64, 128);
    fail_ret += test_batch(ff,  65, 129);
    fail_ret += test_batch(ff, 201, 200);
    fail_ret += test_batch(ff, 217,   2);

    gf2e_free(ff);
  }

  return fail_ret;
}

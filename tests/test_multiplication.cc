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
#include <m4rie/m4rie.h>

using namespace M4RIE;

int test_addmul() {
  int fail_ret = 0;
  for(size_t k=2; k<=10; k++) {
    FiniteField *F = (FiniteField*)(new GFqDom<int>(2,k));
    gf2e *ff = gf2e_init_givgfq(F);
    for(size_t i=0; i<(k*256/(1<<k)); i++) {
      size_t m = random() & 255;
      size_t n = random() & 255;
      size_t l = random() & 255;
      m = m ? m : 1;
      n = n ? n : 1;
      l = l ? n : 1;
      mzed_t *A = mzed_init(ff,m,l);
      mzed_randomize(A);
      mzed_t *B = mzed_init(ff,l,n);
      mzed_randomize(B);

      mzed_t *C0 = mzed_init(ff,m,n);
      mzed_randomize(A);
      mzed_t *C1 = mzed_copy(NULL, C0);
      mzed_t *C2 = mzed_copy(NULL, C0);
      mzed_t *C3 = mzed_copy(NULL, C0);
        
      mzed_addmul_travolta(C0, A, B);
      mzed_addmul_naive(C1, A, B);
      mzed_addmul_strassen(C2, A, B, 64);
      printf("addmul: m: %5d, l: %5d, n: %5d, k: %2d ... ",m,l,n,k);

      m4rie_check( mzed_cmp(C0, C1) == 0);
      m4rie_check( mzed_cmp(C1, C2) == 0);

      if (k == 2) {
        mzed_addmul_karatsuba(C3, A, B);
        m4rie_check( mzed_cmp(C2, C3) == 0);
      }

      if (fail_ret == 0)
        printf("pass\n");
      else
        printf("FAIL\n");

      mzed_free(A);
      mzed_free(B);
      mzed_free(C0);
      mzed_free(C1);
      mzed_free(C2);
      mzed_free(C3);
    }
    gf2e_free(ff);
    delete F;
  }
  return fail_ret;
}

int test_mul() {
  int fail_ret = 0;
  for(size_t k=2; k<=10; k++) {
    FiniteField *F = (FiniteField*)(new GFqDom<int>(2,k));
    gf2e *ff = gf2e_init_givgfq(F);
    for(size_t i=0; i<(k*256/(1<<k)); i++) {
      size_t m = random() & 255;
      size_t n = random() & 255;
      size_t l = random() & 255;
      m = m ? m : 1;
      n = n ? n : 1;
      l = l ? n : 1;
      mzed_t *A = mzed_init(ff,m,l);
      mzed_randomize(A);
      mzed_t *B = mzed_init(ff,l,n);
      mzed_randomize(B);

      
      mzed_t *C0 = mzed_mul_travolta(NULL, A, B);
      mzed_t *C1 = mzed_mul_naive(NULL, A, B);
      mzed_t *C2 = mzed_mul_strassen(NULL, A, B, 64);
      printf("mul: m: %5d, l: %5d, n: %5d, k: %2d ... ",m,l,n,k);

      m4rie_check( mzed_cmp(C0, C1) == 0);
      m4rie_check( mzed_cmp(C1, C2) == 0);

      if (k == 2) {
        mzed_t *C3 = mzed_mul_karatsuba(NULL, A, B);
        m4rie_check( mzed_cmp(C2, C3) == 0);
        mzed_free(C3);
      }

      if (fail_ret == 0)
        printf("pass\n");
      else 
        printf("FAIL\n");

      mzed_free(A);
      mzed_free(B);
      mzed_free(C0);
      mzed_free(C1);
      mzed_free(C2);
    }
    gf2e_free(ff);
    delete F;
  }
  return fail_ret;
}

int main(int argc, char **argv) {
  return test_mul() + test_addmul();
}

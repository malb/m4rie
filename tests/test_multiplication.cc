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

#include <gf2e_cxx/finite_field_givaro.h>
#include <m4rie/m4rie.h>

using namespace M4RIE;

int test_addmul() {
  int pass = 0;
  for(size_t k=2; k<=10; k++) {
    FiniteField *F = (FiniteField*)(new GFqDom<int>(2,k));
    gf2e *ff = gf2e_init_givgfq(F);
    for(size_t i=0; i<(k*512/(1<<k)); i++) {
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
      
      mzed_addmul_travolta(C0, A, B);
      mzed_addmul_naive(C1, A, B);
      mzed_addmul(C2, A, B);
      printf("addmul: m: %5d, l: %5d, n: %5d, k: %2d ... ",m,l,n,k);

      if (!mzed_cmp(C0,C1) && !mzed_cmp(C1,C2)) {
        printf("pass\n");
      } else {
        printf("FAIL\n");
        pass = 0;
      }
      mzed_free(A);
      mzed_free(B);
      mzed_free(C0);
      mzed_free(C1);
      mzed_free(C2);
    }
    gf2e_free(ff);
    delete F;
  }
  return pass;
}

int test_mul() {
  int pass = 0;
  for(size_t k=2; k<=10; k++) {
    FiniteField *F = (FiniteField*)(new GFqDom<int>(2,k));
    gf2e *ff = gf2e_init_givgfq(F);
    for(size_t i=0; i<(k*512/(1<<k)); i++) {
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
      mzed_t *C2 = mzed_mul(NULL, A, B);
      printf("mul: m: %5d, l: %5d, n: %5d, k: %2d ... ",m,l,n,k);

      if (!mzed_cmp(C0,C1) && !mzed_cmp(C1,C2)) {
        printf("pass\n");
      } else {
        printf("FAIL\n");
        pass = 0;
      }
      mzed_free(A);
      mzed_free(B);
      mzed_free(C0);
      mzed_free(C1);
      mzed_free(C2);
    }
    gf2e_free(ff);
    delete F;
  }
  return pass;
}

int main(int argc, char **argv) {
  return test_mul() + test_addmul();
}

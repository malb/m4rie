#ifndef M4RIE_MZD_POLY_H
#define M4RIE_MZD_POLY_H

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

#include <m4ri/m4ri.h>
#include "mzd_ptr.h"
#include "gf2x.h"

/**
 * \brief will be the data type for matrices over \GF2[x] in the future
 *
 * @warning Do not use yet.
 */

typedef struct {
  mzd_t **x;   /**< Coefficients. */
  rci_t nrows; /**< Number of rows. */
  rci_t ncols; /**< Number of columns. */
  deg_t depth;   /**< Degree +1      */
} mzd_poly_t;

static inline mzd_poly_t *_mzd_poly_add(mzd_poly_t *C, const mzd_poly_t *A, const mzd_poly_t *B, unsigned int offset) {
  _mzd_ptr_add(C->x+offset, (const mzd_t**)A->x, (const mzd_t**)B->x, A->depth);
  return C;
}

static inline mzd_poly_t *mzd_poly_add(mzd_poly_t *C, const mzd_poly_t *A, const mzd_poly_t *B) {
  assert(C->depth >= A->depth && A->depth == B->depth);
  return _mzd_poly_add(C, A, B, 0);
}

static inline mzd_poly_t *mzd_poly_init(const deg_t d, const rci_t m, const rci_t n) {
  mzd_poly_t *A = (mzd_poly_t*)m4ri_mm_malloc(sizeof(mzd_poly_t));
  A->x = (mzd_t**)m4ri_mm_malloc(sizeof(mzd_t*)*(d+1));

  A->nrows = m;
  A->ncols = n;
  A->depth = d+1;

  for(int i=0; i<A->depth; i++)
    A->x[i] = mzd_init(m,n);
  return A;
}

static inline void mzd_poly_free(mzd_poly_t *A) {
  for(int i=0; i<A->depth; i++)
   mzd_free(A->x[i]);
#if __M4RI_USE_MM_MALLOC
  _mm_free(A);
#else
  free(A);
#endif
}

static inline mzd_poly_t *_mzd_poly_adapt_depth(mzd_poly_t *A, const deg_t new_depth) {
  if (new_depth < A->depth) {
    for(int i=new_depth; i<A->depth; i++) {
      mzd_free(A->x[i]);
      A->x[i] = NULL;
    }
  } else {
    for(int i=A->depth; i<new_depth; i++) {
      A->x[i] = mzd_init(A->nrows,A->ncols);
    }
  }
  A->depth = new_depth;
  return A;
}

static inline mzd_poly_t *_mzd_poly_addmul_naive(mzd_poly_t *C, const mzd_poly_t *A, const mzd_poly_t *B) {
  if (C == NULL)
    C = mzd_poly_init(A->depth+B->depth-1, A->nrows, B->ncols);

  for(unsigned int i=0; i<A->depth; i++) {
    for(unsigned int j=0; j<B->depth; j++) {
      mzd_addmul(C->x[i+j], A->x[i], B->x[j], 0);
    }
  }
  return C;
}

static inline mzd_poly_t *_mzd_poly_addmul_balanced(mzd_poly_t *C, const mzd_poly_t *A, const mzd_poly_t *B) {
  assert(A->depth == B->depth);

  if (C == NULL)
    C = mzd_poly_init(A->depth+B->depth-1, A->nrows, B->ncols);
  switch(A->depth) {
  case 0:
    m4ri_die("depth 0: seriously?");
  case 1: mzd_addmul(C->x[0], A->x[0], B->x[0], 0); break;
  case 2: _mzd_ptr_addmul_karatsuba2(NULL, C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
  case 3: _mzd_ptr_addmul_karatsuba3(NULL, C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
  case 4: _mzd_ptr_addmul_karatsuba4(NULL, C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
  case 5: _mzd_ptr_addmul_karatsuba5(NULL, C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
  case 6: _mzd_ptr_addmul_karatsuba6(NULL, C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
  case 7: _mzd_ptr_addmul_karatsuba7(NULL, C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
  case 8: _mzd_ptr_addmul_karatsuba8(NULL, C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
  default:
    _mzd_poly_addmul_naive(C, A, B); break;
  }
  return C;
}

/**
 * \brief C += A*B using arithmetic in GF(2^log2(d)) if C has degree d.
 */

mzd_poly_t *_mzd_poly_addmul_ext1(mzd_poly_t *C, mzd_poly_t *A, mzd_poly_t *B);

/**
 * \brief Return -1,0,1 if if A < B, A == B or A > B respectively.
 *
 * \param A Matrix.
 * \param B Matrix.
 *
 * \note This comparison is not well defined (except for !=0) mathematically and relatively
 * arbitrary.
 *
 * \ingroup Comparison
 */

static inline int mzd_poly_cmp(mzd_poly_t *A, mzd_poly_t *B) {
  int r = 0;
  if ((A->depth != B->depth) ) {
    if (A->depth < B->depth)
      return -1;
    else 
      return 1;
  }
  for(int i=0; i<A->depth; i++)
    r |= mzd_cmp(A->x[i],B->x[i]);
  return r;
}

/**
 * \brief Fill matrix A with random elements.
 *
 * \param A Matrix
 *
 * \todo Allow the user to provide a RNG callback.
 *
 * \ingroup Assignment
 */

static inline void mzd_poly_randomize(mzd_poly_t *A) {
  for(int i=0; i<A->depth; i++)
    mzd_randomize(A->x[i]);
}

#endif //M4RIE_MZD_POLY_H

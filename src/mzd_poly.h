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
#include <stdarg.h>

/********************************************************************
 * Internal representation
 *******************************************************************/

static inline void _poly_add(mzd_t **c, const mzd_t **a, const mzd_t **b,const unsigned int length) {
  switch(length) {
  case 16: mzd_add(c[15], a[15], b[15]);
  case 15: mzd_add(c[14], a[14], b[14]);
  case 14: mzd_add(c[13], a[13], b[13]);
  case 13: mzd_add(c[12], a[12], b[12]);
  case 12: mzd_add(c[11], a[11], b[11]);
  case 11: mzd_add(c[10], a[10], b[10]);
  case 10: mzd_add(c[ 9], a[ 9], b[ 9]);
  case  9: mzd_add(c[ 8], a[ 8], b[ 8]);
  case  8: mzd_add(c[ 7], a[ 7], b[ 7]);
  case  7: mzd_add(c[ 6], a[ 6], b[ 6]);
  case  6: mzd_add(c[ 5], a[ 5], b[ 5]);
  case  5: mzd_add(c[ 4], a[ 4], b[ 4]);
  case  4: mzd_add(c[ 3], a[ 3], b[ 3]);
  case  3: mzd_add(c[ 2], a[ 2], b[ 2]);
  case  2: mzd_add(c[ 1], a[ 1], b[ 1]);
  case  1: mzd_add(c[ 0], a[ 0], b[ 0]);
  case  0:
    break;
  default:
    for(int i=0; i<length; i++)
      mzd_add(c[ i], a[ i], b[ i]);
  }
}

void _poly_addmul2(mzd_t **X, const mzd_t **a, const mzd_t **b);
void _poly_addmul4(mzd_t **X, const mzd_t **a, const mzd_t **b);

/*********************************************************************
 * mzd_poly_t will be the data type for matrices over GF(2)[x] in the
 * future
 *
 * DO NOT USE YET.
 *
 *********************************************************************/

typedef int deg_t;

typedef struct {
  mzd_t **x;
  rci_t nrows; /**< Number of rows. */
  rci_t ncols; /**< Number of columns. */
  deg_t depth;   /**< Degree +1      */
} mzd_poly_t;

static inline mzd_poly_t *_mzd_poly_add(mzd_poly_t *C, const mzd_poly_t *A, const mzd_poly_t *B, unsigned int offset) {
  _poly_add(C->x+offset, (const mzd_t**)A->x, (const mzd_t**)B->x, A->depth);
  return C;
}

static inline mzd_poly_t *mzd_poly_add(mzd_poly_t *C, const mzd_poly_t *A, const mzd_poly_t *B) {
  assert(C->depth >= A->depth && A->depth == B->depth);
  return _mzd_poly_add(C, A, B, 0);
}

static inline mzd_poly_t *mzd_poly_init(const deg_t d, const rci_t m, const rci_t n) {
  mzd_poly_t *A = m4ri_mm_malloc(sizeof(mzd_poly_t));
  A->x = (mzd_t**)m4ri_mm_malloc(sizeof(mzd_t*)*d);

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

#endif //M4RIE_MZD_POLY_H

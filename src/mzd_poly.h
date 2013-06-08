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

static inline void _poly_add(mzd_t **c, const mzd_t **a, const mzd_t **b, const unsigned int length) {
  switch(length) {
  case 32: mzd_add(c[31], a[31], b[31]);
  case 31: mzd_add(c[30], a[30], b[30]);
  case 30: mzd_add(c[29], a[29], b[29]);
  case 29: mzd_add(c[28], a[28], b[28]);
  case 28: mzd_add(c[27], a[27], b[27]);
  case 27: mzd_add(c[26], a[26], b[26]);
  case 26: mzd_add(c[25], a[25], b[25]);
  case 25: mzd_add(c[24], a[24], b[24]);
  case 24: mzd_add(c[23], a[23], b[23]);
  case 23: mzd_add(c[22], a[22], b[22]);
  case 22: mzd_add(c[21], a[21], b[21]);
  case 21: mzd_add(c[20], a[20], b[20]);
  case 20: mzd_add(c[19], a[19], b[19]);
  case 19: mzd_add(c[18], a[18], b[18]);
  case 18: mzd_add(c[17], a[17], b[17]);
  case 17: mzd_add(c[16], a[16], b[16]);
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
  case 1:
    mzd_addmul(C->x[0], A->x[0], B->x[0], 0); break;
  case 2:
    _poly_addmul2(C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
  case 3:
    _mzd_poly_addmul_naive(C, A, B); break;
  case 4:
    _poly_addmul4(C->x, (const mzd_t**)A->x, (const mzd_t**)B->x); break;
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

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

#include "travolta.h"

/**
 * \brief Perform Gaussian reduction to reduced row echelon form on a
 * submatrix.
 * 
 * The submatrix has dimension at most k starting at r x c of A. Checks
 * for pivot rows up to row endrow (exclusive). Terminates as soon as
 * finding a pivot column fails.
 *
 * \param A Matrix.
 * \param r First row.
 * \param c First column.
 * \param k Maximal dimension of identity matrix to produce.
 * \param end_row Maximal row index (exclusive) for rows to consider
 * for inclusion.
 */

static inline int _mzed_gauss_submatrix_full(mzed_t *A, size_t r, size_t c, size_t end_row, int k) {
  size_t i,j,l;
  size_t start_row = r;
  int found;
  word tmp;

  gf2e *ff = A->finite_field;

  for (j=c; j<c+k; j++) {
    found = 0;
    for (i=start_row; i< end_row; i++) {
      /* first we need to clear the first columns */
      for (l=0; l<j-c; l++) {
        tmp = mzed_read_elem(A, i, c+l);
        if (tmp) mzed_add_multiple_of_row(A, i, A, r+l, ff->mul[tmp], c+l);
      }
      /* pivot? */
      const word x = mzed_read_elem(A, i, j);
      if (x) {
        mzed_rescale_row(A, i, j, ff->inv[x]);
        mzd_row_swap(A->x, i, start_row);

        /* clear above */
        for (l=r; l<start_row; l++) {
          tmp = mzed_read_elem(A, l, j);
          if (tmp) mzed_add_multiple_of_row(A, l, A, start_row, ff->mul[tmp], j);
        }
        start_row++;
        found = 1;
        break;
      }
    }
    if (found==0) {
      return j - c;
    }
  }
  return j - c;
}


void mzed_make_table(mzed_t *T, mzed_t *A, size_t r, size_t c, gf2e *ff) {
  mzd_set_ui(T->x,0);

  for(size_t i=0; i< TWOPOW(ff->degree); i++) {
    word *X = ff->mul[i];
    mzed_add_multiple_of_row(T, i, A, r, X, c);
  }
}

size_t mzed_echelonize_travolta(mzed_t *A, int full) {
  gf2e* ff = A->finite_field;

  size_t r,c,i;

  size_t k = ff->degree;
  size_t kk = 1;
  size_t kbar = 0;

  mzed_t *T0 = mzed_init(ff, TWOPOW(k), A->ncols);
  mzed_t *T1 = mzed_init(ff, TWOPOW(k), A->ncols);
  
  /* this is dummy, we keep it for compatibility with the M4RI functions */
  size_t *L0 = (size_t *)m4ri_mm_calloc(TWOPOW(k), sizeof(size_t));
  for(i=0;i<TWOPOW(k);i++) L0[i] = i;

  r = 0;
  c = 0;
  while(c < A->ncols) {
    if(c+kk > A->ncols) kk = A->ncols - c;
    kbar = _mzed_gauss_submatrix_full(A, r, c, A->nrows, kk);
    if(kbar > k) {
      mzed_make_table(T0, A, r, c, ff);
      mzed_make_table(T1, A, r+1, c+1, ff);
      if(kbar == kk)
        mzd_process_rows2(A->x, r+1, A->nrows, c*A->width, kbar*A->width, T0->x, L0, T1->x, L0);
      if(full)
        mzd_process_rows2(A->x, 0, r, c*A->width, kbar*A->width, T0->x, L0, T1->x, L0);
    } else if (kbar > 0) {
      mzed_make_table(T0, A, r, c, ff);
      if(kbar == kk)
        mzd_process_rows(A->x, r+1, A->nrows, c*A->width, kbar*A->width, T0->x, L0);
      if(full)
        mzd_process_rows(A->x, 0, r, c*A->width, kbar*A->width, T0->x, L0);
    } else {
      c++;
    }
    r += kbar;
    c += kbar;
  }
  m4ri_mm_free(L0);
  mzed_free(T0);
  mzed_free(T1);
  return r;
}


mzed_t *mzed_addmul_travolta(mzed_t *C, mzed_t *A, mzed_t *B) {
  if (C->nrows != A->nrows || C->ncols != B->ncols || C->finite_field != A->finite_field) {
    m4ri_die("mzd_mul_naive: Provided return matrix has wrong dimensions or wrong base field.\n");
  }
  return _mzed_mul_travolta(C, A, B);
}


mzed_t *_mzed_mul_travolta(mzed_t *C, mzed_t *A, mzed_t *B) {
  int k = A->finite_field->degree;

  mzed_t *T0 = mzed_init(C->finite_field, TWOPOW(k), B->ncols);
  
  for(size_t c=0; c<B->nrows; c++) {
    mzed_make_table(T0, B, c, 0, A->finite_field);
    for(size_t r=0; r<A->nrows; r++) {
      mzd_combine(C->x, r, 0, C->x, r, 0, T0->x, mzed_read_elem(A, r, c), 0);
    }
  }

  mzed_free(T0);

  return C;
}

mzed_t *mzed_mul_travolta(mzed_t *C, mzed_t *A, mzed_t *B) {
  if (C==NULL) {
    C = mzed_init(A->finite_field, A->nrows, B->ncols);
  } else {
    if (C->nrows != A->nrows || C->ncols != B->ncols || C->finite_field != A->finite_field) {
      m4ri_die("mzd_mul_naive: Provided return matrix has wrong dimensions or wrong base field.\n");
    }
    mzd_set_ui(C->x, 0);
  }
  return _mzed_mul_travolta(C, A, B);
}


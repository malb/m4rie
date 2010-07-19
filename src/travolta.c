
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

int _mzed_gauss_submatrix_full(mzed_t *A, size_t r, size_t c, size_t end_row, int k) {
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

  /** cf. mzd_echelonize_m4ri **/
  size_t kk = m4ri_opt_k(A->x->nrows, A->x->ncols, 0);
  if (kk>=7) 
    kk = 7;
  if ( (6*(1<<kk)*A->ncols / 8.0) > CPU_L2_CACHE / 2.0 )
    kk -= 1;
  kk = (6*kk)/k;

  /** enforcing bounds **/
  if (kk == 0)
    kk = 1;
  else if (kk > 6)
    kk = 6;

  size_t kbar = 0;

  mzed_t *T0 = mzed_init(ff, TWOPOW(k), A->ncols);
  mzed_t *T1 = mzed_init(ff, TWOPOW(k), A->ncols);
  mzed_t *T2 = mzed_init(ff, TWOPOW(k), A->ncols);
  mzed_t *T3 = mzed_init(ff, TWOPOW(k), A->ncols);
  mzed_t *T4 = mzed_init(ff, TWOPOW(k), A->ncols);
  mzed_t *T5 = mzed_init(ff, TWOPOW(k), A->ncols);
  
  /* this is dummy, we keep it for compatibility with the M4RI functions */
  size_t *L = (size_t *)m4ri_mm_calloc(TWOPOW(k), sizeof(size_t));
  for(i=0;i<TWOPOW(k);i++) L[i] = i;

  r = 0;
  c = 0;
  while(c < A->ncols) {
    if(c+kk > A->ncols) kk = A->ncols - c;
    kbar = _mzed_gauss_submatrix_full(A, r, c, A->nrows, kk);
    if (kbar == 6)  {
      mzed_make_table(T0, A, r,   c,   ff);
      mzed_make_table(T1, A, r+1, c+1, ff);
      mzed_make_table(T2, A, r+2, c+2, ff);
      mzed_make_table(T3, A, r+3, c+3, ff);
      mzed_make_table(T4, A, r+4, c+4, ff);
      mzed_make_table(T5, A, r+5, c+5, ff);
      if(kbar == kk)
        mzd_process_rows6(A->x, r+6, A->nrows, c*A->width, kbar*A->width, T0->x, L, T1->x, L, T2->x, L, T3->x, L, T4->x, L, T5->x, L);
      if(full)
        mzd_process_rows6(A->x, 0, r, c*A->width, kbar*A->width, T0->x, L, T1->x, L, T2->x, L, T3->x, L, T4->x, L, T5->x, L);
    } else if(kbar == 5) {
      mzed_make_table(T0, A, r, c, ff);
      mzed_make_table(T1, A, r+1, c+1, ff);
      mzed_make_table(T2, A, r+2, c+2, ff);
      mzed_make_table(T3, A, r+3, c+3, ff);
      mzed_make_table(T4, A, r+4, c+4, ff);
      if(kbar == kk)
        mzd_process_rows5(A->x, r+5, A->nrows, c*A->width, kbar*A->width, T0->x, L, T1->x, L, T2->x, L, T3->x, L, T4->x, L);
      if(full)
        mzd_process_rows5(A->x, 0, r, c*A->width, kbar*A->width, T0->x, L, T1->x, L, T2->x, L, T3->x, L, T4->x, L);
    } else if(kbar == 4) {
      mzed_make_table(T0, A, r, c, ff);
      mzed_make_table(T1, A, r+1, c+1, ff);
      mzed_make_table(T2, A, r+2, c+2, ff);
      mzed_make_table(T3, A, r+3, c+3, ff);
      if(kbar == kk)
        mzd_process_rows4(A->x, r+4, A->nrows, c*A->width, kbar*A->width, T0->x, L, T1->x, L, T2->x, L, T3->x, L);
      if(full)
        mzd_process_rows4(A->x, 0, r, c*A->width, kbar*A->width, T0->x, L, T1->x, L, T2->x, L, T3->x, L);
    } else if(kbar == 3) {
      mzed_make_table(T0, A, r, c, ff);
      mzed_make_table(T1, A, r+1, c+1, ff);
      mzed_make_table(T2, A, r+2, c+2, ff);
      if(kbar == kk)
        mzd_process_rows3(A->x, r+3, A->nrows, c*A->width, kbar*A->width, T0->x, L, T1->x, L, T2->x, L);
      if(full)
        mzd_process_rows3(A->x, 0, r, c*A->width, kbar*A->width, T0->x, L, T1->x, L, T2->x, L);
    } else if(kbar == 2) {
      mzed_make_table(T0, A, r, c, ff);
      mzed_make_table(T1, A, r+1, c+1, ff);
      if(kbar == kk)
        mzd_process_rows2(A->x, r+2, A->nrows, c*A->width, 2*A->width, T0->x, L, T1->x, L);
      if(full)
        mzd_process_rows2(A->x, 0, r, c*A->width, 2*A->width, T0->x, L, T1->x, L);
    } else if (kbar == 1) {
      mzed_make_table(T0, A, r, c, ff);
      if(kbar == kk)
        mzd_process_rows(A->x, r+1, A->nrows, c*A->width, kbar*A->width, T0->x, L);
      if(full)
        mzd_process_rows(A->x, 0, r, c*A->width, kbar*A->width, T0->x, L);
    } else {
      c++;
    }
    r += kbar;
    c += kbar;
  }
  m4ri_mm_free(L);
  mzed_free(T0);
  mzed_free(T1);
  mzed_free(T2);
  mzed_free(T3);
  mzed_free(T4);
  mzed_free(T5);
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


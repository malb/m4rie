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

void mzed_make_table(mzed_t *T, mzed_t *A, size_t r, size_t c, gf2e *ff) {
  mzd_set_ui(T->x,0);

  for(size_t i=0; i< TWOPOW(ff->degree); i++) {
    word *X = ff->mul[i];
    mzed_add_multiple_of_row(T, i, A, r, X, c);
  }
}

size_t mzed_echelonize_travolta(mzed_t *A, int full) {
  gf2e* ff = A->finite_field;

  size_t start_row,r,c,i;
  size_t nr = A->nrows;
  size_t nc = A->ncols;

  int k = ff->degree;

  mzed_t *T0 = mzed_init(ff, TWOPOW(k), A->ncols);
  size_t *L0 = (size_t *)m4ri_mm_calloc(TWOPOW(k), sizeof(size_t));
  for(i=0;i<TWOPOW(k);i++)
    L0[i] = i;

  start_row = 0;

  for(c=0; c<nc; c++) {
    for(r=start_row; r<nr; r++) {
      word x = mzed_read_elem(A,r,c);
      if (x) {
        word *x_inverse = ff->mul[ff->inv[x]];
        /* rescale row */
        for(i=c; i<nc; i++) {
          mzed_write_elem(A, r, i, x_inverse[mzed_read_elem(A, r, i)]);
        }
        mzd_row_swap(A->x, r, start_row);

        mzed_make_table(T0, A, start_row, c, ff);

        mzd_process_rows(A->x, start_row+1, A->nrows, c*A->width, A->width, T0->x, L0);
        if(full)
          mzd_process_rows(A->x, 0, start_row, c*A->width, A->width, T0->x, L0);
        start_row++;
        break;
      }
    }
  }
  m4ri_mm_free(L0);
  mzed_free(T0);
  return start_row;
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


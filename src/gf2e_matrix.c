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

#include "gf2e_matrix.h"
#include <stdlib.h>

mzed_t *mzed_init(gf2e* k, size_t m, size_t n) {
  mzed_t *A = (mzed_t *)m4ri_mm_malloc(sizeof(mzed_t));

  A->finite_field = k;
  
  switch(k->degree) {
  case 2: 
    A->width = 2; 
    break;
  case 3: 
  case 4: 
    A->width = 4;
    break;
  case 5: 
  case 6: 
  case 7: 
  case 8: 
    A->width = 8;
    break;
  case 9: 
  case 10: 
    A->width = 16;
    break;
  default:
    m4ri_die("degree %d not supported.\n",k->degree);
  }

  A->nrows = m;
  A->ncols = n;
  A->x = mzd_init(m, A->width*n);
  return A;
}

void mzed_free(mzed_t *A) {
  mzd_free(A->x);
  m4ri_mm_free(A);
}

void mzed_randomize(mzed_t *A) {
  int bitmask = (1<<A->finite_field->degree)-1;
  for(size_t r=0; r<A->nrows; r++) {
    for(size_t c=0; c<A->ncols; c++) {
      mzed_write_elem(A,r,c, random()&bitmask);
    }
  }
}

mzed_t *mzed_add(mzed_t *C, const mzed_t *A, const mzed_t *B) {
  if (A->nrows != B->nrows || A->ncols != B->ncols || A->finite_field != B->finite_field) {
    m4ri_die("mzed_add: rows, columns and fields must match.\n");
  }
  if (C == NULL) {
    C = mzed_init(A->finite_field, A->nrows, A->ncols);
  } else if (C != A) {
    if (C->finite_field != A->finite_field || C->nrows != A->nrows || C->ncols != A->ncols) {
      m4ri_die("mzed_add: rows and columns of returned matrix must match.\n");
    }
  }
  mzd_add(C->x, A->x, B->x);
  return C;
}

mzed_t *mzed_copy(mzed_t *A, const mzed_t *B) {
  if (A == B)
    return A;
  if (A == NULL)
    A = mzed_init(B->finite_field, B->nrows, B->ncols);
  if (A->finite_field != B->finite_field || A->nrows != B->nrows || A->ncols != B->ncols) {
    m4ri_die("mzed_copy: target matrix has wrong dimensions or base field.");
  }
  mzd_copy(A->x, B->x);
  return A;
}


size_t mzed_echelonize_naive(mzed_t *A, int full) {
  size_t start_row,r,c,i,elim_start;
  int x,x_inverse = 0;

  size_t nr = A->nrows;
  size_t nc = A->ncols;

  gf2e *ff = A->finite_field;

  start_row = 0;

  for(c=0; c<nc; c++) {
    for(r=start_row; r<nr; r++) {
      x = mzed_read_elem(A, r, c);
      if (x) {
        x_inverse = ff->inv[x];
        /* rescale row */
        for(i=c; i<nc; i++) {
          x = ff->mul[x_inverse][mzed_read_elem(A, r, i)];
          mzed_write_elem(A, r, i, x);
        }
        mzd_row_swap(A->x, r, start_row);
        if (full)
          elim_start = 0;
        else
          elim_start = start_row + 1;
        for(i=elim_start; i<nr; i++) {
          if (i==start_row) 
            continue;
          x = mzed_read_elem(A,i,c);
          if(!x) continue;
          /* clear row */
          mzed_add_multiple_of_row(A, i, A, start_row, ff->mul[x], c);
        }
        start_row++;
        break;
      }
    }
  }
  return start_row;
}



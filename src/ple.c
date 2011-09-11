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

#include "ple.h"

rci_t mzed_ple_naive(mzed_t *A, mzp_t *P, mzp_t *Q) {
  rci_t col_pos = 0;
  rci_t row_pos = 0;
  word tmp = 0;
  gf2e *ff = A->finite_field;
  rci_t i,j;
  int found = 0;

  while (row_pos < A->nrows && col_pos < A->ncols) {
    found = 0;
    for(j=col_pos; j<A->ncols; j++) {
      for(i=row_pos; i<A->nrows; i++) {
        if( (tmp = mzed_read_elem(A, i,j)) != 0) {
          found = 1;
          break;
        }
      }
      if (found)
        break;
    }
    if (found) {
      P->values[row_pos] = i;
      Q->values[row_pos] = j;
      mzed_row_swap(A, row_pos, i);

      mzed_rescale_row(A, row_pos, j+1, ff->mul[ff->inv[tmp]]);

      if(j+1 < A->ncols) {
        for(rci_t l=row_pos+1; l<A->nrows; l++) {
          if ((tmp = mzed_read_elem(A,l,j))) 
            mzed_add_multiple_of_row(A, l, A, row_pos, ff->mul[tmp], j+1);
        }
      }
      row_pos++;
      col_pos = j + 1;
    } else {  
      break;
    }
  }

  for (rci_t i=0; i < row_pos; i++) {
    mzed_col_swap_in_rows(A, i, Q->values[i], i, A->nrows);
  }
  return row_pos;
}

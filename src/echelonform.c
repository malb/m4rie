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

#include "echelonform.h"
#include "travolta.h"
#include "permutation.h"
#include "trsm.h"
#include "ple.h"

rci_t mzed_echelonize(mzed_t *A, int full) {
  if (A->finite_field->degree > A->nrows/2) {
    return mzed_echelonize_naive(A, full);
  } else {
    return mzed_echelonize_travolta(A, full);
  }
}

rci_t mzd_slice_echelonize_ple(mzd_slice_t *A, int full) {
  mzp_t *P = mzp_init(A->nrows);
  mzp_t *Q = mzp_init(A->ncols);
  rci_t r;

  if(full) {
    r = mzd_slice_pluq(A, P, Q);

    mzd_slice_t *U = mzd_slice_init_window(A, 0, 0, r, r);
    mzd_slice_t *B = mzd_slice_init_window(A, 0, r, r, A->ncols);


    if(r!=A->ncols) 
      mzd_slice_trsm_upper_left(U, B);
    if(r) 
      mzd_slice_set_ui(U, 0);
    for(rci_t i = 0; i < r; ++i)
      mzd_slice_write_elem(A, i, i, 1);

    mzd_slice_free_window(U);
    mzd_slice_free_window(B);

    if(r) {
      mzd_slice_t *A0 = mzd_slice_init_window(A, 0, 0, r, A->ncols);
      mzd_slice_apply_p_right(A0, Q);
      mzd_slice_free_window(A0);
    } else {
      mzd_slice_apply_p_right(A, Q);
    }

  } else {
    r = mzd_slice_ple(A, P, Q);

    for(rci_t i = 0; i < r; ++i) {
      for(int e=0; e < A->depth; e++) {
        for(rci_t j = 0; j <= i; j++) {
          int const length = MIN(m4ri_radix, i - j + 1);
          mzd_clear_bits(A->x[e], i, j, length);
        }
      }
      mzd_slice_write_elem(A, i, Q->values[i], 1);
    }
  }

  if(r != A->nrows) {
    mzd_slice_t *R = mzd_slice_init_window(A, r, 0, A->nrows, A->ncols);
    mzd_slice_set_ui(R, 0);
    mzd_slice_free_window(R);
  }
  
  mzp_free(P);
  mzp_free(Q);

  return r;
}



#include "trsm.h"

void mzed_trsm_lower_left_naive(const mzed_t *L, mzed_t *B) {
  assert(L->finite_field == B->finite_field);
  assert(L->nrows == L->ncols);
  assert(B->nrows == L->ncols);

  gf2e *ff = L->finite_field;
  for(rci_t i=0; i<B->nrows; i++) {
    for(rci_t k=0; k<i; k++) {
      mzed_add_multiple_of_row(B, i, B, k, ff->mul[mzed_read_elem(L, i, k)], 0);
    }
    mzed_rescale_row(B, i, 0, ff->mul[ff->inv[mzed_read_elem(L, i, i)]]);
  }
}



void mzed_trsm_lower_left(const mzed_t *L, mzed_t *B) {
  assert(L->finite_field == B->finite_field);
  assert(L->nrows == L->ncols);
  assert(B->nrows == L->ncols);

  const rci_t c = MZED_TRSM_CUTOFF;

  if (L->nrows < c) {
    mzed_trsm_lower_left_naive(L,B);
    return;
  }

  /**
  \verbatim  
  |\           ______
  | \         |      |
  |  \        |  B0  |
  |L00\       |      |
  |____\      |______|
  |    |\     |      |
  |    | \    |      |
  |    |  \   |  B1  |
  |L10 |L11\  |      |
  |____|____\ |______|
  \endverbatim 
  * \li L00 L10 B0 and B1 are possibly located at uneven locations.
  * \li Their column dimension is lower than 64.
  * \li The first column of L01, L11, B1 are aligned to words.
  */

  mzed_t *B0  = mzed_init_window(B, 0, 0, c, B->ncols);
  mzed_t *B1  = mzed_init_window(B, c, 0, B->nrows, B->ncols);
  const mzed_t *L00 = (const mzed_t*)mzed_init_window((mzed_t*)L, 0, 0, c, c);
  const mzed_t *L10 = (const mzed_t*)mzed_init_window((mzed_t*)L, c, 0, B->nrows, c);
  const mzed_t *L11 = (const mzed_t*)mzed_init_window((mzed_t*)L, c, c, B->nrows, B->nrows);
    
  mzed_trsm_lower_left_naive(L00, B0);
  mzed_addmul(B1, L10, B0);
  mzed_trsm_lower_left(L11, B1);
    
  mzed_free_window(B0);
  mzed_free_window(B1);
    
  mzed_free_window((mzed_t*)L00);
  mzed_free_window((mzed_t*)L10);
  mzed_free_window((mzed_t*)L11);
}


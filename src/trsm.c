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

void mzd_slice_trsm_lower_left_naive(const mzd_slice_t *L, mzd_slice_t *B) {
  assert(L->finite_field == B->finite_field);
  assert(L->nrows == L->ncols);
  assert(B->nrows == L->ncols);

  const mzed_t *Le = mzed_cling(NULL, L);
  mzed_t *Be = mzed_cling(NULL, B);
  mzed_trsm_lower_left_naive(Le, Be);

  mzed_slice(B, Be);
  mzed_free((mzed_t*)Le);
  mzed_free(Be);
}


#include "mzed_t_intro.inl"
#include "trsm.inl"
#include "mzed_t_outro.inl"


#include "mzd_slice_t_intro.inl"
#include "trsm.inl"
#include "mzd_slice_t_outro.inl"


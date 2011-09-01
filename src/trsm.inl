/**
 * \brief inline template for TRSM routines 
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 *
 * \note We want to keep this library in C, hence we cannot use of C++
 * templates.
 */

void matrix_trsm_lower_left(const matrix_t *L, matrix_t *B) {
  assert((L->finite_field == B->finite_field) && (L->nrows == L->ncols) && (B->nrows == L->ncols));

  const rci_t c = MZED_TRSM_CUTOFF;

  if (L->nrows < c)
    return matrix_trsm_lower_left_naive(L,B);

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
  */

  matrix_t *B0  = matrix_init_window(B, 0, 0, c, B->ncols);
  matrix_t *B1  = matrix_init_window(B, c, 0, B->nrows, B->ncols);
  const matrix_t *L00 = (const matrix_t*)matrix_init_window((matrix_t*)L, 0, 0, c, c);
  const matrix_t *L10 = (const matrix_t*)matrix_init_window((matrix_t*)L, c, 0, B->nrows, c);
  const matrix_t *L11 = (const matrix_t*)matrix_init_window((matrix_t*)L, c, c, B->nrows, B->nrows);
    
  matrix_trsm_lower_left_naive(L00, B0);
  matrix_addmul(B1, L10, B0);
  matrix_trsm_lower_left(L11, B1);
    
  matrix_free_window(B0);
  matrix_free_window(B1);
  matrix_free_window((matrix_t*)L00);
  matrix_free_window((matrix_t*)L10);
  matrix_free_window((matrix_t*)L11);
}


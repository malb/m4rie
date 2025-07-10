/**
 * \file mzed.h
 *
 * \brief Dense matrices over \GF2E represented as packed matrices.
 *
 * This file implements the data type mzed_t. That is, matrices over
 * \GF2E in row major representation.

 * For example, let \f$ a = \sum a_i x_i / <f>\f$ and \f$b = \sum b_i x_i / <f>\f$
 * be elements in \f$\mathbb{F}_{2^6}\f$ with minimal polynomial \f$f\f$. Then, the
 * \f$ 1 \times 2\f$ matrix [b a] would be stored as
\verbatim
 [...| 0 0 b5 b4 b3 b2 b1 b0 | 0 0 a5 a4 a3 a2 a1 a0]
\endverbatim
 *
 * Internally M4RI matrices are used to store bits with allows to
 * re-use existing M4RI methods (such as mzd_add) when implementing
 * functions for mzed_t.
 *
 * This data type is preferable when Newton-John tables ought be used
 * or when the matrix is small (\f$ m \times n \times e < L2\f$).
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */


#ifndef M4RIE_MZED_H
#define M4RIE_MZED_H

/******************************************************************************
*
*            M4RIE: Linear Algebra over GF(2^e)
*
*    Copyright (C) 2010,2011 Martin Albrecht <martinralbrecht@googlemail.com>
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
#include <m4rie/gf2e.h>
#include <m4rie/m4ri_functions.h>

/**
 * \brief Dense matrices over \GF2E represented as packed matrices.
 *
 * \ingroup Definitions
 */

typedef struct {
  mzd_t *x; /**< \f$m \times n\f$ matrices over \GF2E are represented as \f$m \times (en)\f$ matrices over \GF2. */
  const gf2e *finite_field; /**< A finite field \GF2E. */
  rci_t nrows; /**< Number of rows. */
  rci_t ncols; /**< Number of columns. */
  wi_t w;   /**< The internal width of elements (must divide 64). */
} mzed_t;


/**
 * \brief Create a new matrix of dimension m x n over ff
 *
 * Use mzed_free() to kill it.
 *
 * \param ff Finite field
 * \param m Number of rows
 * \param n Number of columns
 *
 * \ingroup Constructions
 */

mzed_t *mzed_init(const gf2e *ff, const rci_t m, const rci_t n);

/**
 * \brief Free a matrix created with mzed_init().
 *
 * \param A Matrix
 *
 * \ingroup Constructions
 */

void mzed_free(mzed_t *A);


/**
 * \brief Concatenate B to A and write the result to C.
 *
 * That is,
 \verbatim
 [ A ], [ B ] -> [ A  B ] = C
 \endverbatim
 * The inputs are not modified but a new matrix is created.
 *
 * \param C Matrix, may be NULL for automatic creation
 * \param A Matrix
 * \param B Matrix
 *
 * \note This is sometimes called augment.
 *
 * \ingroup Constructions
 */

static inline mzed_t *mzed_concat(mzed_t *C, const mzed_t *A, const mzed_t *B) {
  if(C==NULL)
    C = mzed_init(A->finite_field, A->nrows, A->ncols + B->ncols);
  mzd_concat(C->x, A->x, B->x);
  return C;
}

/**
 * \brief Stack A on top of B and write the result to C.
 *
 * That is,
 \verbatim
 [ A ], [ B ] -> [ A ] = C
                 [ B ]
 \endverbatim
 * The inputs are not modified but a new matrix is created.
 *
 * \param C Matrix, may be NULL for automatic creation
 * \param A Matrix
 * \param B Matrix
 *
 * \ingroup Constructions
 */

static inline mzed_t *mzed_stack(mzed_t *C, const mzed_t *A, const mzed_t *B) {
  if(C==NULL)
    C = mzed_init(A->finite_field, A->nrows + B->nrows, A->ncols);
  mzd_stack(C->x, A->x, B->x);
  return C;
}


/**
 * \brief Copy a submatrix.
 *
 * Note that the upper bounds are not included.
 *
 * \param S Preallocated space for submatrix, may be NULL for automatic creation.
 * \param M Matrix
 * \param lowr start rows
 * \param lowc start column
 * \param highr stop row (this row is \em not included)
 * \param highc stop column (this column is \em not included)
 *
 * \ingroup Constructions
 */
static inline mzed_t *mzed_submatrix(mzed_t *S, const mzed_t *M, const rci_t lowr, const rci_t lowc, const rci_t highr, const rci_t highc) {
  if(S==NULL)
    S = mzed_init(M->finite_field, highr - lowr, highc - lowc);

  mzd_submatrix(S->x, M->x, lowr, lowc*M->w, highr, highc*M->w);
  return S;
}

/**
 * \brief Create a window/view into the matrix A.
 *
 * A matrix window for A is a meta structure on the matrix A. It is
 * setup to point into the matrix so M \em must \em not be freed while
 * the matrix window is used.
 *
 * This function puts the restriction on the provided parameters that
 * all parameters must be within range for A which is not currently
 * enforced.
 *
 * Use mzed_free_window() to free the window.
 *
 * \param A Matrix
 * \param lowr Starting row (inclusive)
 * \param lowc Starting column (inclusive)
 * \param highr End row (exclusive)
 * \param highc End column (exclusive)
 *
 * \ingroup Constructions
 */

static inline mzed_t *mzed_init_window(const mzed_t *A, const rci_t lowr, const rci_t lowc, const rci_t highr, const rci_t highc) {
  mzed_t *B = (mzed_t *)m4ri_mm_malloc(sizeof(mzed_t));
  B->finite_field = A->finite_field;
  B->w = gf2e_degree_to_w(A->finite_field);
  B->nrows = highr - lowr;
  B->ncols = highc - lowc;
  B->x = mzd_init_window(A->x, lowr, B->w*lowc, highr, B->w*highc);
  return B;
}

/**
 * \brief Free a matrix window created with mzed_init_window().
 *
 * \param A Matrix
 *
 * \ingroup Constructions
 */

static inline void mzed_free_window(mzed_t *A) {
  mzd_free_window(A->x);
  m4ri_mm_free(A);
}

/**
 * \brief \f$ C = A+B \f$.
 *
 * C is also returned. If C is NULL then a new matrix is created which
 * must be freed by mzed_free().
 *
 * \param C Preallocated sum matrix, may be NULL for automatic creation.
 * \param A Matrix
 * \param B Matrix
 *
 * \ingroup Addition
 */

mzed_t *mzed_add(mzed_t *C, const mzed_t *A, const mzed_t *B);

/**
 * \brief \f$ C = A+B \f$.
 *
 * \param C Preallocated sum matrix, may be NULL for automatic creation.
 * \param A Matrix
 * \param B Matrix
 *
 * \ingroup Addition
 */

mzed_t *_mzed_add(mzed_t *C, const mzed_t *A, const mzed_t *B);

/**
 * \brief \f$ C = A+B \f$.
 *
 * \param C Preallocated difference matrix, may be NULL for automatic creation.
 * \param A Matrix
 * \param B Matrix
 *
 * \ingroup Addition
 */

#define mzed_sub mzed_add

/**
 * \brief \f$ C = A+B \f$.
 *
 * \param C Preallocated difference matrix, may be NULL for automatic creation.
 * \param A Matrix
 * \param B Matrix
 *
 * \ingroup Addition
 */

#define _mzed_sub _mzed_add

/**
 * \brief \f$ C = A \cdot B \f$.
 *
 * \param C Preallocated return matrix, may be NULL for automatic creation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \ingroup Multiplication
 */

mzed_t *mzed_mul(mzed_t *C, const mzed_t *A, const mzed_t *B);

/**
 * \brief \f$ C = C + A \cdot B \f$.
 *
 * \param C Preallocated product matrix, may be NULL for automatic creation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \ingroup Multiplication
 */

mzed_t *mzed_addmul(mzed_t *C, const mzed_t *A, const mzed_t *B);

/**
 * \brief \f$ C = A \cdot B \f$.
 *
 * \param C Preallocated product matrix.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \ingroup Multiplication
 */

mzed_t *_mzed_mul(mzed_t *C, const mzed_t *A, const mzed_t *B);

/**
 * \brief \f$ C = C + A \cdot B \f$.
 *
 * \param C Preallocated product matrix.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \ingroup Multiplication
 */

mzed_t *_mzed_addmul(mzed_t *C, const mzed_t *A, const mzed_t *B);


/**
 * \brief \f$ C = C + A \cdot B \f$ using naive cubic multiplication.
 *
 * \param C Preallocated product matrix.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \note There is no reason to call this function except for checking
 * the correctness of other algorithms. It is very slow.
 *
 * \ingroup Multiplication
 */

mzed_t *mzed_addmul_naive(mzed_t *C, const mzed_t *A, const mzed_t *B);

/**
 * \brief \f$ C = A \cdot B \f$ using naive cubic multiplication.
 *
 * \param C Preallocated product matrix, may be NULL for automatic creation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \note There is no reason to call this function except for checking
 * the correctness of other algorithms. It is very slow.
 *
 * \ingroup Multiplication
 */

mzed_t *mzed_mul_naive(mzed_t *C, const mzed_t *A, const mzed_t *B);

/**
 * \brief \f$ C = C + A \cdot B \f$ using naive cubic multiplication.
 *
 * \param C Preallocated product matrix.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \ingroup Multiplication
 */

mzed_t *_mzed_mul_naive(mzed_t *C, const mzed_t *A, const mzed_t *B);

/**
 * \brief \f$ C = a \cdot B \f$.
 *
 * \param C Preallocated product matrix or NULL.
 * \param a finite field element.
 * \param B Input matrix B.
 *
 * \ingroup Multiplication
 */

mzed_t *mzed_mul_scalar(mzed_t *C, const word a, const mzed_t *B);

/**
 * Check whether C, A and B match in sizes and fields for
 * multiplication
 *
 * \param C Output matrix, if NULL a new matrix is created.
 * \param A Input matrix.
 * \param B Input matrix.
 * \param clear Write zeros to C or not.
 */

mzed_t *_mzed_mul_init(mzed_t *C, const mzed_t *A, const mzed_t *B, int clear);

/**
 * \brief Fill matrix A with random elements.
 *
 * \param A Matrix
 *
 * \todo Allow the user to provide a RNG callback.
 *
 * \ingroup Assignment
 */

void mzed_randomize(mzed_t *A);

/**
 * \brief Copy matrix A to B.
 *
 * \param B May be NULL for automatic creation.
 * \param A Source matrix.
 *
 * \ingroup Assignment
 */

mzed_t *mzed_copy(mzed_t *B, const mzed_t *A);

/**
 * \brief Transpose a matrix.
 *
 * \param DST Preallocated return matrix, may be NULL for automatic creation.
 * \param A Matrix
 */
mzed_t *mzed_transpose(mzed_t *DST, const mzed_t *A);

/**
 * \brief Return diagonal matrix with value on the diagonal.
 *
 * If the matrix is not square then the largest possible square
 * submatrix is used.
 *
 * \param A Matrix
 * \param value Finite Field element
 *
 * \ingroup Assignment
 */

void mzed_set_ui(mzed_t *A, word value);


/**
 * \brief Get the element at position (row,col) from the matrix A.
 *
 * \param A Source matrix.
 * \param row Starting row.
 * \param col Starting column.
 *
 * \ingroup Assignment
 */

static inline word mzed_read_elem(const mzed_t *A, const rci_t row, const rci_t col) {
  return __mzd_read_bits(A->x, row, A->w*col, A->w);
}

/**
 * \brief At the element elem to the element at position (row,col) in the matrix A.
 *
 * \param A Target matrix.
 * \param row Starting row.
 * \param col Starting column.
 * \param elem finite field element.
 *
 * \ingroup Assignment
 */

static inline void mzed_add_elem(mzed_t *A, const rci_t row, const rci_t col, const word elem) {
  __mzd_xor_bits(A->x, row, A->w*col, A->w, elem);
}

/**
 * \brief Write the element elem to the position (row,col) in the matrix A.
 *
 * \param A Target matrix.
 * \param row Starting row.
 * \param col Starting column.
 * \param elem finite field element.
 *
 * \ingroup Assignment
 */

static inline void mzed_write_elem(mzed_t *A, const rci_t row, const rci_t col, const word elem) {
  __mzd_clear_bits(A->x, row, A->w*col, A->w);
  __mzd_xor_bits(A->x, row, A->w*col, A->w, elem);
}

/**
 * \brief Return -1,0,1 if if A < B, A == B or A > B respectively.
 *
 * \param A Matrix.
 * \param B Matrix.
 *
 * \note This comparison is not well defined mathematically and
 * relatively arbitrary since elements of \GF2E don't have an
 * ordering.
 *
 * \ingroup Comparison
 */

static inline int mzed_cmp(mzed_t *A, mzed_t *B) {
  return mzd_cmp(A->x,B->x);
}


/**
 * \brief Zero test for matrix.
 *
 * \param A Input matrix.
 *
 * \ingroup Comparison
 */
static inline int mzed_is_zero(const mzed_t *A) {
  return mzd_is_zero(A->x);
}

/**
 * A[ar,c] = A[ar,c] + x*B[br,c] for all c >= startcol.
 *
 * \param A Matrix.
 * \param ar Row index in A.
 * \param B Matrix.
 * \param br Row index in B.
 * \param x Finite field element.
 * \param start_col Column index.
 *
 * \ingroup RowOperations
 */

void mzed_add_multiple_of_row(mzed_t *A, rci_t ar, const mzed_t *B, rci_t br, word x, rci_t start_col);

/**
 * A[ar,c] = A[ar,c] + B[br,c] for all c >= startcol.
 *
 * \param A Matrix.
 * \param ar Row index in A.
 * \param B Matrix.
 * \param br Row index in B.
 * \param start_col Column index.
 *
 * \ingroup RowOperations
 */

static inline void mzed_add_row(mzed_t *A, rci_t ar, const mzed_t *B, rci_t br, rci_t start_col) {
  assert(A->ncols == B->ncols && A->finite_field == B->finite_field);
  assert(start_col < A->ncols);

  const rci_t start = A->w*start_col;
  const wi_t startblock = start/m4ri_radix;
  const word bitmask_begin = __M4RI_RIGHT_BITMASK(m4ri_radix - (start%m4ri_radix));
  const word bitmask_end = A->x->high_bitmask;

  word *_a = mzd_row(A->x, ar);
  const word *_b = mzd_row(B->x, br);
  wi_t j;

  if (A->x->width - startblock > 1) {
    _a[startblock] ^= _b[startblock] & bitmask_begin;
    for(j=startblock+1; j<A->x->width-1; j++)
      _a[j] ^= _b[j];
    _a[j] ^= _b[j] & bitmask_end;
  } else {
    _a[startblock] ^= _b[startblock] & (bitmask_begin & bitmask_end);
  }
}

/**
 * \brief Rescale the row r in A by X starting c.
 *
 * \param A Matrix
 * \param r Row index.
 * \param start_col Column index.
 * \param x Multiplier
 *
 * \ingroup RowOperations
 */

void mzed_rescale_row(mzed_t *A, rci_t r, rci_t start_col, const word x);

/**
 * \brief Swap the two rows rowa and rowb.
 *
 * \param M Matrix
 * \param rowa Row index.
 * \param rowb Row index.
 *
 * \ingroup RowOperations
 */

static inline void mzed_row_swap(mzed_t *M, const rci_t rowa, const rci_t rowb) {
  mzd_row_swap(M->x, rowa, rowb);
}

/**
 * \brief copy row j from A to row i from B.
 *
 * The the number of columns of A must be less than or equal to the number of columns of B.
 *
 * \param B Target matrix.
 * \param i Target row index.
 * \param A Source matrix.
 * \param j Source row index.
 *
 * \ingroup RowOperations
 */

static inline void mzed_copy_row(mzed_t* B, rci_t i, const mzed_t* A, rci_t j) {
  mzd_copy_row(B->x, i, A->x, j);
}

/**
 * \brief Swap the two columns cola and colb.
 *
 * \param M Matrix.
 * \param cola Column index.
 * \param colb Column index.
 *
 * \ingroup RowOperations
 */

static inline void mzed_col_swap(mzed_t *M, const rci_t cola, const rci_t colb) {
  for(rci_t i=0; i<M->w; i++)
    mzd_col_swap(M->x,M->w*cola+i, M->w*colb+i);
}

/**
 * \brief Swap the two columns cola and colb but only between start_row and stop_row.
 *
 * \param A Matrix.
 * \param cola Column index.
 * \param colb Column index.
 * \param start_row Row index.
 * \param stop_row Row index (exclusive).
 *
 * \ingroup RowOperations
 */

static inline void mzed_col_swap_in_rows(mzed_t *A, const rci_t cola, const rci_t colb, const rci_t start_row, rci_t stop_row) {
  for(int e=0; e < A->finite_field->degree; e++) {
    mzd_col_swap_in_rows(A->x, A->w*cola+e, A->w*colb+e, start_row, stop_row);
  };
}

/**
 * \brief Add the rows sourcerow and destrow and stores the total in
 * the row destrow.
 *
 * \param M Matrix
 * \param sourcerow Index of source row
 * \param destrow Index of target row
 *
 * \note this can be done much faster with mzed_combine.
 *
 * \ingroup RowOperations
 */

static inline void mzed_row_add(mzed_t *M, const rci_t sourcerow, const rci_t destrow) {
  mzd_row_add(M->x, sourcerow, destrow);
}

/**
 * \brief Return the first row with all zero entries.
 *
 * If no such row can be found returns nrows.
 *
 * \param A Matrix
 *
 * \ingroup RowOperations
 */

static inline rci_t mzed_first_zero_row(mzed_t *A) {
  return mzd_first_zero_row(A->x);
}


/**
 * \brief Gaussian elimination.
 *
 * Perform Gaussian elimination on the matrix A.  If full=0, then it
 * will do triangular style elimination, and if full=1, it will do
 * Gauss-Jordan style, or full elimination.
 *
 * \param A Matrix
 * \param full Gauss-Jordan style or upper unit-triangular form only.
 *
 * \ingroup Echelon
 */

rci_t mzed_echelonize_naive(mzed_t *A, int full);

/**
 * \brief Print a matrix to stdout.
 *
 * \param M Matrix
 *
 * \ingroup StringConversions
 */

void mzed_print(const mzed_t *M);

#endif //M4RIE_MATRIX_H

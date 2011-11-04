/**
 * \file bitslice.h
 *
 * \brief Matrices using a bitsliced representation.
 *
 * Matrices over \GF2E can be represented as polynomials with matrix
 * coefficients where the matrices are in \GF2.
 *
 * In this file, matrices over \GF2E are implemented as \e slices of
 * matrices over \GF2 where each slice holds the coefficients of one
 * degree when viewing elements of \GF2E as polynomials over \GF2.
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

#ifndef M4RIE_BITSLICE_H
#define M4RIE_BITSLICE_H

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
#include "gf2e_matrix.h"

/**
 * Degree up to which Karatsuba multiplication and slicing/cling is
 * implemented.
 */

#define __M4RIE_MAX_KARATSUBA_DEGREE 8

/**
 * \brief Dense matrices over \GF2E represented as slices of matrices over \GF2.
 *
 * This is one of two fundamental data types of this library, the
 * other being mzed_t. For large matrices (\f$m \times n \times e > L2\f$)
 * it is advisable to use this data type because multiplication
 * is faster in this representation. Hence, compared to mzed_t one
 * saves the time to convert betwen representations and - more
 * importantly - memory.
 *
 * \ingroup Definitions
 */

typedef struct {
  /**
   * A->x[e][i,j] is the \e-th bit of the entry A[i,j].
   */
  mzd_t *x[16];  // We only support 10 but 16 migh help with alignment.

  /**
   * A finite field \GF2E.
   */

  const gf2e *finite_field;

  /**
   * Number of rows.
   */

  rci_t nrows;

  /**
   * Number of columns.
   */

  rci_t ncols;

  /**
   * Number of slices
   *
   * \note This value may be greater than finite_field->degree in some
   * situations
   */

  int depth;

} mzd_slice_t;


/**
 * \brief Create a new matrix of dimension \f$ m \times n\f$ over ff
 *
 * Use mzd_slice_free to free it.
 *
 * \param ff Finite field
 * \param m Number of rows
 * \param n Number of columns
 *
 * \ingroup Constructions
 */

static inline mzd_slice_t *mzd_slice_init(const gf2e *ff, const rci_t m, const rci_t n) {

  /**
   * \todo avoid all these malloc() calls and call it once only
   */
  mzd_slice_t *A;

#if __M4RI_USE_MM_MALLOC
  A = (mzd_slice_t*)_mm_malloc(sizeof(mzd_slice_t), 64);
#elif __M4RI_USE_POSIX_MEMALIGN
  int error = posix_memalign(&A, 64, sizeof(mzd_slice_t));
  if (error) A = NULL;
#else
  A = (mzd_slice_t*)malloc(sizeof(mzd_slice_t));
#endif

  if(__M4RI_UNLIKELY(A == NULL))
    m4ri_die("m4ri_slice_init: could not allocate memory.\n");

  A->finite_field = ff;
  A->nrows = m;
  A->ncols = n;
  A->depth = ff->degree;

  for(int i=0; i<A->depth; i++)
    A->x[i] = mzd_init(m,n);
  return A;
}

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

void mzd_slice_set_ui(mzd_slice_t *A, word value);

/**
 * \brief Extend or truncate the depth of A to depth new_depth.
 *
 * We may think of mzd_slice_t as polynomials over matrices over
 * GF(2). This function then truncates/extends these polynomials to
 * degree new_depth-1. In case of extension, all newly created
 * coefficients are zero, hence the mathematical content of A is not
 * changed. In case of truncation higher degree terms are simply
 * deleted and A's mathematical content modified.
 *
 * \param A Matrix, modifed in place
 * \param new_depth Integer >= A->finite_field->degree
 */

static inline mzd_slice_t *_mzd_slice_adapt_depth(mzd_slice_t *A, const int new_depth) {
  assert(A->finite_field->degree <= new_depth);

  if (new_depth < A->depth) { 
    for(int i=new_depth; i<A->depth; i++) {
      mzd_free(A->x[i]);
      A->x[i] = NULL;
    }
  } else {
    for(int i=A->depth; i<new_depth; i++) {
      A->x[i] = mzd_init(A->nrows,A->ncols);
    }
  }
  A->depth = new_depth;
  return A;
}


/**
 * \brief Free a matrix created with mzd_slice_init.
 * 
 * \param A Matrix
 *
 * \ingroup Constructions
 */

static inline void mzd_slice_free(mzd_slice_t *A) {
  for(int i=0; i<A->depth; i++)
   mzd_free(A->x[i]);
#if __M4RI_USE_MM_MALLOC
  _mm_free(A);
#else
  free(A);
#endif
}

/**
 * \brief Concatenate B to A and write the result to C.
 * 
 * That is,
 *
 \verbatim
 [ A ], [ B ] -> [ A  B ] = C
 \endverbatim
 *
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

static inline mzd_slice_t *mzd_slice_concat(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B) {
  if(C == NULL)
    C = mzd_slice_init(A->finite_field, A->nrows, A->ncols + B->ncols);

  for(int i=0; i<A->depth; i++) {
    mzd_concat(C->x[i], A->x[i], B->x[i]);
  }
  return C;
}

/**
 * \brief Stack A on top of B and write the result to C.
 *
 * That is, 
 *
 \verbatim
 [ A ], [ B ] -> [ A ] = C
                 [ B ]
 \endverbatim
 *
 * The inputs are not modified but a new matrix is created.
 *
 * \param C Matrix, may be NULL for automatic creation
 * \param A Matrix
 * \param B Matrix
 *
 * \ingroup Constructions
 */

static inline mzd_slice_t *mzd_slice_stack(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B) {
  if(C == NULL)
    C = mzd_slice_init(A->finite_field, A->nrows + B->nrows, A->ncols);

  for(int i=0; i<A->depth; i++) {
    mzd_stack(C->x[i], A->x[i], B->x[i]);
  }
  return C;
}

/**
 * \brief Copy a submatrix.
 *
 * Note that the upper bounds are not included.
 *
 * \param S Preallocated space for submatrix, may be NULL for automatic creation.
 * \param A Matrix
 * \param lowr start rows
 * \param lowc start column
 * \param highr stop row (this row is \em not included)
 * \param highc stop column (this column is \em not included)
 *
 * \ingroup Constructions
 */

static inline mzd_slice_t *mzd_slice_submatrix(mzd_slice_t *S, const mzd_slice_t *A,
                                               const size_t lowr, const size_t lowc, const size_t highr, const size_t highc) {
  if(S==NULL)
    S = mzd_slice_init(A->finite_field, highr - lowr, highc - lowc);

  for(int i=0; i<A->depth; i++) {
    mzd_submatrix(S->x[i], A->x[i], lowr, lowc, highr, highc);
  }
  return S;
}

/**
 * \brief Create a window/view into the matrix M.
 *
 * A matrix window for M is a meta structure on the matrix M. It is
 * setup to point into the matrix so M \em must \em not be freed while the
 * matrix window is used.
 *
 * This function puts the restriction on the provided parameters that
 * all parameters must be within range for M which is not currently
 * enforced.
 *
 * Use mzd_slice_free_window to free the window.
 *
 * \param A Matrix
 * \param lowr Starting row (inclusive)
 * \param lowc Starting column (inclusive)
 * \param highr End row (exclusive)
 * \param highc End column (exclusive)
 *
 * \ingroup Constructions
 */

static inline mzd_slice_t *mzd_slice_init_window(const mzd_slice_t *A,
                                                 const size_t lowr, const size_t lowc,
                                                 const size_t highr, const size_t highc) {
  mzd_slice_t *B = (mzd_slice_t *)m4ri_mm_malloc(sizeof(mzd_slice_t));
  B->finite_field = A->finite_field;
  B->depth = A->depth;
  B->nrows = highr - lowr;
  B->ncols = highc - lowc;
  for(int i=0; i<A->depth; i++) {
    B->x[i] = mzd_init_window(A->x[i], lowr, lowc, highr, highc);
  }
  return B;
}

/**
 * \brief Free a matrix window created with mzd_slice_init_window.
 * 
 * \param A Matrix
 *
 * \ingroup Constructions
 */

static inline void mzd_slice_free_window(mzd_slice_t *A) {
  for(int i=0; i<A->depth; i++) {
    mzd_free_window(A->x[i]);
  }
  m4ri_mm_free(A);
}

/**
 * \brief Same as mzd_add but without any checks on the input.
 *
 * \param C Preallocated sum matrix, may be NULL for automatic creation.
 * \param A Matrix
 * \param B Matrix
 *
 * \ingroup Addition
 */

static inline mzd_slice_t *_mzd_slice_add(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B) {
  for(int i=0; i<A->depth; i++) 
    _mzd_add(C->x[i], A->x[i], B->x[i]);
  return C;
}

/**
 * \brief Set C = A+B.
 *
 * C is also returned. If C is NULL then a new matrix is created which
 * must be freed by mzd_slice_free.
 *
 * \param C Preallocated sum matrix, may be NULL for automatic creation.
 * \param A Matrix
 * \param B Matrix
 *
 * \ingroup Addition
 */

static inline mzd_slice_t *mzd_slice_add(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B) {
  if ( (A->finite_field != B->finite_field) | (A->nrows != B->nrows) | (A->ncols != B->ncols) ) 
    m4ri_die("mzd_slice_add: input matrices A (%d x %d) and B (%d x %d) do not match.\n",A->nrows,A->ncols, B->nrows,B->ncols);

  if(C == NULL) 
    mzd_slice_init(A->finite_field, A->nrows, A->ncols);
  else if ( (A->finite_field != C->finite_field) | (A->nrows != C->nrows) | (A->ncols != C->ncols) )
    m4ri_die("mzd_slice_add: input matrix A (%d x %d) and output matrix (%d x %d) do not match.\n",A->nrows,A->ncols, C->nrows, C->ncols);

  return _mzd_slice_add(C,A,B);
}

/**
 * \brief Set C = A+B.
 *
 * C is also returned. If C is NULL then a new matrix is created which
 * must be freed by mzd_slice_free.
 *
 * \param C Preallocated sum matrix, may be NULL for automatic creation.
 * \param A Matrix
 * \param B Matrix
 *
 * \ingroup Addition
 */

#define mzd_slice_sub mzd_slice_add

/**
 * \brief Same as mzd_add but without any checks on the input.
 *
 * \param C Preallocated sum matrix, may be NULL for automatic creation.
 * \param A Matrix
 * \param B Matrix
 *
 * \ingroup Addition
 */

#define _mzd_slice_sub _mzd_slice_add

/**
 * \brief Compute C = A*B over GF(2^e) using naive polynomials multiplication.
 *
 * \param C Preallocated return matrix, may be NULL for automatic creation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 */

mzd_slice_t *_mzd_slice_mul_naive(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B);

/**
 * \brief Compute C = A*B over GF(2^2) using 3 multiplications over GF(2) and 2 temporary GF(2) matrices..
 *
 * \param C Preallocated return matrix, may be NULL for automatic creation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa _mzd_slice_mul_karatsuba
 */
mzd_slice_t *_mzd_slice_mul_karatsuba2(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B);

/**
 * \brief Compute C = A*B over GF(2^3) using 6 multiplications over GF(2) and 3 temporary GF(2) matrices..
 *
 * The formula was taken from Peter L. Montgomery. "Five, Six, and
 * Seven-Term Karatsuba-Like Formulae" in IEEE TRANSACTIONS ON
 * COMPUTERS, VOL. 54, NO. 3, MARCH 2005/
 *
 * \param C Preallocated return matrix, may be NULL for automatic creation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa _mzd_slice_mul_karatsuba
 */

mzd_slice_t *_mzd_slice_mul_karatsuba3(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B);

/**
 * \brief Compute C = A*B over GF(2^4) using 9 multiplications over GF(2) and 3 temporary GF(2) matrices..
 *
 * \param C Preallocated return matrix, may be NULL for automatic creation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa _mzd_slice_mul_karatsuba
 */

mzd_slice_t *_mzd_slice_mul_karatsuba4(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B);

/**
 * \brief Compute C = A*B over GF(2^5) using 13 multiplications over GF(2) and 3 temporary GF(2) matrices..
 *
 * The formula was taken from Peter L. Montgomery. "Five, Six, and
 * Seven-Term Karatsuba-Like Formulae" in IEEE TRANSACTIONS ON
 * COMPUTERS, VOL. 54, NO. 3, MARCH 2005/
 *
 * \param C Preallocated return matrix, may be NULL for automatic creation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa _mzd_slice_mul_karatsuba
 */

mzd_slice_t *_mzd_slice_mul_karatsuba5(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B);

/**
 * \brief Compute C = A*B over GF(2^6) using 17 multiplications over GF(2) and 3 temporary GF(2) matrices.
 *
 * The formula was taken from Peter L. Montgomery. "Five, Six, and
 * Seven-Term Karatsuba-Like Formulae" in IEEE TRANSACTIONS ON
 * COMPUTERS, VOL. 54, NO. 3, MARCH 2005/
 *
 * \param C Preallocated return matrix, may be NULL for automatic creation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa _mzd_slice_mul_karatsuba
 */

mzd_slice_t *_mzd_slice_mul_karatsuba6(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B);

/**
 * \brief Compute C = A*B over GF(2^7) using 22 multiplications over GF(2) and 3 temporary GF(2) matrices.
 *
 * The formula was taken from Peter L. Montgomery. "Five, Six, and
 * Seven-Term Karatsuba-Like Formulae" in IEEE TRANSACTIONS ON
 * COMPUTERS, VOL. 54, NO. 3, MARCH 2005/
 *
 * \param C Preallocated return matrix, may be NULL for automatic creation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa _mzd_slice_mul_karatsuba
 */

mzd_slice_t *_mzd_slice_mul_karatsuba7(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B);

/**
 * \brief Compute C = A*B over GF(2^8) using 27 multiplications over GF(2) and 15 temporary GF(2) matrices.
 *
 * \param C Preallocated return matrix, may be NULL for automatic creation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa _mzd_slice_mul_karatsuba
 */


mzd_slice_t *_mzd_slice_mul_karatsuba8(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B);

/**
 * \brief Compute C += A*B using Karatsuba multiplication of polynomials over GF(2).
 *
 * This functionreduces matrix multiplication over GF(2^e) to matrix
 * multiplication over GF(2).
 * 
 * As an example consider GF(2^2). The minimal polynomial is x^2 + x +
 * 1. The matrix A can be represented as A0*x + A1 and the matrix B
 * can be represented as B0*x + B1. Their product C is
 * \f[
 A0*B0*x^2 + (A0*B1 + A1*B0)*x + A1*B1.
 * \f]
 * Reduction modulo x^2 + x + 1 gives
 * \f[
 (A0*B0 + A0*B1 + A1*B0)*x + A1*B1 + A0*B0.
 * \f]
 * This can be re-written as
 * \f[
 ((A0 + A1)*(B0 + B1) + A1*B1)*x + A1*B1 + A0*B0
 * \f]
 * and thus this multiplication costs 3 matrix multiplications over
 * GF(2) and 4 matrix additions over GF(2).
 *
 * This technique was proposed in Tomas J. Boothby and Robert
 * W. Bradshaw; Bitslicing and the Method of Four Russians Over Larger
 * Finite Fields; 2009; http://arxiv.org/abs/0901.1413
 *
 * \param C Preallocated return matrix, may be NULL for automatic creation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa mzed_mul() mzd_slice_mul mzd_slice_mul_karatsuba
 */

static inline mzd_slice_t *_mzd_slice_mul_karatsuba(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B) {
  switch(A->finite_field->degree) {
  case  2: C = _mzd_slice_mul_karatsuba2(C, A, B); break;
  case  3: C = _mzd_slice_mul_karatsuba3(C, A, B); break;
  case  4: C = _mzd_slice_mul_karatsuba4(C, A, B); break;
  case  5: C = _mzd_slice_mul_karatsuba5(C, A, B); break;
  case  6: C = _mzd_slice_mul_karatsuba6(C, A, B); break;
  case  7: C = _mzd_slice_mul_karatsuba7(C, A, B); break;
  case  8: C = _mzd_slice_mul_karatsuba8(C, A, B); break;
  case  9: C = _mzd_slice_mul_naive(C, A, B); break;
  case 10: C = _mzd_slice_mul_naive(C, A, B); break;
  default:
    m4ri_die("_mzd_slice_mul_karatsuba: only implemented for GF(2^e) with e <= 4");
  }
  return C;
}

/**
 * \brief Compute C = A*B using Karatsuba multiplication of polynomials over GF(2).
 *
 * \param C Preallocated return matrix, may be NULL for automatic creation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa _mzd_slice_mul_karatsuba
 */

static inline mzd_slice_t *mzd_slice_mul_karatsuba(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B) {
  if (A->ncols != B->nrows || A->finite_field != B->finite_field) 
    m4ri_die("mzd_slice_mul_karatsuba: rows, columns and fields must match.\n");
  if (C != NULL) {
    if (C->finite_field != A->finite_field || C->nrows != A->nrows || C->ncols != B->ncols) 
      m4ri_die("mzd_slice_mul_karatsuba: rows and columns of returned matrix must match.\n");
    mzd_slice_set_ui(C,0);
  }
  return _mzd_slice_mul_karatsuba(C, A, B);
}

/**
 * \brief Compute C += A*B using Karatsuba multiplication of polynomials over GF(2).
 *
 * \param C Preallocated return matrix.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa _mzd_slice_mul_karatsuba
 */

static inline mzd_slice_t *mzd_slice_addmul_karatsuba(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B) {
  assert(C != NULL);
  if (A->ncols != B->nrows || A->finite_field != B->finite_field)
    m4ri_die("mzd_slice_addmul_karatsuba: rows, columns and fields must match.\n");
  if (C->finite_field != A->finite_field || C->nrows != A->nrows || C->ncols != B->ncols)
    m4ri_die("mzd_slice_addmul_karatsuba: rows and columns of returned matrix must match.\n");
  return _mzd_slice_mul_karatsuba(C, A, B);
}

/**
 * \brief Compute C = A*B.
 *
 * \param C Preallocated return matrix, may be NULL for automatic creation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa _mzd_slice_mul_karatsuba
 */

static inline mzd_slice_t *mzd_slice_mul(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B) {
  return mzd_slice_mul_karatsuba(C,A,B);
}

/**
 * \brief Compute C += A*B.
 *
 * \param C Preallocated return matrix.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa _mzd_slice_mul_karatsuba
 */

static inline mzd_slice_t *mzd_slice_addmul(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B) {
  return mzd_slice_addmul_karatsuba(C,A,B);
}

/**
 * \brief Fill matrix A with random elements.
 *
 * \param A Matrix
 *
 * \todo Allow the user to provide a RNG callback.
 *
 * \ingroup Assignment
 */

static inline void mzd_slice_randomize(mzd_slice_t *A) {
  switch(A->depth) {
  case 10: mzd_randomize(A->x[9]);
  case  9: mzd_randomize(A->x[8]);
  case  8: mzd_randomize(A->x[7]);
  case  7: mzd_randomize(A->x[6]);
  case  6: mzd_randomize(A->x[5]);
  case  5: mzd_randomize(A->x[4]);
  case  4: mzd_randomize(A->x[3]);
  case  3: mzd_randomize(A->x[2]);
  case  2: mzd_randomize(A->x[1]);
  case  1: mzd_randomize(A->x[0]); break;
  default:
    m4ri_die("impossible");
  }
}

/**
 * \brief Copy matrix A to B.
 *
 * \param B May be NULL for automatic creation.
 * \param A Source matrix.
 *
 * \ingroup Assignment
 */

static inline mzd_slice_t *mzd_slice_copy(mzd_slice_t *B, const mzd_slice_t *A) {
  if(B == NULL)
    B = mzd_slice_init(A->finite_field, A->nrows, A->ncols);

  for(int i=0; i<A->depth; i++) {
    mzd_copy(B->x[i],A->x[i]);
  }
  return B;
}

/**
 * \brief Get the element at position (row,col) from the matrix A.
 *
 * \param A Source matrix.
 * \param row Starting row.
 * \param col Starting column.
 *
 * \todo This function is considerably slower than it needs to be.
 *
 * \ingroup Assignment
 */ 

static inline word mzd_slice_read_elem(const mzd_slice_t *A, const rci_t row, const rci_t col) {
  word ret = 0;
  for(int i=0; i<A->depth; i++) {
    ret |= mzd_read_bit(A->x[i], row, col)<<i;
  }
  return ret;
}

/**
 * \brief At the element elem to the element at position (row,col) in the matrix A.
 *
 * \param A Target matrix.
 * \param row Starting row.
 * \param col Starting column.
 * \param elem finite field element.
 *
 * \todo This function is considerably slower than it needs to be.
 *
 * \ingroup Assignment
 */ 

static inline void mzd_slice_add_elem(mzd_slice_t *A, const rci_t row, const rci_t col, word elem) {
  for(int i=0; i<A->depth; i++) {
    __mzd_xor_bits(A->x[i], row, col, 1, elem&1);
    elem=elem>>1;
  }
}

/**
 * \brief Write the element elem to the position (row,col) in the matrix A.
 *
 * \param A Target matrix.
 * \param row Starting row.
 * \param col Starting column.
 * \param elem finite field element.
 *
 * \todo This function is considerably slower than it needs to be.
 *
 * \ingroup Assignment
 */

static inline void mzd_slice_write_elem(mzd_slice_t *A, const rci_t row, const rci_t col, word elem) {
  for(int i=0; i<A->depth; i++) {
    mzd_write_bit(A->x[i], row, col, elem&1);
    elem=elem>>1;
  }
}

/**
 * \brief Return -1,0,1 if if A < B, A == B or A > B respectively.
 *
 * \param A Matrix.
 * \param B Matrix.
 *
 * \note This comparison is not well defined (except for !=0)
 * mathematically and relatively arbitrary since elements of GF(2^k)
 * don't have an ordering.
 *
 * \ingroup Comparison
 */

static inline int mzd_slice_cmp(mzd_slice_t *A, mzd_slice_t *B) {
  int r = 0;
  if ((A->finite_field != B->finite_field) | (A->depth != B->depth) )
    return -1;
  for(int i=0; i<A->depth; i++)
    r |= mzd_cmp(A->x[i],B->x[i]);
  return r;
}

/**
 * \brief Zero test for matrix.
 *
 * \param A Input matrix.
 *
 * \ingroup Comparison
 */

static inline int mzd_slice_is_zero(const mzd_slice_t *A) {
  for(int i=0; i<A->depth; i++) {
    if (!mzd_is_zero(A->x[i]))
      return 0;
  }
  return 1;
}

/**
 * \brief Swap the two rows rowa and rowb.
 *
 * \param A Matrix
 * \param rowa Row index.
 * \param rowb Row index.
 *
 * \ingroup RowOperations
 *
 * \wordoffset
 */

static inline void mzd_slice_row_swap(mzd_slice_t *A, const rci_t rowa, const rci_t rowb) {
  for(int i=0; i<A->depth; i++) {
    mzd_row_swap(A->x[i], rowa, rowb);
  }
}

/**
 * \brief copy row j from A to row i from B.
 *
 * The offsets of A and B must match and the number of columns of A
 * must be less than or equal to the number of columns of B.
 *
 * \param B Target matrix.
 * \param i Target row index.
 * \param A Source matrix.
 * \param j Source row index.
 *
 * \ingroup RowOperations
 */

static inline void mzd_slice_copy_row(mzd_slice_t* B, size_t i, const mzd_slice_t* A, size_t j) {
  for(int ii=0; ii<A->depth; ii++)
    mzd_copy_row(B->x[ii], i, A->x[ii], j);
}

/**
 * \brief Swap the two columns cola and colb.
 *
 * \param A Matrix.
 * \param cola Column index.
 * \param colb Column index.
 *
 * \ingroup RowOperations
 */

static inline void mzd_slice_col_swap(mzd_slice_t *A, const rci_t cola, const rci_t colb) {
  for(int i=0; i<A->depth; i++)
    mzd_col_swap(A->x[i], cola, colb);
}

/**
 * \brief Swap the two columns cola and colb but only between start_row and stop_row.
 *
 * \param A Matrix.
 * \param cola Column index.
 * \param colb Column index.
 * \param start_row Row index.
 * \param stop_row Row index (exclusive).
 */

static inline void mzd_slice_col_swap_in_rows(mzd_slice_t *A, const rci_t cola, const rci_t colb, const rci_t start_row, rci_t stop_row) {
  for(int e=0; e < A->finite_field->degree; e++) {
    mzd_col_swap_in_rows(A->x[e], cola, colb, start_row, stop_row);
  };
}

/**
 * \brief Add the rows sourcerow and destrow and stores the total in
 * the row destrow.
 *
 * \param A Matrix
 * \param sourcerow Index of source row
 * \param destrow Index of target row
 *
 * \note this can be done much faster with mzd_combine.
 *
 * \ingroup RowOperations
 */

static inline void mzd_slice_row_add(mzd_slice_t *A, const rci_t sourcerow, const rci_t destrow) {
  for(int i=0; i<A->depth; i++)
    mzd_row_add(A->x[i], sourcerow, destrow);
}

/**
 * \brief Clear the given row, but only begins at the column coloffset.
 *
 * \param A Matrix
 * \param row Index of row
 * \param coloffset Column offset
 *
 * \ingroup RowOperations
 */

static inline void mzd_slice_row_clear_offset(mzd_slice_t *A, const rci_t row, const rci_t coloffset) {
  for(int i=0; i<A->depth; i++)
    mzd_row_clear_offset(A->x[i], row, coloffset);
}

/**
 * \brief Print a matrix to stdout.
 *
 * \param A Matrix
 *
 * \ingroup StringConversions
 */

void mzd_slice_print(const mzd_slice_t *A);

/**
 * \brief Move the submatrix L of rank r2 starting at column n1 to the left to column r1.
 *
 * \param A Matrix
 * \param r1 Integer < n1
 * \param n1 Integer > r1
 * \param r2 Integer <= A->ncols - n1
 */

static inline void _mzd_slice_compress_l(mzd_slice_t *A, const rci_t r1, const rci_t n1, const rci_t r2) {
  switch(A->finite_field->degree) {
  case 10: _mzd_compress_l(A->x[9], r1, n1, r2);
  case  9: _mzd_compress_l(A->x[8], r1, n1, r2);
  case  8: _mzd_compress_l(A->x[7], r1, n1, r2);
  case  7: _mzd_compress_l(A->x[6], r1, n1, r2);
  case  6: _mzd_compress_l(A->x[5], r1, n1, r2);
  case  5: _mzd_compress_l(A->x[4], r1, n1, r2);
  case  4: _mzd_compress_l(A->x[3], r1, n1, r2);
  case  3: _mzd_compress_l(A->x[2], r1, n1, r2);
  case  2: _mzd_compress_l(A->x[1], r1, n1, r2);
  case  1: _mzd_compress_l(A->x[0], r1, n1, r2); break;
  default:
    m4ri_die("impossible");
  };
}

#endif //M4RIE_BITSLICE_H

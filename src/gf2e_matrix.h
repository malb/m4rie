/**
 * \file gf2e_matrix.h
 * \brief Dense matrices over GF(2^n) represented by M4RI matrices.
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

#ifndef GF2E_MATRIX_H
#define GF2E_MATRIX_H

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

#include <m4ri/m4ri.h>
#include "finite_field.h"
#include "m4ri_functions.h"

/**
 * \brief Dense matrices over GF(2^n). 
 * 
 * The most fundamental data type in this library.
 */
 
typedef struct {

  /**
   * The internal representation of our matrices.
   */
  mzd_t *x;

  /**
   * A finite field GF(2^n).
   */

  gf2e *finite_field;

  /**
   * Number of rows.
   */

  size_t nrows;

  /**
   * Number of columns.
   */

  size_t ncols;

  /**
   * The internal width of elements (must divide 64)
   */

  size_t width;

} mzed_t;


/**
 * \brief Create a new matrix of dimension m x n over ff
 *
 * Use mzed_free to kill it.
 *
 * \param ff Finite field
 * \param m Number of rows
 * \param n Number of columns
 *
 */

mzed_t *mzed_init(gf2e *ff, const size_t m, const size_t n);

/**
 * \brief Free a matrix created with mzed_init.
 * 
 * \param A Matrix
 */

void mzed_free(mzed_t *A);

/**
 * \brief Fill matrix A with random elements.
 *
 * \param A Matrix
 *
 * \todo Allow the user to provide a RNG callback.
 *
 * \wordoffset
 */

void mzed_randomize(mzed_t *A);

/**
 * \brief Set C = A+B.
 *
 * C is also returned. If C is NULL then a new matrix is created which
 * must be freed by mzed_free.
 *
 * \param C Preallocated sum matrix, may be NULL for automatic creation.
 * \param A Matrix
 * \param B Matrix
 */

mzed_t *mzed_add(mzed_t *C, const mzed_t *A, const mzed_t *B);

/**
 * \brief Same as mzed_add but without any checks on the input.
 *
 * \param C Preallocated sum matrix, may be NULL for automatic creation.
 * \param A Matrix
 * \param B Matrix
 *
 * \wordoffset
 */

mzed_t *_mzed_add(mzed_t *C, const mzed_t *A, const mzed_t *B);

/**
 * \brief Same as mzed_add.
 *
 * \param C Preallocated difference matrix, may be NULL for automatic creation.
 * \param A Matrix
 * \param B Matrix
 *
 * \wordoffset
 */

#define mzed_sub mzed_add

/**
 * \brief Same as mzed_sub but without any checks on the input.
 *
 * \param C Preallocated difference matrix, may be NULL for automatic creation.
 * \param A Matrix
 * \param B Matrix
 *
 * \wordoffset
 */

#define _mzed_sub _mzed_add


/**
 * \brief Copy matrix  A to DST.
 *
 * \param DST May be NULL for automatic creation.
 * \param A Source matrix.
 *
 * \wordoffset
 */

mzed_t *mzed_copy(mzed_t *o, const mzed_t *i);

/**
 * \brief Return -1,0,1 if if A < B, A == B or A > B respectively.
 *
 * \param A Matrix.
 * \param B Matrix.
 *
 * \note This comparison is not well defined mathematically and
 * relatively arbitrary since elements of GF(2) don't have an
 * ordering.
 *
 * \wordoffset
 */

static inline int mzed_cmp(mzed_t *l, mzed_t *r) {
  return mzd_cmp(l->x,r->x);
}

/**
 * Get the element at position (row,col) from the matrix A.
 *
 * \param A Source matrix.
 * \param row Starting row.
 * \param col Starting column.
 */ 

static inline int mzed_read_elem(const mzed_t *A, const size_t row, const size_t col) {
  return (int)__mzd_read_bits(A->x, row, A->width*col, A->width);
}

/**
 * At the element elem to the element at position (row,col) in the matrix A.
 *
 * \param A Target matrix.
 * \param row Starting row.
 * \param col Starting column.
 * \param elem finite field element.
 */ 


static inline void mzed_add_elem(mzed_t *a, const size_t row, const size_t col, const int elem) {
  __mzd_xor_bits(a->x, row, a->width*col, a->width, elem);
}

/**
 * Write the element elem to the position (row,col) in the matrix A.
 *
 * \param A Target matrix.
 * \param row Starting row.
 * \param col Starting column.
 * \param elem finite field element.
 */ 

static inline void mzed_write_elem(mzed_t *a, const size_t row, const size_t col, const int elem) {
  __mzd_clear_bits(a->x, row, a->width*col, a->width);
  __mzd_xor_bits(a->x, row, a->width*col, a->width, elem);
}

/**
 *  Compute A[ar,c] = A[ar,c] + X*B[br,c] for all c >= startcol.
 *
 * \param A Matrix.
 * \param ar Row index in A.
 * \param B Matrix.
 * \param br Row index in B.
 * \param X Lookup table for multiplication with some finite field element x.
 * \param start_col Column index.
 */

static inline void mzed_add_multiple_of_row(mzed_t *A, size_t ar, mzed_t *B, size_t br, word *X, size_t start_col) {
  assert(A->ncols == B->ncols && A->finite_field == B->finite_field);
  assert(A->x->offset == 0 && B->x->offset == 0);
  if(A->width == 4) {
    size_t startblock = start_col/RADIX;
    mzd_t *from_x = B->x;
    mzd_t *to_x = A->x;
    word *_f = from_x->rows[br];
    word *_t = to_x->rows[ar];
    size_t j;
    register word __t, __f;
    for(j=startblock; j<to_x->width -1; j++) {
      __f = _f[j], __t = _t[j];
      __t ^= (X[((__f)& 0x000000000000000FULL)])<<0;   __f >>= 4;
      __t ^= (X[((__f)& 0x000000000000000FULL)])<<4;   __f >>= 4;
      __t ^= (X[((__f)& 0x000000000000000FULL)])<<8;   __f >>= 4;
      __t ^= (X[((__f)& 0x000000000000000FULL)])<<12;  __f >>= 4;
      __t ^= (X[((__f)& 0x000000000000000FULL)])<<16;  __f >>= 4;
      __t ^= (X[((__f)& 0x000000000000000FULL)])<<20;  __f >>= 4;
      __t ^= (X[((__f)& 0x000000000000000FULL)])<<24;  __f >>= 4;
      __t ^= (X[((__f)& 0x000000000000000FULL)])<<28;  __f >>= 4;
      __t ^= (X[((__f)& 0x000000000000000FULL)])<<32;  __f >>= 4;
      __t ^= (X[((__f)& 0x000000000000000FULL)])<<36;  __f >>= 4;
      __t ^= (X[((__f)& 0x000000000000000FULL)])<<40;  __f >>= 4;
      __t ^= (X[((__f)& 0x000000000000000FULL)])<<44;  __f >>= 4;
      __t ^= (X[((__f)& 0x000000000000000FULL)])<<48;  __f >>= 4;
      __t ^= (X[((__f)& 0x000000000000000FULL)])<<52;  __f >>= 4;
      __t ^= (X[((__f)& 0x000000000000000FULL)])<<56;  __f >>= 4;
      __t ^= (X[((__f)& 0x000000000000000FULL)])<<60;  __f >>= 4;
      _t[j] = __t;
    }

    switch(to_x->ncols % RADIX) {
    case  0: _t[j] ^= ((word)X[(int)((_f[j] & 0x000000000000000FULL)>> 0)])<< 0;
    case 60: _t[j] ^= ((word)X[(int)((_f[j] & 0x00000000000000F0ULL)>> 4)])<< 4;
    case 56: _t[j] ^= ((word)X[(int)((_f[j] & 0x0000000000000F00ULL)>> 8)])<< 8;
    case 52: _t[j] ^= ((word)X[(int)((_f[j] & 0x000000000000F000ULL)>>12)])<<12;
    case 48: _t[j] ^= ((word)X[(int)((_f[j] & 0x00000000000F0000ULL)>>16)])<<16;
    case 44: _t[j] ^= ((word)X[(int)((_f[j] & 0x0000000000F00000ULL)>>20)])<<20;
    case 40: _t[j] ^= ((word)X[(int)((_f[j] & 0x000000000F000000ULL)>>24)])<<24;
    case 36: _t[j] ^= ((word)X[(int)((_f[j] & 0x00000000F0000000ULL)>>28)])<<28;
    case 32: _t[j] ^= ((word)X[(int)((_f[j] & 0x0000000F00000000ULL)>>32)])<<32;
    case 28: _t[j] ^= ((word)X[(int)((_f[j] & 0x000000F000000000ULL)>>36)])<<36;
    case 24: _t[j] ^= ((word)X[(int)((_f[j] & 0x00000F0000000000ULL)>>40)])<<40;
    case 20: _t[j] ^= ((word)X[(int)((_f[j] & 0x0000F00000000000ULL)>>44)])<<44;
    case 16: _t[j] ^= ((word)X[(int)((_f[j] & 0x000F000000000000ULL)>>48)])<<48;
    case 12: _t[j] ^= ((word)X[(int)((_f[j] & 0x00F0000000000000ULL)>>52)])<<52;
    case  8: _t[j] ^= ((word)X[(int)((_f[j] & 0x0F00000000000000ULL)>>56)])<<56;
    case  4: _t[j] ^= ((word)X[(int)((_f[j] & 0xF000000000000000ULL)>>60)])<<60;
    };

  } else if (A->width == 8) {
    size_t startblock = start_col/RADIX;
    mzd_t *from_x = B->x;
    mzd_t *to_x = A->x;
    word *_f = from_x->rows[br];
    word *_t = to_x->rows[ar];
    size_t j;
    register word __t, __f;

    /* for(j=startblock; j<to_x->width -1; j++) { */
    /*   tmp = _f[j]; */
    /*   _t[j] ^= ((word)X[(int)((tmp & 0xFF00000000000000ULL)>>56)])<<56; */
    /*   _t[j] ^= ((word)X[(int)((tmp & 0x00FF000000000000ULL)>>48)])<<48; */
    /*   _t[j] ^= ((word)X[(int)((tmp & 0x0000FF0000000000ULL)>>40)])<<40; */
    /*   _t[j] ^= ((word)X[(int)((tmp & 0x000000FF00000000ULL)>>32)])<<32; */
    /*   _t[j] ^= ((word)X[(int)((tmp & 0x00000000FF000000ULL)>>24)])<<24; */
    /*   _t[j] ^= ((word)X[(int)((tmp & 0x0000000000FF0000ULL)>>16)])<<16; */
    /*   _t[j] ^= ((word)X[(int)((tmp & 0x000000000000FF00ULL)>> 8)])<< 8; */
    /*   _t[j] ^= ((word)X[(int)((tmp & 0x00000000000000FFULL)>> 0)])<< 0; */
    /* } */

    /**
     * TODO: revert this and benchmark which is faster
     */
    for(j=startblock; j<to_x->width -1; j++) {
      __f = _f[j], __t = _t[j];
      __t ^= (X[((__f)& 0x00000000000000FFULL)])<<0;  __f >>= 8;
      __t ^= (X[((__f)& 0x00000000000000FFULL)])<<8;  __f >>= 8;
      __t ^= (X[((__f)& 0x00000000000000FFULL)])<<16; __f >>= 8;
      __t ^= (X[((__f)& 0x00000000000000FFULL)])<<24; __f >>= 8;
      __t ^= (X[((__f)& 0x00000000000000FFULL)])<<32; __f >>= 8;
      __t ^= (X[((__f)& 0x00000000000000FFULL)])<<40; __f >>= 8;
      __t ^= (X[((__f)& 0x00000000000000FFULL)])<<48; __f >>= 8;
      __t ^= (X[((__f)& 0x00000000000000FFULL)])<<56; __f >>= 8;
      _t[j] = __t;
    }
    
    switch(to_x->ncols % RADIX) {
    case  0: _t[j] ^= ((word)X[(int)((_f[j] & 0x00000000000000FFULL)>> 0)])<< 0;
    case 56: _t[j] ^= ((word)X[(int)((_f[j] & 0x000000000000FF00ULL)>> 8)])<< 8;
    case 48: _t[j] ^= ((word)X[(int)((_f[j] & 0x0000000000FF0000ULL)>>16)])<<16;
    case 40: _t[j] ^= ((word)X[(int)((_f[j] & 0x00000000FF000000ULL)>>24)])<<24;
    case 32: _t[j] ^= ((word)X[(int)((_f[j] & 0x000000FF00000000ULL)>>32)])<<32;
    case 24: _t[j] ^= ((word)X[(int)((_f[j] & 0x0000FF0000000000ULL)>>40)])<<40;
    case 16: _t[j] ^= ((word)X[(int)((_f[j] & 0x00FF000000000000ULL)>>48)])<<48;
    case  8: _t[j] ^= ((word)X[(int)((_f[j] & 0xFF00000000000000ULL)>>56)])<<56;
    };

  } else if (A->width == 16) {
    size_t startblock = start_col/RADIX;
    mzd_t *from_x = B->x;
    mzd_t *to_x = A->x;
    word *_f = from_x->rows[br];
    word *_t = to_x->rows[ar];
    size_t j;
    register word __t, __f;

    for(j=startblock; j+4<to_x->width; j+=4) {
      __f = _f[j], __t = _t[j];
      __t ^= (X[((__f)& 0x000000000000FFFFULL)])<<0;  __f >>= 16;
      __t ^= (X[((__f)& 0x000000000000FFFFULL)])<<16; __f >>= 16;
      __t ^= (X[((__f)& 0x000000000000FFFFULL)])<<32; __f >>= 16;
      __t ^= (X[((__f)& 0x000000000000FFFFULL)])<<48; __f >>= 16;
      _t[j] = __t;

      __f = _f[j+1], __t = _t[j+1];
      __t ^= (X[((__f)& 0x000000000000FFFFULL)])<<0;  __f >>= 16;
      __t ^= (X[((__f)& 0x000000000000FFFFULL)])<<16; __f >>= 16;
      __t ^= (X[((__f)& 0x000000000000FFFFULL)])<<32; __f >>= 16;
      __t ^= (X[((__f)& 0x000000000000FFFFULL)])<<48; __f >>= 16;
      _t[j+1] = __t;


      __f = _f[j+2], __t = _t[j+2];
      __t ^= (X[((__f)& 0x000000000000FFFFULL)])<<0;  __f >>= 16;
      __t ^= (X[((__f)& 0x000000000000FFFFULL)])<<16; __f >>= 16;
      __t ^= (X[((__f)& 0x000000000000FFFFULL)])<<32; __f >>= 16;
      __t ^= (X[((__f)& 0x000000000000FFFFULL)])<<48; __f >>= 16;
      _t[j+2] = __t;

      __f = _f[j+3], __t = _t[j+3];
      __t ^= (X[((__f)& 0x000000000000FFFFULL)])<<0;  __f >>= 16;
      __t ^= (X[((__f)& 0x000000000000FFFFULL)])<<16; __f >>= 16;
      __t ^= (X[((__f)& 0x000000000000FFFFULL)])<<32; __f >>= 16;
      __t ^= (X[((__f)& 0x000000000000FFFFULL)])<<48; __f >>= 16;
      _t[j+3] = __t;
    }
    for( ; j<to_x->width-1; j++) {
      __f = _f[j], __t = _t[j];
      __t ^= (X[((__f)& 0x000000000000FFFFULL)])<<0;  __f >>= 16;
      __t ^= (X[((__f)& 0x000000000000FFFFULL)])<<16; __f >>= 16;
      __t ^= (X[((__f)& 0x000000000000FFFFULL)])<<32; __f >>= 16;
      __t ^= (X[((__f)& 0x000000000000FFFFULL)])<<48; __f >>= 16;
      _t[j] = __t;
    }

    switch(to_x->ncols % RADIX) {
    case  0: _t[j] ^= ((word)X[(int)((_f[j] & 0x000000000000FFFFULL)>> 0)])<< 0;
    case 48: _t[j] ^= ((word)X[(int)((_f[j] & 0x00000000FFFF0000ULL)>>16)])<<16;
    case 32: _t[j] ^= ((word)X[(int)((_f[j] & 0x0000FFFF00000000ULL)>>32)])<<32;
    case 16: _t[j] ^= ((word)X[(int)((_f[j] & 0xFFFF000000000000ULL)>>48)])<<48;
    };

  }  else {
    for(size_t j=start_col; j<B->ncols; j++) {
      mzed_add_elem(A, ar, j, X[mzed_read_elem(B, br, j)]);
    }
  }
}

/**
 * \brief Gaussian elimination.
 * 
 * This will do Gaussian elimination on the matrix A.  If full=0, then
 * it will do triangular style elimination, and if full=1, it will do
 * Gauss-Jordan style, or full elimination.
 *
 * \param A Matrix
 * \param full Gauss-Jordan style or upper triangular form only.
 *
 * \wordoffset
 * 
 */

size_t mzed_echelonize_naive(mzed_t *A, int full);


static inline void mzed_rescale_row(mzed_t *A, size_t r, size_t c, word x) {
  gf2e *ff = A->finite_field;
  /* rescale row */
  word *X = ff->mul[x];
  for(size_t l=c; l<A->ncols; l++) {
    mzed_write_elem(A, r, l, X[mzed_read_elem(A, r, l)]);
  }
}

/**************** TODO: *****************/

/**
 * \brief Create a window/view into the matrix M.
 *
 * A matrix window for M is a meta structure on the matrix M. It is
 * setup to point into the matrix so M \em must \em not be freed while the
 * matrix window is used.
 *
 * This function puts the restriction on the provided parameters that
 * all parameters must be within range for M which is not enforced
 * currently .
 *
 * Use mzed_free_window to free the window.
 *
 * \param M Matrix
 * \param lowr Starting row (inclusive)
 * \param lowc Starting column (inclusive)
 * \param highr End row (exclusive)
 * \param highc End column (exclusive)
 *
 */

mzed_t *mzed_init_window(const mzed_t *M, const size_t lowr, const size_t lowc, const size_t highr, const size_t highc);

/**
 * \brief Free a matrix window created with mzed_init_window.
 * 
 * \param A Matrix
 */

void mzed_free_window(mzed_t *A);
 
/**
 * \brief Swap the two rows rowa and rowb.
 * 
 * \param M Matrix
 * \param rowa Row index.
 * \param rowb Row index.
 */

void mzed_row_swap(mzed_t *M, const size_t rowa, const size_t rowb);

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
 */

void mzed_copy_row(mzed_t* B, size_t i, const mzed_t* A, size_t j);

/**
 * \brief Swap the two columns cola and colb.
 * 
 * \param M Matrix.
 * \param cola Column index.
 * \param colb Column index.
 */
 
void mzed_col_swap(mzed_t *M, const size_t cola, const size_t colb);

/**
 * \brief Print a matrix to stdout. 
 *
 * The output will contain colons between every 4-th column.
 *
 * \param M Matrix
 */

void mzed_print(const mzed_t *M);

/**
 * \brief Add the rows sourcerow and destrow and stores the total in
 * the row destrow.
 *
 * \param M Matrix
 * \param sourcerow Index of source row
 * \param destrow Index of target row
 *
 * \note this can be done much faster with mzed_combine.
 */

void mzed_row_add(mzed_t *M, const size_t sourcerow, const size_t destrow);

/**
 * \brief Transpose a matrix.
 *
 * This function uses the fact that:
\verbatim
   [ A B ]T    [AT CT]
   [ C D ]  =  [BT DT] 
 \endverbatim 
 * and thus rearranges the blocks recursively. 
 *
 * \param DST Preallocated return matrix, may be NULL for automatic creation.
 * \param A Matrix
 */

mzed_t *mzed_transpose(mzed_t *DST, const mzed_t *A );

/**
 * \brief Naive cubic matrix multiplication and addition
 *
 * That is, compute C such that C == C + AB.
 *
 * \param C Preallocated product matrix.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \note Normally, if you will multiply several times by b, it is
 * smarter to calculate bT yourself, and keep it, and then use the
 * function called _mzed_mul_naive
 */

mzed_t *mzed_addmul_naive(mzed_t *C, const mzed_t *A, const mzed_t *B);

/**
 * \brief Naive cubic matrix multiplication with the pre-transposed B.
 *
 * That is, compute C such that C == AB^t.
 *
 * \param C Preallocated product matrix.
 * \param A Input matrix A.
 * \param B Pre-transposed input matrix B.
 * \param clear Whether to clear C before accumulating AB
 */

mzed_t *_mzed_mul_naive(mzed_t *C, const mzed_t *A, const mzed_t *B, const int clear);

/**
 * \brief Matrix multiplication optimized for v*A where v is a vector.
 *
 * \param C Preallocated product matrix.
 * \param v Input matrix v.
 * \param A Input matrix A.
 * \param clear If set clear C first, otherwise add result to C.
 *
 */
mzed_t *_mzed_mul_va(mzed_t *C, const mzed_t *v, const mzed_t *A, const int clear);

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
 * \wordoffset
 */

mzed_t *mzed_concat(mzed_t *C, const mzed_t *A, const mzed_t *B);

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
 * \wordoffset
 */

mzed_t *mzed_stack(mzed_t *C, const mzed_t *A, const mzed_t *B);

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
 */
mzed_t *mzed_submatrix(mzed_t *S, const mzed_t *M, const size_t lowr, const size_t lowc, const size_t highr, const size_t highc);

/**
 * \brief Invert the matrix target using Gaussian elimination. 
 *
 * To avoid recomputing the identity matrix over and over again, I may
 * be passed in as identity parameter.
 *
 * \param INV Preallocated space for inversion matrix, may be NULL for automatic creation.
 * \param A Matrix to be reduced.
 * \param I Identity matrix.
 *
 * \wordoffset
 */

mzed_t *mzed_invert_naive(mzed_t *INV, mzed_t *A, const mzed_t *I);

/**
 * \brief Zero test for matrix.
 *
 * \param A Input matrix.
 *
 */
int mzed_is_zero(mzed_t *A);

/**
 * \brief Clear the given row, but only begins at the column coloffset.
 *
 * \param M Matrix
 * \param row Index of row
 * \param coloffset Column offset
 */

void mzed_row_clear_offset(mzed_t *M, const size_t row, const size_t coloffset);

/**
 * \brief Find the next nonzero entry in M starting at start_row and start_col. 
 *
 * This function walks down rows in the inner loop and columns in the
 * outer loop. If a nonzero entry is found this function returns 1 and
 * zero otherwise.
 *
 * If and only if a nonzero entry is found r and c are updated.
 *
 * \param M Matrix
 * \param start_row Index of row where to start search
 * \param start_col Index of column where to start search
 * \param r Row index updated if pivot is found
 * \param c Column index updated if pivot is found
 */

int mzed_find_pivot(mzed_t *M, size_t start_row, size_t start_col, size_t *r, size_t *c);


/**
 * \brief Return the number of nonzero entries divided by nrows *
 * ncols
 *
 * If res = 0 then 100 samples per row are made, if res > 0 the
 * function takes res sized steps within each row (res = 1 uses every
 * word).
 *
 * \param A Matrix
 * \param res Resolution of sampling
 */

double mzed_density(mzed_t *A, int res);

/**
 * \brief Return the number of nonzero entries divided by nrows *
 * ncols considering only the submatrix starting at (r,c).
 *
 * If res = 0 then 100 samples per row are made, if res > 0 the
 * function takes res sized steps within each row (res = 1 uses every
 * word).
 *
 * \param A Matrix
 * \param res Resolution of sampling
 * \param r Row to start counting
 * \param c Column to start counting
 */

double _mzed_density(mzed_t *A, int res, size_t r, size_t c);


/**
 * \brief Return the first row with all zero entries.
 *
 * If no such row can be found returns nrows.
 *
 * \param A Matrix
 */

size_t mzed_first_zero_row(mzed_t *A);


#endif //MATRIX_H

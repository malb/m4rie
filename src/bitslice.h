/**
 * \file bitslice.h
 * \brief Bitsliced Extension Matrices
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

#ifndef BITSLICE_H
#define BITSLICE_H

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
#include "gf2e_matrix.h"

/**
 * Dense matrices over GF(2^k) represented as slices of matrices over GF(2).
 */

typedef struct {
  /**
   * A->x[e][i,j] is the e^th bit of the entry A[i,j].
   */
  mzd_t *x[16];  // We only support 10 but 16 migh help with alignment.
  gf2e *finite_field;
  rci_t nrows;
  rci_t ncols;
  int depth; // the number of slices
} mzd_slice_t;

static inline mzd_slice_t *mzd_slice_init(gf2e *ff, const rci_t m, const rci_t n) {

  /**
   * @TODO: avoid all these malloc() calls and call it once only
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

static inline mzd_slice_t *_mzd_slice_adapt_depth(mzd_slice_t *A, const int new_depth) {
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
 * \brief Pack a bitslice matrix into a classical represenation.
 *
 * \param A Matrix over GF(2^k) or NULL
 * \param Z Bitslice matrix over GF(2^k)
 */

mzed_t *mzed_cling(mzed_t *A, const mzd_slice_t *Z);

/**
 * \brief Unpack the matrix Z into bitslice representation.
 *
 * \param A Bitslice matrix or NULL
 * \param Z Input matrix
 */

mzd_slice_t *mzed_slice(mzd_slice_t *A, const mzed_t *Z);

static inline mzd_slice_t *mzd_slice_concat(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B) {
  if(C == NULL)
    C = mzd_slice_init(A->finite_field, A->nrows, A->ncols + B->ncols);

  for(int i=0; i<A->depth; i++) {
    mzd_concat(C->x[i], A->x[i], B->x[i]);
  }
  return C;
}

static inline mzd_slice_t *mzd_slice_stack(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B) {
  if(C == NULL)
    C = mzd_slice_init(A->finite_field, A->nrows + B->nrows, A->ncols);

  for(int i=0; i<A->depth; i++) {
    mzd_stack(C->x[i], A->x[i], B->x[i]);
  }
  return C;
}

static inline mzd_slice_t *mzd_slice_submatrix(mzd_slice_t *S, const mzd_slice_t *A, 
                                               const size_t lowr, const size_t lowc, const size_t highr, const size_t highc) {
  if(S==NULL)
    S = mzd_slice_init(A->finite_field, highr - lowr, highc - lowc);

  for(int i=0; i<A->depth; i++) {
    mzd_submatrix(S->x[i], A->x[i], lowr, lowc, highr, highc);
  }
  return S;
}

static inline mzd_slice_t *mzd_slice_init_window(const mzd_slice_t *A, 
                                                 const size_t lowr, const size_t lowc, const size_t highr, const size_t highc) {
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

static inline void mzd_slice_free_window(mzd_slice_t *A) {
  for(int i=0; i<A->depth; i++) {
    mzd_free_window(A->x[i]);
  }
  m4ri_mm_free(A);
}

static inline mzd_slice_t *_mzd_slice_add(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B) {
  for(int i=0; i<A->depth; i++) 
    _mzd_add(C->x[i], A->x[i], B->x[i]);
  return C;
}

static inline mzd_slice_t *mzd_slice_add(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B) {
  if ( (A->finite_field != B->finite_field) | (A->nrows != B->nrows) | (A->ncols != B->ncols) ) 
    m4ri_die("mzd_slice_add: input matrices A (%d x %d) and B (%d x %d) do not match.\n",A->nrows,A->ncols, B->nrows,B->ncols);

  if(C == NULL) 
    mzd_slice_init(A->finite_field, A->nrows, A->ncols);
  else if ( (A->finite_field != C->finite_field) | (A->nrows != C->nrows) | (A->ncols != C->ncols) )
    m4ri_die("mzd_slice_add: input matrix A (%d x %d) and output matrix (%d x %d) do not match.\n",A->nrows,A->ncols, C->nrows, C->ncols);

  return _mzd_slice_add(C,A,B);
}

#define mzd_slice_sub mzd_slice_add

#define _mzd_slice_sub _mzd_slice_add

static inline void mzd_slice_randomize(mzd_slice_t *A) {
  for(int i=0; i<-A->depth; i++) {
    mzd_randomize(A->x[i]);
  }
}
 
static inline mzd_slice_t *mzd_slice_copy(mzd_slice_t *B, const mzd_slice_t *A) {
  if(B == NULL)
    B = mzd_slice_init(A->finite_field, A->nrows, A->ncols);
  
  for(int i=0; i<A->depth; i++) {
    mzd_copy(B->x[i],A->x[i]);
  }
  return B;
}

void mzd_slice_set_ui(mzd_slice_t *A, word value);

static inline word mzd_slice_read_elem(const mzd_slice_t *A, const rci_t row, const rci_t col) {
  word ret = 0;
  for(int i=0; i<A->depth; i++) {
    ret |= mzd_read_bit(A->x[i], row, col)<<i;
  }
  return ret;
}

static inline void mzd_slice_add_elem(mzd_slice_t *A, const rci_t row, const rci_t col, word elem) {
  for(int i=0; i<A->depth; i++) {
    __mzd_xor_bits(A->x[i], row, col, 1, elem&1);
    elem=elem>>1;
  }
}

static inline void mzd_slice_write_elem(mzd_slice_t *A, const rci_t row, const rci_t col, word elem) {
  for(int i=0; i<A->depth; i++) {
    mzd_write_bit(A->x[i], row, col, elem&1);
    elem=elem>>1;
  }
}

static inline int mzd_slice_cmp(mzd_slice_t *A, mzd_slice_t *B) {
  int r = 0;
  if ((A->finite_field != B->finite_field) | (A->depth != B->depth) )
    return -1;
  for(int i=0; i<A->depth; i++)
    r |= mzd_cmp(A->x[i],B->x[i]);
  return r;
}

static inline int mzd_slice_is_zero(mzd_slice_t *A) {
  for(int i=0; i<A->depth; i++) {
    if (mzd_is_zero(A->x[i]))
      return 1;
  }
  return 0;
}

static inline void mzd_slice_rescale_row(mzd_slice_t *A, rci_t r, rci_t c, word *X) {
  mzd_slice_t *A_w = mzd_slice_init_window(A, r, 0, r+1, A->ncols);
  mzed_t *A_we = mzed_cling(NULL, A_w);

  mzed_rescale_row(A_we, r, c, X);

  mzed_slice(A_w, A_we);
  mzed_free(A_we);
  mzd_slice_free_window(A_w);
}

static inline void mzd_slice_row_swap(mzd_slice_t *A, const rci_t rowa, const rci_t rowb) {
  for(int i=0; i<A->depth; i++) {
    mzd_row_swap(A->x[i], rowa, rowb);
  }
}

static inline void mzd_slice_copy_row(mzd_slice_t* B, size_t i, const mzd_slice_t* A, size_t j) {
  for(int ii=0; ii<A->depth; ii++)
    mzd_copy_row(B->x[ii], i, A->x[ii], j);
}
 
static inline void mzd_slice_col_swap(mzd_slice_t *A, const rci_t cola, const rci_t colb) {
  for(int i=0; i<A->depth; i++)
    mzd_col_swap(A->x[i], cola, colb);
}

static inline void mzd_slice_row_add(mzd_slice_t *A, const rci_t sourcerow, const rci_t destrow) {
  for(int i=0; i<A->depth; i++)
    mzd_row_add(A->x[i], sourcerow, destrow);
}

static inline void mzd_slice_row_clear_offset(mzd_slice_t *A, const rci_t row, const rci_t coloffset) {
  for(int i=0; i<A->depth; i++) 
    mzd_row_clear_offset(A->x[i], row, coloffset);
}

static inline void mzd_slice_print(const mzd_slice_t *A) {
  mzed_t *B = mzed_cling(NULL, A);
  mzed_print(B);
  mzed_free(B);
}

/**
 * \brief Unpack the matrix Z over GF(2^2) into bitslice representation.
 *
 * Elements in GF(2^2) can be represented as x*a + y where a is a root
 * of x^2 + x + 1. A0 contains the coefficients for x while A1
 * contains the coefficients for y.
 *
 * \param A Zero bitslice matrix over GF(2^2)
 * \param Z Matrix over GF(2^2)
 */

mzd_slice_t *_mzed_slice2(mzd_slice_t *A, const mzed_t *Z);

/**
 * \brief Unpack the matrix Z over GF(2^2) into bitslice representation.
 *
 * Elements in GF(2^2) can be represented as x*a + y where a is a root
 * of x^2 + x + 1. A0 contains the coefficients for x while A1
 * contains the coefficients for y.
 *
 * \param A Bitslice matrix over GF(2^2) or NULL
 * \param Z Matrix over GF(2^2)
 */

static inline mzd_slice_t *mzed_slice2(mzd_slice_t *A, const mzed_t *Z) {
  if (A == NULL)
    A = mzd_slice_init(Z->finite_field, Z->nrows, Z->ncols);
  else
    mzd_slice_set_ui(A, 0);
 
  return _mzed_slice2(A, Z);
}

/**
 */

mzd_slice_t *_mzed_slice4(mzd_slice_t *A, const mzed_t *Z);

/**
 * \brief Pack a bitslice matrix into a classical represenation over GF(2^2).
 *
 * Elements in GF(2^2) can be represented as c_1*a + c_0 where a is a
 * root of x^2 + x + 1. A1 contains the coefficients for c_1 while A0
 * contains the coefficients for c_0.
 *
 * \param A Matrix over GF(2^2), must be zero
 * \param Z Bitslice matrix over GF(2^2)
 */

mzed_t *_mzed_cling2(mzed_t *A, const mzd_slice_t *Z);

/**
 * \brief Pack a bitslice matrix into a classical represenation over GF(2^2).
 *
 * Elements in GF(2^2) can be represented as c_1*a + c_0 where a is a
 * root of x^2 + x + 1. A1 contains the coefficients for c_1 while A0
 * contains the coefficients for c_0.
 *
 * \param A Matrix over GF(2^2) or NULL
 * \param Z Bitslice matrix over GF(2^2)
 */

static inline mzed_t* mzed_cling2(mzed_t *A, const mzd_slice_t *Z) {
  if (A == NULL) 
    A = mzed_init(Z->finite_field, Z->nrows, Z->ncols);
  else
    mzed_set_ui(A, 0);

  _mzed_cling2(A, Z);
  return A;
}

mzed_t *_mzed_cling4(mzed_t *A, const mzd_slice_t *Z);

static inline mzed_t* mzed_cling4(mzed_t *A, const mzd_slice_t *Z) {
  if (A == NULL) 
    A = mzed_init(Z->finite_field, Z->nrows, Z->ncols);
  else
    mzed_set_ui(A, 0);

  _mzed_cling2(A, Z);
  return A;
}

mzd_slice_t *_mzd_slice_mul_karatsuba2(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B);
mzd_slice_t *_mzd_slice_mul_karatsuba3(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B);
mzd_slice_t *_mzd_slice_mul_karatsuba4(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B);

/**
 * \brief Compute C += A*B using Karatsuba multiplication of polynomials over GF(2).
 *
 * Matrices over GF(2^e) can be represented as polynomials with matrix
 * coefficients where the matrices are in GF(2). This function uses
 * this fact to reduce matrix multiplication over GF(2^e) to matrix
 * multiplication over GF(2).
 *
 * As an example consider GF(2^2), the minimal polynomial is x^2 + x +
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
 * \note This function is only implemented for GF(2^e) for e<=4 so far.
 *
 * \sa mzed_mul()
 *
 * \wordoffset
 */

static inline mzd_slice_t *_mzd_slice_mul_karatsuba(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B) {
  switch(A->finite_field->degree) {
  case  2:
    C = _mzd_slice_mul_karatsuba2(C, A, B); break;
  case  3:
    C = _mzd_slice_mul_karatsuba3(C, A, B); break;
  case  4:
    C = _mzd_slice_mul_karatsuba4(C, A, B); break;
  case  5:
  case  6:
  case  7:
  case  8:
  case  9:
  case 10:
  default:
    m4ri_die("mzed_mul_karatsuba: only implemented for GF(2^e) with e <= 4");
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
 *
 * \wordoffset
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
 *
 * \wordoffset
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
 *
 * \wordoffset
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
 *
 * \wordoffset
 */

static inline mzd_slice_t *mzd_slice_addmul(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B) {
  return mzd_slice_addmul_karatsuba(C,A,B);
}

/**
 * \brief Compute C += A*B using Karatsuba multiplication of polynomials over GF(2).
 *
 * \param C Preallocated return matrix, may be NULL for automatic creation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa _mzd_slice_mul_karatsuba
 *
 * \wordoffset
 */

static inline mzed_t *_mzed_mul_karatsuba(mzed_t *C, const mzed_t *A, const mzed_t *B) {
  mzd_slice_t *As,*Bs,*Cs;
  if(C)
    Cs = mzed_slice(NULL,C);
  else
    Cs = NULL;
  As = mzed_slice(NULL,A);
  Bs = mzed_slice(NULL,B);

  _mzd_slice_mul_karatsuba(Cs, As, Bs);

  C = mzed_cling(C, Cs);

  mzd_slice_free(As);
  mzd_slice_free(Bs);
  mzd_slice_free(Cs);
  return C;
}

/**
 * \brief Compute C = A*B.
 *
 * \param C Preallocated return matrix, may be NULL for automatic creation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa _mzd_slice_mul_karatsuba
 *
 * \wordoffset
 */

static inline mzed_t *mzed_mul_karatsuba(mzed_t *C, const mzed_t *A, const mzed_t *B) {
  if (A->ncols != B->nrows || A->finite_field != B->finite_field) 
    m4ri_die("mzed_mul_karatsuba: rows, columns and fields must match.\n");
  if (C != NULL) {
    if (C->finite_field != A->finite_field || C->nrows != A->nrows || C->ncols != B->ncols) 
      m4ri_die("mzed_mul_karatsuba: rows and columns of returned matrix must match.\n");
    mzed_set_ui(C,0);
  }
  return _mzed_mul_karatsuba(C, A, B);
}

/**
 * \brief Compute C += A*B.
 *
 * \param C Preallocated return matrix.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa _mzd_slice_mul_karatsuba
 *
 * \wordoffset
 */

static inline mzed_t *mzed_addmul_karatsuba(mzed_t *C, const mzed_t *A, const mzed_t *B) {
  assert(C != NULL);
  if (A->ncols != B->nrows || A->finite_field != B->finite_field) 
    m4ri_die("mzed_addmul_karatsuba: rows, columns and fields must match.\n");
  if (C->finite_field != A->finite_field || C->nrows != A->nrows || C->ncols != B->ncols) 
    m4ri_die("mzed_addmul_karatsuba: rows and columns of returned matrix must match.\n");
  return _mzed_mul_karatsuba(C, A, B);
}

#endif //BITSLICE_H

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
 * \brief Unpack the matrix A over GF(2^2) into the two matrices A0 and A1 over GF(2).
 *
 * Elements in GF(2^2) can be represented as x*a + y where a is a root
 * of x^2 + x + 1. A0 contains the coefficients for x while A1
 * contains the coefficients for y.
 *
 * \param A1 Matrix over GF(2), must be zero 
 * \param A0 Matrix over GF(2), must be zero
 * \param A Matrix over GF(2^2)
 */

void _mzed_slice2(mzd_t *A1, mzd_t *A0, const mzed_t *A);


/**
 * \brief Unpack the matrix A over GF(2^2) into the two matrices A0 and A1 over GF(2).
 *
 * Elements in GF(2^2) can be represented as x*a + y where a is a root
 * of x^2 + x + 1. A0 contains the coefficients for x while A1
 * contains the coefficients for y.
 *
 * \param A1 Matrix over GF(2) 
 * \param A0 Matrix over GF(2)
 * \param A Matrix over GF(2^2)
 */

static inline void mzed_slice2(mzd_t *A1, mzd_t *A0, const mzed_t *A) {
  mzd_set_ui(A0, 0);
  mzd_set_ui(A1, 0);
  _mzed_slice2(A0, A1, A);
}

/**
 * \brief Pack the matrices A0 and A1 over GF(2) to A over GF(2^2).
 *
 * Elements in GF(2^2) can be represented as c_1*a + c_0 where a is a
 * root of x^2 + x + 1. A1 contains the coefficients for c_1 while A0
 * contains the coefficients for c_0.
 *
 * \param A Matrix over GF(2^2), must be zero
 * \param A1 Matrix over GF(2)
 * \param A0 Matrix over GF(2) 
 */

void _mzed_cling2(mzed_t *A, const mzd_t *A1, const mzd_t *A0);

/**
 * \brief Pack the matrices A0 and A1 over GF(2) to A over GF(2^2).
 *
 * Elements in GF(2^2) can be represented as c_1*a + c_0 where a is a root
 * of x^2 + x + 1. A0 contains the coefficients for c_0 while A1
 * contains the coefficients for c_1.
 *
 * \param A Matrix over GF(2^2)
 * \param A1 Matrix over GF(2) 
 * \param A0 Matrix over GF(2)
 */

static inline void mzed_cling2(mzed_t *A, const mzd_t *A1, const mzd_t *A0) {
  mzed_set_ui(A, 0);
  _mzed_cling2(A, A0, A1);
}

/**
 * \brief Compute C == A*B using Karatsuba multiplication of polynomials over GF(2).
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
 * A0*B0*x^2 + (A0*B1 + A1*B0)*x + A1*B1.
 * \f]
 * Reduction modulo x^2 + x + 1 gives
 * \f[
 * (A0*B0 + A0*B1 + A1*B0)*x + A1*B1 + A0*B0.
 * \f]
 * This can be re-written as
 * \f[
 * ((A0 + A1)*(B0 + B1) + A1*B1)*x + A1*B1 + A0*B0
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
 * \note This function is only implemented for GF(2^2) so far.
 *
 * \sa mzed_mul()
 *
 * \wordoffset
 */

mzed_t *mzed_mul_karatsuba(mzed_t *C, const mzed_t *A, const mzed_t *B);

/**
 * \brief Compute C == A*B over GF(2^2) using Karatsuba multiplication.
 *
 * \sa mzed_mul_karatsuba()
 *
 * \param C Preallocated return matrix, may be NULL for automatic creation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \wordoffset
 */

mzed_t *_mzed_mul_karatsuba2(mzed_t *C, const mzed_t *A, const mzed_t *B);

#endif //BITSLICE_H

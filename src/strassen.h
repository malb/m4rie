#ifndef GF2E_STRASSEN_H
#define GF2E_STRASSEN_H 

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

/**
 * \brief C such that C == AB.
 *
 * This function uses Strassen-Winograd multiplication (Bodrato
 * variant) recursively until it reaches the cutoff, where it switches
 * to Travolta table based multiplication or naive multiplication.
 *
 * \param C Preallocated product matrix, may be NULL for allocation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \ingroup Multiplication
 */

mzed_t *mzed_mul_strassen(mzed_t *C, const mzed_t *A, const mzed_t *B, int cutoff);

/**
 * \brief C such that C == C + AB.
 *
 * This function uses Strassen-Winograd multiplication (Bodrato
 * variant) recursively until it reaches the cutoff, where it switches
 * to Travolta table based multiplication or naive multiplication.
 *
 * \param C Preallocated product matrix, may be NULL for allocation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \ingroup Multiplication
 */

mzed_t *mzed_addmul_strassen(mzed_t *C, const mzed_t *A, const mzed_t *B, int cutoff);

/**
 * \brief C such that C == AB.
 *
 * This function uses Strassen-Winograd multiplication (Bodrato
 * variant) recursively until it reaches the cutoff, where it switches
 * to Travolta table based multiplication or naive multiplication.
 *
 * \param C Preallocated product matrix.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \ingroup Multiplication
 *
 * \wordoffset
 */

mzed_t *_mzed_mul_strassen_even(mzed_t *C, const mzed_t *A, const mzed_t *B, int cutoff);


/**
 * \brief C such that C == C + AB.
 *
 * This function uses Strassen-Winograd multiplication (Bodrato
 * variant) recursively until it reaches the cutoff, where it switches
 * to Travolta table based multiplication or naive multiplication.
 *
 * \param C Preallocated product matrix.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \ingroup Multiplication
 *
 * \wordoffset
 */

mzed_t *_mzed_addmul_strassen_even(mzed_t *C, const mzed_t *A, const mzed_t *B, int cutoff);

/**
 * \brief C such that C == AB.
 *
 * This function uses Strassen-Winograd multiplication (Bodrato
 * variant) recursively until it reaches the cutoff, where it switches
 * to Travolta table based multiplication or naive multiplication.
 *
 * \param C Preallocated product matrix.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \ingroup Multiplication
 */

static inline mzed_t *_mzed_mul_strassen(mzed_t *C, const mzed_t *A, const mzed_t *B, int cutoff) {
  return _mzed_mul_strassen_even(C, A, B, cutoff);
}

/**
 * \brief C such that C == C + AB.
 *
 * This function uses Strassen-Winograd multiplication (Bodrato
 * variant) recursively until it reaches the cutoff, where it switches
 * to Travolta table based multiplication or naive multiplication.
 *
 * \param C Preallocated product matrix.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \ingroup Multiplication
 */

static inline mzed_t *_mzed_addmul_strassen(mzed_t *C, const mzed_t *A, const mzed_t *B, int cutoff) {
  return _mzed_addmul_strassen_even(C, A, B, cutoff);
}


rci_t _mzed_strassen_cutoff(const mzed_t *C, const mzed_t *A, const mzed_t *B);

#endif //GF2E_STRASSEN_H

#ifndef TRAVOLTA_H
#define TRAVOLTA_H

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

#include "finite_field.h"
#include "gf2e_matrix.h"

void mzed_make_table(const mzed_t *A, size_t r, size_t c, mzed_t *T, size_t *L, gf2e *ff);

/**
 * \brief Compute C such that C == AB using Travolta tables.
 *
 * \param C Preallocated return matrix, may be NULL for automatic creation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa mzed_mul_travolta1 mzed_mul
 *
 * \wordoffset
 */

mzed_t *mzed_mul_travolta(mzed_t *C, const mzed_t *A, const mzed_t *B);

/**
 * \brief Compute C such that C == C + AB using Travolta tables.
 *
 * \param C Preallocated product matrix, may be NULL for automatic creation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa _mzed_mul_travolta mzed_mul
 *
 * \wordoffset
 */

mzed_t *mzed_addmul_travolta(mzed_t *C, const mzed_t *A, const mzed_t *B);

/**
 * \brief Compute C such that C == C + AB using Travolta tables.
 *
 * This is a simple implementation.
 *
 * \param C Preallocated product matrix.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa mzed_mul_travolta mzed_mul
 *
 * \wordoffset
 */

mzed_t *_mzed_mul_travolta0(mzed_t *C, const mzed_t *A, const mzed_t *B);

/**
 * \brief Compute C such that C == C + AB using Travolta tables.
 *
 * This is an optimised implementation.
 *
 * \param C Preallocated product matrix.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa mzed_mul_travolta0 mzed_mul
 *
 * \wordoffset
 */

mzed_t *_mzed_mul_travolta(mzed_t *C, const mzed_t *A, const mzed_t *B);

/**
 * \brief Reduce matrix A to row echelon form using Gauss-Travolta
 * elimination.
 *
 * \param A Matrix to be reduced.
 * \param full If set to true, the reduced row echelon form will be
 * computed.
 *
 * \wordoffset
 */

size_t mzed_echelonize_travolta(mzed_t *A, int full);

/**
 * \brief Invert the matrix A using Gauss-Travolta elimination. 
 *
 * \param B Preallocated space for inversion matrix, may be NULL for
 * automatic creation.
 * \param A Matrix to be inverted.
 *
 * \wordoffset
 */

mzed_t *mzed_invert_travolta(mzed_t *B, const mzed_t *A);

/**
 * \brief The function looks up 6 entries from position i,startcol in
 * each row and adds the appropriate row from T to the row i.
 *
 * This process is iterated for i from startrow to stoprow
 * (exclusive).
 *
 * \param M Matrix to operate on
 * \param startrow top row which is operated on
 * \param endrow bottom row which is operated on
 * \param startcol Starting column for addition
 * \param T contains the correct row to be added
 * \param L Contains row number to be added
 *
 * \wordoffset
 */

static inline void mzed_process_rows(mzed_t *M, size_t startrow, size_t endrow, size_t startcol, mzed_t *T, size_t *L) {
  mzd_process_rows(M->x, startrow, endrow, startcol*M->w, M->w, T->x, L);
}

/**
 * \brief Same as mzed_process_rows but works with two Travolta tables
 * in parallel.
 *
 * \param M Matrix to operate on
 * \param startrow top row which is operated on
 * \param endrow bottom row which is operated on
 * \param startcol Starting column for addition
 * \param T0 contains the correct row to be added
 * \param L0 Contains row number to be added
 * \param T1 contains the correct row to be added
 * \param L1 Contains row number to be added
 *
 * \wordoffset
 */

static inline void mzed_process_rows2(mzed_t *M, size_t startrow, size_t endrow, size_t startcol, mzed_t *T0, size_t *L0, mzed_t *T1, size_t *L1) {
  mzd_process_rows2(M->x, startrow, endrow, startcol*M->w, 2*M->w, T0->x, L0, T1->x, L1);
}

/**
 * \brief Same as mzed_process_rows but works with three Travolta
 * tables in parallel.
 *
 * \param M Matrix to operate on
 * \param startrow top row which is operated on
 * \param endrow bottom row which is operated on
 * \param startcol Starting column for addition
 * \param T0 contains the correct row to be added
 * \param L0 Contains row number to be added
 * \param T1 contains the correct row to be added
 * \param L1 Contains row number to be added
 * \param T2 contains the correct row to be added
 * \param L2 Contains row number to be added
 *
 * \wordoffset
 */

static inline void mzed_process_rows3(mzed_t *M, size_t startrow, size_t endrow, size_t startcol, 
                                      mzed_t *T0, size_t *L0, mzed_t *T1, size_t *L1,
                                      mzed_t *T2, size_t *L2) {
  mzd_process_rows3(M->x, startrow, endrow, startcol*M->w, 3*M->w, T0->x, L0, T1->x, L1, T2->x, L2);
}

/**
 * \brief Same as mzed_process_rows but works with four Travolta
 * tables in parallel.
 *
 * \param M Matrix to operate on
 * \param startrow top row which is operated on
 * \param endrow bottom row which is operated on
 * \param startcol Starting column for addition
 * \param T0 contains the correct row to be added
 * \param L0 Contains row number to be added
 * \param T1 contains the correct row to be added
 * \param L1 Contains row number to be added
 * \param T2 contains the correct row to be added
 * \param L2 Contains row number to be added
 * \param T3 contains the correct row to be added
 * \param L3 Contains row number to be added
 *
 * \wordoffset
 */

static inline void mzed_process_rows4(mzed_t *M, size_t startrow, size_t endrow, size_t startcol,
                                      mzed_t *T0, size_t *L0, mzed_t *T1, size_t *L1,
                                      mzed_t *T2, size_t *L2, mzed_t *T3, size_t *L3) {
  mzd_process_rows4(M->x, startrow, endrow, startcol*M->w, 4*M->w, T0->x, L0, T1->x, L1, T2->x, L2, T3->x, L3);
}


/**
 * \brief Same as mzed_process_rows but works with five Travolta
 * tables in parallel.
 *
 * \param M Matrix to operate on
 * \param startrow top row which is operated on
 * \param endrow bottom row which is operated on
 * \param startcol Starting column for addition
 * \param T0 contains the correct row to be added
 * \param L0 Contains row number to be added
 * \param T1 contains the correct row to be added
 * \param L1 Contains row number to be added
 * \param T2 contains the correct row to be added
 * \param L2 Contains row number to be added
 * \param T3 contains the correct row to be added
 * \param L3 Contains row number to be added
 * \param T4 contains the correct row to be added
 * \param L4 Contains row number to be added
 *
 * \wordoffset
 */

static inline void mzed_process_rows5(mzed_t *M, size_t startrow, size_t endrow, size_t startcol,
                                      mzed_t *T0, size_t *L0, mzed_t *T1, size_t *L1,
                                      mzed_t *T2, size_t *L2, mzed_t *T3, size_t *L3,
                                      mzed_t* T4, size_t *L4) {
  mzd_process_rows5(M->x, startrow, endrow, startcol*M->w, 5*M->w, T0->x, L0, T1->x, L1, T2->x, L2, T3->x, L3, T4->x, L4);
}


/**
 * \brief Same as mzed_process_rows but works with six Travolta tables
 * in parallel.
 *
 * \param M Matrix to operate on
 * \param startrow top row which is operated on
 * \param endrow bottom row which is operated on
 * \param startcol Starting column for addition
 * \param T0 contains the correct row to be added
 * \param L0 Contains row number to be added
 * \param T1 contains the correct row to be added
 * \param L1 Contains row number to be added
 * \param T2 contains the correct row to be added
 * \param L2 Contains row number to be added
 * \param T3 contains the correct row to be added
 * \param L3 Contains row number to be added
 * \param T4 contains the correct row to be added
 * \param L4 Contains row number to be added
 * \param T5 contains the correct row to be added
 * \param L5 Contains row number to be added
 *
 * \wordoffset
 */

static inline void mzed_process_rows6(mzed_t *M, size_t startrow, size_t endrow, size_t startcol,
                                      mzed_t *T0, size_t *L0, mzed_t *T1, size_t *L1,
                                      mzed_t *T2, size_t *L2, mzed_t *T3, size_t *L3,
                                      mzed_t* T4, size_t *L4, mzed_t *T5, size_t *L5) {
  mzd_process_rows6(M->x, startrow, endrow, startcol*M->w, 6*M->w, T0->x, L0, T1->x, L1, T2->x, L2, T3->x, L3, T4->x, L4, T5->x, L5);
}

#endif //TRAVOLTA_H

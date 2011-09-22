#ifndef M4RIE_TRAVOLTA_H
#define M4RIE_TRAVOLTA_H

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
#include "bitslice.h"

typedef struct {
  rci_t *L;
  mzed_t *T;
  mzed_t *M;
} tt_mzed_t;

tt_mzed_t *tt_mzed_init(const gf2e *ff, const rci_t ncols);

void tt_mzed_free(tt_mzed_t *);

/**
 * \brief Construct Travolta table T for row r of A, and element A[r,c].
 *
 * \param T Preallocated Travolta table or NULL.
 * \param A Matrix.
 * \param r Row index.
 * \param c Column index.
 */

tt_mzed_t * mzed_make_table(tt_mzed_t *T, const mzed_t *A, const rci_t r, const rci_t c);

/**
 * \brief Compute C such that C == AB using Travolta tables.
 *
 * \param C Preallocated return matrix, may be NULL for automatic creation.
 * \param A Input matrix A.
 * \param B Input matrix B.
 *
 * \sa mzed_mul _mzed_mul_travolta0
 * 
 * \ingroup Multiplication
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
 * \ingroup Multiplication
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
 * \ingroup Multiplication
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
 * \sa mzed_mul
 * 
 * \ingroup Multiplication
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
 * \ingroup Echelon
 */

rci_t mzed_echelonize_travolta(mzed_t *A, int full);

/**
 * \brief Invert the matrix A using Gauss-Travolta elimination. 
 *
 * \param B Preallocated space for inversion matrix, may be NULL for
 * automatic creation.
 * \param A Matrix to be inverted.
 * 
 * \ingroup Echelon
 */

mzed_t *mzed_invert_travolta(mzed_t *B, const mzed_t *A);

void mzed_trsm_lower_left_travolta(const mzed_t *L, mzed_t *B);

void mzd_slice_trsm_lower_left_travolta(const mzd_slice_t *L, mzd_slice_t *B);

/**
 * \brief L*E = P*A
 *
 * \ingroup PLE
 */

rci_t mzed_ple_travolta(mzed_t *A, mzp_t *P, mzp_t *Q);

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
 * \param T Travolta table
 * 
 * \ingroup RowOperations
 */

static inline void mzed_process_rows(mzed_t *M, const rci_t startrow, const rci_t endrow, rci_t startcol, const tt_mzed_t *T) {
  mzd_process_rows(M->x, startrow, endrow, startcol*M->w, M->w, T->T->x, T->L);
}

/**
 * \brief Same as mzed_process_rows but works with two Travolta tables
 * in parallel.
 *
 * \param M Matrix to operate on
 * \param startrow top row which is operated on
 * \param endrow bottom row which is operated on
 * \param startcol Starting column for addition
 * \param T0 Travolta table
 * \param T1 Travolta table
 * 
 * \ingroup RowOperations
 */

static inline void mzed_process_rows2(mzed_t *M, const rci_t startrow, const rci_t endrow, const rci_t startcol, 
                                      const tt_mzed_t *T0, const tt_mzed_t *T1) {
  mzd_process_rows2(M->x, startrow, endrow, startcol*M->w, 2*M->w, T0->T->x, T0->L, T1->T->x, T1->L);
}

/**
 * \brief Same as mzed_process_rows but works with three Travolta
 * tables in parallel.
 *
 * \param M Matrix to operate on
 * \param startrow top row which is operated on
 * \param endrow bottom row which is operated on
 * \param startcol Starting column for addition
 * \param T0 Travolta table
 * \param T1 Travolta table
 * \param T2 Travolta table
 * 
 * \ingroup RowOperations
 */

static inline void mzed_process_rows3(mzed_t *M, const rci_t startrow, const rci_t endrow, const rci_t startcol, 
                                      const tt_mzed_t *T0, const tt_mzed_t *T1, const tt_mzed_t *T2) {
  mzd_process_rows3(M->x, startrow, endrow, startcol*M->w, 3*M->w, T0->T->x, T0->L, T1->T->x, T1->L, T2->T->x, T2->L);
}

/**
 * \brief Same as mzed_process_rows but works with four Travolta
 * tables in parallel.
 *
 * \param M Matrix to operate on
 * \param startrow top row which is operated on
 * \param endrow bottom row which is operated on
 * \param startcol Starting column for addition
 * \param T0 Travolta table
 * \param T1 Travolta table
 * \param T2 Travolta table
 * \param T3 Travolta table
 * 
 * \ingroup RowOperations
 */

static inline void mzed_process_rows4(mzed_t *M, const rci_t startrow, const rci_t endrow, const rci_t startcol,
                                      const tt_mzed_t *T0, const tt_mzed_t *T1, const tt_mzed_t *T2, const tt_mzed_t *T3) {
  mzd_process_rows4(M->x, startrow, endrow, startcol*M->w, 4*M->w, T0->T->x, T0->L, T1->T->x, T1->L, T2->T->x, T2->L, T3->T->x, T3->L);
}


/**
 * \brief Same as mzed_process_rows but works with five Travolta
 * tables in parallel.
 *
 * \param M Matrix to operate on
 * \param startrow top row which is operated on
 * \param endrow bottom row which is operated on
 * \param startcol Starting column for addition
 * \param T0 Travolta table
 * \param T1 Travolta table
 * \param T2 Travolta table
 * \param T3 Travolta table
 * \param T4 Travolta table
 * 
 * \ingroup RowOperations
 */

static inline void mzed_process_rows5(mzed_t *M, const rci_t startrow, const rci_t endrow, const rci_t startcol,
                                      const tt_mzed_t *T0, const tt_mzed_t *T1, const tt_mzed_t *T2, const tt_mzed_t *T3, const tt_mzed_t *T4) {
  mzd_process_rows5(M->x, startrow, endrow, startcol*M->w, 5*M->w, T0->T->x, T0->L, T1->T->x, T1->L, T2->T->x, T2->L, T3->T->x, T3->L, T4->T->x, T4->L);
}


/**
 * \brief Same as mzed_process_rows but works with six Travolta tables
 * in parallel.
 *
 * \param M Matrix to operate on
 * \param startrow top row which is operated on
 * \param endrow bottom row which is operated on
 * \param startcol Starting column for addition
 * \param T0 Travolta table
 * \param T1 Travolta table
 * \param T2 Travolta table
 * \param T3 Travolta table
 * \param T4 Travolta table
 * \param T5 Travolta table
 * 
 * \ingroup RowOperations
 */

static inline void mzed_process_rows6(mzed_t *M, const rci_t startrow, const rci_t endrow, const rci_t startcol,
                                      const tt_mzed_t *T0, const tt_mzed_t *T1, const tt_mzed_t *T2, 
                                      const tt_mzed_t *T3, const tt_mzed_t *T4, const tt_mzed_t *T5) {
  mzd_process_rows6(M->x, startrow, endrow, startcol*M->w, 6*M->w, T0->T->x, T0->L, T1->T->x, T1->L, T2->T->x, T2->L, T3->T->x, T3->L, T4->T->x, T4->L, T5->T->x, T5->L);
}


#endif //M4RIE_TRAVOLTA_H

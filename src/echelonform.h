#ifndef M4RIE_ECHELONFORM_H
#define M4RIE_ECHELONFORM_H

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

#include "gf2e_matrix.h"
#include "bitslice.h"
#include "conversion.h"

/**
 * \brief Compute row echelon forms.
 * 
 * Compute the (reduced) row echelon form of the matrix A.  If full=0,
 * then return the reduced REF.
 *
 * \param A Matrix
 * \param full REF or RREF.
 *
 * \ingroup Echelon
 */

rci_t mzed_echelonize(mzed_t *A, int full);

/**
 * \brief Compute row echelon forms.
 * 
 * Compute the (reduced) row echelon form of the matrix A.  If full=0,
 * then return the reduced REF. This function reduces echelon forms to
 * PLE (or PLUQ) decomposition.
 *
 * \param A Matrix
 * \param full REF or RREF.
 *
 * \ingroup Echelon
 */

rci_t mzd_slice_echelonize_ple(mzd_slice_t *A, int full);

/**
 * \brief Compute row echelon forms.
 * 
 * Compute the (reduced) row echelon form of the matrix A.  If full=0,
 * then return the reduced REF.
 *
 * \param A Matrix
 * \param full REF or RREF.
 *
 * \ingroup Echelon
 */

#define mzd_slice_echelonize mzd_slice_echelonize_ple

/**
 * \brief Compute row echelon forms.
 * 
 * Compute the (reduced) row echelon form of the matrix A.  If full=0,
 * then return the reduced REF. This function reduces echelon forms to
 * PLE (or PLUQ) decomposition.
 *
 * \param A Matrix
 * \param full REF or RREF.
 *
 * \ingroup Echelon
 */

static inline rci_t mzed_echelonize_ple(mzed_t *A, int full) {
  mzd_slice_t *a = mzed_slice(NULL, A);
  rci_t r = mzd_slice_echelonize_ple(a, full);
  mzed_cling(A, a);
  mzd_slice_free(a);
  return r;
}


#endif //M4RIE_ECHELONFORM_H

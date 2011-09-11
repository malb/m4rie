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

size_t mzed_echelonize(mzed_t *A, int full);


#endif //M4RIE_ECHELONFORM_H

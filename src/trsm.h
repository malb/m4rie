/**
 * \file trsm.h
 * \brief Triangular System Solving with Matrices (TRSM).
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

#ifndef TRSM_H
#define TRSM_H

/******************************************************************************
*
*            M4RIE: Linear Algebra over GF(2^e)
*
*    Copyright (C) 2011 Martin Albrecht <martinralbrecht@googlemail.com>
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

#include <gf2e_matrix.h>
#include <bitslice.h>

#define MZED_TRSM_CUTOFF 64

void mzed_trsm_lower_left_naive(const mzed_t *L, mzed_t *B);

void mzed_trsm_lower_left(const mzed_t *L, mzed_t *B);


#endif //TRSM_H

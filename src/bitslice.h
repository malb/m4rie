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

void mzed_unpack2(mzd_t *A0, mzd_t *A1, const mzed_t *A);
void mzed_pack2(mzed_t *A, const mzd_t *A0, const mzd_t *A1);

mzed_t *mzed_mul_karatsuba(mzed_t *C, const mzed_t *A, const mzed_t *B);

mzed_t *_mzed_mul_karatsuba2(mzed_t *C, const mzed_t *A, const mzed_t *B);

#endif //BITSLICE_H

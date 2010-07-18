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

void mzed_make_table(mzed_t *T, mzed_t *A, size_t r, size_t c, gf2e *ff);
size_t mzed_echelonize_travolta(mzed_t *A, int full);

mzed_t *mzed_mul_travolta(mzed_t *C, mzed_t *A, mzed_t *B);
mzed_t *mzed_addmul_travolta(mzed_t *C, mzed_t *A, mzed_t *B);

mzed_t *_mzed_mul_travolta(mzed_t *C, mzed_t *A, mzed_t *B);


#endif //TRAVOLTA_H

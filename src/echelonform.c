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

#include "echelonform.h"
#include "travolta.h"

size_t mzed_echelonize(mzed_t *A, int full) {
  if (A->finite_field->degree > A->nrows/2) {
    return mzed_echelonize_naive(A, full);
  } else {
    return mzed_echelonize_travolta(A, full);
  }
}

/**
 * \file ple.h
 * \brief L*E = P*A decomposition.
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

#ifndef M4RIE_PLE_H
#define M4RIE_PLE_H

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

#include <m4ri/m4ri.h>
#include "gf2e_matrix.h"
#include "bitslice.h"

rci_t mzed_ple_naive(mzed_t *A, mzp_t *P, mzp_t *Q);

rci_t _mzd_slice_ple(mzd_slice_t *A, mzp_t *P, mzp_t *Q, rci_t cutoff);

static inline rci_t mzd_slice_ple(mzd_slice_t *A, mzp_t *P, mzp_t *Q) {
  return _mzd_slice_ple(A, P, Q, 0);
}

rci_t _mzed_ple(mzed_t *A, mzp_t *P, mzp_t *Q, rci_t cutoff);

static inline rci_t mzed_ple(mzed_t *A, mzp_t *P, mzp_t *Q) {
  return _mzed_ple(A, P, Q, 0);
}


#define __M4RIE_PLE_CUTOFF 256*256

#endif //M4RIE_PLE_H

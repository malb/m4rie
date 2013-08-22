/**
 * \file blm.h
 *
 * \brief Bilinear Maps on Matrices over GF(2).
 *
 * This is used to realise mzd_poly_t multiplication.
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

#ifndef M4RIE_BLM_H
#define M4RIE_BLM_H

/******************************************************************************
*
*            M4RIE: Linear Algebra over GF(2^e)
*
*    Copyright (C) 2013 Martin Albrecht <martinralbrecht@googlemail.com>
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
#include "m4rie/gf2e.h"

/**
 * \brief Bilinear Maps on Matrices over GF(2).
 *
 * Encodes the bilinear map H*((F*A) x (G*B)) where A,B are vectors of mzd_t, "*" is matrix-vector
 * multiplication and "x" is pointwise multiplication.
 */

typedef struct {
  mzd_t *H;
  mzd_t *F;
  mzd_t *G;
} blm_t;

/**
 * costs[i] = cost of multiplying two polynomials of length i over \GF2.
 */

extern const int costs[17];


/**
 * \brief Apply binlinear map f.
 */

void _mzd_ptr_apply_blm(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B, const blm_t *f);

/**
 * Compute H, F, G such that vec(c) = H*(F*vec(a) x G*vec(b))
 * with poly(c) = poly(a)*poly(b), 
 * deg(poly(a)) = a_ncols -1, deg(poly(b)) = b_ncols -1 and "x" being pointwise multiplication
 *
 * This is realised by a multi-modular computation modulo the primes up to degree deg (which may not
 * be irreducible polynomials, but merely co-prime).
 */

blm_t *blm_init_multimod(const deg_t f_ncols, const deg_t g_ncols, const deg_t deg, const int *primes);

/**
 * Compute H 
 */

blm_t *_blm_finish_polymult(blm_t *f);

void blm_free(blm_t *f);


#endif //M4RIE_BLM_H

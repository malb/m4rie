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
 * \brief We consider at most polynomials of degree M4RIE_MAX_DEGREE in CRT.
 */

#define M4RIE_CRT_LEN (M4RIE_MAX_DEGREE + 1)

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
 * Return the multiplication cost of the multiplication scheme p
 */

static inline int blm_cost_crt(const int p[M4RIE_CRT_LEN]) {
  int cost = costs[p[0]];
  for(deg_t d=1; d<M4RIE_CRT_LEN; d++)
    cost += costs[d] * p[d];
  return cost;
}

/**
 * Find a list of co-prime polynomials p_i such that deg(prod(p_i)) >= f_len*g_len-1.
 *
 * We store the number of polynomials of degree d in p[d]. We store the degree w of (x-infinity)^w
 *  in p[0].
 */

int *crt_init(const deg_t f_len, const deg_t g_len);

/**
 * Compute H, F, G such that vec(c) = H*(F*vec(a) x G*vec(b))
 * with poly(c) = poly(a)*poly(b), 
 * deg(poly(a)) = a_ncols -1, deg(poly(b)) = b_ncols -1 and "x" being pointwise multiplication
 *
 * This is realised by a multi-modular computation modulo the primes up to degree deg (which may not
 * be irreducible polynomials, but merely co-prime).
 */

blm_t *blm_init_crt(const deg_t f_ncols, const deg_t g_ncols, const int *p);

/**
 * Given F and G compute H.
 *
 * \param f Bilinear Map with F and G already computed.
 */

blm_t *_blm_finish_polymult(blm_t *f);

/**
 * Free bilinear map f.
 */

void blm_free(blm_t *f);


void _mzd_ptr_apply_blm_djb(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B, const blm_t *f);

#endif //M4RIE_BLM_H

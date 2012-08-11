/**
 * \file gf2e.h
 *
 * \brief \GF2E
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

#ifndef M4RIE_GF2E_H
#define M4RIE_GF2E_H

/******************************************************************************
*
*            M4RIE: Linear Algebra over GF(2^e)
*
*    Copyright (C) 2010,2011 Martin Albrecht <martinralbrecht@googlemail.com>
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
#include <m4rie/gf2x.h>

#define M4RIE_MAX_DEGREE 10

/**
 * \brief \GF2E
 */

typedef struct {
  unsigned int degree; /**< The degree \e. */
  word minpoly;   /**<  Irreducible polynomial of degree \e. */

  word *pow_gen;   /**< pow_gen[i] holds \f$a^i / <f>\f$ for \f$a\f$ a generator of this field.  */

  word **mul;   /**<
                 * mul[a][b] holds \f$ a \cdot b\f$.
                 * \warning this entry will disappear in future releases. */
} gf2e;

/**
 * Create finite field from minimal polynomial
 *
 * \param minpoly Polynomial represented as series of bits.
 */

gf2e *gf2e_init(const word minpoly);

/**
 * Generate gf2e::pow_gen.
 *
 * \param ff Finite field.
 */

static inline void gf2e_make_pow_gen(gf2e *ff) {
  unsigned int n = 2*ff->degree-1;
  word *m = (word*)m4ri_mm_malloc( n * sizeof(word));
  for(unsigned int i=0; i<n; i++) {
    m[i] = 1<<i;
    for(unsigned int j=i; j>=ff->degree; j--) {
      if (m[i] & 1<<j)
        m[i] ^= ff->minpoly<<(j - ff->degree);
    }
  }
  ff->pow_gen = m;
}

/**
 * Free ff
 *
 * \param ff Finite field.
 */

void gf2e_free(gf2e *ff);

/**
 * \brief a^(-1) % minpoly
 */

static inline word gf2e_inv(const gf2e *ff, word a) {
  return gf2x_invmod(a, ff->minpoly, ff->degree);
}

/**
 * \brief a*b in \GF2E
 */
static inline word gf2e_mul(const gf2e *ff, const word a, const word b) {
  return ff->mul[a][b];
}

/**
 * Return the width used for storing elements of ff
 *
 * \param ff Finite field.
 */

static inline size_t gf2e_degree_to_w(const gf2e *ff) {
  switch(ff->degree) {
  case 2:
    return 2;
  case  3:
  case  4:
    return 4;
  case  5:
  case  6:
  case  7:
  case  8:
    return 8;
  case  9:
  case 10:
  case 11:
  case 12:
  case 13:
  case 14:
  case 15:
  case 16:
    return 16;
  default:
    m4ri_die("degree %d not supported.\n",ff->degree);
  }
  return 0;
}

/**
 * Compute all multiples by a of vectors fitting into 16 bits.
 *
 * \param ff Finite field.
 * \param a Finite field element.
 */

static inline word *gf2e_t16_init(const gf2e *ff, const word a) {
  word *mul = (word*)m4ri_mm_calloc(1<<16, sizeof(word));

  const unsigned int w = gf2e_degree_to_w(ff);
  const word mask_w = (1<<w)-1;

  /**
   * @todo: this is a bit of overkill, we could do better
   */
  for(word i=0; i<1<<16; i++) {
    switch(w) {
    case 2:
      mul[i]  = gf2e_mul(ff, a, ((i>>0)&mask_w))<<0 | gf2e_mul(ff, a, ((i>> 2)&mask_w))<< 2 | gf2e_mul(ff, a, ((i>> 4)&mask_w))<< 4 | gf2e_mul(ff, a, ((i>> 6)&mask_w))<< 6;
      mul[i] |= gf2e_mul(ff, a, ((i>>8)&mask_w))<<8 | gf2e_mul(ff, a, ((i>>10)&mask_w))<<10 | gf2e_mul(ff, a, ((i>>12)&mask_w))<<12 | gf2e_mul(ff, a, ((i>>14)&mask_w))<<14;
      break;
    case 4:
      mul[i]  = gf2e_mul(ff, a, (i&mask_w)) | gf2e_mul(ff, a, ((i>>4)&mask_w))<<4 | gf2e_mul(ff, a, ((i>>8)&mask_w))<<8 | gf2e_mul(ff, a, ((i>>12)&mask_w))<<12;
      break;
    case 8:
      mul[i]  = gf2e_mul(ff, a, (i&mask_w)) | gf2e_mul(ff, a, ((i>>8)&mask_w))<<8;
      break;
    case 16:
      mul[i]  = gf2e_mul(ff, a, (i&mask_w));
      break;
    };
  }
  return mul;
}

/**
 * \brief Free multiplication table.
 *
 * \param mul Multiplication table
 */

static inline void gf2e_t16_free(word *mul) {
  m4ri_mm_free(mul);
}

extern const word* irreducible_polynomials[17];

#endif //M4RIE_GF2E_H

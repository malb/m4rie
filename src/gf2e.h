/**
 * \file gf2e.h
 *
 * \brief \GF2E
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

#ifndef M4RIE_GF2E_H
#define M4RIE_GF2E_H

#include <m4ri/m4ri.h>

#define M4RIE_MAX_DEGREE 10

/**
 * \brief Multiplication in GF(2)[x] where a,b have a degree smaller than d <= 16.
 */

static inline word gf2x_mul(const word a, const word b, unsigned int d) {
  /* double check that integer multiplication is indeed not horribly slow here */
  word res = 0;

  switch(d) {
  case 16: res ^= ((a&(1<<15))>>15) * (b<<15);
  case 15: res ^= ((a&(1<<14))>>14) * (b<<14);
  case 14: res ^= ((a&(1<<13))>>13) * (b<<13);
  case 13: res ^= ((a&(1<<12))>>12) * (b<<12);
  case 12: res ^= ((a&(1<<11))>>11) * (b<<11);
  case 11: res ^= ((a&(1<<10))>>10) * (b<<10);
  case 10: res ^= ((a&(1<< 9))>> 9) * (b<< 9);
  case  9: res ^= ((a&(1<< 8))>> 8) * (b<< 8);
  case  8: res ^= ((a&(1<< 7))>> 7) * (b<< 7);
  case  7: res ^= ((a&(1<< 6))>> 6) * (b<< 6);
  case  6: res ^= ((a&(1<< 5))>> 5) * (b<< 5);
  case  5: res ^= ((a&(1<< 4))>> 4) * (b<< 4);
  case  4: res ^= ((a&(1<< 3))>> 3) * (b<< 3);
  case  3: res ^= ((a&(1<< 2))>> 2) * (b<< 2);
  case  2: res ^= ((a&(1<< 1))>> 1) * (b<< 1);
  case  1: res ^= ((a&(1<< 0))>> 0) * (b<< 0);
    break;
  default:
    m4ri_die("degree %d not supported.",d);
  }
  return res;
}

/**
 * \brief Degree of elements in GF(2)[x].
 */

static inline unsigned int gf2x_deg(word a) {
  unsigned int degree = 0;
  if( (a & 0xffffffff00000000ULL) != 0) { degree += 32; a>>=32; }
  if( (a &         0xffff0000ULL) != 0) { degree += 16; a>>=16; }
  if( (a &             0xff00ULL) != 0) { degree +=  8; a>>= 8; }
  if( (a &               0xf0ULL) != 0) { degree +=  4; a>>= 4; }
  if( (a &                0xcULL) != 0) { degree +=  2; a>>= 2; }
  if( (a &                0x2ULL) != 0) { degree +=  1; a>>= 1; }
  return degree;
}

/**
 * \brief Division in GF(2)[x].
 */

static inline word gf2x_div(word a, word b) {
  word res = 0;
  while(a >= b) {
    unsigned int diff = gf2x_deg(a) - gf2x_deg(b);
    res |= __M4RI_TWOPOW(diff);
    a ^= b<<diff;
  }
  return res;
}

/**
 * \brief Remainders in GF(2)[x].
 */

static inline word gf2x_mod(word a, word b) {
  word res = 0;
  while(a >= b) {
    unsigned int diff = gf2x_deg(a) - gf2x_deg(b);
    res |= __M4RI_TWOPOW(diff);
    a ^= b<<diff;
  }
  return a;
}

/**
 * \brief a^(-1) % b.
 */

static inline word gf2x_invmod(word a, word b, unsigned int d) {
  word x = 0;
  word lastx = 1;
  word y = 1;
  word lasty = 0;

  word tmp = 0;

  while (b != 0) {
    word quotient = gf2x_div(a,b);
    tmp = b; b = gf2x_mod(a, b); a = tmp;
    tmp = x; x = lastx ^ gf2x_mul(quotient, x, d); lastx = tmp;
    tmp = y; y = lasty ^ gf2x_mul(quotient, y, d); lasty = tmp;
  }
  return lastx;
}

/**
 * \brief \GF2E
 */

typedef struct {
  unsigned int degree; /**< The degree \e. */
  word minpoly;   /**<  Irreducible polynomial of degree \e. */

  word *inv; /**< inv[a] holds \f$a^{-1}\f$. */
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
 * Free ff
 *
 * \param ff Finite field.
 */

void gf2e_free(gf2e *ff);

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
 * Compute all multiples by a of vectors fitting into 16 bits.
 *
 * \param ff Finite field.
 * \param a Finite field element.
 */

static inline word *gf2e_t16_init(const gf2e *ff, const word a) {
  word *mul = (word*)m4ri_mm_calloc(1<<16, sizeof(word));

  const unsigned int w = gf2e_degree_to_w(ff);
  const word mask_w = (1<<w)-1;
  const word *x = ff->mul[a];

  /**
   * @todo: this is a bit of overkill, we could do better
   */
  for(word i=0; i<1<<16; i++) {
    switch(w) {
    case 2:
      mul[i]  = x[(i&mask_w)] | x[((i>>2)&mask_w)]<<2 | x[((i>>4)&mask_w)]<<4 | x[((i>>6)&mask_w)]<<6;
      mul[i] |= x[((i>>8)&mask_w)]<<8 | x[((i>>10)&mask_w)]<<10 | x[((i>>12)&mask_w)]<<12 | x[((i>>14)&mask_w)]<<14;
      break;
    case 4:
      mul[i]  = x[(i&mask_w)] | x[((i>>4)&mask_w)]<<4 | x[((i>>8)&mask_w)]<<8 | x[((i>>12)&mask_w)]<<12;
      break;
    case 8:
      mul[i]  = x[(i&mask_w)] | x[((i>>8)&mask_w)]<<8;
      break;
    case 16:
      mul[i]  = x[(i&mask_w)];
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

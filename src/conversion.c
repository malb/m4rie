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

#include "conversion.h"

static const word x80008000 = 0x8000800080008000ULL;
static const word x80808080 = 0x8080808080808080ULL;
static const word x88888888 = 0x8888888888888888ULL;
static const word xaaaaaaaa = 0xaaaaaaaaaaaaaaaaULL;
static const word xcccccccc = 0xccccccccccccccccULL;
static const word xc0c0c0c0 = 0xc0c0c0c0c0c0c0c0ULL;
static const word xf0f0f0f0 = 0xf0f0f0f0f0f0f0f0ULL;
static const word xff00ff00 = 0xff00ff00ff00ff00ULL;
static const word xffff0000 = 0xffff0000ffff0000ULL;
static const word xffffffff = 0xffffffff00000000ULL;
static const word x__left04 = 0xf000000000000000ULL;
static const word x__left08 = 0xff00000000000000ULL;
static const word x__left16 = 0xffff000000000000ULL;
static const word x__left32 = 0xffffffff00000000ULL;

static inline word word_slice_64_02_l(word a) {
  a = (a & xcccccccc) | (a & xcccccccc>> 2)<< 1;
  a = (a & xf0f0f0f0) | (a & xf0f0f0f0>> 4)<< 2;
  a = (a & xff00ff00) | (a & xff00ff00>> 8)<< 4;
  a = (a & xffff0000) | (a & xffff0000>>16)<< 8;
  a = (a & xffffffff) | (a & xffffffff>>32)<<16;
  return a;
}

static inline word word_slice_64_04_l(word a) {
  a = (a & xf0f0f0f0) | (a & xf0f0f0f0>> 4)<< 3;
  a = (a & xff00ff00) | (a & xff00ff00>> 8)<< 6;
  a = (a & xffff0000) | (a & xffff0000>>16)<<12;
  a = (a & xffffffff) | (a & xffffffff>>32)<<24;
  return a;
}

static inline word word_slice_64_08_l(word a) {
  a = (a & xff00ff00) | (a & xff00ff00>> 8)<< 7;
  a = (a & xffff0000) | (a & xffff0000>>16)<<14;
  a = (a & xffffffff) | (a & xffffffff>>32)<<28;
  return a;
}

static inline word word_slice_64_16_l(word a) {
  a = (a & xffff0000) | (a & xffff0000>>16)<<15;
  a = (a & xffffffff) | (a & xffffffff>>32)<<30;
  return a;
}

static inline word word_cling_64_02_l(word a) {
  a = (a & xffff0000 & x__left32) | (a & (xffff0000>>16) & x__left32)>>16;
  a = (a & xff00ff00) | (a & xff00ff00>> 8)>> 8;
  a = (a & xf0f0f0f0) | (a & xf0f0f0f0>> 4)>> 4;
  a = (a & xcccccccc) | (a & xcccccccc>> 2)>> 2;
  a = (a & xaaaaaaaa) | (a & xaaaaaaaa>> 1)>> 1;
  return a;
}

static inline word word_cling_64_04_l(word a) {
  a = (a & xff00ff00 & x__left16) | (a & (xff00ff00>> 8) & x__left16)>>24;
  a = (a & xf0f0f0f0) | (a & xf0f0f0f0>> 4)>>12;
  a = (a & xcccccccc) | (a & xcccccccc>> 2)>> 6;
  a = (a & xaaaaaaaa) | (a & xaaaaaaaa>> 1)>> 3;
  return a;
}

static inline word word_cling_64_08_l(word a) {
  a = (a & xf0f0f0f0 & x__left08) | (a & xf0f0f0f0>> 4  & x__left08)>>28;
  a = (a & xcccccccc) | (a & xcccccccc>> 2)>>14;
  a = (a & xaaaaaaaa) | (a & xaaaaaaaa>> 1)>> 7;
  return a;
}

static inline word word_cling_64_16_l(word a) {
  a = (a & xcccccccc & x__left04) | (a & xcccccccc>> 2  & x__left04)>>30;
  a = (a & xaaaaaaaa) | (a & xaaaaaaaa>> 1)>> 15;
  return a;
}

mzd_slice_t *mzed_slice(mzd_slice_t *A, const mzed_t *Z) {
  if (A == NULL) {
    assert(Z->x->offset == 0);
    A = mzd_slice_init(Z->finite_field, Z->nrows, Z->ncols);
  } else {
    assert((Z->x->offset | A->x[0]->offset) == 0);
    mzd_slice_set_ui(A, 0);
  }

  switch(Z->finite_field->degree) {
  case  2: return _mzed_slice2(A,Z);

  case  3: return _mzed_slice4(A,Z);
  case  4: return _mzed_slice4(A,Z);

  case  5: return _mzed_slice8(A,Z);
  case  6: return _mzed_slice8(A,Z);
  case  7: return _mzed_slice8(A,Z);
  case  8: return _mzed_slice8(A,Z);

  case  9: return _mzed_slice16(A,Z);
  case 10: return _mzed_slice16(A,Z);
  case 11: return _mzed_slice16(A,Z);
  case 12: return _mzed_slice16(A,Z);
  case 13: return _mzed_slice16(A,Z);
  case 14: return _mzed_slice16(A,Z);
  case 15: return _mzed_slice16(A,Z);
  case 16: return _mzed_slice16(A,Z);
  default:
    m4ri_die("slicing not implemented for this degree");
  }
  return A;
}

mzed_t *mzed_cling(mzed_t *A, const mzd_slice_t *Z) {
  if (A == NULL) {
    assert(Z->x[0]->offset == 0);
    A = mzed_init(Z->finite_field, Z->nrows, Z->ncols);
  }
  else {
    assert((A->x->offset | Z->x[0]->offset) == 0);
    mzed_set_ui(A, 0);
  }

  switch(Z->finite_field->degree) {
  case  2: return _mzed_cling2(A,Z);

  case  3: return _mzed_cling4(A,Z);
  case  4: return _mzed_cling4(A,Z);

  case  5: return _mzed_cling8(A,Z);
  case  6: return _mzed_cling8(A,Z);
  case  7: return _mzed_cling8(A,Z);
  case  8: return _mzed_cling8(A,Z);

  case  9: return _mzed_cling16(A,Z);
  case 10: return _mzed_cling16(A,Z);
  case 11: return _mzed_cling16(A,Z);
  case 12: return _mzed_cling16(A,Z);
  case 13: return _mzed_cling16(A,Z);
  case 14: return _mzed_cling16(A,Z);
  case 15: return _mzed_cling16(A,Z);
  case 16: return _mzed_cling16(A,Z);
  default:
    m4ri_die("clinging not implemented for this degree");
  }
  return A;
}

mzd_slice_t *_mzed_slice2(mzd_slice_t *T, const mzed_t *F) {
  assert(T && (T->depth >= 2));
  size_t j, j2 = 0;

  const word bitmask_end = __M4RI_LEFT_BITMASK((T->x[0]->offset + T->ncols) % m4ri_radix);
  register word r0,r1,r2,r3;

  if (mzed_is_zero(F))
    return T;

  for(size_t i=0; i<T->nrows; i++) {
    word *t0 = T->x[0]->rows[i];
    word *t1 = T->x[1]->rows[i];
    const word *f  = F->x->rows[i];

    /* bulk of work */
    for(j=0, j2=0; j+2 < F->x->width; j+=2,j2++) {
      r0 = f[j+0], r1 = f[j+1];
      r2 = word_slice_64_02_l(r0<<1 & xaaaaaaaa);
      r3 = word_slice_64_02_l(r1<<1 & xaaaaaaaa);
      t0[j2] = r3 | (r2>>32);

      r2 = word_slice_64_02_l(r0<<0 & xaaaaaaaa);
      r3 = word_slice_64_02_l(r1<<0 & xaaaaaaaa);
      t1[j2] = r3 | (r2>>32);
    }

    switch(F->x->width - j) {
    case 2:
      r0 = f[j+0]; r1 = f[j+1];

      r2 = word_slice_64_02_l(r0<<1 & xaaaaaaaa);
      r3 = word_slice_64_02_l(r1<<1 & xaaaaaaaa);
      t0[j2] &= ~bitmask_end;
      t0[j2] |= (r3 | (r2>>32)) & bitmask_end;

      r2 = word_slice_64_02_l(r0<<0 & xaaaaaaaa);
      r3 = word_slice_64_02_l(r1<<0 & xaaaaaaaa);
      t1[j2] &= ~bitmask_end;
      t1[j2] |= (r3 | (r2>>32)) & bitmask_end;
      break;
    case 1:
      r0 = f[j+0];

      r2 = word_slice_64_02_l(r0<<1 & xaaaaaaaa);
      t0[j2] &= ~bitmask_end;
      t0[j2] |= (r2>>32) & bitmask_end;

      r2 = word_slice_64_02_l(r0<<0 & xaaaaaaaa);
      t1[j2] &= ~bitmask_end;
      t1[j2] |= (r2>>32) & bitmask_end;
      break;
    default:
      m4ri_die("impossible");
    }
  }

  return T;
}

mzed_t *_mzed_cling2(mzed_t *T, const mzd_slice_t *F) {
  size_t j,j2 = 0;
  register word tmp;

  const word bitmask_end = __M4RI_LEFT_BITMASK((T->x->offset + T->x->ncols) % m4ri_radix);

  if (mzd_slice_is_zero(F))
    return T;

  for(size_t i=0; i<T->nrows; i++) {
    const word *f0 = F->x[0]->rows[i];
    const word *f1 = F->x[1]->rows[i];
    word *t  = T->x->rows[i];

    for(j=0, j2=0; j+2 < T->x->width; j+=2, j2++) {
      t[j+0] = (word_cling_64_02_l(f0[j2]<<32)>>1) | (word_cling_64_02_l(f1[j2]<<32)>>0);
      t[j+1] = (word_cling_64_02_l(f0[j2]<< 0)>>1) | (word_cling_64_02_l(f1[j2]<< 0)>>0);
    }
    switch(T->x->width - j) {
    case 2:
      tmp    = (word_cling_64_02_l(f0[j2]<< 0)>>1) | (word_cling_64_02_l(f1[j2]<< 0)>>0);
      t[j+0] = (word_cling_64_02_l(f0[j2]<<32)>>1) | (word_cling_64_02_l(f1[j2]<<32)>>0);
      t[j+1] = (t[j+1] & ~bitmask_end) | (tmp & bitmask_end);
      break;
    case 1:
      tmp    = (word_cling_64_02_l(f0[j2]<<32)>>1) | (word_cling_64_02_l(f1[j2]<<32)>>0);
      t[j+0] = (t[j+0] & ~bitmask_end) | (tmp & bitmask_end);
      break;
    }
  }
  return T;
}

mzd_slice_t *_mzed_slice4(mzd_slice_t *T, const mzed_t *F) {
  assert(T && (T->depth == 3 || T->depth == 4) && T->x[0]->offset == 0);
  size_t j, j2 = 0;
  register word r0,r1,r2,r3 = 0;

  const word bitmask_end = __M4RI_LEFT_BITMASK((T->x[0]->offset + T->ncols) % m4ri_radix);

  if (mzed_is_zero(F))
    return T;

  if (T->depth == 3) {
    for(size_t i=0; i<T->nrows; i++) {
      word *t0 = T->x[0]->rows[i];
      word *t1 = T->x[1]->rows[i];
      word *t2 = T->x[2]->rows[i];
      const word const *f  = F->x->rows[i];

      /* bulk of work */
      for(j=0, j2=0; j+4 < F->x->width; j+=4,j2++) {
        t0[j2] = word_slice_64_04_l(f[j+0]<<3 & x88888888)>>48 | word_slice_64_04_l(f[j+1]<<3 & x88888888)>>32 \
          |      word_slice_64_04_l(f[j+2]<<3 & x88888888)>>16 | word_slice_64_04_l(f[j+3]<<3 & x88888888)>> 0;
        t1[j2] = word_slice_64_04_l(f[j+0]<<2 & x88888888)>>48 | word_slice_64_04_l(f[j+1]<<2 & x88888888)>>32 \
          |      word_slice_64_04_l(f[j+2]<<2 & x88888888)>>16 | word_slice_64_04_l(f[j+3]<<2 & x88888888)>> 0;
        t2[j2] = word_slice_64_04_l(f[j+0]<<1 & x88888888)>>48 | word_slice_64_04_l(f[j+1]<<1 & x88888888)>>32 \
          |      word_slice_64_04_l(f[j+2]<<1 & x88888888)>>16 | word_slice_64_04_l(f[j+3]<<1 & x88888888)>> 0;
      }
      r0 = r1 = r2 = 0;
      switch(F->x->width - j) {
      case 4:
        r0 |= word_slice_64_04_l(f[j+3]<<3 & x88888888)>> 0;
        r1 |= word_slice_64_04_l(f[j+3]<<2 & x88888888)>> 0;
        r2 |= word_slice_64_04_l(f[j+3]<<1 & x88888888)>> 0;
      case 3:
        r0 |= word_slice_64_04_l(f[j+2]<<3 & x88888888)>>16;
        r1 |= word_slice_64_04_l(f[j+2]<<2 & x88888888)>>16;
        r2 |= word_slice_64_04_l(f[j+2]<<1 & x88888888)>>16;
      case 2:
        r0 |= word_slice_64_04_l(f[j+1]<<3 & x88888888)>>32;
        r1 |= word_slice_64_04_l(f[j+1]<<2 & x88888888)>>32;
        r2 |= word_slice_64_04_l(f[j+1]<<1 & x88888888)>>32;
      case 1:
        r0 |= word_slice_64_04_l(f[j+0]<<3 & x88888888)>>48;
        r1 |= word_slice_64_04_l(f[j+0]<<2 & x88888888)>>48;
        r2 |= word_slice_64_04_l(f[j+0]<<1 & x88888888)>>48;
        break;
      default:
        m4ri_die("impossible");
      }
      t0[j2] |= r0 & bitmask_end;
      t1[j2] |= r1 & bitmask_end;
      t2[j2] |= r2 & bitmask_end;
    }
  } else {
    for(size_t i=0; i<T->nrows; i++) {
      word *t0 = T->x[0]->rows[i];
      word *t1 = T->x[1]->rows[i];
      word *t2 = T->x[2]->rows[i];
      word *t3 = T->x[3]->rows[i];
      const word const *f  = F->x->rows[i];

      /* bulk of work */
      for(j=0, j2=0; j+4 < F->x->width; j+=4,j2++) {
        t0[j2] = word_slice_64_04_l(f[j+0]<<3 & x88888888)>>48 | word_slice_64_04_l(f[j+1]<<3 & x88888888)>>32 \
          |      word_slice_64_04_l(f[j+2]<<3 & x88888888)>>16 | word_slice_64_04_l(f[j+3]<<3 & x88888888)>> 0;
        t1[j2] = word_slice_64_04_l(f[j+0]<<2 & x88888888)>>48 | word_slice_64_04_l(f[j+1]<<2 & x88888888)>>32 \
          |      word_slice_64_04_l(f[j+2]<<2 & x88888888)>>16 | word_slice_64_04_l(f[j+3]<<2 & x88888888)>> 0;
        t2[j2] = word_slice_64_04_l(f[j+0]<<1 & x88888888)>>48 | word_slice_64_04_l(f[j+1]<<1 & x88888888)>>32 \
          |      word_slice_64_04_l(f[j+2]<<1 & x88888888)>>16 | word_slice_64_04_l(f[j+3]<<1 & x88888888)>> 0;
        t3[j2] = word_slice_64_04_l(f[j+0]<<0 & x88888888)>>48 | word_slice_64_04_l(f[j+1]<<0 & x88888888)>>32 \
          |      word_slice_64_04_l(f[j+2]<<0 & x88888888)>>16 | word_slice_64_04_l(f[j+3]<<0 & x88888888)>> 0;
      }
      r0 = r1 = r2 = r3 = 0;
      switch(F->x->width - j) {
      case 4:
        r0 |= word_slice_64_04_l(f[j+3]<<3 & x88888888)>> 0;
        r1 |= word_slice_64_04_l(f[j+3]<<2 & x88888888)>> 0;
        r2 |= word_slice_64_04_l(f[j+3]<<1 & x88888888)>> 0;
        r3 |= word_slice_64_04_l(f[j+3]<<0 & x88888888)>> 0;
      case 3:
        r0 |= word_slice_64_04_l(f[j+2]<<3 & x88888888)>>16;
        r1 |= word_slice_64_04_l(f[j+2]<<2 & x88888888)>>16;
        r2 |= word_slice_64_04_l(f[j+2]<<1 & x88888888)>>16;
        r3 |= word_slice_64_04_l(f[j+2]<<0 & x88888888)>>16;
      case 2:
        r0 |= word_slice_64_04_l(f[j+1]<<3 & x88888888)>>32;
        r1 |= word_slice_64_04_l(f[j+1]<<2 & x88888888)>>32;
        r2 |= word_slice_64_04_l(f[j+1]<<1 & x88888888)>>32;
        r3 |= word_slice_64_04_l(f[j+1]<<0 & x88888888)>>32;
      case 1:
        r0 |= word_slice_64_04_l(f[j+0]<<3 & x88888888)>>48;
        r1 |= word_slice_64_04_l(f[j+0]<<2 & x88888888)>>48;
        r2 |= word_slice_64_04_l(f[j+0]<<1 & x88888888)>>48;
        r3 |= word_slice_64_04_l(f[j+0]<<0 & x88888888)>>48;
        break;
      default:
        m4ri_die("impossible");
      }
      t0[j2] |= r0 & bitmask_end;
      t1[j2] |= r1 & bitmask_end;
      t2[j2] |= r2 & bitmask_end;
      t3[j2] |= r3 & bitmask_end;
    }
  }
  return T;
}

mzed_t *_mzed_cling4(mzed_t *T, const mzd_slice_t *F) {
  size_t j,j2 = 0;

  const word bitmask_end = __M4RI_LEFT_BITMASK((T->x->offset + T->x->ncols) % m4ri_radix);

  if (mzd_slice_is_zero(F))
    return T;

  if (F->finite_field->degree == 4) {
    for(rci_t i=0; i<T->nrows; i++) {
      const word *f0 = F->x[0]->rows[i];
      const word *f1 = F->x[1]->rows[i];
      const word *f2 = F->x[2]->rows[i];
      const word *f3 = F->x[3]->rows[i];
      word *t  = T->x->rows[i];

      for(j=0, j2=0; j+4 < T->x->width; j+=4, j2++) {
        t[j+0] = (word_cling_64_04_l(f0[j2]<<48)>>3) | (word_cling_64_04_l(f1[j2]<<48)>>2) | (word_cling_64_04_l(f2[j2]<<48)>>1) | (word_cling_64_04_l(f3[j2]<<48)>>0);
        t[j+1] = (word_cling_64_04_l(f0[j2]<<32)>>3) | (word_cling_64_04_l(f1[j2]<<32)>>2) | (word_cling_64_04_l(f2[j2]<<32)>>1) | (word_cling_64_04_l(f3[j2]<<32)>>0);
        t[j+2] = (word_cling_64_04_l(f0[j2]<<16)>>3) | (word_cling_64_04_l(f1[j2]<<16)>>2) | (word_cling_64_04_l(f2[j2]<<16)>>1) | (word_cling_64_04_l(f3[j2]<<16)>>0);
        t[j+3] = (word_cling_64_04_l(f0[j2]<< 0)>>3) | (word_cling_64_04_l(f1[j2]<< 0)>>2) | (word_cling_64_04_l(f2[j2]<< 0)>>1) | (word_cling_64_04_l(f3[j2]<< 0)>>0);
      }

      register word tmp=0;
      switch(T->x->width - j) {
      case 4:
        t[j+0] = (word_cling_64_04_l(f0[j2]<<48)>>3) | (word_cling_64_04_l(f1[j2]<<48)>>2) | (word_cling_64_04_l(f2[j2]<<48)>>1) | (word_cling_64_04_l(f3[j2]<<48)>>0);
        t[j+1] = (word_cling_64_04_l(f0[j2]<<32)>>3) | (word_cling_64_04_l(f1[j2]<<32)>>2) | (word_cling_64_04_l(f2[j2]<<32)>>1) | (word_cling_64_04_l(f3[j2]<<32)>>0);
        t[j+2] = (word_cling_64_04_l(f0[j2]<<16)>>3) | (word_cling_64_04_l(f1[j2]<<16)>>2) | (word_cling_64_04_l(f2[j2]<<16)>>1) | (word_cling_64_04_l(f3[j2]<<16)>>0);
        tmp    = (word_cling_64_04_l(f0[j2]<< 0)>>3) | (word_cling_64_04_l(f1[j2]<< 0)>>2) | (word_cling_64_04_l(f2[j2]<< 0)>>1) | (word_cling_64_04_l(f3[j2]<< 0)>>0);
        t[j+3] = (t[j+3] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 3:
        t[j+0] = (word_cling_64_04_l(f0[j2]<<48)>>3) | (word_cling_64_04_l(f1[j2]<<48)>>2) | (word_cling_64_04_l(f2[j2]<<48)>>1) | (word_cling_64_04_l(f3[j2]<<48)>>0);
        t[j+1] = (word_cling_64_04_l(f0[j2]<<32)>>3) | (word_cling_64_04_l(f1[j2]<<32)>>2) | (word_cling_64_04_l(f2[j2]<<32)>>1) | (word_cling_64_04_l(f3[j2]<<32)>>0);
        tmp    = (word_cling_64_04_l(f0[j2]<<16)>>3) | (word_cling_64_04_l(f1[j2]<<16)>>2) | (word_cling_64_04_l(f2[j2]<<16)>>1) | (word_cling_64_04_l(f3[j2]<<16)>>0);
        t[j+2] = (t[j+2] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 2:
        t[j+0] = (word_cling_64_04_l(f0[j2]<<48)>>3) | (word_cling_64_04_l(f1[j2]<<48)>>2) | (word_cling_64_04_l(f2[j2]<<48)>>1) | (word_cling_64_04_l(f3[j2]<<48)>>0);
        tmp    = (word_cling_64_04_l(f0[j2]<<32)>>3) | (word_cling_64_04_l(f1[j2]<<32)>>2) | (word_cling_64_04_l(f2[j2]<<32)>>1) | (word_cling_64_04_l(f3[j2]<<32)>>0);
        t[j+1] = (t[j+1] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 1:
        tmp =    (word_cling_64_04_l(f0[j2]<<48)>>3) | (word_cling_64_04_l(f1[j2]<<48)>>2) | (word_cling_64_04_l(f2[j2]<<48)>>1) | (word_cling_64_04_l(f3[j2]<<48)>>0);
        t[j+0] = (t[j+0] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      default:
        m4ri_die("impossible");
      }
    }
  } else { //degree == 3
    for(rci_t i=0; i<T->nrows; i++) {
      const word *f0 = F->x[0]->rows[i];
      const word *f1 = F->x[1]->rows[i];
      const word *f2 = F->x[2]->rows[i];
      word *t  = T->x->rows[i];

      for(j=0, j2=0; j+4 < T->x->width; j+=4, j2++) {
        t[j+0] = (word_cling_64_04_l(f0[j2]<<48)>>3) | (word_cling_64_04_l(f1[j2]<<48)>>2) | (word_cling_64_04_l(f2[j2]<<48)>>1);
        t[j+1] = (word_cling_64_04_l(f0[j2]<<32)>>3) | (word_cling_64_04_l(f1[j2]<<32)>>2) | (word_cling_64_04_l(f2[j2]<<32)>>1);
        t[j+2] = (word_cling_64_04_l(f0[j2]<<16)>>3) | (word_cling_64_04_l(f1[j2]<<16)>>2) | (word_cling_64_04_l(f2[j2]<<16)>>1);
        t[j+3] = (word_cling_64_04_l(f0[j2]<< 0)>>3) | (word_cling_64_04_l(f1[j2]<< 0)>>2) | (word_cling_64_04_l(f2[j2]<< 0)>>1);
      }

      register word tmp=0;
      switch(T->x->width - j) {
      case 4:
        t[j+0] = (word_cling_64_04_l(f0[j2]<<48)>>3) | (word_cling_64_04_l(f1[j2]<<48)>>2) | (word_cling_64_04_l(f2[j2]<<48)>>1);
        t[j+1] = (word_cling_64_04_l(f0[j2]<<32)>>3) | (word_cling_64_04_l(f1[j2]<<32)>>2) | (word_cling_64_04_l(f2[j2]<<32)>>1);
        t[j+2] = (word_cling_64_04_l(f0[j2]<<16)>>3) | (word_cling_64_04_l(f1[j2]<<16)>>2) | (word_cling_64_04_l(f2[j2]<<16)>>1);
        tmp  =   (word_cling_64_04_l(f0[j2]<< 0)>>3) | (word_cling_64_04_l(f1[j2]<< 0)>>2) | (word_cling_64_04_l(f2[j2]<< 0)>>1);
        t[j+3] = (t[j+3] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 3:
        t[j+0] = (word_cling_64_04_l(f0[j2]<<48)>>3) | (word_cling_64_04_l(f1[j2]<<48)>>2) | (word_cling_64_04_l(f2[j2]<<48)>>1);
        t[j+1] = (word_cling_64_04_l(f0[j2]<<32)>>3) | (word_cling_64_04_l(f1[j2]<<32)>>2) | (word_cling_64_04_l(f2[j2]<<32)>>1);
        tmp =    (word_cling_64_04_l(f0[j2]<<16)>>3) | (word_cling_64_04_l(f1[j2]<<16)>>2) | (word_cling_64_04_l(f2[j2]<<16)>>1);
        t[j+2] = (t[j+2] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 2:
        t[j+0] = (word_cling_64_04_l(f0[j2]<<48)>>3) | (word_cling_64_04_l(f1[j2]<<48)>>2) | (word_cling_64_04_l(f2[j2]<<48)>>1);
        tmp =    (word_cling_64_04_l(f0[j2]<<32)>>3) | (word_cling_64_04_l(f1[j2]<<32)>>2) | (word_cling_64_04_l(f2[j2]<<32)>>1);
        t[j+1] = (t[j+1] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 1:
        tmp =    (word_cling_64_04_l(f0[j2]<<48)>>3) | (word_cling_64_04_l(f1[j2]<<48)>>2) | (word_cling_64_04_l(f2[j2]<<48)>>1);
        t[j+0] = (t[j+0] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      default:
        m4ri_die("impossible");
      }
    }
  }
  return T;
}

mzd_slice_t *_mzed_slice8(mzd_slice_t *T, const mzed_t *F) {
  assert(T && (4 < T->depth && T->depth <= 8) && T->x[0]->offset == 0);
  size_t j, j2 = 0;
  register word r0,r1,r2,r3,r4,r5,r6,r7 = 0;

  const word bitmask_end = __M4RI_LEFT_BITMASK((T->x[0]->offset + T->ncols) % m4ri_radix);

  if (mzed_is_zero(F))
    return T;

  switch(T->depth) {
  case 8: {
      for(size_t i=0; i<T->nrows; i++) {
        word *t0 = T->x[0]->rows[i];
        word *t1 = T->x[1]->rows[i];
        word *t2 = T->x[2]->rows[i];
        word *t3 = T->x[3]->rows[i];
        word *t4 = T->x[4]->rows[i];
        word *t5 = T->x[5]->rows[i];
        word *t6 = T->x[6]->rows[i];
        word *t7 = T->x[7]->rows[i];
        const word const *f  = F->x->rows[i];

        /* bulk of work */
        for(j=0, j2=0; j+8 < F->x->width; j+=8,j2++) {
          t0[j2] |= word_slice_64_08_l(f[j+0]<<7 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<7 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<7 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<7 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<7 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<7 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<7 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<7 & x80808080)>> 0;

          t1[j2] |= word_slice_64_08_l(f[j+0]<<6 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<6 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<6 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<6 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<6 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<6 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<6 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<6 & x80808080)>> 0;

          t2[j2] |= word_slice_64_08_l(f[j+0]<<5 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<5 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<5 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<5 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<5 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<5 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<5 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<5 & x80808080)>> 0;

          t3[j2] |= word_slice_64_08_l(f[j+0]<<4 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<4 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<4 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<4 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<4 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<4 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<4 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<4 & x80808080)>> 0;

          t4[j2] |= word_slice_64_08_l(f[j+0]<<3 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<3 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<3 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<3 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<3 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<3 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<3 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<3 & x80808080)>> 0;

          t5[j2] |= word_slice_64_08_l(f[j+0]<<2 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<2 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<2 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<2 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<2 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<2 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<2 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<2 & x80808080)>> 0;

          t6[j2] |= word_slice_64_08_l(f[j+0]<<1 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<1 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<1 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<1 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<1 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<1 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<1 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<1 & x80808080)>> 0;

          t7[j2] |= word_slice_64_08_l(f[j+0]<<0 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<0 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<0 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<0 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<0 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<0 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<0 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<0 & x80808080)>> 0;
        }
        r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = 0;
        switch(F->x->width - j) {
        case 8:
          r0 |= word_slice_64_08_l(f[j+7]<<7 & x80808080)>> 0;
          r1 |= word_slice_64_08_l(f[j+7]<<6 & x80808080)>> 0;
          r2 |= word_slice_64_08_l(f[j+7]<<5 & x80808080)>> 0;
          r3 |= word_slice_64_08_l(f[j+7]<<4 & x80808080)>> 0;
          r4 |= word_slice_64_08_l(f[j+7]<<3 & x80808080)>> 0;
          r5 |= word_slice_64_08_l(f[j+7]<<2 & x80808080)>> 0;
          r6 |= word_slice_64_08_l(f[j+7]<<1 & x80808080)>> 0;
          r7 |= word_slice_64_08_l(f[j+7]<<0 & x80808080)>> 0;
        case 7:
          r0 |= word_slice_64_08_l(f[j+6]<<7 & x80808080)>> 8;
          r1 |= word_slice_64_08_l(f[j+6]<<6 & x80808080)>> 8;
          r2 |= word_slice_64_08_l(f[j+6]<<5 & x80808080)>> 8;
          r3 |= word_slice_64_08_l(f[j+6]<<4 & x80808080)>> 8;
          r4 |= word_slice_64_08_l(f[j+6]<<3 & x80808080)>> 8;
          r5 |= word_slice_64_08_l(f[j+6]<<2 & x80808080)>> 8;
          r6 |= word_slice_64_08_l(f[j+6]<<1 & x80808080)>> 8;
          r7 |= word_slice_64_08_l(f[j+6]<<0 & x80808080)>> 8;
        case 6:
          r0 |= word_slice_64_08_l(f[j+5]<<7 & x80808080)>>16;
          r1 |= word_slice_64_08_l(f[j+5]<<6 & x80808080)>>16;
          r2 |= word_slice_64_08_l(f[j+5]<<5 & x80808080)>>16;
          r3 |= word_slice_64_08_l(f[j+5]<<4 & x80808080)>>16;
          r4 |= word_slice_64_08_l(f[j+5]<<3 & x80808080)>>16;
          r5 |= word_slice_64_08_l(f[j+5]<<2 & x80808080)>>16;
          r6 |= word_slice_64_08_l(f[j+5]<<1 & x80808080)>>16;
          r7 |= word_slice_64_08_l(f[j+5]<<0 & x80808080)>>16;
        case 5:
          r0 |= word_slice_64_08_l(f[j+4]<<7 & x80808080)>>24;
          r1 |= word_slice_64_08_l(f[j+4]<<6 & x80808080)>>24;
          r2 |= word_slice_64_08_l(f[j+4]<<5 & x80808080)>>24;
          r3 |= word_slice_64_08_l(f[j+4]<<4 & x80808080)>>24;
          r4 |= word_slice_64_08_l(f[j+4]<<3 & x80808080)>>24;
          r5 |= word_slice_64_08_l(f[j+4]<<2 & x80808080)>>24;
          r6 |= word_slice_64_08_l(f[j+4]<<1 & x80808080)>>24;
          r7 |= word_slice_64_08_l(f[j+4]<<0 & x80808080)>>24;
        case 4:
          r0 |= word_slice_64_08_l(f[j+3]<<7 & x80808080)>>32;
          r1 |= word_slice_64_08_l(f[j+3]<<6 & x80808080)>>32;
          r2 |= word_slice_64_08_l(f[j+3]<<5 & x80808080)>>32;
          r3 |= word_slice_64_08_l(f[j+3]<<4 & x80808080)>>32;
          r4 |= word_slice_64_08_l(f[j+3]<<3 & x80808080)>>32;
          r5 |= word_slice_64_08_l(f[j+3]<<2 & x80808080)>>32;
          r6 |= word_slice_64_08_l(f[j+3]<<1 & x80808080)>>32;
          r7 |= word_slice_64_08_l(f[j+3]<<0 & x80808080)>>32;
        case 3:
          r0 |= word_slice_64_08_l(f[j+2]<<7 & x80808080)>>40;
          r1 |= word_slice_64_08_l(f[j+2]<<6 & x80808080)>>40;
          r2 |= word_slice_64_08_l(f[j+2]<<5 & x80808080)>>40;
          r3 |= word_slice_64_08_l(f[j+2]<<4 & x80808080)>>40;
          r4 |= word_slice_64_08_l(f[j+2]<<3 & x80808080)>>40;
          r5 |= word_slice_64_08_l(f[j+2]<<2 & x80808080)>>40;
          r6 |= word_slice_64_08_l(f[j+2]<<1 & x80808080)>>40;
          r7 |= word_slice_64_08_l(f[j+2]<<0 & x80808080)>>40;
        case 2:
          r0 |= word_slice_64_08_l(f[j+1]<<7 & x80808080)>>48;
          r1 |= word_slice_64_08_l(f[j+1]<<6 & x80808080)>>48;
          r2 |= word_slice_64_08_l(f[j+1]<<5 & x80808080)>>48;
          r3 |= word_slice_64_08_l(f[j+1]<<4 & x80808080)>>48;
          r4 |= word_slice_64_08_l(f[j+1]<<3 & x80808080)>>48;
          r5 |= word_slice_64_08_l(f[j+1]<<2 & x80808080)>>48;
          r6 |= word_slice_64_08_l(f[j+1]<<1 & x80808080)>>48;
          r7 |= word_slice_64_08_l(f[j+1]<<0 & x80808080)>>48;
        case 1:
          r0 |= word_slice_64_08_l(f[j+0]<<7 & x80808080)>>56;
          r1 |= word_slice_64_08_l(f[j+0]<<6 & x80808080)>>56;
          r2 |= word_slice_64_08_l(f[j+0]<<5 & x80808080)>>56;
          r3 |= word_slice_64_08_l(f[j+0]<<4 & x80808080)>>56;
          r4 |= word_slice_64_08_l(f[j+0]<<3 & x80808080)>>56;
          r5 |= word_slice_64_08_l(f[j+0]<<2 & x80808080)>>56;
          r6 |= word_slice_64_08_l(f[j+0]<<1 & x80808080)>>56;
          r7 |= word_slice_64_08_l(f[j+0]<<0 & x80808080)>>56;
          break;
        default:
          m4ri_die("impossible");
        }
        t0[j2] |= r0 & bitmask_end;
        t1[j2] |= r1 & bitmask_end;
        t2[j2] |= r2 & bitmask_end;
        t3[j2] |= r3 & bitmask_end;
        t4[j2] |= r4 & bitmask_end;
        t5[j2] |= r5 & bitmask_end;
        t6[j2] |= r6 & bitmask_end;
        t7[j2] |= r7 & bitmask_end;
      }
    }
    break;

  case 7: {
      for(size_t i=0; i<T->nrows; i++) {
        word *t0 = T->x[0]->rows[i];
        word *t1 = T->x[1]->rows[i];
        word *t2 = T->x[2]->rows[i];
        word *t3 = T->x[3]->rows[i];
        word *t4 = T->x[4]->rows[i];
        word *t5 = T->x[5]->rows[i];
        word *t6 = T->x[6]->rows[i];
        const word const *f  = F->x->rows[i];

        /* bulk of work */
        for(j=0, j2=0; j+8 < F->x->width; j+=8,j2++) {
          t0[j2] |= word_slice_64_08_l(f[j+0]<<7 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<7 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<7 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<7 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<7 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<7 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<7 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<7 & x80808080)>> 0;

          t1[j2] |= word_slice_64_08_l(f[j+0]<<6 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<6 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<6 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<6 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<6 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<6 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<6 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<6 & x80808080)>> 0;

          t2[j2] |= word_slice_64_08_l(f[j+0]<<5 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<5 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<5 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<5 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<5 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<5 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<5 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<5 & x80808080)>> 0;

          t3[j2] |= word_slice_64_08_l(f[j+0]<<4 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<4 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<4 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<4 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<4 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<4 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<4 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<4 & x80808080)>> 0;

          t4[j2] |= word_slice_64_08_l(f[j+0]<<3 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<3 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<3 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<3 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<3 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<3 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<3 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<3 & x80808080)>> 0;

          t5[j2] |= word_slice_64_08_l(f[j+0]<<2 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<2 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<2 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<2 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<2 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<2 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<2 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<2 & x80808080)>> 0;

          t6[j2] |= word_slice_64_08_l(f[j+0]<<1 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<1 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<1 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<1 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<1 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<1 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<1 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<1 & x80808080)>> 0;
        }
        r0 = r1 = r2 = r3 = r4 = r5 = r6 = 0;
        switch(F->x->width - j) {
        case 8:
          r0 |= word_slice_64_08_l(f[j+7]<<7 & x80808080)>> 0;
          r1 |= word_slice_64_08_l(f[j+7]<<6 & x80808080)>> 0;
          r2 |= word_slice_64_08_l(f[j+7]<<5 & x80808080)>> 0;
          r3 |= word_slice_64_08_l(f[j+7]<<4 & x80808080)>> 0;
          r4 |= word_slice_64_08_l(f[j+7]<<3 & x80808080)>> 0;
          r5 |= word_slice_64_08_l(f[j+7]<<2 & x80808080)>> 0;
          r6 |= word_slice_64_08_l(f[j+7]<<1 & x80808080)>> 0;
        case 7:
          r0 |= word_slice_64_08_l(f[j+6]<<7 & x80808080)>> 8;
          r1 |= word_slice_64_08_l(f[j+6]<<6 & x80808080)>> 8;
          r2 |= word_slice_64_08_l(f[j+6]<<5 & x80808080)>> 8;
          r3 |= word_slice_64_08_l(f[j+6]<<4 & x80808080)>> 8;
          r4 |= word_slice_64_08_l(f[j+6]<<3 & x80808080)>> 8;
          r5 |= word_slice_64_08_l(f[j+6]<<2 & x80808080)>> 8;
          r6 |= word_slice_64_08_l(f[j+6]<<1 & x80808080)>> 8;
        case 6:
          r0 |= word_slice_64_08_l(f[j+5]<<7 & x80808080)>>16;
          r1 |= word_slice_64_08_l(f[j+5]<<6 & x80808080)>>16;
          r2 |= word_slice_64_08_l(f[j+5]<<5 & x80808080)>>16;
          r3 |= word_slice_64_08_l(f[j+5]<<4 & x80808080)>>16;
          r4 |= word_slice_64_08_l(f[j+5]<<3 & x80808080)>>16;
          r5 |= word_slice_64_08_l(f[j+5]<<2 & x80808080)>>16;
          r6 |= word_slice_64_08_l(f[j+5]<<1 & x80808080)>>16;
        case 5:
          r0 |= word_slice_64_08_l(f[j+4]<<7 & x80808080)>>24;
          r1 |= word_slice_64_08_l(f[j+4]<<6 & x80808080)>>24;
          r2 |= word_slice_64_08_l(f[j+4]<<5 & x80808080)>>24;
          r3 |= word_slice_64_08_l(f[j+4]<<4 & x80808080)>>24;
          r4 |= word_slice_64_08_l(f[j+4]<<3 & x80808080)>>24;
          r5 |= word_slice_64_08_l(f[j+4]<<2 & x80808080)>>24;
          r6 |= word_slice_64_08_l(f[j+4]<<1 & x80808080)>>24;
        case 4:
          r0 |= word_slice_64_08_l(f[j+3]<<7 & x80808080)>>32;
          r1 |= word_slice_64_08_l(f[j+3]<<6 & x80808080)>>32;
          r2 |= word_slice_64_08_l(f[j+3]<<5 & x80808080)>>32;
          r3 |= word_slice_64_08_l(f[j+3]<<4 & x80808080)>>32;
          r4 |= word_slice_64_08_l(f[j+3]<<3 & x80808080)>>32;
          r5 |= word_slice_64_08_l(f[j+3]<<2 & x80808080)>>32;
          r6 |= word_slice_64_08_l(f[j+3]<<1 & x80808080)>>32;
        case 3:
          r0 |= word_slice_64_08_l(f[j+2]<<7 & x80808080)>>40;
          r1 |= word_slice_64_08_l(f[j+2]<<6 & x80808080)>>40;
          r2 |= word_slice_64_08_l(f[j+2]<<5 & x80808080)>>40;
          r3 |= word_slice_64_08_l(f[j+2]<<4 & x80808080)>>40;
          r4 |= word_slice_64_08_l(f[j+2]<<3 & x80808080)>>40;
          r5 |= word_slice_64_08_l(f[j+2]<<2 & x80808080)>>40;
          r6 |= word_slice_64_08_l(f[j+2]<<1 & x80808080)>>40;
        case 2:
          r0 |= word_slice_64_08_l(f[j+1]<<7 & x80808080)>>48;
          r1 |= word_slice_64_08_l(f[j+1]<<6 & x80808080)>>48;
          r2 |= word_slice_64_08_l(f[j+1]<<5 & x80808080)>>48;
          r3 |= word_slice_64_08_l(f[j+1]<<4 & x80808080)>>48;
          r4 |= word_slice_64_08_l(f[j+1]<<3 & x80808080)>>48;
          r5 |= word_slice_64_08_l(f[j+1]<<2 & x80808080)>>48;
          r6 |= word_slice_64_08_l(f[j+1]<<1 & x80808080)>>48;
        case 1:
          r0 |= word_slice_64_08_l(f[j+0]<<7 & x80808080)>>56;
          r1 |= word_slice_64_08_l(f[j+0]<<6 & x80808080)>>56;
          r2 |= word_slice_64_08_l(f[j+0]<<5 & x80808080)>>56;
          r3 |= word_slice_64_08_l(f[j+0]<<4 & x80808080)>>56;
          r4 |= word_slice_64_08_l(f[j+0]<<3 & x80808080)>>56;
          r5 |= word_slice_64_08_l(f[j+0]<<2 & x80808080)>>56;
          r6 |= word_slice_64_08_l(f[j+0]<<1 & x80808080)>>56;
          break;
        default:
          m4ri_die("impossible");
        }
        t0[j2] |= r0 & bitmask_end;
        t1[j2] |= r1 & bitmask_end;
        t2[j2] |= r2 & bitmask_end;
        t3[j2] |= r3 & bitmask_end;
        t4[j2] |= r4 & bitmask_end;
        t5[j2] |= r5 & bitmask_end;
        t6[j2] |= r6 & bitmask_end;
      }
    }
    break;

  case 6: {
      for(size_t i=0; i<T->nrows; i++) {
        word *t0 = T->x[0]->rows[i];
        word *t1 = T->x[1]->rows[i];
        word *t2 = T->x[2]->rows[i];
        word *t3 = T->x[3]->rows[i];
        word *t4 = T->x[4]->rows[i];
        word *t5 = T->x[5]->rows[i];
        const word const *f  = F->x->rows[i];

        /* bulk of work */
        for(j=0, j2=0; j+8 < F->x->width; j+=8,j2++) {
          t0[j2] |= word_slice_64_08_l(f[j+0]<<7 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<7 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<7 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<7 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<7 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<7 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<7 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<7 & x80808080)>> 0;

          t1[j2] |= word_slice_64_08_l(f[j+0]<<6 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<6 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<6 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<6 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<6 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<6 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<6 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<6 & x80808080)>> 0;

          t2[j2] |= word_slice_64_08_l(f[j+0]<<5 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<5 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<5 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<5 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<5 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<5 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<5 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<5 & x80808080)>> 0;

          t3[j2] |= word_slice_64_08_l(f[j+0]<<4 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<4 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<4 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<4 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<4 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<4 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<4 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<4 & x80808080)>> 0;

          t4[j2] |= word_slice_64_08_l(f[j+0]<<3 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<3 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<3 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<3 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<3 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<3 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<3 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<3 & x80808080)>> 0;

          t5[j2] |= word_slice_64_08_l(f[j+0]<<2 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<2 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<2 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<2 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<2 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<2 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<2 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<2 & x80808080)>> 0;
        }
        r0 = r1 = r2 = r3 = r4 = r5 = 0;
        switch(F->x->width - j) {
        case 8:
          r0 |= word_slice_64_08_l(f[j+7]<<7 & x80808080)>> 0;
          r1 |= word_slice_64_08_l(f[j+7]<<6 & x80808080)>> 0;
          r2 |= word_slice_64_08_l(f[j+7]<<5 & x80808080)>> 0;
          r3 |= word_slice_64_08_l(f[j+7]<<4 & x80808080)>> 0;
          r4 |= word_slice_64_08_l(f[j+7]<<3 & x80808080)>> 0;
          r5 |= word_slice_64_08_l(f[j+7]<<2 & x80808080)>> 0;
        case 7:
          r0 |= word_slice_64_08_l(f[j+6]<<7 & x80808080)>> 8;
          r1 |= word_slice_64_08_l(f[j+6]<<6 & x80808080)>> 8;
          r2 |= word_slice_64_08_l(f[j+6]<<5 & x80808080)>> 8;
          r3 |= word_slice_64_08_l(f[j+6]<<4 & x80808080)>> 8;
          r4 |= word_slice_64_08_l(f[j+6]<<3 & x80808080)>> 8;
          r5 |= word_slice_64_08_l(f[j+6]<<2 & x80808080)>> 8;
        case 6:
          r0 |= word_slice_64_08_l(f[j+5]<<7 & x80808080)>>16;
          r1 |= word_slice_64_08_l(f[j+5]<<6 & x80808080)>>16;
          r2 |= word_slice_64_08_l(f[j+5]<<5 & x80808080)>>16;
          r3 |= word_slice_64_08_l(f[j+5]<<4 & x80808080)>>16;
          r4 |= word_slice_64_08_l(f[j+5]<<3 & x80808080)>>16;
          r5 |= word_slice_64_08_l(f[j+5]<<2 & x80808080)>>16;
        case 5:
          r0 |= word_slice_64_08_l(f[j+4]<<7 & x80808080)>>24;
          r1 |= word_slice_64_08_l(f[j+4]<<6 & x80808080)>>24;
          r2 |= word_slice_64_08_l(f[j+4]<<5 & x80808080)>>24;
          r3 |= word_slice_64_08_l(f[j+4]<<4 & x80808080)>>24;
          r4 |= word_slice_64_08_l(f[j+4]<<3 & x80808080)>>24;
          r5 |= word_slice_64_08_l(f[j+4]<<2 & x80808080)>>24;
        case 4:
          r0 |= word_slice_64_08_l(f[j+3]<<7 & x80808080)>>32;
          r1 |= word_slice_64_08_l(f[j+3]<<6 & x80808080)>>32;
          r2 |= word_slice_64_08_l(f[j+3]<<5 & x80808080)>>32;
          r3 |= word_slice_64_08_l(f[j+3]<<4 & x80808080)>>32;
          r4 |= word_slice_64_08_l(f[j+3]<<3 & x80808080)>>32;
          r5 |= word_slice_64_08_l(f[j+3]<<2 & x80808080)>>32;
        case 3:
          r0 |= word_slice_64_08_l(f[j+2]<<7 & x80808080)>>40;
          r1 |= word_slice_64_08_l(f[j+2]<<6 & x80808080)>>40;
          r2 |= word_slice_64_08_l(f[j+2]<<5 & x80808080)>>40;
          r3 |= word_slice_64_08_l(f[j+2]<<4 & x80808080)>>40;
          r4 |= word_slice_64_08_l(f[j+2]<<3 & x80808080)>>40;
          r5 |= word_slice_64_08_l(f[j+2]<<2 & x80808080)>>40;
        case 2:
          r0 |= word_slice_64_08_l(f[j+1]<<7 & x80808080)>>48;
          r1 |= word_slice_64_08_l(f[j+1]<<6 & x80808080)>>48;
          r2 |= word_slice_64_08_l(f[j+1]<<5 & x80808080)>>48;
          r3 |= word_slice_64_08_l(f[j+1]<<4 & x80808080)>>48;
          r4 |= word_slice_64_08_l(f[j+1]<<3 & x80808080)>>48;
          r5 |= word_slice_64_08_l(f[j+1]<<2 & x80808080)>>48;
        case 1:
          r0 |= word_slice_64_08_l(f[j+0]<<7 & x80808080)>>56;
          r1 |= word_slice_64_08_l(f[j+0]<<6 & x80808080)>>56;
          r2 |= word_slice_64_08_l(f[j+0]<<5 & x80808080)>>56;
          r3 |= word_slice_64_08_l(f[j+0]<<4 & x80808080)>>56;
          r4 |= word_slice_64_08_l(f[j+0]<<3 & x80808080)>>56;
          r5 |= word_slice_64_08_l(f[j+0]<<2 & x80808080)>>56;
          break;
        default:
          m4ri_die("impossible");
        }
        t0[j2] |= r0 & bitmask_end;
        t1[j2] |= r1 & bitmask_end;
        t2[j2] |= r2 & bitmask_end;
        t3[j2] |= r3 & bitmask_end;
        t4[j2] |= r4 & bitmask_end;
        t5[j2] |= r5 & bitmask_end;
      }
    }
    break;

  case 5: {
      for(size_t i=0; i<T->nrows; i++) {
        word *t0 = T->x[0]->rows[i];
        word *t1 = T->x[1]->rows[i];
        word *t2 = T->x[2]->rows[i];
        word *t3 = T->x[3]->rows[i];
        word *t4 = T->x[4]->rows[i];
        const word const *f  = F->x->rows[i];

        /* bulk of work */
        for(j=0, j2=0; j+8 < F->x->width; j+=8,j2++) {
          t0[j2] |= word_slice_64_08_l(f[j+0]<<7 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<7 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<7 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<7 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<7 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<7 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<7 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<7 & x80808080)>> 0;

          t1[j2] |= word_slice_64_08_l(f[j+0]<<6 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<6 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<6 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<6 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<6 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<6 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<6 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<6 & x80808080)>> 0;

          t2[j2] |= word_slice_64_08_l(f[j+0]<<5 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<5 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<5 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<5 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<5 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<5 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<5 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<5 & x80808080)>> 0;

          t3[j2] |= word_slice_64_08_l(f[j+0]<<4 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<4 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<4 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<4 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<4 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<4 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<4 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<4 & x80808080)>> 0;

          t4[j2] |= word_slice_64_08_l(f[j+0]<<3 & x80808080)>>56 | word_slice_64_08_l(f[j+1]<<3 & x80808080)>>48 \
            |       word_slice_64_08_l(f[j+2]<<3 & x80808080)>>40 | word_slice_64_08_l(f[j+3]<<3 & x80808080)>>32 \
            |       word_slice_64_08_l(f[j+4]<<3 & x80808080)>>24 | word_slice_64_08_l(f[j+5]<<3 & x80808080)>>16 \
            |       word_slice_64_08_l(f[j+6]<<3 & x80808080)>> 8 | word_slice_64_08_l(f[j+7]<<3 & x80808080)>> 0;
        }
        r0 = r1 = r2 = r3 = r4 = 0;
        switch(F->x->width - j) {
        case 8:
          r0 |= word_slice_64_08_l(f[j+7]<<7 & x80808080)>> 0;
          r1 |= word_slice_64_08_l(f[j+7]<<6 & x80808080)>> 0;
          r2 |= word_slice_64_08_l(f[j+7]<<5 & x80808080)>> 0;
          r3 |= word_slice_64_08_l(f[j+7]<<4 & x80808080)>> 0;
          r4 |= word_slice_64_08_l(f[j+7]<<3 & x80808080)>> 0;
        case 7:
          r0 |= word_slice_64_08_l(f[j+6]<<7 & x80808080)>> 8;
          r1 |= word_slice_64_08_l(f[j+6]<<6 & x80808080)>> 8;
          r2 |= word_slice_64_08_l(f[j+6]<<5 & x80808080)>> 8;
          r3 |= word_slice_64_08_l(f[j+6]<<4 & x80808080)>> 8;
          r4 |= word_slice_64_08_l(f[j+6]<<3 & x80808080)>> 8;
        case 6:
          r0 |= word_slice_64_08_l(f[j+5]<<7 & x80808080)>>16;
          r1 |= word_slice_64_08_l(f[j+5]<<6 & x80808080)>>16;
          r2 |= word_slice_64_08_l(f[j+5]<<5 & x80808080)>>16;
          r3 |= word_slice_64_08_l(f[j+5]<<4 & x80808080)>>16;
          r4 |= word_slice_64_08_l(f[j+5]<<3 & x80808080)>>16;
        case 5:
          r0 |= word_slice_64_08_l(f[j+4]<<7 & x80808080)>>24;
          r1 |= word_slice_64_08_l(f[j+4]<<6 & x80808080)>>24;
          r2 |= word_slice_64_08_l(f[j+4]<<5 & x80808080)>>24;
          r3 |= word_slice_64_08_l(f[j+4]<<4 & x80808080)>>24;
          r4 |= word_slice_64_08_l(f[j+4]<<3 & x80808080)>>24;
        case 4:
          r0 |= word_slice_64_08_l(f[j+3]<<7 & x80808080)>>32;
          r1 |= word_slice_64_08_l(f[j+3]<<6 & x80808080)>>32;
          r2 |= word_slice_64_08_l(f[j+3]<<5 & x80808080)>>32;
          r3 |= word_slice_64_08_l(f[j+3]<<4 & x80808080)>>32;
          r4 |= word_slice_64_08_l(f[j+3]<<3 & x80808080)>>32;
        case 3:
          r0 |= word_slice_64_08_l(f[j+2]<<7 & x80808080)>>40;
          r1 |= word_slice_64_08_l(f[j+2]<<6 & x80808080)>>40;
          r2 |= word_slice_64_08_l(f[j+2]<<5 & x80808080)>>40;
          r3 |= word_slice_64_08_l(f[j+2]<<4 & x80808080)>>40;
          r4 |= word_slice_64_08_l(f[j+2]<<3 & x80808080)>>40;
        case 2:
          r0 |= word_slice_64_08_l(f[j+1]<<7 & x80808080)>>48;
          r1 |= word_slice_64_08_l(f[j+1]<<6 & x80808080)>>48;
          r2 |= word_slice_64_08_l(f[j+1]<<5 & x80808080)>>48;
          r3 |= word_slice_64_08_l(f[j+1]<<4 & x80808080)>>48;
          r4 |= word_slice_64_08_l(f[j+1]<<3 & x80808080)>>48;
        case 1:
          r0 |= word_slice_64_08_l(f[j+0]<<7 & x80808080)>>56;
          r1 |= word_slice_64_08_l(f[j+0]<<6 & x80808080)>>56;
          r2 |= word_slice_64_08_l(f[j+0]<<5 & x80808080)>>56;
          r3 |= word_slice_64_08_l(f[j+0]<<4 & x80808080)>>56;
          r4 |= word_slice_64_08_l(f[j+0]<<3 & x80808080)>>56;
          break;
        default:
          m4ri_die("impossible");
        }
        t0[j2] |= r0 & bitmask_end;
        t1[j2] |= r1 & bitmask_end;
        t2[j2] |= r2 & bitmask_end;
        t3[j2] |= r3 & bitmask_end;
        t4[j2] |= r4 & bitmask_end;
      }
    }
    break;

  default:
    m4ri_die("impossible\n");
  }
  return T;
}

mzed_t *_mzed_cling8(mzed_t *T, const mzd_slice_t *F) {
  size_t j,j2 = 0;

  const word bitmask_end = __M4RI_LEFT_BITMASK((T->x->offset + T->x->ncols) % m4ri_radix);

  if (mzd_slice_is_zero(F))
    return T;

  switch (F->finite_field->degree) {
  case 8: {
    for(rci_t i=0; i<T->nrows; i++) {
      const word *f0 = F->x[0]->rows[i];
      const word *f1 = F->x[1]->rows[i];
      const word *f2 = F->x[2]->rows[i];
      const word *f3 = F->x[3]->rows[i];
      const word *f4 = F->x[4]->rows[i];
      const word *f5 = F->x[5]->rows[i];
      const word *f6 = F->x[6]->rows[i];
      const word *f7 = F->x[7]->rows[i];
      word *t  = T->x->rows[i];

      for(j=0, j2=0; j+8 < T->x->width; j+=8, j2++) {
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3) | (word_cling_64_08_l(f5[j2]<<56)>>2) | (word_cling_64_08_l(f6[j2]<<56)>>1) | (word_cling_64_08_l(f7[j2]<<56)>>0);
        t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3) | (word_cling_64_08_l(f5[j2]<<48)>>2) | (word_cling_64_08_l(f6[j2]<<48)>>1) | (word_cling_64_08_l(f7[j2]<<48)>>0);
        t[j+2] = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<40)>>3) | (word_cling_64_08_l(f5[j2]<<40)>>2) | (word_cling_64_08_l(f6[j2]<<40)>>1) | (word_cling_64_08_l(f7[j2]<<40)>>0);
        t[j+3] = (word_cling_64_08_l(f0[j2]<<32)>>7) | (word_cling_64_08_l(f1[j2]<<32)>>6) | (word_cling_64_08_l(f2[j2]<<32)>>5) | (word_cling_64_08_l(f3[j2]<<32)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<32)>>3) | (word_cling_64_08_l(f5[j2]<<32)>>2) | (word_cling_64_08_l(f6[j2]<<32)>>1) | (word_cling_64_08_l(f7[j2]<<32)>>0);
        t[j+4] = (word_cling_64_08_l(f0[j2]<<24)>>7) | (word_cling_64_08_l(f1[j2]<<24)>>6) | (word_cling_64_08_l(f2[j2]<<24)>>5) | (word_cling_64_08_l(f3[j2]<<24)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<24)>>3) | (word_cling_64_08_l(f5[j2]<<24)>>2) | (word_cling_64_08_l(f6[j2]<<24)>>1) | (word_cling_64_08_l(f7[j2]<<24)>>0);
        t[j+5] = (word_cling_64_08_l(f0[j2]<<16)>>7) | (word_cling_64_08_l(f1[j2]<<16)>>6) | (word_cling_64_08_l(f2[j2]<<16)>>5) | (word_cling_64_08_l(f3[j2]<<16)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<16)>>3) | (word_cling_64_08_l(f5[j2]<<16)>>2) | (word_cling_64_08_l(f6[j2]<<16)>>1) | (word_cling_64_08_l(f7[j2]<<16)>>0);
        t[j+6] = (word_cling_64_08_l(f0[j2]<< 8)>>7) | (word_cling_64_08_l(f1[j2]<< 8)>>6) | (word_cling_64_08_l(f2[j2]<< 8)>>5) | (word_cling_64_08_l(f3[j2]<< 8)>>4) \
          |      (word_cling_64_08_l(f4[j2]<< 8)>>3) | (word_cling_64_08_l(f5[j2]<< 8)>>2) | (word_cling_64_08_l(f6[j2]<< 8)>>1) | (word_cling_64_08_l(f7[j2]<< 8)>>0);
        t[j+7] = (word_cling_64_08_l(f0[j2]<< 0)>>7) | (word_cling_64_08_l(f1[j2]<< 0)>>6) | (word_cling_64_08_l(f2[j2]<< 0)>>5) | (word_cling_64_08_l(f3[j2]<< 0)>>4) \
          |      (word_cling_64_08_l(f4[j2]<< 0)>>3) | (word_cling_64_08_l(f5[j2]<< 0)>>2) | (word_cling_64_08_l(f6[j2]<< 0)>>1) | (word_cling_64_08_l(f7[j2]<< 0)>>0);
      }

      register word tmp = t[T->x->width-1];
      switch(T->x->width - j) {
      case 8:
        t[j+7] = (word_cling_64_08_l(f0[j2]<< 0)>>7) | (word_cling_64_08_l(f1[j2]<< 0)>>6) | (word_cling_64_08_l(f2[j2]<< 0)>>5) | (word_cling_64_08_l(f3[j2]<< 0)>>4) \
          |      (word_cling_64_08_l(f4[j2]<< 0)>>3) | (word_cling_64_08_l(f5[j2]<< 0)>>2) | (word_cling_64_08_l(f6[j2]<< 0)>>1) | (word_cling_64_08_l(f7[j2]<< 0)>>0);
      case 7:
        t[j+6] = (word_cling_64_08_l(f0[j2]<< 8)>>7) | (word_cling_64_08_l(f1[j2]<< 8)>>6) | (word_cling_64_08_l(f2[j2]<< 8)>>5) | (word_cling_64_08_l(f3[j2]<< 8)>>4) \
          |      (word_cling_64_08_l(f4[j2]<< 8)>>3) | (word_cling_64_08_l(f5[j2]<< 8)>>2) | (word_cling_64_08_l(f6[j2]<< 8)>>1) | (word_cling_64_08_l(f7[j2]<< 8)>>0);
      case 6:
        t[j+5] = (word_cling_64_08_l(f0[j2]<<16)>>7) | (word_cling_64_08_l(f1[j2]<<16)>>6) | (word_cling_64_08_l(f2[j2]<<16)>>5) | (word_cling_64_08_l(f3[j2]<<16)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<16)>>3) | (word_cling_64_08_l(f5[j2]<<16)>>2) | (word_cling_64_08_l(f6[j2]<<16)>>1) | (word_cling_64_08_l(f7[j2]<<16)>>0);
      case 5:
        t[j+4] = (word_cling_64_08_l(f0[j2]<<24)>>7) | (word_cling_64_08_l(f1[j2]<<24)>>6) | (word_cling_64_08_l(f2[j2]<<24)>>5) | (word_cling_64_08_l(f3[j2]<<24)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<24)>>3) | (word_cling_64_08_l(f5[j2]<<24)>>2) | (word_cling_64_08_l(f6[j2]<<24)>>1) | (word_cling_64_08_l(f7[j2]<<24)>>0);
      case 4:
        t[j+3] = (word_cling_64_08_l(f0[j2]<<32)>>7) | (word_cling_64_08_l(f1[j2]<<32)>>6) | (word_cling_64_08_l(f2[j2]<<32)>>5) | (word_cling_64_08_l(f3[j2]<<32)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<32)>>3) | (word_cling_64_08_l(f5[j2]<<32)>>2) | (word_cling_64_08_l(f6[j2]<<32)>>1) | (word_cling_64_08_l(f7[j2]<<32)>>0);
      case 3:
        t[j+2] = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<40)>>3) | (word_cling_64_08_l(f5[j2]<<40)>>2) | (word_cling_64_08_l(f6[j2]<<40)>>1) | (word_cling_64_08_l(f7[j2]<<40)>>0);
      case 2:
        t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3) | (word_cling_64_08_l(f5[j2]<<48)>>2) | (word_cling_64_08_l(f6[j2]<<48)>>1) | (word_cling_64_08_l(f7[j2]<<48)>>0);
      case 1:
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3) | (word_cling_64_08_l(f5[j2]<<56)>>2) | (word_cling_64_08_l(f6[j2]<<56)>>1) | (word_cling_64_08_l(f7[j2]<<56)>>0);
        break;
      default:
        m4ri_die("impossible");
      }
      t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
    } // for loop
  }
    break;

  case 7: {
    for(rci_t i=0; i<T->nrows; i++) {
      const word *f0 = F->x[0]->rows[i];
      const word *f1 = F->x[1]->rows[i];
      const word *f2 = F->x[2]->rows[i];
      const word *f3 = F->x[3]->rows[i];
      const word *f4 = F->x[4]->rows[i];
      const word *f5 = F->x[5]->rows[i];
      const word *f6 = F->x[6]->rows[i];
      word *t  = T->x->rows[i];

      for(j=0, j2=0; j+8 < T->x->width; j+=8, j2++) {
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3) | (word_cling_64_08_l(f5[j2]<<56)>>2) | (word_cling_64_08_l(f6[j2]<<56)>>1);
        t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3) | (word_cling_64_08_l(f5[j2]<<48)>>2) | (word_cling_64_08_l(f6[j2]<<48)>>1);
        t[j+2] = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<40)>>3) | (word_cling_64_08_l(f5[j2]<<40)>>2) | (word_cling_64_08_l(f6[j2]<<40)>>1);
        t[j+3] = (word_cling_64_08_l(f0[j2]<<32)>>7) | (word_cling_64_08_l(f1[j2]<<32)>>6) | (word_cling_64_08_l(f2[j2]<<32)>>5) | (word_cling_64_08_l(f3[j2]<<32)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<32)>>3) | (word_cling_64_08_l(f5[j2]<<32)>>2) | (word_cling_64_08_l(f6[j2]<<32)>>1);
        t[j+4] = (word_cling_64_08_l(f0[j2]<<24)>>7) | (word_cling_64_08_l(f1[j2]<<24)>>6) | (word_cling_64_08_l(f2[j2]<<24)>>5) | (word_cling_64_08_l(f3[j2]<<24)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<24)>>3) | (word_cling_64_08_l(f5[j2]<<24)>>2) | (word_cling_64_08_l(f6[j2]<<24)>>1);
        t[j+5] = (word_cling_64_08_l(f0[j2]<<16)>>7) | (word_cling_64_08_l(f1[j2]<<16)>>6) | (word_cling_64_08_l(f2[j2]<<16)>>5) | (word_cling_64_08_l(f3[j2]<<16)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<16)>>3) | (word_cling_64_08_l(f5[j2]<<16)>>2) | (word_cling_64_08_l(f6[j2]<<16)>>1);
        t[j+6] = (word_cling_64_08_l(f0[j2]<< 8)>>7) | (word_cling_64_08_l(f1[j2]<< 8)>>6) | (word_cling_64_08_l(f2[j2]<< 8)>>5) | (word_cling_64_08_l(f3[j2]<< 8)>>4) \
          |      (word_cling_64_08_l(f4[j2]<< 8)>>3) | (word_cling_64_08_l(f5[j2]<< 8)>>2) | (word_cling_64_08_l(f6[j2]<< 8)>>1);
        t[j+7] = (word_cling_64_08_l(f0[j2]<< 0)>>7) | (word_cling_64_08_l(f1[j2]<< 0)>>6) | (word_cling_64_08_l(f2[j2]<< 0)>>5) | (word_cling_64_08_l(f3[j2]<< 0)>>4) \
          |      (word_cling_64_08_l(f4[j2]<< 0)>>3) | (word_cling_64_08_l(f5[j2]<< 0)>>2) | (word_cling_64_08_l(f6[j2]<< 0)>>1);
      }

      register word tmp= t[T->x->width-1];
      switch(T->x->width - j) {
      case 8:
        t[j+7] = (word_cling_64_08_l(f0[j2]<< 0)>>7) | (word_cling_64_08_l(f1[j2]<< 0)>>6) | (word_cling_64_08_l(f2[j2]<< 0)>>5) | (word_cling_64_08_l(f3[j2]<< 0)>>4) \
          |      (word_cling_64_08_l(f4[j2]<< 0)>>3) | (word_cling_64_08_l(f5[j2]<< 0)>>2) | (word_cling_64_08_l(f6[j2]<< 0)>>1);
      case 7:
        t[j+6] = (word_cling_64_08_l(f0[j2]<< 8)>>7) | (word_cling_64_08_l(f1[j2]<< 8)>>6) | (word_cling_64_08_l(f2[j2]<< 8)>>5) | (word_cling_64_08_l(f3[j2]<< 8)>>4) \
          |      (word_cling_64_08_l(f4[j2]<< 8)>>3) | (word_cling_64_08_l(f5[j2]<< 8)>>2) | (word_cling_64_08_l(f6[j2]<< 8)>>1);
      case 6:
        t[j+5] = (word_cling_64_08_l(f0[j2]<<16)>>7) | (word_cling_64_08_l(f1[j2]<<16)>>6) | (word_cling_64_08_l(f2[j2]<<16)>>5) | (word_cling_64_08_l(f3[j2]<<16)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<16)>>3) | (word_cling_64_08_l(f5[j2]<<16)>>2) | (word_cling_64_08_l(f6[j2]<<16)>>1);
      case 5:
        t[j+4] = (word_cling_64_08_l(f0[j2]<<24)>>7) | (word_cling_64_08_l(f1[j2]<<24)>>6) | (word_cling_64_08_l(f2[j2]<<24)>>5) | (word_cling_64_08_l(f3[j2]<<24)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<24)>>3) | (word_cling_64_08_l(f5[j2]<<24)>>2) | (word_cling_64_08_l(f6[j2]<<24)>>1);
      case 4:
        t[j+3] = (word_cling_64_08_l(f0[j2]<<32)>>7) | (word_cling_64_08_l(f1[j2]<<32)>>6) | (word_cling_64_08_l(f2[j2]<<32)>>5) | (word_cling_64_08_l(f3[j2]<<32)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<32)>>3) | (word_cling_64_08_l(f5[j2]<<32)>>2) | (word_cling_64_08_l(f6[j2]<<32)>>1);
      case 3:
        t[j+2] = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<40)>>3) | (word_cling_64_08_l(f5[j2]<<40)>>2) | (word_cling_64_08_l(f6[j2]<<40)>>1);
      case 2:
        t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3) | (word_cling_64_08_l(f5[j2]<<48)>>2) | (word_cling_64_08_l(f6[j2]<<48)>>1);
      case 1:
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3) | (word_cling_64_08_l(f5[j2]<<56)>>2) | (word_cling_64_08_l(f6[j2]<<56)>>1);
        break;
      default:
        m4ri_die("impossible");
      }
      t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
    }
  }
    break;

  case 6: {
    for(rci_t i=0; i<T->nrows; i++) {
      const word *f0 = F->x[0]->rows[i];
      const word *f1 = F->x[1]->rows[i];
      const word *f2 = F->x[2]->rows[i];
      const word *f3 = F->x[3]->rows[i];
      const word *f4 = F->x[4]->rows[i];
      const word *f5 = F->x[5]->rows[i];
      word *t  = T->x->rows[i];

      for(j=0, j2=0; j+8 < T->x->width; j+=8, j2++) {
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3) | (word_cling_64_08_l(f5[j2]<<56)>>2);
        t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3) | (word_cling_64_08_l(f5[j2]<<48)>>2);
        t[j+2] = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<40)>>3) | (word_cling_64_08_l(f5[j2]<<40)>>2);
        t[j+3] = (word_cling_64_08_l(f0[j2]<<32)>>7) | (word_cling_64_08_l(f1[j2]<<32)>>6) | (word_cling_64_08_l(f2[j2]<<32)>>5) | (word_cling_64_08_l(f3[j2]<<32)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<32)>>3) | (word_cling_64_08_l(f5[j2]<<32)>>2);
        t[j+4] = (word_cling_64_08_l(f0[j2]<<24)>>7) | (word_cling_64_08_l(f1[j2]<<24)>>6) | (word_cling_64_08_l(f2[j2]<<24)>>5) | (word_cling_64_08_l(f3[j2]<<24)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<24)>>3) | (word_cling_64_08_l(f5[j2]<<24)>>2);
        t[j+5] = (word_cling_64_08_l(f0[j2]<<16)>>7) | (word_cling_64_08_l(f1[j2]<<16)>>6) | (word_cling_64_08_l(f2[j2]<<16)>>5) | (word_cling_64_08_l(f3[j2]<<16)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<16)>>3) | (word_cling_64_08_l(f5[j2]<<16)>>2);
        t[j+6] = (word_cling_64_08_l(f0[j2]<< 8)>>7) | (word_cling_64_08_l(f1[j2]<< 8)>>6) | (word_cling_64_08_l(f2[j2]<< 8)>>5) | (word_cling_64_08_l(f3[j2]<< 8)>>4) \
          |      (word_cling_64_08_l(f4[j2]<< 8)>>3) | (word_cling_64_08_l(f5[j2]<< 8)>>2);
        t[j+7] = (word_cling_64_08_l(f0[j2]<< 0)>>7) | (word_cling_64_08_l(f1[j2]<< 0)>>6) | (word_cling_64_08_l(f2[j2]<< 0)>>5) | (word_cling_64_08_l(f3[j2]<< 0)>>4) \
          |      (word_cling_64_08_l(f4[j2]<< 0)>>3) | (word_cling_64_08_l(f5[j2]<< 0)>>2);
      }

      register word tmp = t[T->x->width-1];
      switch(T->x->width - j) {
      case 8:
        t[j+7] = (word_cling_64_08_l(f0[j2]<< 0)>>7) | (word_cling_64_08_l(f1[j2]<< 0)>>6) | (word_cling_64_08_l(f2[j2]<< 0)>>5) | (word_cling_64_08_l(f3[j2]<< 0)>>4) \
          |      (word_cling_64_08_l(f4[j2]<< 0)>>3) | (word_cling_64_08_l(f5[j2]<< 0)>>2);
      case 7:
        t[j+6] = (word_cling_64_08_l(f0[j2]<< 8)>>7) | (word_cling_64_08_l(f1[j2]<< 8)>>6) | (word_cling_64_08_l(f2[j2]<< 8)>>5) | (word_cling_64_08_l(f3[j2]<< 8)>>4) \
          |      (word_cling_64_08_l(f4[j2]<< 8)>>3) | (word_cling_64_08_l(f5[j2]<< 8)>>2);
      case 6:
        t[j+5] = (word_cling_64_08_l(f0[j2]<<16)>>7) | (word_cling_64_08_l(f1[j2]<<16)>>6) | (word_cling_64_08_l(f2[j2]<<16)>>5) | (word_cling_64_08_l(f3[j2]<<16)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<16)>>3) | (word_cling_64_08_l(f5[j2]<<16)>>2);
      case 5:
        t[j+4] = (word_cling_64_08_l(f0[j2]<<24)>>7) | (word_cling_64_08_l(f1[j2]<<24)>>6) | (word_cling_64_08_l(f2[j2]<<24)>>5) | (word_cling_64_08_l(f3[j2]<<24)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<24)>>3) | (word_cling_64_08_l(f5[j2]<<24)>>2);
      case 4:
        t[j+3] = (word_cling_64_08_l(f0[j2]<<32)>>7) | (word_cling_64_08_l(f1[j2]<<32)>>6) | (word_cling_64_08_l(f2[j2]<<32)>>5) | (word_cling_64_08_l(f3[j2]<<32)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<32)>>3) | (word_cling_64_08_l(f5[j2]<<32)>>2);
      case 3:
        t[j+2] = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<40)>>3) | (word_cling_64_08_l(f5[j2]<<40)>>2);
      case 2:
        t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3) | (word_cling_64_08_l(f5[j2]<<48)>>2);
      case 1:
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3) | (word_cling_64_08_l(f5[j2]<<56)>>2);
        break;
      default:
        m4ri_die("impossible");
      }
      t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
    }
  }
    break;

  case 5: {
    for(rci_t i=0; i<T->nrows; i++) {
      const word *f0 = F->x[0]->rows[i];
      const word *f1 = F->x[1]->rows[i];
      const word *f2 = F->x[2]->rows[i];
      const word *f3 = F->x[3]->rows[i];
      const word *f4 = F->x[4]->rows[i];
      word *t  = T->x->rows[i];

      for(j=0, j2=0; j+8 < T->x->width; j+=8, j2++) {
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) | (word_cling_64_08_l(f4[j2]<<56)>>3);
        t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) | (word_cling_64_08_l(f4[j2]<<48)>>3);
        t[j+2] = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) | (word_cling_64_08_l(f4[j2]<<40)>>3);
        t[j+3] = (word_cling_64_08_l(f0[j2]<<32)>>7) | (word_cling_64_08_l(f1[j2]<<32)>>6) | (word_cling_64_08_l(f2[j2]<<32)>>5) | (word_cling_64_08_l(f3[j2]<<32)>>4) | (word_cling_64_08_l(f4[j2]<<32)>>3);
        t[j+4] = (word_cling_64_08_l(f0[j2]<<24)>>7) | (word_cling_64_08_l(f1[j2]<<24)>>6) | (word_cling_64_08_l(f2[j2]<<24)>>5) | (word_cling_64_08_l(f3[j2]<<24)>>4) | (word_cling_64_08_l(f4[j2]<<24)>>3);
        t[j+5] = (word_cling_64_08_l(f0[j2]<<16)>>7) | (word_cling_64_08_l(f1[j2]<<16)>>6) | (word_cling_64_08_l(f2[j2]<<16)>>5) | (word_cling_64_08_l(f3[j2]<<16)>>4) | (word_cling_64_08_l(f4[j2]<<16)>>3);
        t[j+6] = (word_cling_64_08_l(f0[j2]<< 8)>>7) | (word_cling_64_08_l(f1[j2]<< 8)>>6) | (word_cling_64_08_l(f2[j2]<< 8)>>5) | (word_cling_64_08_l(f3[j2]<< 8)>>4) | (word_cling_64_08_l(f4[j2]<< 8)>>3);
        t[j+7] = (word_cling_64_08_l(f0[j2]<< 0)>>7) | (word_cling_64_08_l(f1[j2]<< 0)>>6) | (word_cling_64_08_l(f2[j2]<< 0)>>5) | (word_cling_64_08_l(f3[j2]<< 0)>>4) | (word_cling_64_08_l(f4[j2]<< 0)>>3);
      }

      register word tmp = t[T->x->width - 1];
      switch(T->x->width - j) {
      case 8: t[j+7] = (word_cling_64_08_l(f0[j2]<< 0)>>7) | (word_cling_64_08_l(f1[j2]<< 0)>>6) | (word_cling_64_08_l(f2[j2]<< 0)>>5) | (word_cling_64_08_l(f3[j2]<< 0)>>4) | (word_cling_64_08_l(f4[j2]<< 0)>>3);
      case 7: t[j+6] = (word_cling_64_08_l(f0[j2]<< 8)>>7) | (word_cling_64_08_l(f1[j2]<< 8)>>6) | (word_cling_64_08_l(f2[j2]<< 8)>>5) | (word_cling_64_08_l(f3[j2]<< 8)>>4) | (word_cling_64_08_l(f4[j2]<< 8)>>3);
      case 6: t[j+5] = (word_cling_64_08_l(f0[j2]<<16)>>7) | (word_cling_64_08_l(f1[j2]<<16)>>6) | (word_cling_64_08_l(f2[j2]<<16)>>5) | (word_cling_64_08_l(f3[j2]<<16)>>4) | (word_cling_64_08_l(f4[j2]<<16)>>3);
      case 5: t[j+4] = (word_cling_64_08_l(f0[j2]<<24)>>7) | (word_cling_64_08_l(f1[j2]<<24)>>6) | (word_cling_64_08_l(f2[j2]<<24)>>5) | (word_cling_64_08_l(f3[j2]<<24)>>4) | (word_cling_64_08_l(f4[j2]<<24)>>3);
      case 4: t[j+3] = (word_cling_64_08_l(f0[j2]<<32)>>7) | (word_cling_64_08_l(f1[j2]<<32)>>6) | (word_cling_64_08_l(f2[j2]<<32)>>5) | (word_cling_64_08_l(f3[j2]<<32)>>4) | (word_cling_64_08_l(f4[j2]<<32)>>3);
      case 3: t[j+2] = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) | (word_cling_64_08_l(f4[j2]<<40)>>3);
      case 2: t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) | (word_cling_64_08_l(f4[j2]<<48)>>3);
      case 1: t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) | (word_cling_64_08_l(f4[j2]<<56)>>3);
        break;
      default:
        m4ri_die("impossible");
      }
      t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
    }
  }
    break;
  default:
    m4ri_die("impossible");
  }
  return T;
}

/* we define these things to keep code compact below. */

#define word_slice_64_16_l_combine_bulk(T, Ti, F, Fi, shift)  \
   T[Ti] |= word_slice_64_16_l(F[Fi+ 0]<<shift & x80008000)>>60 | word_slice_64_16_l(F[Fi+ 1]<<shift & x80008000)>>56 \
    |       word_slice_64_16_l(F[Fi+ 2]<<shift & x80008000)>>52 | word_slice_64_16_l(F[Fi+ 3]<<shift & x80008000)>>48 \
    |       word_slice_64_16_l(F[Fi+ 4]<<shift & x80008000)>>44 | word_slice_64_16_l(F[Fi+ 5]<<shift & x80008000)>>40 \
    |       word_slice_64_16_l(F[Fi+ 6]<<shift & x80008000)>>36 | word_slice_64_16_l(F[Fi+ 7]<<shift & x80008000)>>32 \
    |       word_slice_64_16_l(F[Fi+ 8]<<shift & x80008000)>>28 | word_slice_64_16_l(F[Fi+ 9]<<shift & x80008000)>>24 \
    |       word_slice_64_16_l(F[Fi+10]<<shift & x80008000)>>20 | word_slice_64_16_l(F[Fi+11]<<shift & x80008000)>>16 \
    |       word_slice_64_16_l(F[Fi+12]<<shift & x80008000)>>12 | word_slice_64_16_l(F[Fi+13]<<shift & x80008000)>> 8 \
    |       word_slice_64_16_l(F[Fi+14]<<shift & x80008000)>> 4 | word_slice_64_16_l(F[Fi+15]<<shift & x80008000)>> 0;

#define word_slice_64_16_l_slice_rest(F, Fi, shift)         \
  r0 |= word_slice_64_16_l(F[Fi]<<15 & x80008000)>> shift;         \
  r1 |= word_slice_64_16_l(F[Fi]<<14 & x80008000)>> shift;         \
  r2 |= word_slice_64_16_l(F[Fi]<<13 & x80008000)>> shift;         \
  r3 |= word_slice_64_16_l(F[Fi]<<12 & x80008000)>> shift;         \
  r4 |= word_slice_64_16_l(F[Fi]<<11 & x80008000)>> shift;         \
  r5 |= word_slice_64_16_l(F[Fi]<<10 & x80008000)>> shift;         \
  r6 |= word_slice_64_16_l(F[Fi]<< 9 & x80008000)>> shift;         \
  r7 |= word_slice_64_16_l(F[Fi]<< 8 & x80008000)>> shift;

mzd_slice_t *_mzed_slice16(mzd_slice_t *T, const mzed_t *F) {
  assert(T && (8 < T->depth && T->depth <= 16) && T->x[0]->offset == 0);
  size_t j, j2 = 0;
  register word r0,r1,r2,r3,r4,r5,r6,r7 = 0;

  const word bitmask_end = __M4RI_LEFT_BITMASK((T->x[0]->offset + T->ncols) % m4ri_radix);

  if (mzed_is_zero(F))
    return T;

  /* we do multiple runs over T to make the code more compact, we start by doing the first eight
     bits */

  for(size_t i=0; i<T->nrows; i++) {
    word *t0 = T->x[0]->rows[i];
    word *t1 = T->x[1]->rows[i];
    word *t2 = T->x[2]->rows[i];
    word *t3 = T->x[3]->rows[i];
    word *t4 = T->x[4]->rows[i];
    word *t5 = T->x[5]->rows[i];
    word *t6 = T->x[6]->rows[i];
    word *t7 = T->x[7]->rows[i];
    const word const *f  = F->x->rows[i];

    /* bulk of work */
    for(j=0, j2=0; j+16 < F->x->width; j+=16,j2++) {
      word_slice_64_16_l_combine_bulk(t0, j2, f, j, 15);
      word_slice_64_16_l_combine_bulk(t1, j2, f, j, 14);
      word_slice_64_16_l_combine_bulk(t2, j2, f, j, 13);
      word_slice_64_16_l_combine_bulk(t3, j2, f, j, 12);
      word_slice_64_16_l_combine_bulk(t4, j2, f, j, 11);
      word_slice_64_16_l_combine_bulk(t5, j2, f, j, 10);
      word_slice_64_16_l_combine_bulk(t6, j2, f, j,  9);
      word_slice_64_16_l_combine_bulk(t7, j2, f, j,  8);
    }
    r0 = r1 = r2 = r3 = r4 = r5 = r6 = r7 = 0;
    switch(F->x->width - j) {
    case 16: word_slice_64_16_l_slice_rest(f, j+15,  0);
    case 15: word_slice_64_16_l_slice_rest(f, j+14,  4);
    case 14: word_slice_64_16_l_slice_rest(f, j+13,  8);
    case 13: word_slice_64_16_l_slice_rest(f, j+12, 12);
    case 12: word_slice_64_16_l_slice_rest(f, j+11, 16);
    case 11: word_slice_64_16_l_slice_rest(f, j+10, 20);
    case 10: word_slice_64_16_l_slice_rest(f, j+ 9, 24);
    case  9: word_slice_64_16_l_slice_rest(f, j+ 8, 28);
    case  8: word_slice_64_16_l_slice_rest(f, j+ 7, 32);
    case  7: word_slice_64_16_l_slice_rest(f, j+ 6, 36);
    case  6: word_slice_64_16_l_slice_rest(f, j+ 5, 40);
    case  5: word_slice_64_16_l_slice_rest(f, j+ 4, 44);
    case  4: word_slice_64_16_l_slice_rest(f, j+ 3, 48);
    case  3: word_slice_64_16_l_slice_rest(f, j+ 2, 52);
    case  2: word_slice_64_16_l_slice_rest(f, j+ 1, 56);
    case  1: word_slice_64_16_l_slice_rest(f, j+ 0, 60);
      break;
    default:
      m4ri_die("impossible");
    }
    t0[j2] |= r0 & bitmask_end;
    t1[j2] |= r1 & bitmask_end;
    t2[j2] |= r2 & bitmask_end;
    t3[j2] |= r3 & bitmask_end;
    t4[j2] |= r4 & bitmask_end;
    t5[j2] |= r5 & bitmask_end;
    t6[j2] |= r6 & bitmask_end;
    t7[j2] |= r7 & bitmask_end;
  }
  if(T->depth >= 12) {
    for(size_t i=0; i<T->nrows; i++) {
      word *t0 = T->x[ 8]->rows[i];
      word *t1 = T->x[ 9]->rows[i];
      word *t2 = T->x[10]->rows[i];
      word *t3 = T->x[11]->rows[i];
      const word const *f  = F->x->rows[i];

      /* bulk of work */
      for(j=0, j2=0; j+16 < F->x->width; j+=16,j2++) {
        word_slice_64_16_l_combine_bulk(t0, j2, f, j,  7);
        word_slice_64_16_l_combine_bulk(t1, j2, f, j,  6);
        word_slice_64_16_l_combine_bulk(t2, j2, f, j,  5);
        word_slice_64_16_l_combine_bulk(t3, j2, f, j,  4);
      }
      r0 = r1 = r2 = r3 = 0;
      switch(F->x->width - j) {
      case 16: r0 |= word_slice_64_16_l(f[j+15]<< 7 & x80008000)>>  0; r1 |= word_slice_64_16_l(f[j+15]<< 6 & x80008000)>>  0; r2 |= word_slice_64_16_l(f[j+15]<< 5 & x80008000)>>  0; r3 |= word_slice_64_16_l(f[j+15]<< 4 & x80008000)>>  0;
      case 15: r0 |= word_slice_64_16_l(f[j+14]<< 7 & x80008000)>>  4; r1 |= word_slice_64_16_l(f[j+14]<< 6 & x80008000)>>  4; r2 |= word_slice_64_16_l(f[j+14]<< 5 & x80008000)>>  4; r3 |= word_slice_64_16_l(f[j+14]<< 4 & x80008000)>>  4;
      case 14: r0 |= word_slice_64_16_l(f[j+13]<< 7 & x80008000)>>  8; r1 |= word_slice_64_16_l(f[j+13]<< 6 & x80008000)>>  8; r2 |= word_slice_64_16_l(f[j+13]<< 5 & x80008000)>>  8; r3 |= word_slice_64_16_l(f[j+13]<< 4 & x80008000)>>  8;
      case 13: r0 |= word_slice_64_16_l(f[j+12]<< 7 & x80008000)>> 12; r1 |= word_slice_64_16_l(f[j+12]<< 6 & x80008000)>> 12; r2 |= word_slice_64_16_l(f[j+12]<< 5 & x80008000)>> 12; r3 |= word_slice_64_16_l(f[j+12]<< 4 & x80008000)>> 12;
      case 12: r0 |= word_slice_64_16_l(f[j+11]<< 7 & x80008000)>> 16; r1 |= word_slice_64_16_l(f[j+11]<< 6 & x80008000)>> 16; r2 |= word_slice_64_16_l(f[j+11]<< 5 & x80008000)>> 16; r3 |= word_slice_64_16_l(f[j+11]<< 4 & x80008000)>> 16;
      case 11: r0 |= word_slice_64_16_l(f[j+10]<< 7 & x80008000)>> 20; r1 |= word_slice_64_16_l(f[j+10]<< 6 & x80008000)>> 20; r2 |= word_slice_64_16_l(f[j+10]<< 5 & x80008000)>> 20; r3 |= word_slice_64_16_l(f[j+10]<< 4 & x80008000)>> 20;
      case 10: r0 |= word_slice_64_16_l(f[j+ 9]<< 7 & x80008000)>> 24; r1 |= word_slice_64_16_l(f[j+ 9]<< 6 & x80008000)>> 24; r2 |= word_slice_64_16_l(f[j+ 9]<< 5 & x80008000)>> 24; r3 |= word_slice_64_16_l(f[j+ 9]<< 4 & x80008000)>> 24;
      case  9: r0 |= word_slice_64_16_l(f[j+ 8]<< 7 & x80008000)>> 28; r1 |= word_slice_64_16_l(f[j+ 8]<< 6 & x80008000)>> 28; r2 |= word_slice_64_16_l(f[j+ 8]<< 5 & x80008000)>> 28; r3 |= word_slice_64_16_l(f[j+ 8]<< 4 & x80008000)>> 28;
      case  8: r0 |= word_slice_64_16_l(f[j+ 7]<< 7 & x80008000)>> 32; r1 |= word_slice_64_16_l(f[j+ 7]<< 6 & x80008000)>> 32; r2 |= word_slice_64_16_l(f[j+ 7]<< 5 & x80008000)>> 32; r3 |= word_slice_64_16_l(f[j+ 7]<< 4 & x80008000)>> 32;
      case  7: r0 |= word_slice_64_16_l(f[j+ 6]<< 7 & x80008000)>> 36; r1 |= word_slice_64_16_l(f[j+ 6]<< 6 & x80008000)>> 36; r2 |= word_slice_64_16_l(f[j+ 6]<< 5 & x80008000)>> 36; r3 |= word_slice_64_16_l(f[j+ 6]<< 4 & x80008000)>> 36;
      case  6: r0 |= word_slice_64_16_l(f[j+ 5]<< 7 & x80008000)>> 40; r1 |= word_slice_64_16_l(f[j+ 5]<< 6 & x80008000)>> 40; r2 |= word_slice_64_16_l(f[j+ 5]<< 5 & x80008000)>> 40; r3 |= word_slice_64_16_l(f[j+ 5]<< 4 & x80008000)>> 40;
      case  5: r0 |= word_slice_64_16_l(f[j+ 4]<< 7 & x80008000)>> 44; r1 |= word_slice_64_16_l(f[j+ 4]<< 6 & x80008000)>> 44; r2 |= word_slice_64_16_l(f[j+ 4]<< 5 & x80008000)>> 44; r3 |= word_slice_64_16_l(f[j+ 4]<< 4 & x80008000)>> 44;
      case  4: r0 |= word_slice_64_16_l(f[j+ 3]<< 7 & x80008000)>> 48; r1 |= word_slice_64_16_l(f[j+ 3]<< 6 & x80008000)>> 48; r2 |= word_slice_64_16_l(f[j+ 3]<< 5 & x80008000)>> 48; r3 |= word_slice_64_16_l(f[j+ 3]<< 4 & x80008000)>> 48;
      case  3: r0 |= word_slice_64_16_l(f[j+ 2]<< 7 & x80008000)>> 52; r1 |= word_slice_64_16_l(f[j+ 2]<< 6 & x80008000)>> 52; r2 |= word_slice_64_16_l(f[j+ 2]<< 5 & x80008000)>> 52; r3 |= word_slice_64_16_l(f[j+ 2]<< 4 & x80008000)>> 52;
      case  2: r0 |= word_slice_64_16_l(f[j+ 1]<< 7 & x80008000)>> 56; r1 |= word_slice_64_16_l(f[j+ 1]<< 6 & x80008000)>> 56; r2 |= word_slice_64_16_l(f[j+ 1]<< 5 & x80008000)>> 56; r3 |= word_slice_64_16_l(f[j+ 1]<< 4 & x80008000)>> 56;
      case  1: r0 |= word_slice_64_16_l(f[j+ 0]<< 7 & x80008000)>> 60; r1 |= word_slice_64_16_l(f[j+ 0]<< 6 & x80008000)>> 60; r2 |= word_slice_64_16_l(f[j+ 0]<< 5 & x80008000)>> 60; r3 |= word_slice_64_16_l(f[j+ 0]<< 4 & x80008000)>> 60;
        break;
      default:
        m4ri_die("impossible");
      }
    t0[j2] |= r0 & bitmask_end;
    t1[j2] |= r1 & bitmask_end;
    t2[j2] |= r2 & bitmask_end;
    t3[j2] |= r3 & bitmask_end;
    }

    switch(T->depth) {
    case 16: {
      for(size_t i=0; i<T->nrows; i++) {
        word *t0 = T->x[12]->rows[i];
        word *t1 = T->x[13]->rows[i];
        word *t2 = T->x[14]->rows[i];
        word *t3 = T->x[15]->rows[i];
        const word const *f  = F->x->rows[i];

        /* bulk of work */
        for(j=0, j2=0; j+16 < F->x->width; j+=16,j2++) {
          word_slice_64_16_l_combine_bulk(t0, j2, f, j,  3);
          word_slice_64_16_l_combine_bulk(t1, j2, f, j,  2);
          word_slice_64_16_l_combine_bulk(t2, j2, f, j,  1);
          word_slice_64_16_l_combine_bulk(t3, j2, f, j,  0);
        }
        r0 = r1 = r2 = r3 = 0;
        switch(F->x->width - j) {
        case 16: r0 |= word_slice_64_16_l(f[j+15]<< 3 & x80008000)>>  0; r1 |= word_slice_64_16_l(f[j+15]<< 2 & x80008000)>>  0; r2 |= word_slice_64_16_l(f[j+15]<< 1 & x80008000)>>  0; r3 |= word_slice_64_16_l(f[j+15]<< 0 & x80008000)>>  0;
        case 15: r0 |= word_slice_64_16_l(f[j+14]<< 3 & x80008000)>>  4; r1 |= word_slice_64_16_l(f[j+14]<< 2 & x80008000)>>  4; r2 |= word_slice_64_16_l(f[j+14]<< 1 & x80008000)>>  4; r3 |= word_slice_64_16_l(f[j+14]<< 0 & x80008000)>>  4;
        case 14: r0 |= word_slice_64_16_l(f[j+13]<< 3 & x80008000)>>  8; r1 |= word_slice_64_16_l(f[j+13]<< 2 & x80008000)>>  8; r2 |= word_slice_64_16_l(f[j+13]<< 1 & x80008000)>>  8; r3 |= word_slice_64_16_l(f[j+13]<< 0 & x80008000)>>  8;
        case 13: r0 |= word_slice_64_16_l(f[j+12]<< 3 & x80008000)>> 12; r1 |= word_slice_64_16_l(f[j+12]<< 2 & x80008000)>> 12; r2 |= word_slice_64_16_l(f[j+12]<< 1 & x80008000)>> 12; r3 |= word_slice_64_16_l(f[j+12]<< 0 & x80008000)>> 12;
        case 12: r0 |= word_slice_64_16_l(f[j+11]<< 3 & x80008000)>> 16; r1 |= word_slice_64_16_l(f[j+11]<< 2 & x80008000)>> 16; r2 |= word_slice_64_16_l(f[j+11]<< 1 & x80008000)>> 16; r3 |= word_slice_64_16_l(f[j+11]<< 0 & x80008000)>> 16;
        case 11: r0 |= word_slice_64_16_l(f[j+10]<< 3 & x80008000)>> 20; r1 |= word_slice_64_16_l(f[j+10]<< 2 & x80008000)>> 20; r2 |= word_slice_64_16_l(f[j+10]<< 1 & x80008000)>> 20; r3 |= word_slice_64_16_l(f[j+10]<< 0 & x80008000)>> 20;
        case 10: r0 |= word_slice_64_16_l(f[j+ 9]<< 3 & x80008000)>> 24; r1 |= word_slice_64_16_l(f[j+ 9]<< 2 & x80008000)>> 24; r2 |= word_slice_64_16_l(f[j+ 9]<< 1 & x80008000)>> 24; r3 |= word_slice_64_16_l(f[j+ 9]<< 0 & x80008000)>> 24;
        case  9: r0 |= word_slice_64_16_l(f[j+ 8]<< 3 & x80008000)>> 28; r1 |= word_slice_64_16_l(f[j+ 8]<< 2 & x80008000)>> 28; r2 |= word_slice_64_16_l(f[j+ 8]<< 1 & x80008000)>> 28; r3 |= word_slice_64_16_l(f[j+ 8]<< 0 & x80008000)>> 28;
        case  8: r0 |= word_slice_64_16_l(f[j+ 7]<< 3 & x80008000)>> 32; r1 |= word_slice_64_16_l(f[j+ 7]<< 2 & x80008000)>> 32; r2 |= word_slice_64_16_l(f[j+ 7]<< 1 & x80008000)>> 32; r3 |= word_slice_64_16_l(f[j+ 7]<< 0 & x80008000)>> 32;
        case  7: r0 |= word_slice_64_16_l(f[j+ 6]<< 3 & x80008000)>> 36; r1 |= word_slice_64_16_l(f[j+ 6]<< 2 & x80008000)>> 36; r2 |= word_slice_64_16_l(f[j+ 6]<< 1 & x80008000)>> 36; r3 |= word_slice_64_16_l(f[j+ 6]<< 0 & x80008000)>> 36;
        case  6: r0 |= word_slice_64_16_l(f[j+ 5]<< 3 & x80008000)>> 40; r1 |= word_slice_64_16_l(f[j+ 5]<< 2 & x80008000)>> 40; r2 |= word_slice_64_16_l(f[j+ 5]<< 1 & x80008000)>> 40; r3 |= word_slice_64_16_l(f[j+ 5]<< 0 & x80008000)>> 40;
        case  5: r0 |= word_slice_64_16_l(f[j+ 4]<< 3 & x80008000)>> 44; r1 |= word_slice_64_16_l(f[j+ 4]<< 2 & x80008000)>> 44; r2 |= word_slice_64_16_l(f[j+ 4]<< 1 & x80008000)>> 44; r3 |= word_slice_64_16_l(f[j+ 4]<< 0 & x80008000)>> 44;
        case  4: r0 |= word_slice_64_16_l(f[j+ 3]<< 3 & x80008000)>> 48; r1 |= word_slice_64_16_l(f[j+ 3]<< 2 & x80008000)>> 48; r2 |= word_slice_64_16_l(f[j+ 3]<< 1 & x80008000)>> 48; r3 |= word_slice_64_16_l(f[j+ 3]<< 0 & x80008000)>> 48;
        case  3: r0 |= word_slice_64_16_l(f[j+ 2]<< 3 & x80008000)>> 52; r1 |= word_slice_64_16_l(f[j+ 2]<< 2 & x80008000)>> 52; r2 |= word_slice_64_16_l(f[j+ 2]<< 1 & x80008000)>> 52; r3 |= word_slice_64_16_l(f[j+ 2]<< 0 & x80008000)>> 52;
        case  2: r0 |= word_slice_64_16_l(f[j+ 1]<< 3 & x80008000)>> 56; r1 |= word_slice_64_16_l(f[j+ 1]<< 2 & x80008000)>> 56; r2 |= word_slice_64_16_l(f[j+ 1]<< 1 & x80008000)>> 56; r3 |= word_slice_64_16_l(f[j+ 1]<< 0 & x80008000)>> 56;
        case  1: r0 |= word_slice_64_16_l(f[j+ 0]<< 3 & x80008000)>> 60; r1 |= word_slice_64_16_l(f[j+ 0]<< 2 & x80008000)>> 60; r2 |= word_slice_64_16_l(f[j+ 0]<< 1 & x80008000)>> 60; r3 |= word_slice_64_16_l(f[j+ 0]<< 0 & x80008000)>> 60;
          break;
        default:
          m4ri_die("impossible");
        }
        t0[j2] |= r0 & bitmask_end;
        t1[j2] |= r1 & bitmask_end;
        t2[j2] |= r2 & bitmask_end;
        t3[j2] |= r3 & bitmask_end;
      }
    } break;
    case 15: {
      for(size_t i=0; i<T->nrows; i++) {
        word *t0 = T->x[12]->rows[i];
        word *t1 = T->x[13]->rows[i];
        word *t2 = T->x[14]->rows[i];
        const word const *f  = F->x->rows[i];

        /* bulk of work */
        for(j=0, j2=0; j+16 < F->x->width; j+=16,j2++) {
          word_slice_64_16_l_combine_bulk(t0, j2, f, j,  3);
          word_slice_64_16_l_combine_bulk(t1, j2, f, j,  2);
          word_slice_64_16_l_combine_bulk(t2, j2, f, j,  1);
        }
        r0 = r1 = r2 = 0;
        switch(F->x->width - j) {
        case 16: r0 |= word_slice_64_16_l(f[j+15]<< 3 & x80008000)>>  0; r1 |= word_slice_64_16_l(f[j+15]<< 2 & x80008000)>>  0; r2 |= word_slice_64_16_l(f[j+15]<< 1 & x80008000)>>  0;
        case 15: r0 |= word_slice_64_16_l(f[j+14]<< 3 & x80008000)>>  4; r1 |= word_slice_64_16_l(f[j+14]<< 2 & x80008000)>>  4; r2 |= word_slice_64_16_l(f[j+14]<< 1 & x80008000)>>  4;
        case 14: r0 |= word_slice_64_16_l(f[j+13]<< 3 & x80008000)>>  8; r1 |= word_slice_64_16_l(f[j+13]<< 2 & x80008000)>>  8; r2 |= word_slice_64_16_l(f[j+13]<< 1 & x80008000)>>  8;
        case 13: r0 |= word_slice_64_16_l(f[j+12]<< 3 & x80008000)>> 12; r1 |= word_slice_64_16_l(f[j+12]<< 2 & x80008000)>> 12; r2 |= word_slice_64_16_l(f[j+12]<< 1 & x80008000)>> 12;
        case 12: r0 |= word_slice_64_16_l(f[j+11]<< 3 & x80008000)>> 16; r1 |= word_slice_64_16_l(f[j+11]<< 2 & x80008000)>> 16; r2 |= word_slice_64_16_l(f[j+11]<< 1 & x80008000)>> 16;
        case 11: r0 |= word_slice_64_16_l(f[j+10]<< 3 & x80008000)>> 20; r1 |= word_slice_64_16_l(f[j+10]<< 2 & x80008000)>> 20; r2 |= word_slice_64_16_l(f[j+10]<< 1 & x80008000)>> 20;
        case 10: r0 |= word_slice_64_16_l(f[j+ 9]<< 3 & x80008000)>> 24; r1 |= word_slice_64_16_l(f[j+ 9]<< 2 & x80008000)>> 24; r2 |= word_slice_64_16_l(f[j+ 9]<< 1 & x80008000)>> 24;
        case  9: r0 |= word_slice_64_16_l(f[j+ 8]<< 3 & x80008000)>> 28; r1 |= word_slice_64_16_l(f[j+ 8]<< 2 & x80008000)>> 28; r2 |= word_slice_64_16_l(f[j+ 8]<< 1 & x80008000)>> 28;
        case  8: r0 |= word_slice_64_16_l(f[j+ 7]<< 3 & x80008000)>> 32; r1 |= word_slice_64_16_l(f[j+ 7]<< 2 & x80008000)>> 32; r2 |= word_slice_64_16_l(f[j+ 7]<< 1 & x80008000)>> 32;
        case  7: r0 |= word_slice_64_16_l(f[j+ 6]<< 3 & x80008000)>> 36; r1 |= word_slice_64_16_l(f[j+ 6]<< 2 & x80008000)>> 36; r2 |= word_slice_64_16_l(f[j+ 6]<< 1 & x80008000)>> 36;
        case  6: r0 |= word_slice_64_16_l(f[j+ 5]<< 3 & x80008000)>> 40; r1 |= word_slice_64_16_l(f[j+ 5]<< 2 & x80008000)>> 40; r2 |= word_slice_64_16_l(f[j+ 5]<< 1 & x80008000)>> 40;
        case  5: r0 |= word_slice_64_16_l(f[j+ 4]<< 3 & x80008000)>> 44; r1 |= word_slice_64_16_l(f[j+ 4]<< 2 & x80008000)>> 44; r2 |= word_slice_64_16_l(f[j+ 4]<< 1 & x80008000)>> 44;
        case  4: r0 |= word_slice_64_16_l(f[j+ 3]<< 3 & x80008000)>> 48; r1 |= word_slice_64_16_l(f[j+ 3]<< 2 & x80008000)>> 48; r2 |= word_slice_64_16_l(f[j+ 3]<< 1 & x80008000)>> 48;
        case  3: r0 |= word_slice_64_16_l(f[j+ 2]<< 3 & x80008000)>> 52; r1 |= word_slice_64_16_l(f[j+ 2]<< 2 & x80008000)>> 52; r2 |= word_slice_64_16_l(f[j+ 2]<< 1 & x80008000)>> 52;
        case  2: r0 |= word_slice_64_16_l(f[j+ 1]<< 3 & x80008000)>> 56; r1 |= word_slice_64_16_l(f[j+ 1]<< 2 & x80008000)>> 56; r2 |= word_slice_64_16_l(f[j+ 1]<< 1 & x80008000)>> 56;
        case  1: r0 |= word_slice_64_16_l(f[j+ 0]<< 3 & x80008000)>> 60; r1 |= word_slice_64_16_l(f[j+ 0]<< 2 & x80008000)>> 60; r2 |= word_slice_64_16_l(f[j+ 0]<< 1 & x80008000)>> 60;
          break;
        default:
          m4ri_die("impossible");
        }
        t0[j2] |= r0 & bitmask_end;
        t1[j2] |= r1 & bitmask_end;
        t2[j2] |= r2 & bitmask_end;
      }
    } break;
    case 14: {
      for(size_t i=0; i<T->nrows; i++) {
        word *t0 = T->x[12]->rows[i];
        word *t1 = T->x[13]->rows[i];
        const word const *f  = F->x->rows[i];

        /* bulk of work */
        for(j=0, j2=0; j+16 < F->x->width; j+=16,j2++) {
          word_slice_64_16_l_combine_bulk(t0, j2, f, j,  3);
          word_slice_64_16_l_combine_bulk(t1, j2, f, j,  2);
        }
        r0 = r1 = 0;
        switch(F->x->width - j) {
        case 16: r0 |= word_slice_64_16_l(f[j+15]<< 3 & x80008000)>>  0; r1 |= word_slice_64_16_l(f[j+15]<< 2 & x80008000)>>  0;
        case 15: r0 |= word_slice_64_16_l(f[j+14]<< 3 & x80008000)>>  4; r1 |= word_slice_64_16_l(f[j+14]<< 2 & x80008000)>>  4;
        case 14: r0 |= word_slice_64_16_l(f[j+13]<< 3 & x80008000)>>  8; r1 |= word_slice_64_16_l(f[j+13]<< 2 & x80008000)>>  8;
        case 13: r0 |= word_slice_64_16_l(f[j+12]<< 3 & x80008000)>> 12; r1 |= word_slice_64_16_l(f[j+12]<< 2 & x80008000)>> 12;
        case 12: r0 |= word_slice_64_16_l(f[j+11]<< 3 & x80008000)>> 16; r1 |= word_slice_64_16_l(f[j+11]<< 2 & x80008000)>> 16;
        case 11: r0 |= word_slice_64_16_l(f[j+10]<< 3 & x80008000)>> 20; r1 |= word_slice_64_16_l(f[j+10]<< 2 & x80008000)>> 20;
        case 10: r0 |= word_slice_64_16_l(f[j+ 9]<< 3 & x80008000)>> 24; r1 |= word_slice_64_16_l(f[j+ 9]<< 2 & x80008000)>> 24;
        case  9: r0 |= word_slice_64_16_l(f[j+ 8]<< 3 & x80008000)>> 28; r1 |= word_slice_64_16_l(f[j+ 8]<< 2 & x80008000)>> 28;
        case  8: r0 |= word_slice_64_16_l(f[j+ 7]<< 3 & x80008000)>> 32; r1 |= word_slice_64_16_l(f[j+ 7]<< 2 & x80008000)>> 32;
        case  7: r0 |= word_slice_64_16_l(f[j+ 6]<< 3 & x80008000)>> 36; r1 |= word_slice_64_16_l(f[j+ 6]<< 2 & x80008000)>> 36;
        case  6: r0 |= word_slice_64_16_l(f[j+ 5]<< 3 & x80008000)>> 40; r1 |= word_slice_64_16_l(f[j+ 5]<< 2 & x80008000)>> 40;
        case  5: r0 |= word_slice_64_16_l(f[j+ 4]<< 3 & x80008000)>> 44; r1 |= word_slice_64_16_l(f[j+ 4]<< 2 & x80008000)>> 44;
        case  4: r0 |= word_slice_64_16_l(f[j+ 3]<< 3 & x80008000)>> 48; r1 |= word_slice_64_16_l(f[j+ 3]<< 2 & x80008000)>> 48;
        case  3: r0 |= word_slice_64_16_l(f[j+ 2]<< 3 & x80008000)>> 52; r1 |= word_slice_64_16_l(f[j+ 2]<< 2 & x80008000)>> 52;
        case  2: r0 |= word_slice_64_16_l(f[j+ 1]<< 3 & x80008000)>> 56; r1 |= word_slice_64_16_l(f[j+ 1]<< 2 & x80008000)>> 56;
        case  1: r0 |= word_slice_64_16_l(f[j+ 0]<< 3 & x80008000)>> 60; r1 |= word_slice_64_16_l(f[j+ 0]<< 2 & x80008000)>> 60;
          break;
        default:
          m4ri_die("impossible");
        }
        t0[j2] |= r0 & bitmask_end;
        t1[j2] |= r1 & bitmask_end;
      }
    } break;
    case 13: {
      for(size_t i=0; i<T->nrows; i++) {
        word *t0 = T->x[12]->rows[i];
        const word const *f  = F->x->rows[i];

        /* bulk of work */
        for(j=0, j2=0; j+16 < F->x->width; j+=16,j2++) {
          word_slice_64_16_l_combine_bulk(t0, j2, f, j,  3);
        }
        r0 = 0;
        switch(F->x->width - j) {
        case 16: r0 |= word_slice_64_16_l(f[j+15]<< 3 & x80008000)>>  0;
        case 15: r0 |= word_slice_64_16_l(f[j+14]<< 3 & x80008000)>>  4;
        case 14: r0 |= word_slice_64_16_l(f[j+13]<< 3 & x80008000)>>  8;
        case 13: r0 |= word_slice_64_16_l(f[j+12]<< 3 & x80008000)>> 12;
        case 12: r0 |= word_slice_64_16_l(f[j+11]<< 3 & x80008000)>> 16;
        case 11: r0 |= word_slice_64_16_l(f[j+10]<< 3 & x80008000)>> 20;
        case 10: r0 |= word_slice_64_16_l(f[j+ 9]<< 3 & x80008000)>> 24;
        case  9: r0 |= word_slice_64_16_l(f[j+ 8]<< 3 & x80008000)>> 28;
        case  8: r0 |= word_slice_64_16_l(f[j+ 7]<< 3 & x80008000)>> 32;
        case  7: r0 |= word_slice_64_16_l(f[j+ 6]<< 3 & x80008000)>> 36;
        case  6: r0 |= word_slice_64_16_l(f[j+ 5]<< 3 & x80008000)>> 40;
        case  5: r0 |= word_slice_64_16_l(f[j+ 4]<< 3 & x80008000)>> 44;
        case  4: r0 |= word_slice_64_16_l(f[j+ 3]<< 3 & x80008000)>> 48;
        case  3: r0 |= word_slice_64_16_l(f[j+ 2]<< 3 & x80008000)>> 52;
        case  2: r0 |= word_slice_64_16_l(f[j+ 1]<< 3 & x80008000)>> 56;
        case  1: r0 |= word_slice_64_16_l(f[j+ 0]<< 3 & x80008000)>> 60;
          break;
        default:
          m4ri_die("impossible");
        }
        t0[j2] |= r0 & bitmask_end;
      }
    } break;
    }
  } else {
    switch(T->depth) {
    case 11: {
      for(size_t i=0; i<T->nrows; i++) {
        word *t0 = T->x[ 8]->rows[i];
        word *t1 = T->x[ 9]->rows[i];
        word *t2 = T->x[10]->rows[i];
        const word const *f  = F->x->rows[i];

        /* bulk of work */
        for(j=0, j2=0; j+16 < F->x->width; j+=16,j2++) {
          word_slice_64_16_l_combine_bulk(t0, j2, f, j,  7);
          word_slice_64_16_l_combine_bulk(t1, j2, f, j,  6);
          word_slice_64_16_l_combine_bulk(t2, j2, f, j,  5);
        }
        r0 = r1 = r2 = 0;
        switch(F->x->width - j) {
        case 16: r0 |= word_slice_64_16_l(f[j+15]<< 7 & x80008000)>>  0; r1 |= word_slice_64_16_l(f[j+15]<< 6 & x80008000)>>  0; r2 |= word_slice_64_16_l(f[j+15]<< 5 & x80008000)>>  0;
        case 15: r0 |= word_slice_64_16_l(f[j+14]<< 7 & x80008000)>>  4; r1 |= word_slice_64_16_l(f[j+14]<< 6 & x80008000)>>  4; r2 |= word_slice_64_16_l(f[j+14]<< 5 & x80008000)>>  4;
        case 14: r0 |= word_slice_64_16_l(f[j+13]<< 7 & x80008000)>>  8; r1 |= word_slice_64_16_l(f[j+13]<< 6 & x80008000)>>  8; r2 |= word_slice_64_16_l(f[j+13]<< 5 & x80008000)>>  8;
        case 13: r0 |= word_slice_64_16_l(f[j+12]<< 7 & x80008000)>> 12; r1 |= word_slice_64_16_l(f[j+12]<< 6 & x80008000)>> 12; r2 |= word_slice_64_16_l(f[j+12]<< 5 & x80008000)>> 12;
        case 12: r0 |= word_slice_64_16_l(f[j+11]<< 7 & x80008000)>> 16; r1 |= word_slice_64_16_l(f[j+11]<< 6 & x80008000)>> 16; r2 |= word_slice_64_16_l(f[j+11]<< 5 & x80008000)>> 16;
        case 11: r0 |= word_slice_64_16_l(f[j+10]<< 7 & x80008000)>> 20; r1 |= word_slice_64_16_l(f[j+10]<< 6 & x80008000)>> 20; r2 |= word_slice_64_16_l(f[j+10]<< 5 & x80008000)>> 20;
        case 10: r0 |= word_slice_64_16_l(f[j+ 9]<< 7 & x80008000)>> 24; r1 |= word_slice_64_16_l(f[j+ 9]<< 6 & x80008000)>> 24; r2 |= word_slice_64_16_l(f[j+ 9]<< 5 & x80008000)>> 24;
        case  9: r0 |= word_slice_64_16_l(f[j+ 8]<< 7 & x80008000)>> 28; r1 |= word_slice_64_16_l(f[j+ 8]<< 6 & x80008000)>> 28; r2 |= word_slice_64_16_l(f[j+ 8]<< 5 & x80008000)>> 28;
        case  8: r0 |= word_slice_64_16_l(f[j+ 7]<< 7 & x80008000)>> 32; r1 |= word_slice_64_16_l(f[j+ 7]<< 6 & x80008000)>> 32; r2 |= word_slice_64_16_l(f[j+ 7]<< 5 & x80008000)>> 32;
        case  7: r0 |= word_slice_64_16_l(f[j+ 6]<< 7 & x80008000)>> 36; r1 |= word_slice_64_16_l(f[j+ 6]<< 6 & x80008000)>> 36; r2 |= word_slice_64_16_l(f[j+ 6]<< 5 & x80008000)>> 36;
        case  6: r0 |= word_slice_64_16_l(f[j+ 5]<< 7 & x80008000)>> 40; r1 |= word_slice_64_16_l(f[j+ 5]<< 6 & x80008000)>> 40; r2 |= word_slice_64_16_l(f[j+ 5]<< 5 & x80008000)>> 40;
        case  5: r0 |= word_slice_64_16_l(f[j+ 4]<< 7 & x80008000)>> 44; r1 |= word_slice_64_16_l(f[j+ 4]<< 6 & x80008000)>> 44; r2 |= word_slice_64_16_l(f[j+ 4]<< 5 & x80008000)>> 44;
        case  4: r0 |= word_slice_64_16_l(f[j+ 3]<< 7 & x80008000)>> 48; r1 |= word_slice_64_16_l(f[j+ 3]<< 6 & x80008000)>> 48; r2 |= word_slice_64_16_l(f[j+ 3]<< 5 & x80008000)>> 48;
        case  3: r0 |= word_slice_64_16_l(f[j+ 2]<< 7 & x80008000)>> 52; r1 |= word_slice_64_16_l(f[j+ 2]<< 6 & x80008000)>> 52; r2 |= word_slice_64_16_l(f[j+ 2]<< 5 & x80008000)>> 52;
        case  2: r0 |= word_slice_64_16_l(f[j+ 1]<< 7 & x80008000)>> 56; r1 |= word_slice_64_16_l(f[j+ 1]<< 6 & x80008000)>> 56; r2 |= word_slice_64_16_l(f[j+ 1]<< 5 & x80008000)>> 56;
        case  1: r0 |= word_slice_64_16_l(f[j+ 0]<< 7 & x80008000)>> 60; r1 |= word_slice_64_16_l(f[j+ 0]<< 6 & x80008000)>> 60; r2 |= word_slice_64_16_l(f[j+ 0]<< 5 & x80008000)>> 60;
          break;
        default:
          m4ri_die("impossible");
        }
        t0[j2] |= r0 & bitmask_end;
        t1[j2] |= r1 & bitmask_end;
        t2[j2] |= r2 & bitmask_end;
      }
    } break;
    case 10: {
      for(size_t i=0; i<T->nrows; i++) {
        word *t0 = T->x[ 8]->rows[i];
        word *t1 = T->x[ 9]->rows[i];
        const word const *f  = F->x->rows[i];

        /* bulk of work */
        for(j=0, j2=0; j+16 < F->x->width; j+=16,j2++) {
          word_slice_64_16_l_combine_bulk(t0, j2, f, j,  7);
          word_slice_64_16_l_combine_bulk(t1, j2, f, j,  6);
        }
        r0 = r1 = 0;
        switch(F->x->width - j) {
        case 16: r0 |= word_slice_64_16_l(f[j+15]<< 7 & x80008000)>>  0; r1 |= word_slice_64_16_l(f[j+15]<< 6 & x80008000)>>  0;
        case 15: r0 |= word_slice_64_16_l(f[j+14]<< 7 & x80008000)>>  4; r1 |= word_slice_64_16_l(f[j+14]<< 6 & x80008000)>>  4;
        case 14: r0 |= word_slice_64_16_l(f[j+13]<< 7 & x80008000)>>  8; r1 |= word_slice_64_16_l(f[j+13]<< 6 & x80008000)>>  8;
        case 13: r0 |= word_slice_64_16_l(f[j+12]<< 7 & x80008000)>> 12; r1 |= word_slice_64_16_l(f[j+12]<< 6 & x80008000)>> 12;
        case 12: r0 |= word_slice_64_16_l(f[j+11]<< 7 & x80008000)>> 16; r1 |= word_slice_64_16_l(f[j+11]<< 6 & x80008000)>> 16;
        case 11: r0 |= word_slice_64_16_l(f[j+10]<< 7 & x80008000)>> 20; r1 |= word_slice_64_16_l(f[j+10]<< 6 & x80008000)>> 20;
        case 10: r0 |= word_slice_64_16_l(f[j+ 9]<< 7 & x80008000)>> 24; r1 |= word_slice_64_16_l(f[j+ 9]<< 6 & x80008000)>> 24;
        case  9: r0 |= word_slice_64_16_l(f[j+ 8]<< 7 & x80008000)>> 28; r1 |= word_slice_64_16_l(f[j+ 8]<< 6 & x80008000)>> 28;
        case  8: r0 |= word_slice_64_16_l(f[j+ 7]<< 7 & x80008000)>> 32; r1 |= word_slice_64_16_l(f[j+ 7]<< 6 & x80008000)>> 32;
        case  7: r0 |= word_slice_64_16_l(f[j+ 6]<< 7 & x80008000)>> 36; r1 |= word_slice_64_16_l(f[j+ 6]<< 6 & x80008000)>> 36;
        case  6: r0 |= word_slice_64_16_l(f[j+ 5]<< 7 & x80008000)>> 40; r1 |= word_slice_64_16_l(f[j+ 5]<< 6 & x80008000)>> 40;
        case  5: r0 |= word_slice_64_16_l(f[j+ 4]<< 7 & x80008000)>> 44; r1 |= word_slice_64_16_l(f[j+ 4]<< 6 & x80008000)>> 44;
        case  4: r0 |= word_slice_64_16_l(f[j+ 3]<< 7 & x80008000)>> 48; r1 |= word_slice_64_16_l(f[j+ 3]<< 6 & x80008000)>> 48;
        case  3: r0 |= word_slice_64_16_l(f[j+ 2]<< 7 & x80008000)>> 52; r1 |= word_slice_64_16_l(f[j+ 2]<< 6 & x80008000)>> 52;
        case  2: r0 |= word_slice_64_16_l(f[j+ 1]<< 7 & x80008000)>> 56; r1 |= word_slice_64_16_l(f[j+ 1]<< 6 & x80008000)>> 56;
        case  1: r0 |= word_slice_64_16_l(f[j+ 0]<< 7 & x80008000)>> 60; r1 |= word_slice_64_16_l(f[j+ 0]<< 6 & x80008000)>> 60;
          break;
        default:
          m4ri_die("impossible");
        }
        t0[j2] |= r0 & bitmask_end;
        t1[j2] |= r1 & bitmask_end;
      }
    } break;
    case  9: {
      for(size_t i=0; i<T->nrows; i++) {
        word *t0 = T->x[ 8]->rows[i];
        const word const *f  = F->x->rows[i];

        /* bulk of work */
        for(j=0, j2=0; j+16 < F->x->width; j+=16,j2++) {
          word_slice_64_16_l_combine_bulk(t0, j2, f, j,  7);
        }
        r0 = 0;
        switch(F->x->width - j) {
        case 16: r0 |= word_slice_64_16_l(f[j+15]<< 7 & x80008000)>>  0;
        case 15: r0 |= word_slice_64_16_l(f[j+14]<< 7 & x80008000)>>  4;
        case 14: r0 |= word_slice_64_16_l(f[j+13]<< 7 & x80008000)>>  8;
        case 13: r0 |= word_slice_64_16_l(f[j+12]<< 7 & x80008000)>> 12;
        case 12: r0 |= word_slice_64_16_l(f[j+11]<< 7 & x80008000)>> 16;
        case 11: r0 |= word_slice_64_16_l(f[j+10]<< 7 & x80008000)>> 20;
        case 10: r0 |= word_slice_64_16_l(f[j+ 9]<< 7 & x80008000)>> 24;
        case  9: r0 |= word_slice_64_16_l(f[j+ 8]<< 7 & x80008000)>> 28;
        case  8: r0 |= word_slice_64_16_l(f[j+ 7]<< 7 & x80008000)>> 32;
        case  7: r0 |= word_slice_64_16_l(f[j+ 6]<< 7 & x80008000)>> 36;
        case  6: r0 |= word_slice_64_16_l(f[j+ 5]<< 7 & x80008000)>> 40;
        case  5: r0 |= word_slice_64_16_l(f[j+ 4]<< 7 & x80008000)>> 44;
        case  4: r0 |= word_slice_64_16_l(f[j+ 3]<< 7 & x80008000)>> 48;
        case  3: r0 |= word_slice_64_16_l(f[j+ 2]<< 7 & x80008000)>> 52;
        case  2: r0 |= word_slice_64_16_l(f[j+ 1]<< 7 & x80008000)>> 56;
        case  1: r0 |= word_slice_64_16_l(f[j+ 0]<< 7 & x80008000)>> 60;
          break;
        default:
          m4ri_die("impossible");
        }
        t0[j2] |= r0 & bitmask_end;
      }
    } break;
    default:
      m4ri_die("impossible");
    }
  }
  return T;
}

mzed_t *_mzed_cling16(mzed_t *T, const mzd_slice_t *F) {
  wi_t j,j2 = 0;

  const word bitmask_end = __M4RI_LEFT_BITMASK((T->x->offset + T->x->ncols) % m4ri_radix);

  if (mzd_slice_is_zero(F))
    return T;

  for(rci_t i=0; i<T->nrows; i++) {
    const word *f00 = F->x[ 0]->rows[i];
    const word *f01 = F->x[ 1]->rows[i];
    const word *f02 = F->x[ 2]->rows[i];
    const word *f03 = F->x[ 3]->rows[i];
    const word *f04 = F->x[ 4]->rows[i];
    const word *f05 = F->x[ 5]->rows[i];
    const word *f06 = F->x[ 6]->rows[i];
    const word *f07 = F->x[ 7]->rows[i];
    word *t  = T->x->rows[i];

    for(j=0, j2=0; j+16 < T->x->width; j+=16, j2++) {
      t[j+ 0] = (word_cling_64_16_l(f00[j2]<<60)>>15) | (word_cling_64_16_l(f01[j2]<<60)>>14) | (word_cling_64_16_l(f02[j2]<<60)>>13) | (word_cling_64_16_l(f03[j2]<<60)>>12) \
        |       (word_cling_64_16_l(f04[j2]<<60)>>11) | (word_cling_64_16_l(f05[j2]<<60)>>10) | (word_cling_64_16_l(f06[j2]<<60)>> 9) | (word_cling_64_16_l(f07[j2]<<60)>> 8);
      t[j+ 1] = (word_cling_64_16_l(f00[j2]<<56)>>15) | (word_cling_64_16_l(f01[j2]<<56)>>14) | (word_cling_64_16_l(f02[j2]<<56)>>13) | (word_cling_64_16_l(f03[j2]<<56)>>12) \
        |       (word_cling_64_16_l(f04[j2]<<56)>>11) | (word_cling_64_16_l(f05[j2]<<56)>>10) | (word_cling_64_16_l(f06[j2]<<56)>> 9) | (word_cling_64_16_l(f07[j2]<<56)>> 8);
      t[j+ 2] = (word_cling_64_16_l(f00[j2]<<52)>>15) | (word_cling_64_16_l(f01[j2]<<52)>>14) | (word_cling_64_16_l(f02[j2]<<52)>>13) | (word_cling_64_16_l(f03[j2]<<52)>>12) \
        |       (word_cling_64_16_l(f04[j2]<<52)>>11) | (word_cling_64_16_l(f05[j2]<<52)>>10) | (word_cling_64_16_l(f06[j2]<<52)>> 9) | (word_cling_64_16_l(f07[j2]<<52)>> 8);
      t[j+ 3] = (word_cling_64_16_l(f00[j2]<<48)>>15) | (word_cling_64_16_l(f01[j2]<<48)>>14) | (word_cling_64_16_l(f02[j2]<<48)>>13) | (word_cling_64_16_l(f03[j2]<<48)>>12) \
        |       (word_cling_64_16_l(f04[j2]<<48)>>11) | (word_cling_64_16_l(f05[j2]<<48)>>10) | (word_cling_64_16_l(f06[j2]<<48)>> 9) | (word_cling_64_16_l(f07[j2]<<48)>> 8);
      t[j+ 4] = (word_cling_64_16_l(f00[j2]<<44)>>15) | (word_cling_64_16_l(f01[j2]<<44)>>14) | (word_cling_64_16_l(f02[j2]<<44)>>13) | (word_cling_64_16_l(f03[j2]<<44)>>12) \
        |       (word_cling_64_16_l(f04[j2]<<44)>>11) | (word_cling_64_16_l(f05[j2]<<44)>>10) | (word_cling_64_16_l(f06[j2]<<44)>> 9) | (word_cling_64_16_l(f07[j2]<<44)>> 8);
      t[j+ 5] = (word_cling_64_16_l(f00[j2]<<40)>>15) | (word_cling_64_16_l(f01[j2]<<40)>>14) | (word_cling_64_16_l(f02[j2]<<40)>>13) | (word_cling_64_16_l(f03[j2]<<40)>>12) \
        |       (word_cling_64_16_l(f04[j2]<<40)>>11) | (word_cling_64_16_l(f05[j2]<<40)>>10) | (word_cling_64_16_l(f06[j2]<<40)>> 9) | (word_cling_64_16_l(f07[j2]<<40)>> 8);
      t[j+ 6] = (word_cling_64_16_l(f00[j2]<<36)>>15) | (word_cling_64_16_l(f01[j2]<<36)>>14) | (word_cling_64_16_l(f02[j2]<<36)>>13) | (word_cling_64_16_l(f03[j2]<<36)>>12) \
        |       (word_cling_64_16_l(f04[j2]<<36)>>11) | (word_cling_64_16_l(f05[j2]<<36)>>10) | (word_cling_64_16_l(f06[j2]<<36)>> 9) | (word_cling_64_16_l(f07[j2]<<36)>> 8);
      t[j+ 7] = (word_cling_64_16_l(f00[j2]<<32)>>15) | (word_cling_64_16_l(f01[j2]<<32)>>14) | (word_cling_64_16_l(f02[j2]<<32)>>13) | (word_cling_64_16_l(f03[j2]<<32)>>12) \
        |       (word_cling_64_16_l(f04[j2]<<32)>>11) | (word_cling_64_16_l(f05[j2]<<32)>>10) | (word_cling_64_16_l(f06[j2]<<32)>> 9) | (word_cling_64_16_l(f07[j2]<<32)>> 8);
      t[j+ 8] = (word_cling_64_16_l(f00[j2]<<28)>>15) | (word_cling_64_16_l(f01[j2]<<28)>>14) | (word_cling_64_16_l(f02[j2]<<28)>>13) | (word_cling_64_16_l(f03[j2]<<28)>>12) \
        |       (word_cling_64_16_l(f04[j2]<<28)>>11) | (word_cling_64_16_l(f05[j2]<<28)>>10) | (word_cling_64_16_l(f06[j2]<<28)>> 9) | (word_cling_64_16_l(f07[j2]<<28)>> 8);
      t[j+ 9] = (word_cling_64_16_l(f00[j2]<<24)>>15) | (word_cling_64_16_l(f01[j2]<<24)>>14) | (word_cling_64_16_l(f02[j2]<<24)>>13) | (word_cling_64_16_l(f03[j2]<<24)>>12) \
        |       (word_cling_64_16_l(f04[j2]<<24)>>11) | (word_cling_64_16_l(f05[j2]<<24)>>10) | (word_cling_64_16_l(f06[j2]<<24)>> 9) | (word_cling_64_16_l(f07[j2]<<24)>> 8);
      t[j+10] = (word_cling_64_16_l(f00[j2]<<20)>>15) | (word_cling_64_16_l(f01[j2]<<20)>>14) | (word_cling_64_16_l(f02[j2]<<20)>>13) | (word_cling_64_16_l(f03[j2]<<20)>>12) \
        |       (word_cling_64_16_l(f04[j2]<<20)>>11) | (word_cling_64_16_l(f05[j2]<<20)>>10) | (word_cling_64_16_l(f06[j2]<<20)>> 9) | (word_cling_64_16_l(f07[j2]<<20)>> 8);
      t[j+11] = (word_cling_64_16_l(f00[j2]<<16)>>15) | (word_cling_64_16_l(f01[j2]<<16)>>14) | (word_cling_64_16_l(f02[j2]<<16)>>13) | (word_cling_64_16_l(f03[j2]<<16)>>12) \
        |       (word_cling_64_16_l(f04[j2]<<16)>>11) | (word_cling_64_16_l(f05[j2]<<16)>>10) | (word_cling_64_16_l(f06[j2]<<16)>> 9) | (word_cling_64_16_l(f07[j2]<<16)>> 8);
      t[j+12] = (word_cling_64_16_l(f00[j2]<<12)>>15) | (word_cling_64_16_l(f01[j2]<<12)>>14) | (word_cling_64_16_l(f02[j2]<<12)>>13) | (word_cling_64_16_l(f03[j2]<<12)>>12) \
        |       (word_cling_64_16_l(f04[j2]<<12)>>11) | (word_cling_64_16_l(f05[j2]<<12)>>10) | (word_cling_64_16_l(f06[j2]<<12)>> 9) | (word_cling_64_16_l(f07[j2]<<12)>> 8);
      t[j+13] = (word_cling_64_16_l(f00[j2]<< 8)>>15) | (word_cling_64_16_l(f01[j2]<< 8)>>14) | (word_cling_64_16_l(f02[j2]<< 8)>>13) | (word_cling_64_16_l(f03[j2]<< 8)>>12) \
        |       (word_cling_64_16_l(f04[j2]<< 8)>>11) | (word_cling_64_16_l(f05[j2]<< 8)>>10) | (word_cling_64_16_l(f06[j2]<< 8)>> 9) | (word_cling_64_16_l(f07[j2]<< 8)>> 8);
      t[j+14] = (word_cling_64_16_l(f00[j2]<< 4)>>15) | (word_cling_64_16_l(f01[j2]<< 4)>>14) | (word_cling_64_16_l(f02[j2]<< 4)>>13) | (word_cling_64_16_l(f03[j2]<< 4)>>12) \
        |       (word_cling_64_16_l(f04[j2]<< 4)>>11) | (word_cling_64_16_l(f05[j2]<< 4)>>10) | (word_cling_64_16_l(f06[j2]<< 4)>> 9) | (word_cling_64_16_l(f07[j2]<< 4)>> 8);
      t[j+15] = (word_cling_64_16_l(f00[j2]<< 0)>>15) | (word_cling_64_16_l(f01[j2]<< 0)>>14) | (word_cling_64_16_l(f02[j2]<< 0)>>13) | (word_cling_64_16_l(f03[j2]<< 0)>>12) \
        |       (word_cling_64_16_l(f04[j2]<< 0)>>11) | (word_cling_64_16_l(f05[j2]<< 0)>>10) | (word_cling_64_16_l(f06[j2]<< 0)>> 9) | (word_cling_64_16_l(f07[j2]<< 0)>> 8);
    }

    register word tmp = t[T->x->width-1];
    switch(T->x->width - j) {
    case 16: t[j+15] = (word_cling_64_16_l(f00[j2]<< 0)>>15) | (word_cling_64_16_l(f01[j2]<< 0)>>14) | (word_cling_64_16_l(f02[j2]<< 0)>>13) | (word_cling_64_16_l(f03[j2]<< 0)>>12) | \
                       (word_cling_64_16_l(f04[j2]<< 0)>>11) | (word_cling_64_16_l(f05[j2]<< 0)>>10) | (word_cling_64_16_l(f06[j2]<< 0)>> 9) | (word_cling_64_16_l(f07[j2]<< 0)>> 8);
    case 15: t[j+14] = (word_cling_64_16_l(f00[j2]<< 4)>>15) | (word_cling_64_16_l(f01[j2]<< 4)>>14) | (word_cling_64_16_l(f02[j2]<< 4)>>13) | (word_cling_64_16_l(f03[j2]<< 4)>>12) | \
                       (word_cling_64_16_l(f04[j2]<< 4)>>11) | (word_cling_64_16_l(f05[j2]<< 4)>>10) | (word_cling_64_16_l(f06[j2]<< 4)>> 9) | (word_cling_64_16_l(f07[j2]<< 4)>> 8);
    case 14: t[j+13] = (word_cling_64_16_l(f00[j2]<< 8)>>15) | (word_cling_64_16_l(f01[j2]<< 8)>>14) | (word_cling_64_16_l(f02[j2]<< 8)>>13) | (word_cling_64_16_l(f03[j2]<< 8)>>12) | \
                       (word_cling_64_16_l(f04[j2]<< 8)>>11) | (word_cling_64_16_l(f05[j2]<< 8)>>10) | (word_cling_64_16_l(f06[j2]<< 8)>> 9) | (word_cling_64_16_l(f07[j2]<< 8)>> 8);
    case 13: t[j+12] = (word_cling_64_16_l(f00[j2]<<12)>>15) | (word_cling_64_16_l(f01[j2]<<12)>>14) | (word_cling_64_16_l(f02[j2]<<12)>>13) | (word_cling_64_16_l(f03[j2]<<12)>>12) | \
                       (word_cling_64_16_l(f04[j2]<<12)>>11) | (word_cling_64_16_l(f05[j2]<<12)>>10) | (word_cling_64_16_l(f06[j2]<<12)>> 9) | (word_cling_64_16_l(f07[j2]<<12)>> 8);
    case 12: t[j+11] = (word_cling_64_16_l(f00[j2]<<16)>>15) | (word_cling_64_16_l(f01[j2]<<16)>>14) | (word_cling_64_16_l(f02[j2]<<16)>>13) | (word_cling_64_16_l(f03[j2]<<16)>>12) | \
                       (word_cling_64_16_l(f04[j2]<<16)>>11) | (word_cling_64_16_l(f05[j2]<<16)>>10) | (word_cling_64_16_l(f06[j2]<<16)>> 9) | (word_cling_64_16_l(f07[j2]<<16)>> 8);
    case 11: t[j+10] = (word_cling_64_16_l(f00[j2]<<20)>>15) | (word_cling_64_16_l(f01[j2]<<20)>>14) | (word_cling_64_16_l(f02[j2]<<20)>>13) | (word_cling_64_16_l(f03[j2]<<20)>>12) | \
                       (word_cling_64_16_l(f04[j2]<<20)>>11) | (word_cling_64_16_l(f05[j2]<<20)>>10) | (word_cling_64_16_l(f06[j2]<<20)>> 9) | (word_cling_64_16_l(f07[j2]<<20)>> 8);
    case 10: t[j+ 9] = (word_cling_64_16_l(f00[j2]<<24)>>15) | (word_cling_64_16_l(f01[j2]<<24)>>14) | (word_cling_64_16_l(f02[j2]<<24)>>13) | (word_cling_64_16_l(f03[j2]<<24)>>12) | \
                       (word_cling_64_16_l(f04[j2]<<24)>>11) | (word_cling_64_16_l(f05[j2]<<24)>>10) | (word_cling_64_16_l(f06[j2]<<24)>> 9) | (word_cling_64_16_l(f07[j2]<<24)>> 8);
    case  9: t[j+ 8] = (word_cling_64_16_l(f00[j2]<<28)>>15) | (word_cling_64_16_l(f01[j2]<<28)>>14) | (word_cling_64_16_l(f02[j2]<<28)>>13) | (word_cling_64_16_l(f03[j2]<<28)>>12) | \
                       (word_cling_64_16_l(f04[j2]<<28)>>11) | (word_cling_64_16_l(f05[j2]<<28)>>10) | (word_cling_64_16_l(f06[j2]<<28)>> 9) | (word_cling_64_16_l(f07[j2]<<28)>> 8);
    case  8: t[j+ 7] = (word_cling_64_16_l(f00[j2]<<32)>>15) | (word_cling_64_16_l(f01[j2]<<32)>>14) | (word_cling_64_16_l(f02[j2]<<32)>>13) | (word_cling_64_16_l(f03[j2]<<32)>>12) | \
                       (word_cling_64_16_l(f04[j2]<<32)>>11) | (word_cling_64_16_l(f05[j2]<<32)>>10) | (word_cling_64_16_l(f06[j2]<<32)>> 9) | (word_cling_64_16_l(f07[j2]<<32)>> 8);
    case  7: t[j+ 6] = (word_cling_64_16_l(f00[j2]<<36)>>15) | (word_cling_64_16_l(f01[j2]<<36)>>14) | (word_cling_64_16_l(f02[j2]<<36)>>13) | (word_cling_64_16_l(f03[j2]<<36)>>12) | \
                       (word_cling_64_16_l(f04[j2]<<36)>>11) | (word_cling_64_16_l(f05[j2]<<36)>>10) | (word_cling_64_16_l(f06[j2]<<36)>> 9) | (word_cling_64_16_l(f07[j2]<<36)>> 8);
    case  6: t[j+ 5] = (word_cling_64_16_l(f00[j2]<<40)>>15) | (word_cling_64_16_l(f01[j2]<<40)>>14) | (word_cling_64_16_l(f02[j2]<<40)>>13) | (word_cling_64_16_l(f03[j2]<<40)>>12) | \
                       (word_cling_64_16_l(f04[j2]<<40)>>11) | (word_cling_64_16_l(f05[j2]<<40)>>10) | (word_cling_64_16_l(f06[j2]<<40)>> 9) | (word_cling_64_16_l(f07[j2]<<40)>> 8);
    case  5: t[j+ 4] = (word_cling_64_16_l(f00[j2]<<44)>>15) | (word_cling_64_16_l(f01[j2]<<44)>>14) | (word_cling_64_16_l(f02[j2]<<44)>>13) | (word_cling_64_16_l(f03[j2]<<44)>>12) | \
                       (word_cling_64_16_l(f04[j2]<<44)>>11) | (word_cling_64_16_l(f05[j2]<<44)>>10) | (word_cling_64_16_l(f06[j2]<<44)>> 9) | (word_cling_64_16_l(f07[j2]<<44)>> 8);
    case  4: t[j+ 3] = (word_cling_64_16_l(f00[j2]<<48)>>15) | (word_cling_64_16_l(f01[j2]<<48)>>14) | (word_cling_64_16_l(f02[j2]<<48)>>13) | (word_cling_64_16_l(f03[j2]<<48)>>12) | \
                       (word_cling_64_16_l(f04[j2]<<48)>>11) | (word_cling_64_16_l(f05[j2]<<48)>>10) | (word_cling_64_16_l(f06[j2]<<48)>> 9) | (word_cling_64_16_l(f07[j2]<<48)>> 8);
    case  3: t[j+ 2] = (word_cling_64_16_l(f00[j2]<<52)>>15) | (word_cling_64_16_l(f01[j2]<<52)>>14) | (word_cling_64_16_l(f02[j2]<<52)>>13) | (word_cling_64_16_l(f03[j2]<<52)>>12) | \
                       (word_cling_64_16_l(f04[j2]<<52)>>11) | (word_cling_64_16_l(f05[j2]<<52)>>10) | (word_cling_64_16_l(f06[j2]<<52)>> 9) | (word_cling_64_16_l(f07[j2]<<52)>> 8);
    case  2: t[j+ 1] = (word_cling_64_16_l(f00[j2]<<56)>>15) | (word_cling_64_16_l(f01[j2]<<56)>>14) | (word_cling_64_16_l(f02[j2]<<56)>>13) | (word_cling_64_16_l(f03[j2]<<56)>>12) | \
                       (word_cling_64_16_l(f04[j2]<<56)>>11) | (word_cling_64_16_l(f05[j2]<<56)>>10) | (word_cling_64_16_l(f06[j2]<<56)>> 9) | (word_cling_64_16_l(f07[j2]<<56)>> 8);
    case  1: t[j+ 0] = (word_cling_64_16_l(f00[j2]<<60)>>15) | (word_cling_64_16_l(f01[j2]<<60)>>14) | (word_cling_64_16_l(f02[j2]<<60)>>13) | (word_cling_64_16_l(f03[j2]<<60)>>12) | \
                       (word_cling_64_16_l(f04[j2]<<60)>>11) | (word_cling_64_16_l(f05[j2]<<60)>>10) | (word_cling_64_16_l(f06[j2]<<60)>> 9) | (word_cling_64_16_l(f07[j2]<<60)>> 8);
      break;
    default:
      m4ri_die("impossible");
    }
    t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
  }

  if(T->finite_field->degree < 12) {
    switch(T->finite_field->degree) {
    case 9: {
      for(rci_t i=0; i<T->nrows; i++) {
        const word *f00 = F->x[ 8]->rows[i];
        word *t  = T->x->rows[i];

        for(j=0, j2=0; j+16 < T->x->width; j+=16, j2++) {
          t[j+ 0] |= (word_cling_64_16_l(f00[j2]<<60)>>7);
          t[j+ 1] |= (word_cling_64_16_l(f00[j2]<<56)>>7);
          t[j+ 2] |= (word_cling_64_16_l(f00[j2]<<52)>>7);
          t[j+ 3] |= (word_cling_64_16_l(f00[j2]<<48)>>7);
          t[j+ 4] |= (word_cling_64_16_l(f00[j2]<<44)>>7);
          t[j+ 5] |= (word_cling_64_16_l(f00[j2]<<40)>>7);
          t[j+ 6] |= (word_cling_64_16_l(f00[j2]<<36)>>7);
          t[j+ 7] |= (word_cling_64_16_l(f00[j2]<<32)>>7);
          t[j+ 8] |= (word_cling_64_16_l(f00[j2]<<28)>>7);
          t[j+ 9] |= (word_cling_64_16_l(f00[j2]<<24)>>7);
          t[j+10] |= (word_cling_64_16_l(f00[j2]<<20)>>7);
          t[j+11] |= (word_cling_64_16_l(f00[j2]<<16)>>7);
          t[j+12] |= (word_cling_64_16_l(f00[j2]<<12)>>7);
          t[j+13] |= (word_cling_64_16_l(f00[j2]<< 8)>>7);
          t[j+14] |= (word_cling_64_16_l(f00[j2]<< 4)>>7);
          t[j+15] |= (word_cling_64_16_l(f00[j2]<< 0)>>7);
        }

        register word tmp = t[T->x->width-1];
        switch(T->x->width - j) {
        case 16: t[j+15] |= (word_cling_64_16_l(f00[j2]<< 0)>>7);
        case 15: t[j+14] |= (word_cling_64_16_l(f00[j2]<< 4)>>7);
        case 14: t[j+13] |= (word_cling_64_16_l(f00[j2]<< 8)>>7);
        case 13: t[j+12] |= (word_cling_64_16_l(f00[j2]<<12)>>7);
        case 12: t[j+11] |= (word_cling_64_16_l(f00[j2]<<16)>>7);
        case 11: t[j+10] |= (word_cling_64_16_l(f00[j2]<<20)>>7);
        case 10: t[j+ 9] |= (word_cling_64_16_l(f00[j2]<<24)>>7);
        case  9: t[j+ 8] |= (word_cling_64_16_l(f00[j2]<<28)>>7);
        case  8: t[j+ 7] |= (word_cling_64_16_l(f00[j2]<<32)>>7);
        case  7: t[j+ 6] |= (word_cling_64_16_l(f00[j2]<<36)>>7);
        case  6: t[j+ 5] |= (word_cling_64_16_l(f00[j2]<<40)>>7);
        case  5: t[j+ 4] |= (word_cling_64_16_l(f00[j2]<<44)>>7);
        case  4: t[j+ 3] |= (word_cling_64_16_l(f00[j2]<<48)>>7);
        case  3: t[j+ 2] |= (word_cling_64_16_l(f00[j2]<<52)>>7);
        case  2: t[j+ 1] |= (word_cling_64_16_l(f00[j2]<<56)>>7);
        case  1: t[j+ 0] |= (word_cling_64_16_l(f00[j2]<<60)>>7);
          break;
        default:
          m4ri_die("impossible");
        }
        t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
      }
    }
      break;
    case 10: {
      for(rci_t i=0; i<T->nrows; i++) {
        const word *f00 = F->x[ 8]->rows[i];
        const word *f01 = F->x[ 9]->rows[i];
        word *t  = T->x->rows[i];

        for(j=0, j2=0; j+16 < T->x->width; j+=16, j2++) {
          t[j+ 0] |= (word_cling_64_16_l(f00[j2]<<60)>>7) | (word_cling_64_16_l(f01[j2]<<60)>>6);
          t[j+ 1] |= (word_cling_64_16_l(f00[j2]<<56)>>7) | (word_cling_64_16_l(f01[j2]<<56)>>6);
          t[j+ 2] |= (word_cling_64_16_l(f00[j2]<<52)>>7) | (word_cling_64_16_l(f01[j2]<<52)>>6);
          t[j+ 3] |= (word_cling_64_16_l(f00[j2]<<48)>>7) | (word_cling_64_16_l(f01[j2]<<48)>>6);
          t[j+ 4] |= (word_cling_64_16_l(f00[j2]<<44)>>7) | (word_cling_64_16_l(f01[j2]<<44)>>6);
          t[j+ 5] |= (word_cling_64_16_l(f00[j2]<<40)>>7) | (word_cling_64_16_l(f01[j2]<<40)>>6);
          t[j+ 6] |= (word_cling_64_16_l(f00[j2]<<36)>>7) | (word_cling_64_16_l(f01[j2]<<36)>>6);
          t[j+ 7] |= (word_cling_64_16_l(f00[j2]<<32)>>7) | (word_cling_64_16_l(f01[j2]<<32)>>6);
          t[j+ 8] |= (word_cling_64_16_l(f00[j2]<<28)>>7) | (word_cling_64_16_l(f01[j2]<<28)>>6);
          t[j+ 9] |= (word_cling_64_16_l(f00[j2]<<24)>>7) | (word_cling_64_16_l(f01[j2]<<24)>>6);
          t[j+10] |= (word_cling_64_16_l(f00[j2]<<20)>>7) | (word_cling_64_16_l(f01[j2]<<20)>>6);
          t[j+11] |= (word_cling_64_16_l(f00[j2]<<16)>>7) | (word_cling_64_16_l(f01[j2]<<16)>>6);
          t[j+12] |= (word_cling_64_16_l(f00[j2]<<12)>>7) | (word_cling_64_16_l(f01[j2]<<12)>>6);
          t[j+13] |= (word_cling_64_16_l(f00[j2]<< 8)>>7) | (word_cling_64_16_l(f01[j2]<< 8)>>6);
          t[j+14] |= (word_cling_64_16_l(f00[j2]<< 4)>>7) | (word_cling_64_16_l(f01[j2]<< 4)>>6);
          t[j+15] |= (word_cling_64_16_l(f00[j2]<< 0)>>7) | (word_cling_64_16_l(f01[j2]<< 0)>>6);
        }

        register word tmp = t[T->x->width-1];
        switch(T->x->width - j) {
        case 16: t[j+15] |= (word_cling_64_16_l(f00[j2]<< 0)>>7) | (word_cling_64_16_l(f01[j2]<< 0)>>6);
        case 15: t[j+14] |= (word_cling_64_16_l(f00[j2]<< 4)>>7) | (word_cling_64_16_l(f01[j2]<< 4)>>6);
        case 14: t[j+13] |= (word_cling_64_16_l(f00[j2]<< 8)>>7) | (word_cling_64_16_l(f01[j2]<< 8)>>6);
        case 13: t[j+12] |= (word_cling_64_16_l(f00[j2]<<12)>>7) | (word_cling_64_16_l(f01[j2]<<12)>>6);
        case 12: t[j+11] |= (word_cling_64_16_l(f00[j2]<<16)>>7) | (word_cling_64_16_l(f01[j2]<<16)>>6);
        case 11: t[j+10] |= (word_cling_64_16_l(f00[j2]<<20)>>7) | (word_cling_64_16_l(f01[j2]<<20)>>6);
        case 10: t[j+ 9] |= (word_cling_64_16_l(f00[j2]<<24)>>7) | (word_cling_64_16_l(f01[j2]<<24)>>6);
        case  9: t[j+ 8] |= (word_cling_64_16_l(f00[j2]<<28)>>7) | (word_cling_64_16_l(f01[j2]<<28)>>6);
        case  8: t[j+ 7] |= (word_cling_64_16_l(f00[j2]<<32)>>7) | (word_cling_64_16_l(f01[j2]<<32)>>6);
        case  7: t[j+ 6] |= (word_cling_64_16_l(f00[j2]<<36)>>7) | (word_cling_64_16_l(f01[j2]<<36)>>6);
        case  6: t[j+ 5] |= (word_cling_64_16_l(f00[j2]<<40)>>7) | (word_cling_64_16_l(f01[j2]<<40)>>6);
        case  5: t[j+ 4] |= (word_cling_64_16_l(f00[j2]<<44)>>7) | (word_cling_64_16_l(f01[j2]<<44)>>6);
        case  4: t[j+ 3] |= (word_cling_64_16_l(f00[j2]<<48)>>7) | (word_cling_64_16_l(f01[j2]<<48)>>6);
        case  3: t[j+ 2] |= (word_cling_64_16_l(f00[j2]<<52)>>7) | (word_cling_64_16_l(f01[j2]<<52)>>6);
        case  2: t[j+ 1] |= (word_cling_64_16_l(f00[j2]<<56)>>7) | (word_cling_64_16_l(f01[j2]<<56)>>6);
        case  1: t[j+ 0] |= (word_cling_64_16_l(f00[j2]<<60)>>7) | (word_cling_64_16_l(f01[j2]<<60)>>6);
          break;
        default:
          m4ri_die("impossible");
        }
        t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
      }
    }
      break;
    case 11: {
      for(rci_t i=0; i<T->nrows; i++) {
        const word *f00 = F->x[ 8]->rows[i];
        const word *f01 = F->x[ 9]->rows[i];
        const word *f02 = F->x[10]->rows[i];
        word *t  = T->x->rows[i];

        for(j=0, j2=0; j+16 < T->x->width; j+=16, j2++) {
          t[j+ 0] |= (word_cling_64_16_l(f00[j2]<<60)>>7) | (word_cling_64_16_l(f01[j2]<<60)>>6) | (word_cling_64_16_l(f02[j2]<<60)>>5);
          t[j+ 1] |= (word_cling_64_16_l(f00[j2]<<56)>>7) | (word_cling_64_16_l(f01[j2]<<56)>>6) | (word_cling_64_16_l(f02[j2]<<56)>>5);
          t[j+ 2] |= (word_cling_64_16_l(f00[j2]<<52)>>7) | (word_cling_64_16_l(f01[j2]<<52)>>6) | (word_cling_64_16_l(f02[j2]<<52)>>5);
          t[j+ 3] |= (word_cling_64_16_l(f00[j2]<<48)>>7) | (word_cling_64_16_l(f01[j2]<<48)>>6) | (word_cling_64_16_l(f02[j2]<<48)>>5);
          t[j+ 4] |= (word_cling_64_16_l(f00[j2]<<44)>>7) | (word_cling_64_16_l(f01[j2]<<44)>>6) | (word_cling_64_16_l(f02[j2]<<44)>>5);
          t[j+ 5] |= (word_cling_64_16_l(f00[j2]<<40)>>7) | (word_cling_64_16_l(f01[j2]<<40)>>6) | (word_cling_64_16_l(f02[j2]<<40)>>5);
          t[j+ 6] |= (word_cling_64_16_l(f00[j2]<<36)>>7) | (word_cling_64_16_l(f01[j2]<<36)>>6) | (word_cling_64_16_l(f02[j2]<<36)>>5);
          t[j+ 7] |= (word_cling_64_16_l(f00[j2]<<32)>>7) | (word_cling_64_16_l(f01[j2]<<32)>>6) | (word_cling_64_16_l(f02[j2]<<32)>>5);
          t[j+ 8] |= (word_cling_64_16_l(f00[j2]<<28)>>7) | (word_cling_64_16_l(f01[j2]<<28)>>6) | (word_cling_64_16_l(f02[j2]<<28)>>5);
          t[j+ 9] |= (word_cling_64_16_l(f00[j2]<<24)>>7) | (word_cling_64_16_l(f01[j2]<<24)>>6) | (word_cling_64_16_l(f02[j2]<<24)>>5);
          t[j+10] |= (word_cling_64_16_l(f00[j2]<<20)>>7) | (word_cling_64_16_l(f01[j2]<<20)>>6) | (word_cling_64_16_l(f02[j2]<<20)>>5);
          t[j+11] |= (word_cling_64_16_l(f00[j2]<<16)>>7) | (word_cling_64_16_l(f01[j2]<<16)>>6) | (word_cling_64_16_l(f02[j2]<<16)>>5);
          t[j+12] |= (word_cling_64_16_l(f00[j2]<<12)>>7) | (word_cling_64_16_l(f01[j2]<<12)>>6) | (word_cling_64_16_l(f02[j2]<<12)>>5);
          t[j+13] |= (word_cling_64_16_l(f00[j2]<< 8)>>7) | (word_cling_64_16_l(f01[j2]<< 8)>>6) | (word_cling_64_16_l(f02[j2]<< 8)>>5);
          t[j+14] |= (word_cling_64_16_l(f00[j2]<< 4)>>7) | (word_cling_64_16_l(f01[j2]<< 4)>>6) | (word_cling_64_16_l(f02[j2]<< 4)>>5);
          t[j+15] |= (word_cling_64_16_l(f00[j2]<< 0)>>7) | (word_cling_64_16_l(f01[j2]<< 0)>>6) | (word_cling_64_16_l(f02[j2]<< 0)>>5);
        }

        register word tmp = t[T->x->width-1];
        switch(T->x->width - j) {
        case 16: t[j+15] |= (word_cling_64_16_l(f00[j2]<< 0)>>7) | (word_cling_64_16_l(f01[j2]<< 0)>>6) | (word_cling_64_16_l(f02[j2]<< 0)>>5);
        case 15: t[j+14] |= (word_cling_64_16_l(f00[j2]<< 4)>>7) | (word_cling_64_16_l(f01[j2]<< 4)>>6) | (word_cling_64_16_l(f02[j2]<< 4)>>5);
        case 14: t[j+13] |= (word_cling_64_16_l(f00[j2]<< 8)>>7) | (word_cling_64_16_l(f01[j2]<< 8)>>6) | (word_cling_64_16_l(f02[j2]<< 8)>>5);
        case 13: t[j+12] |= (word_cling_64_16_l(f00[j2]<<12)>>7) | (word_cling_64_16_l(f01[j2]<<12)>>6) | (word_cling_64_16_l(f02[j2]<<12)>>5);
        case 12: t[j+11] |= (word_cling_64_16_l(f00[j2]<<16)>>7) | (word_cling_64_16_l(f01[j2]<<16)>>6) | (word_cling_64_16_l(f02[j2]<<16)>>5);
        case 11: t[j+10] |= (word_cling_64_16_l(f00[j2]<<20)>>7) | (word_cling_64_16_l(f01[j2]<<20)>>6) | (word_cling_64_16_l(f02[j2]<<20)>>5);
        case 10: t[j+ 9] |= (word_cling_64_16_l(f00[j2]<<24)>>7) | (word_cling_64_16_l(f01[j2]<<24)>>6) | (word_cling_64_16_l(f02[j2]<<24)>>5);
        case  9: t[j+ 8] |= (word_cling_64_16_l(f00[j2]<<28)>>7) | (word_cling_64_16_l(f01[j2]<<28)>>6) | (word_cling_64_16_l(f02[j2]<<28)>>5);
        case  8: t[j+ 7] |= (word_cling_64_16_l(f00[j2]<<32)>>7) | (word_cling_64_16_l(f01[j2]<<32)>>6) | (word_cling_64_16_l(f02[j2]<<32)>>5);
        case  7: t[j+ 6] |= (word_cling_64_16_l(f00[j2]<<36)>>7) | (word_cling_64_16_l(f01[j2]<<36)>>6) | (word_cling_64_16_l(f02[j2]<<36)>>5);
        case  6: t[j+ 5] |= (word_cling_64_16_l(f00[j2]<<40)>>7) | (word_cling_64_16_l(f01[j2]<<40)>>6) | (word_cling_64_16_l(f02[j2]<<40)>>5);
        case  5: t[j+ 4] |= (word_cling_64_16_l(f00[j2]<<44)>>7) | (word_cling_64_16_l(f01[j2]<<44)>>6) | (word_cling_64_16_l(f02[j2]<<44)>>5);
        case  4: t[j+ 3] |= (word_cling_64_16_l(f00[j2]<<48)>>7) | (word_cling_64_16_l(f01[j2]<<48)>>6) | (word_cling_64_16_l(f02[j2]<<48)>>5);
        case  3: t[j+ 2] |= (word_cling_64_16_l(f00[j2]<<52)>>7) | (word_cling_64_16_l(f01[j2]<<52)>>6) | (word_cling_64_16_l(f02[j2]<<52)>>5);
        case  2: t[j+ 1] |= (word_cling_64_16_l(f00[j2]<<56)>>7) | (word_cling_64_16_l(f01[j2]<<56)>>6) | (word_cling_64_16_l(f02[j2]<<56)>>5);
        case  1: t[j+ 0] |= (word_cling_64_16_l(f00[j2]<<60)>>7) | (word_cling_64_16_l(f01[j2]<<60)>>6) | (word_cling_64_16_l(f02[j2]<<60)>>5);
          break;
        default:
          m4ri_die("impossible");
        }
        t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
      }
    }
      break;
    }
  } else {
    for(rci_t i=0; i<T->nrows; i++) {
      const word *f00 = F->x[ 8]->rows[i];
      const word *f01 = F->x[ 9]->rows[i];
      const word *f02 = F->x[10]->rows[i];
      const word *f03 = F->x[11]->rows[i];
      word *t  = T->x->rows[i];

      for(j=0, j2=0; j+16 < T->x->width; j+=16, j2++) {
        t[j+ 0] |= (word_cling_64_16_l(f00[j2]<<60)>>7) | (word_cling_64_16_l(f01[j2]<<60)>>6) | (word_cling_64_16_l(f02[j2]<<60)>>5) | (word_cling_64_16_l(f03[j2]<<60)>>4);
        t[j+ 1] |= (word_cling_64_16_l(f00[j2]<<56)>>7) | (word_cling_64_16_l(f01[j2]<<56)>>6) | (word_cling_64_16_l(f02[j2]<<56)>>5) | (word_cling_64_16_l(f03[j2]<<56)>>4);
        t[j+ 2] |= (word_cling_64_16_l(f00[j2]<<52)>>7) | (word_cling_64_16_l(f01[j2]<<52)>>6) | (word_cling_64_16_l(f02[j2]<<52)>>5) | (word_cling_64_16_l(f03[j2]<<52)>>4);
        t[j+ 3] |= (word_cling_64_16_l(f00[j2]<<48)>>7) | (word_cling_64_16_l(f01[j2]<<48)>>6) | (word_cling_64_16_l(f02[j2]<<48)>>5) | (word_cling_64_16_l(f03[j2]<<48)>>4);
        t[j+ 4] |= (word_cling_64_16_l(f00[j2]<<44)>>7) | (word_cling_64_16_l(f01[j2]<<44)>>6) | (word_cling_64_16_l(f02[j2]<<44)>>5) | (word_cling_64_16_l(f03[j2]<<44)>>4);
        t[j+ 5] |= (word_cling_64_16_l(f00[j2]<<40)>>7) | (word_cling_64_16_l(f01[j2]<<40)>>6) | (word_cling_64_16_l(f02[j2]<<40)>>5) | (word_cling_64_16_l(f03[j2]<<40)>>4);
        t[j+ 6] |= (word_cling_64_16_l(f00[j2]<<36)>>7) | (word_cling_64_16_l(f01[j2]<<36)>>6) | (word_cling_64_16_l(f02[j2]<<36)>>5) | (word_cling_64_16_l(f03[j2]<<36)>>4);
        t[j+ 7] |= (word_cling_64_16_l(f00[j2]<<32)>>7) | (word_cling_64_16_l(f01[j2]<<32)>>6) | (word_cling_64_16_l(f02[j2]<<32)>>5) | (word_cling_64_16_l(f03[j2]<<32)>>4);
        t[j+ 8] |= (word_cling_64_16_l(f00[j2]<<28)>>7) | (word_cling_64_16_l(f01[j2]<<28)>>6) | (word_cling_64_16_l(f02[j2]<<28)>>5) | (word_cling_64_16_l(f03[j2]<<28)>>4);
        t[j+ 9] |= (word_cling_64_16_l(f00[j2]<<24)>>7) | (word_cling_64_16_l(f01[j2]<<24)>>6) | (word_cling_64_16_l(f02[j2]<<24)>>5) | (word_cling_64_16_l(f03[j2]<<24)>>4);
        t[j+10] |= (word_cling_64_16_l(f00[j2]<<20)>>7) | (word_cling_64_16_l(f01[j2]<<20)>>6) | (word_cling_64_16_l(f02[j2]<<20)>>5) | (word_cling_64_16_l(f03[j2]<<20)>>4);
        t[j+11] |= (word_cling_64_16_l(f00[j2]<<16)>>7) | (word_cling_64_16_l(f01[j2]<<16)>>6) | (word_cling_64_16_l(f02[j2]<<16)>>5) | (word_cling_64_16_l(f03[j2]<<16)>>4);
        t[j+12] |= (word_cling_64_16_l(f00[j2]<<12)>>7) | (word_cling_64_16_l(f01[j2]<<12)>>6) | (word_cling_64_16_l(f02[j2]<<12)>>5) | (word_cling_64_16_l(f03[j2]<<12)>>4);
        t[j+13] |= (word_cling_64_16_l(f00[j2]<< 8)>>7) | (word_cling_64_16_l(f01[j2]<< 8)>>6) | (word_cling_64_16_l(f02[j2]<< 8)>>5) | (word_cling_64_16_l(f03[j2]<< 8)>>4);
        t[j+14] |= (word_cling_64_16_l(f00[j2]<< 4)>>7) | (word_cling_64_16_l(f01[j2]<< 4)>>6) | (word_cling_64_16_l(f02[j2]<< 4)>>5) | (word_cling_64_16_l(f03[j2]<< 4)>>4);
        t[j+15] |= (word_cling_64_16_l(f00[j2]<< 0)>>7) | (word_cling_64_16_l(f01[j2]<< 0)>>6) | (word_cling_64_16_l(f02[j2]<< 0)>>5) | (word_cling_64_16_l(f03[j2]<< 0)>>4);
      }

      register word tmp = t[T->x->width-1];
      switch(T->x->width - j) {
      case 16: t[j+15] |= (word_cling_64_16_l(f00[j2]<< 0)>>7) | (word_cling_64_16_l(f01[j2]<< 0)>>6) | (word_cling_64_16_l(f02[j2]<< 0)>>5) | (word_cling_64_16_l(f03[j2]<< 0)>>4);
      case 15: t[j+14] |= (word_cling_64_16_l(f00[j2]<< 4)>>7) | (word_cling_64_16_l(f01[j2]<< 4)>>6) | (word_cling_64_16_l(f02[j2]<< 4)>>5) | (word_cling_64_16_l(f03[j2]<< 4)>>4);
      case 14: t[j+13] |= (word_cling_64_16_l(f00[j2]<< 8)>>7) | (word_cling_64_16_l(f01[j2]<< 8)>>6) | (word_cling_64_16_l(f02[j2]<< 8)>>5) | (word_cling_64_16_l(f03[j2]<< 8)>>4);
      case 13: t[j+12] |= (word_cling_64_16_l(f00[j2]<<12)>>7) | (word_cling_64_16_l(f01[j2]<<12)>>6) | (word_cling_64_16_l(f02[j2]<<12)>>5) | (word_cling_64_16_l(f03[j2]<<12)>>4);
      case 12: t[j+11] |= (word_cling_64_16_l(f00[j2]<<16)>>7) | (word_cling_64_16_l(f01[j2]<<16)>>6) | (word_cling_64_16_l(f02[j2]<<16)>>5) | (word_cling_64_16_l(f03[j2]<<16)>>4);
      case 11: t[j+10] |= (word_cling_64_16_l(f00[j2]<<20)>>7) | (word_cling_64_16_l(f01[j2]<<20)>>6) | (word_cling_64_16_l(f02[j2]<<20)>>5) | (word_cling_64_16_l(f03[j2]<<20)>>4);
      case 10: t[j+ 9] |= (word_cling_64_16_l(f00[j2]<<24)>>7) | (word_cling_64_16_l(f01[j2]<<24)>>6) | (word_cling_64_16_l(f02[j2]<<24)>>5) | (word_cling_64_16_l(f03[j2]<<24)>>4);
      case  9: t[j+ 8] |= (word_cling_64_16_l(f00[j2]<<28)>>7) | (word_cling_64_16_l(f01[j2]<<28)>>6) | (word_cling_64_16_l(f02[j2]<<28)>>5) | (word_cling_64_16_l(f03[j2]<<28)>>4);
      case  8: t[j+ 7] |= (word_cling_64_16_l(f00[j2]<<32)>>7) | (word_cling_64_16_l(f01[j2]<<32)>>6) | (word_cling_64_16_l(f02[j2]<<32)>>5) | (word_cling_64_16_l(f03[j2]<<32)>>4);
      case  7: t[j+ 6] |= (word_cling_64_16_l(f00[j2]<<36)>>7) | (word_cling_64_16_l(f01[j2]<<36)>>6) | (word_cling_64_16_l(f02[j2]<<36)>>5) | (word_cling_64_16_l(f03[j2]<<36)>>4);
      case  6: t[j+ 5] |= (word_cling_64_16_l(f00[j2]<<40)>>7) | (word_cling_64_16_l(f01[j2]<<40)>>6) | (word_cling_64_16_l(f02[j2]<<40)>>5) | (word_cling_64_16_l(f03[j2]<<40)>>4);
      case  5: t[j+ 4] |= (word_cling_64_16_l(f00[j2]<<44)>>7) | (word_cling_64_16_l(f01[j2]<<44)>>6) | (word_cling_64_16_l(f02[j2]<<44)>>5) | (word_cling_64_16_l(f03[j2]<<44)>>4);
      case  4: t[j+ 3] |= (word_cling_64_16_l(f00[j2]<<48)>>7) | (word_cling_64_16_l(f01[j2]<<48)>>6) | (word_cling_64_16_l(f02[j2]<<48)>>5) | (word_cling_64_16_l(f03[j2]<<48)>>4);
      case  3: t[j+ 2] |= (word_cling_64_16_l(f00[j2]<<52)>>7) | (word_cling_64_16_l(f01[j2]<<52)>>6) | (word_cling_64_16_l(f02[j2]<<52)>>5) | (word_cling_64_16_l(f03[j2]<<52)>>4);
      case  2: t[j+ 1] |= (word_cling_64_16_l(f00[j2]<<56)>>7) | (word_cling_64_16_l(f01[j2]<<56)>>6) | (word_cling_64_16_l(f02[j2]<<56)>>5) | (word_cling_64_16_l(f03[j2]<<56)>>4);
      case  1: t[j+ 0] |= (word_cling_64_16_l(f00[j2]<<60)>>7) | (word_cling_64_16_l(f01[j2]<<60)>>6) | (word_cling_64_16_l(f02[j2]<<60)>>5) | (word_cling_64_16_l(f03[j2]<<60)>>4);
        break;
      default:
        m4ri_die("impossible");
      }
      t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);

      switch(T->finite_field->degree) {
      case 13: {
        for(rci_t i=0; i<T->nrows; i++) {
          const word *f00 = F->x[12]->rows[i];
          word *t  = T->x->rows[i];

          for(j=0, j2=0; j+16 < T->x->width; j+=16, j2++) {
            t[j+ 0] |= (word_cling_64_16_l(f00[j2]<<60)>>3);
            t[j+ 1] |= (word_cling_64_16_l(f00[j2]<<56)>>3);
            t[j+ 2] |= (word_cling_64_16_l(f00[j2]<<52)>>3);
            t[j+ 3] |= (word_cling_64_16_l(f00[j2]<<48)>>3);
            t[j+ 4] |= (word_cling_64_16_l(f00[j2]<<44)>>3);
            t[j+ 5] |= (word_cling_64_16_l(f00[j2]<<40)>>3);
            t[j+ 6] |= (word_cling_64_16_l(f00[j2]<<36)>>3);
            t[j+ 7] |= (word_cling_64_16_l(f00[j2]<<32)>>3);
            t[j+ 8] |= (word_cling_64_16_l(f00[j2]<<28)>>3);
            t[j+ 9] |= (word_cling_64_16_l(f00[j2]<<24)>>3);
            t[j+10] |= (word_cling_64_16_l(f00[j2]<<20)>>3);
            t[j+11] |= (word_cling_64_16_l(f00[j2]<<16)>>3);
            t[j+12] |= (word_cling_64_16_l(f00[j2]<<12)>>3);
            t[j+13] |= (word_cling_64_16_l(f00[j2]<< 8)>>3);
            t[j+14] |= (word_cling_64_16_l(f00[j2]<< 4)>>3);
            t[j+15] |= (word_cling_64_16_l(f00[j2]<< 0)>>3);
          }

          register word tmp = t[T->x->width-1];
          switch(T->x->width - j) {
          case 16: t[j+15] |= (word_cling_64_16_l(f00[j2]<< 0)>>3);
          case 15: t[j+14] |= (word_cling_64_16_l(f00[j2]<< 4)>>3);
          case 14: t[j+13] |= (word_cling_64_16_l(f00[j2]<< 8)>>3);
          case 13: t[j+12] |= (word_cling_64_16_l(f00[j2]<<12)>>3);
          case 12: t[j+11] |= (word_cling_64_16_l(f00[j2]<<16)>>3);
          case 11: t[j+10] |= (word_cling_64_16_l(f00[j2]<<20)>>3);
          case 10: t[j+ 9] |= (word_cling_64_16_l(f00[j2]<<24)>>3);
          case  9: t[j+ 8] |= (word_cling_64_16_l(f00[j2]<<28)>>3);
          case  8: t[j+ 7] |= (word_cling_64_16_l(f00[j2]<<32)>>3);
          case  7: t[j+ 6] |= (word_cling_64_16_l(f00[j2]<<36)>>3);
          case  6: t[j+ 5] |= (word_cling_64_16_l(f00[j2]<<40)>>3);
          case  5: t[j+ 4] |= (word_cling_64_16_l(f00[j2]<<44)>>3);
          case  4: t[j+ 3] |= (word_cling_64_16_l(f00[j2]<<48)>>3);
          case  3: t[j+ 2] |= (word_cling_64_16_l(f00[j2]<<52)>>3);
          case  2: t[j+ 1] |= (word_cling_64_16_l(f00[j2]<<56)>>3);
          case  1: t[j+ 0] |= (word_cling_64_16_l(f00[j2]<<60)>>3);
            break;
          default:
            m4ri_die("impossible");
          }
          t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
        }
      }
        break;
      case 14: {
        for(rci_t i=0; i<T->nrows; i++) {
          const word *f00 = F->x[12]->rows[i];
          const word *f01 = F->x[13]->rows[i];
          word *t  = T->x->rows[i];

          for(j=0, j2=0; j+16 < T->x->width; j+=16, j2++) {
            t[j+ 0] |= (word_cling_64_16_l(f00[j2]<<60)>>3) | (word_cling_64_16_l(f01[j2]<<60)>>2);
            t[j+ 1] |= (word_cling_64_16_l(f00[j2]<<56)>>3) | (word_cling_64_16_l(f01[j2]<<56)>>2);
            t[j+ 2] |= (word_cling_64_16_l(f00[j2]<<52)>>3) | (word_cling_64_16_l(f01[j2]<<52)>>2);
            t[j+ 3] |= (word_cling_64_16_l(f00[j2]<<48)>>3) | (word_cling_64_16_l(f01[j2]<<48)>>2);
            t[j+ 4] |= (word_cling_64_16_l(f00[j2]<<44)>>3) | (word_cling_64_16_l(f01[j2]<<44)>>2);
            t[j+ 5] |= (word_cling_64_16_l(f00[j2]<<40)>>3) | (word_cling_64_16_l(f01[j2]<<40)>>2);
            t[j+ 6] |= (word_cling_64_16_l(f00[j2]<<36)>>3) | (word_cling_64_16_l(f01[j2]<<36)>>2);
            t[j+ 7] |= (word_cling_64_16_l(f00[j2]<<32)>>3) | (word_cling_64_16_l(f01[j2]<<32)>>2);
            t[j+ 8] |= (word_cling_64_16_l(f00[j2]<<28)>>3) | (word_cling_64_16_l(f01[j2]<<28)>>2);
            t[j+ 9] |= (word_cling_64_16_l(f00[j2]<<24)>>3) | (word_cling_64_16_l(f01[j2]<<24)>>2);
            t[j+10] |= (word_cling_64_16_l(f00[j2]<<20)>>3) | (word_cling_64_16_l(f01[j2]<<20)>>2);
            t[j+11] |= (word_cling_64_16_l(f00[j2]<<16)>>3) | (word_cling_64_16_l(f01[j2]<<16)>>2);
            t[j+12] |= (word_cling_64_16_l(f00[j2]<<12)>>3) | (word_cling_64_16_l(f01[j2]<<12)>>2);
            t[j+13] |= (word_cling_64_16_l(f00[j2]<< 8)>>3) | (word_cling_64_16_l(f01[j2]<< 8)>>2);
            t[j+14] |= (word_cling_64_16_l(f00[j2]<< 4)>>3) | (word_cling_64_16_l(f01[j2]<< 4)>>2);
            t[j+15] |= (word_cling_64_16_l(f00[j2]<< 0)>>3) | (word_cling_64_16_l(f01[j2]<< 0)>>2);
          }

          register word tmp = t[T->x->width-1];
          switch(T->x->width - j) {
          case 16: t[j+15] |= (word_cling_64_16_l(f00[j2]<< 0)>>3) | (word_cling_64_16_l(f01[j2]<< 0)>>2);
          case 15: t[j+14] |= (word_cling_64_16_l(f00[j2]<< 4)>>3) | (word_cling_64_16_l(f01[j2]<< 4)>>2);
          case 14: t[j+13] |= (word_cling_64_16_l(f00[j2]<< 8)>>3) | (word_cling_64_16_l(f01[j2]<< 8)>>2);
          case 13: t[j+12] |= (word_cling_64_16_l(f00[j2]<<12)>>3) | (word_cling_64_16_l(f01[j2]<<12)>>2);
          case 12: t[j+11] |= (word_cling_64_16_l(f00[j2]<<16)>>3) | (word_cling_64_16_l(f01[j2]<<16)>>2);
          case 11: t[j+10] |= (word_cling_64_16_l(f00[j2]<<20)>>3) | (word_cling_64_16_l(f01[j2]<<20)>>2);
          case 10: t[j+ 9] |= (word_cling_64_16_l(f00[j2]<<24)>>3) | (word_cling_64_16_l(f01[j2]<<24)>>2);
          case  9: t[j+ 8] |= (word_cling_64_16_l(f00[j2]<<28)>>3) | (word_cling_64_16_l(f01[j2]<<28)>>2);
          case  8: t[j+ 7] |= (word_cling_64_16_l(f00[j2]<<32)>>3) | (word_cling_64_16_l(f01[j2]<<32)>>2);
          case  7: t[j+ 6] |= (word_cling_64_16_l(f00[j2]<<36)>>3) | (word_cling_64_16_l(f01[j2]<<36)>>2);
          case  6: t[j+ 5] |= (word_cling_64_16_l(f00[j2]<<40)>>3) | (word_cling_64_16_l(f01[j2]<<40)>>2);
          case  5: t[j+ 4] |= (word_cling_64_16_l(f00[j2]<<44)>>3) | (word_cling_64_16_l(f01[j2]<<44)>>2);
          case  4: t[j+ 3] |= (word_cling_64_16_l(f00[j2]<<48)>>3) | (word_cling_64_16_l(f01[j2]<<48)>>2);
          case  3: t[j+ 2] |= (word_cling_64_16_l(f00[j2]<<52)>>3) | (word_cling_64_16_l(f01[j2]<<52)>>2);
          case  2: t[j+ 1] |= (word_cling_64_16_l(f00[j2]<<56)>>3) | (word_cling_64_16_l(f01[j2]<<56)>>2);
          case  1: t[j+ 0] |= (word_cling_64_16_l(f00[j2]<<60)>>3) | (word_cling_64_16_l(f01[j2]<<60)>>2);
            break;
          default:
            m4ri_die("impossible");
          }
          t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
        }
      }
        break;
      case 15: {
        for(rci_t i=0; i<T->nrows; i++) {
          const word *f00 = F->x[12]->rows[i];
          const word *f01 = F->x[13]->rows[i];
          const word *f02 = F->x[14]->rows[i];
          word *t  = T->x->rows[i];

          for(j=0, j2=0; j+16 < T->x->width; j+=16, j2++) {
            t[j+ 0] |= (word_cling_64_16_l(f00[j2]<<60)>>3) | (word_cling_64_16_l(f01[j2]<<60)>>2) | (word_cling_64_16_l(f02[j2]<<60)>>1);
            t[j+ 1] |= (word_cling_64_16_l(f00[j2]<<56)>>3) | (word_cling_64_16_l(f01[j2]<<56)>>2) | (word_cling_64_16_l(f02[j2]<<56)>>1);
            t[j+ 2] |= (word_cling_64_16_l(f00[j2]<<52)>>3) | (word_cling_64_16_l(f01[j2]<<52)>>2) | (word_cling_64_16_l(f02[j2]<<52)>>1);
            t[j+ 3] |= (word_cling_64_16_l(f00[j2]<<48)>>3) | (word_cling_64_16_l(f01[j2]<<48)>>2) | (word_cling_64_16_l(f02[j2]<<48)>>1);
            t[j+ 4] |= (word_cling_64_16_l(f00[j2]<<44)>>3) | (word_cling_64_16_l(f01[j2]<<44)>>2) | (word_cling_64_16_l(f02[j2]<<44)>>1);
            t[j+ 5] |= (word_cling_64_16_l(f00[j2]<<40)>>3) | (word_cling_64_16_l(f01[j2]<<40)>>2) | (word_cling_64_16_l(f02[j2]<<40)>>1);
            t[j+ 6] |= (word_cling_64_16_l(f00[j2]<<36)>>3) | (word_cling_64_16_l(f01[j2]<<36)>>2) | (word_cling_64_16_l(f02[j2]<<36)>>1);
            t[j+ 7] |= (word_cling_64_16_l(f00[j2]<<32)>>3) | (word_cling_64_16_l(f01[j2]<<32)>>2) | (word_cling_64_16_l(f02[j2]<<32)>>1);
            t[j+ 8] |= (word_cling_64_16_l(f00[j2]<<28)>>3) | (word_cling_64_16_l(f01[j2]<<28)>>2) | (word_cling_64_16_l(f02[j2]<<28)>>1);
            t[j+ 9] |= (word_cling_64_16_l(f00[j2]<<24)>>3) | (word_cling_64_16_l(f01[j2]<<24)>>2) | (word_cling_64_16_l(f02[j2]<<24)>>1);
            t[j+10] |= (word_cling_64_16_l(f00[j2]<<20)>>3) | (word_cling_64_16_l(f01[j2]<<20)>>2) | (word_cling_64_16_l(f02[j2]<<20)>>1);
            t[j+11] |= (word_cling_64_16_l(f00[j2]<<16)>>3) | (word_cling_64_16_l(f01[j2]<<16)>>2) | (word_cling_64_16_l(f02[j2]<<16)>>1);
            t[j+12] |= (word_cling_64_16_l(f00[j2]<<12)>>3) | (word_cling_64_16_l(f01[j2]<<12)>>2) | (word_cling_64_16_l(f02[j2]<<12)>>1);
            t[j+13] |= (word_cling_64_16_l(f00[j2]<< 8)>>3) | (word_cling_64_16_l(f01[j2]<< 8)>>2) | (word_cling_64_16_l(f02[j2]<< 8)>>1);
            t[j+14] |= (word_cling_64_16_l(f00[j2]<< 4)>>3) | (word_cling_64_16_l(f01[j2]<< 4)>>2) | (word_cling_64_16_l(f02[j2]<< 4)>>1);
            t[j+15] |= (word_cling_64_16_l(f00[j2]<< 0)>>3) | (word_cling_64_16_l(f01[j2]<< 0)>>2) | (word_cling_64_16_l(f02[j2]<< 0)>>1);
          }

          register word tmp = t[T->x->width-1];
          switch(T->x->width - j) {
          case 16: t[j+15] |= (word_cling_64_16_l(f00[j2]<< 0)>>3) | (word_cling_64_16_l(f01[j2]<< 0)>>2) | (word_cling_64_16_l(f02[j2]<< 0)>>1);
          case 15: t[j+14] |= (word_cling_64_16_l(f00[j2]<< 4)>>3) | (word_cling_64_16_l(f01[j2]<< 4)>>2) | (word_cling_64_16_l(f02[j2]<< 4)>>1);
          case 14: t[j+13] |= (word_cling_64_16_l(f00[j2]<< 8)>>3) | (word_cling_64_16_l(f01[j2]<< 8)>>2) | (word_cling_64_16_l(f02[j2]<< 8)>>1);
          case 13: t[j+12] |= (word_cling_64_16_l(f00[j2]<<12)>>3) | (word_cling_64_16_l(f01[j2]<<12)>>2) | (word_cling_64_16_l(f02[j2]<<12)>>1);
          case 12: t[j+11] |= (word_cling_64_16_l(f00[j2]<<16)>>3) | (word_cling_64_16_l(f01[j2]<<16)>>2) | (word_cling_64_16_l(f02[j2]<<16)>>1);
          case 11: t[j+10] |= (word_cling_64_16_l(f00[j2]<<20)>>3) | (word_cling_64_16_l(f01[j2]<<20)>>2) | (word_cling_64_16_l(f02[j2]<<20)>>1);
          case 10: t[j+ 9] |= (word_cling_64_16_l(f00[j2]<<24)>>3) | (word_cling_64_16_l(f01[j2]<<24)>>2) | (word_cling_64_16_l(f02[j2]<<24)>>1);
          case  9: t[j+ 8] |= (word_cling_64_16_l(f00[j2]<<28)>>3) | (word_cling_64_16_l(f01[j2]<<28)>>2) | (word_cling_64_16_l(f02[j2]<<28)>>1);
          case  8: t[j+ 7] |= (word_cling_64_16_l(f00[j2]<<32)>>3) | (word_cling_64_16_l(f01[j2]<<32)>>2) | (word_cling_64_16_l(f02[j2]<<32)>>1);
          case  7: t[j+ 6] |= (word_cling_64_16_l(f00[j2]<<36)>>3) | (word_cling_64_16_l(f01[j2]<<36)>>2) | (word_cling_64_16_l(f02[j2]<<36)>>1);
          case  6: t[j+ 5] |= (word_cling_64_16_l(f00[j2]<<40)>>3) | (word_cling_64_16_l(f01[j2]<<40)>>2) | (word_cling_64_16_l(f02[j2]<<40)>>1);
          case  5: t[j+ 4] |= (word_cling_64_16_l(f00[j2]<<44)>>3) | (word_cling_64_16_l(f01[j2]<<44)>>2) | (word_cling_64_16_l(f02[j2]<<44)>>1);
          case  4: t[j+ 3] |= (word_cling_64_16_l(f00[j2]<<48)>>3) | (word_cling_64_16_l(f01[j2]<<48)>>2) | (word_cling_64_16_l(f02[j2]<<48)>>1);
          case  3: t[j+ 2] |= (word_cling_64_16_l(f00[j2]<<52)>>3) | (word_cling_64_16_l(f01[j2]<<52)>>2) | (word_cling_64_16_l(f02[j2]<<52)>>1);
          case  2: t[j+ 1] |= (word_cling_64_16_l(f00[j2]<<56)>>3) | (word_cling_64_16_l(f01[j2]<<56)>>2) | (word_cling_64_16_l(f02[j2]<<56)>>1);
          case  1: t[j+ 0] |= (word_cling_64_16_l(f00[j2]<<60)>>3) | (word_cling_64_16_l(f01[j2]<<60)>>2) | (word_cling_64_16_l(f02[j2]<<60)>>1);
            break;
          default:
            m4ri_die("impossible");
          }
          t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
        }
      }
        break;
      case 16: {
        for(rci_t i=0; i<T->nrows; i++) {
          const word *f00 = F->x[12]->rows[i];
          const word *f01 = F->x[13]->rows[i];
          const word *f02 = F->x[14]->rows[i];
          const word *f03 = F->x[15]->rows[i];
          word *t  = T->x->rows[i];

          for(j=0, j2=0; j+16 < T->x->width; j+=16, j2++) {
            t[j+ 0] |= (word_cling_64_16_l(f00[j2]<<60)>>3) | (word_cling_64_16_l(f01[j2]<<60)>>2) | (word_cling_64_16_l(f02[j2]<<60)>>1) | (word_cling_64_16_l(f03[j2]<<60)>>0);
            t[j+ 1] |= (word_cling_64_16_l(f00[j2]<<56)>>3) | (word_cling_64_16_l(f01[j2]<<56)>>2) | (word_cling_64_16_l(f02[j2]<<56)>>1) | (word_cling_64_16_l(f03[j2]<<56)>>0);
            t[j+ 2] |= (word_cling_64_16_l(f00[j2]<<52)>>3) | (word_cling_64_16_l(f01[j2]<<52)>>2) | (word_cling_64_16_l(f02[j2]<<52)>>1) | (word_cling_64_16_l(f03[j2]<<52)>>0);
            t[j+ 3] |= (word_cling_64_16_l(f00[j2]<<48)>>3) | (word_cling_64_16_l(f01[j2]<<48)>>2) | (word_cling_64_16_l(f02[j2]<<48)>>1) | (word_cling_64_16_l(f03[j2]<<48)>>0);
            t[j+ 4] |= (word_cling_64_16_l(f00[j2]<<44)>>3) | (word_cling_64_16_l(f01[j2]<<44)>>2) | (word_cling_64_16_l(f02[j2]<<44)>>1) | (word_cling_64_16_l(f03[j2]<<44)>>0);
            t[j+ 5] |= (word_cling_64_16_l(f00[j2]<<40)>>3) | (word_cling_64_16_l(f01[j2]<<40)>>2) | (word_cling_64_16_l(f02[j2]<<40)>>1) | (word_cling_64_16_l(f03[j2]<<40)>>0);
            t[j+ 6] |= (word_cling_64_16_l(f00[j2]<<36)>>3) | (word_cling_64_16_l(f01[j2]<<36)>>2) | (word_cling_64_16_l(f02[j2]<<36)>>1) | (word_cling_64_16_l(f03[j2]<<36)>>0);
            t[j+ 7] |= (word_cling_64_16_l(f00[j2]<<32)>>3) | (word_cling_64_16_l(f01[j2]<<32)>>2) | (word_cling_64_16_l(f02[j2]<<32)>>1) | (word_cling_64_16_l(f03[j2]<<32)>>0);
            t[j+ 8] |= (word_cling_64_16_l(f00[j2]<<28)>>3) | (word_cling_64_16_l(f01[j2]<<28)>>2) | (word_cling_64_16_l(f02[j2]<<28)>>1) | (word_cling_64_16_l(f03[j2]<<28)>>0);
            t[j+ 9] |= (word_cling_64_16_l(f00[j2]<<24)>>3) | (word_cling_64_16_l(f01[j2]<<24)>>2) | (word_cling_64_16_l(f02[j2]<<24)>>1) | (word_cling_64_16_l(f03[j2]<<24)>>0);
            t[j+10] |= (word_cling_64_16_l(f00[j2]<<20)>>3) | (word_cling_64_16_l(f01[j2]<<20)>>2) | (word_cling_64_16_l(f02[j2]<<20)>>1) | (word_cling_64_16_l(f03[j2]<<20)>>0);
            t[j+11] |= (word_cling_64_16_l(f00[j2]<<16)>>3) | (word_cling_64_16_l(f01[j2]<<16)>>2) | (word_cling_64_16_l(f02[j2]<<16)>>1) | (word_cling_64_16_l(f03[j2]<<16)>>0);
            t[j+12] |= (word_cling_64_16_l(f00[j2]<<12)>>3) | (word_cling_64_16_l(f01[j2]<<12)>>2) | (word_cling_64_16_l(f02[j2]<<12)>>1) | (word_cling_64_16_l(f03[j2]<<12)>>0);
            t[j+13] |= (word_cling_64_16_l(f00[j2]<< 8)>>3) | (word_cling_64_16_l(f01[j2]<< 8)>>2) | (word_cling_64_16_l(f02[j2]<< 8)>>1) | (word_cling_64_16_l(f03[j2]<< 8)>>0);
            t[j+14] |= (word_cling_64_16_l(f00[j2]<< 4)>>3) | (word_cling_64_16_l(f01[j2]<< 4)>>2) | (word_cling_64_16_l(f02[j2]<< 4)>>1) | (word_cling_64_16_l(f03[j2]<< 4)>>0);
            t[j+15] |= (word_cling_64_16_l(f00[j2]<< 0)>>3) | (word_cling_64_16_l(f01[j2]<< 0)>>2) | (word_cling_64_16_l(f02[j2]<< 0)>>1) | (word_cling_64_16_l(f03[j2]<< 0)>>0);
          }

          register word tmp = t[T->x->width-1];
          switch(T->x->width - j) {
          case 16: t[j+15] |= (word_cling_64_16_l(f00[j2]<< 0)>>3) | (word_cling_64_16_l(f01[j2]<< 0)>>2) | (word_cling_64_16_l(f02[j2]<< 0)>>1) | (word_cling_64_16_l(f03[j2]<< 0)>>0);
          case 15: t[j+14] |= (word_cling_64_16_l(f00[j2]<< 4)>>3) | (word_cling_64_16_l(f01[j2]<< 4)>>2) | (word_cling_64_16_l(f02[j2]<< 4)>>1) | (word_cling_64_16_l(f03[j2]<< 4)>>0);
          case 14: t[j+13] |= (word_cling_64_16_l(f00[j2]<< 8)>>3) | (word_cling_64_16_l(f01[j2]<< 8)>>2) | (word_cling_64_16_l(f02[j2]<< 8)>>1) | (word_cling_64_16_l(f03[j2]<< 8)>>0);
          case 13: t[j+12] |= (word_cling_64_16_l(f00[j2]<<12)>>3) | (word_cling_64_16_l(f01[j2]<<12)>>2) | (word_cling_64_16_l(f02[j2]<<12)>>1) | (word_cling_64_16_l(f03[j2]<<12)>>0);
          case 12: t[j+11] |= (word_cling_64_16_l(f00[j2]<<16)>>3) | (word_cling_64_16_l(f01[j2]<<16)>>2) | (word_cling_64_16_l(f02[j2]<<16)>>1) | (word_cling_64_16_l(f03[j2]<<16)>>0);
          case 11: t[j+10] |= (word_cling_64_16_l(f00[j2]<<20)>>3) | (word_cling_64_16_l(f01[j2]<<20)>>2) | (word_cling_64_16_l(f02[j2]<<20)>>1) | (word_cling_64_16_l(f03[j2]<<20)>>0);
          case 10: t[j+ 9] |= (word_cling_64_16_l(f00[j2]<<24)>>3) | (word_cling_64_16_l(f01[j2]<<24)>>2) | (word_cling_64_16_l(f02[j2]<<24)>>1) | (word_cling_64_16_l(f03[j2]<<24)>>0);
          case  9: t[j+ 8] |= (word_cling_64_16_l(f00[j2]<<28)>>3) | (word_cling_64_16_l(f01[j2]<<28)>>2) | (word_cling_64_16_l(f02[j2]<<28)>>1) | (word_cling_64_16_l(f03[j2]<<28)>>0);
          case  8: t[j+ 7] |= (word_cling_64_16_l(f00[j2]<<32)>>3) | (word_cling_64_16_l(f01[j2]<<32)>>2) | (word_cling_64_16_l(f02[j2]<<32)>>1) | (word_cling_64_16_l(f03[j2]<<32)>>0);
          case  7: t[j+ 6] |= (word_cling_64_16_l(f00[j2]<<36)>>3) | (word_cling_64_16_l(f01[j2]<<36)>>2) | (word_cling_64_16_l(f02[j2]<<36)>>1) | (word_cling_64_16_l(f03[j2]<<36)>>0);
          case  6: t[j+ 5] |= (word_cling_64_16_l(f00[j2]<<40)>>3) | (word_cling_64_16_l(f01[j2]<<40)>>2) | (word_cling_64_16_l(f02[j2]<<40)>>1) | (word_cling_64_16_l(f03[j2]<<40)>>0);
          case  5: t[j+ 4] |= (word_cling_64_16_l(f00[j2]<<44)>>3) | (word_cling_64_16_l(f01[j2]<<44)>>2) | (word_cling_64_16_l(f02[j2]<<44)>>1) | (word_cling_64_16_l(f03[j2]<<44)>>0);
          case  4: t[j+ 3] |= (word_cling_64_16_l(f00[j2]<<48)>>3) | (word_cling_64_16_l(f01[j2]<<48)>>2) | (word_cling_64_16_l(f02[j2]<<48)>>1) | (word_cling_64_16_l(f03[j2]<<48)>>0);
          case  3: t[j+ 2] |= (word_cling_64_16_l(f00[j2]<<52)>>3) | (word_cling_64_16_l(f01[j2]<<52)>>2) | (word_cling_64_16_l(f02[j2]<<52)>>1) | (word_cling_64_16_l(f03[j2]<<52)>>0);
          case  2: t[j+ 1] |= (word_cling_64_16_l(f00[j2]<<56)>>3) | (word_cling_64_16_l(f01[j2]<<56)>>2) | (word_cling_64_16_l(f02[j2]<<56)>>1) | (word_cling_64_16_l(f03[j2]<<56)>>0);
          case  1: t[j+ 0] |= (word_cling_64_16_l(f00[j2]<<60)>>3) | (word_cling_64_16_l(f01[j2]<<60)>>2) | (word_cling_64_16_l(f02[j2]<<60)>>1) | (word_cling_64_16_l(f03[j2]<<60)>>0);
            break;
          default:
            m4ri_die("impossible");
          }
          t[T->x->width-1] = (t[T->x->width-1] & bitmask_end) | (tmp & ~bitmask_end);
        }
      }
        break;
      }
    }
  }
  return T;
}

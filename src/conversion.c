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

static const word x80808080 = 0x8080808080808080ULL;
static const word x88888888 = 0x8888888888888888ULL;
static const word xaaaaaaaa = 0xaaaaaaaaaaaaaaaaULL;
static const word xcccccccc = 0xccccccccccccccccULL;
static const word xf0f0f0f0 = 0xf0f0f0f0f0f0f0f0ULL;
static const word xff00ff00 = 0xff00ff00ff00ff00ULL;
static const word xffff0000 = 0xffff0000ffff0000ULL;
static const word xffffffff = 0xffffffff00000000ULL;
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
  case  9:
  case 10:
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
  case  9:
  case 10:
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

      register word tmp=0;
      switch(T->x->width - j) {
      case 8:
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
        tmp    = (word_cling_64_08_l(f0[j2]<< 0)>>7) | (word_cling_64_08_l(f1[j2]<< 0)>>6) | (word_cling_64_08_l(f2[j2]<< 0)>>5) | (word_cling_64_08_l(f3[j2]<< 0)>>4) \
          |      (word_cling_64_08_l(f4[j2]<< 0)>>3) | (word_cling_64_08_l(f5[j2]<< 0)>>2) | (word_cling_64_08_l(f6[j2]<< 0)>>1) | (word_cling_64_08_l(f7[j2]<< 0)>>0);
        t[j+7] = (t[j+7] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 7:
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
        tmp    = (word_cling_64_08_l(f0[j2]<< 8)>>7) | (word_cling_64_08_l(f1[j2]<< 8)>>6) | (word_cling_64_08_l(f2[j2]<< 8)>>5) | (word_cling_64_08_l(f3[j2]<< 8)>>4) \
          |      (word_cling_64_08_l(f4[j2]<< 8)>>3) | (word_cling_64_08_l(f5[j2]<< 8)>>2) | (word_cling_64_08_l(f6[j2]<< 8)>>1) | (word_cling_64_08_l(f7[j2]<< 8)>>0);
        t[j+6] = (t[j+6] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 6:
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
        tmp    = (word_cling_64_08_l(f0[j2]<<16)>>7) | (word_cling_64_08_l(f1[j2]<<16)>>6) | (word_cling_64_08_l(f2[j2]<<16)>>5) | (word_cling_64_08_l(f3[j2]<<16)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<16)>>3) | (word_cling_64_08_l(f5[j2]<<16)>>2) | (word_cling_64_08_l(f6[j2]<<16)>>1) | (word_cling_64_08_l(f7[j2]<<16)>>0);
        t[j+5] = (t[j+5] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 5:
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3) | (word_cling_64_08_l(f5[j2]<<56)>>2) | (word_cling_64_08_l(f6[j2]<<56)>>1) | (word_cling_64_08_l(f7[j2]<<56)>>0);
        t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3) | (word_cling_64_08_l(f5[j2]<<48)>>2) | (word_cling_64_08_l(f6[j2]<<48)>>1) | (word_cling_64_08_l(f7[j2]<<48)>>0);
        t[j+2] = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<40)>>3) | (word_cling_64_08_l(f5[j2]<<40)>>2) | (word_cling_64_08_l(f6[j2]<<40)>>1) | (word_cling_64_08_l(f7[j2]<<40)>>0);
        t[j+3] = (word_cling_64_08_l(f0[j2]<<32)>>7) | (word_cling_64_08_l(f1[j2]<<32)>>6) | (word_cling_64_08_l(f2[j2]<<32)>>5) | (word_cling_64_08_l(f3[j2]<<32)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<32)>>3) | (word_cling_64_08_l(f5[j2]<<32)>>2) | (word_cling_64_08_l(f6[j2]<<32)>>1) | (word_cling_64_08_l(f7[j2]<<32)>>0);
        tmp    = (word_cling_64_08_l(f0[j2]<<24)>>7) | (word_cling_64_08_l(f1[j2]<<24)>>6) | (word_cling_64_08_l(f2[j2]<<24)>>5) | (word_cling_64_08_l(f3[j2]<<24)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<24)>>3) | (word_cling_64_08_l(f5[j2]<<24)>>2) | (word_cling_64_08_l(f6[j2]<<24)>>1) | (word_cling_64_08_l(f7[j2]<<24)>>0);
        t[j+4] = (t[j+4] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 4:
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3) | (word_cling_64_08_l(f5[j2]<<56)>>2) | (word_cling_64_08_l(f6[j2]<<56)>>1) | (word_cling_64_08_l(f7[j2]<<56)>>0);
        t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3) | (word_cling_64_08_l(f5[j2]<<48)>>2) | (word_cling_64_08_l(f6[j2]<<48)>>1) | (word_cling_64_08_l(f7[j2]<<48)>>0);
        t[j+2] = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<40)>>3) | (word_cling_64_08_l(f5[j2]<<40)>>2) | (word_cling_64_08_l(f6[j2]<<40)>>1) | (word_cling_64_08_l(f7[j2]<<40)>>0);
        tmp    = (word_cling_64_08_l(f0[j2]<<32)>>7) | (word_cling_64_08_l(f1[j2]<<32)>>6) | (word_cling_64_08_l(f2[j2]<<32)>>5) | (word_cling_64_08_l(f3[j2]<<32)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<32)>>3) | (word_cling_64_08_l(f5[j2]<<32)>>2) | (word_cling_64_08_l(f6[j2]<<32)>>1) | (word_cling_64_08_l(f7[j2]<<32)>>0);
        t[j+3] = (t[j+3] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 3:
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3) | (word_cling_64_08_l(f5[j2]<<56)>>2) | (word_cling_64_08_l(f6[j2]<<56)>>1) | (word_cling_64_08_l(f7[j2]<<56)>>0);
        t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3) | (word_cling_64_08_l(f5[j2]<<48)>>2) | (word_cling_64_08_l(f6[j2]<<48)>>1) | (word_cling_64_08_l(f7[j2]<<48)>>0);
        tmp    = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<40)>>3) | (word_cling_64_08_l(f5[j2]<<40)>>2) | (word_cling_64_08_l(f6[j2]<<40)>>1) | (word_cling_64_08_l(f7[j2]<<40)>>0);
        t[j+2] = (t[j+2] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 2:
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3) | (word_cling_64_08_l(f5[j2]<<56)>>2) | (word_cling_64_08_l(f6[j2]<<56)>>1) | (word_cling_64_08_l(f7[j2]<<56)>>0);
        tmp    = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3) | (word_cling_64_08_l(f5[j2]<<48)>>2) | (word_cling_64_08_l(f6[j2]<<48)>>1) | (word_cling_64_08_l(f7[j2]<<48)>>0);
        t[j+1] = (t[j+1] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 1:
        tmp    = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3) | (word_cling_64_08_l(f5[j2]<<56)>>2) | (word_cling_64_08_l(f6[j2]<<56)>>1) | (word_cling_64_08_l(f7[j2]<<56)>>0);
        t[j+0] = (t[j+0] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      default:
        m4ri_die("impossible");
      } //switch
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

      register word tmp=0;
      switch(T->x->width - j) {
      case 8:
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
        tmp    = (word_cling_64_08_l(f0[j2]<< 0)>>7) | (word_cling_64_08_l(f1[j2]<< 0)>>6) | (word_cling_64_08_l(f2[j2]<< 0)>>5) | (word_cling_64_08_l(f3[j2]<< 0)>>4) \
          |      (word_cling_64_08_l(f4[j2]<< 0)>>3) | (word_cling_64_08_l(f5[j2]<< 0)>>2) | (word_cling_64_08_l(f6[j2]<< 0)>>1);
        t[j+7] = (t[j+7] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 7:
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
        tmp    = (word_cling_64_08_l(f0[j2]<< 8)>>7) | (word_cling_64_08_l(f1[j2]<< 8)>>6) | (word_cling_64_08_l(f2[j2]<< 8)>>5) | (word_cling_64_08_l(f3[j2]<< 8)>>4) \
          |      (word_cling_64_08_l(f4[j2]<< 8)>>3) | (word_cling_64_08_l(f5[j2]<< 8)>>2) | (word_cling_64_08_l(f6[j2]<< 8)>>1);
        t[j+6] = (t[j+6] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 6:
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
        tmp    = (word_cling_64_08_l(f0[j2]<<16)>>7) | (word_cling_64_08_l(f1[j2]<<16)>>6) | (word_cling_64_08_l(f2[j2]<<16)>>5) | (word_cling_64_08_l(f3[j2]<<16)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<16)>>3) | (word_cling_64_08_l(f5[j2]<<16)>>2) | (word_cling_64_08_l(f6[j2]<<16)>>1);
        t[j+5] = (t[j+5] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 5:
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3) | (word_cling_64_08_l(f5[j2]<<56)>>2) | (word_cling_64_08_l(f6[j2]<<56)>>1);
        t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3) | (word_cling_64_08_l(f5[j2]<<48)>>2) | (word_cling_64_08_l(f6[j2]<<48)>>1);
        t[j+2] = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<40)>>3) | (word_cling_64_08_l(f5[j2]<<40)>>2) | (word_cling_64_08_l(f6[j2]<<40)>>1);
        t[j+3] = (word_cling_64_08_l(f0[j2]<<32)>>7) | (word_cling_64_08_l(f1[j2]<<32)>>6) | (word_cling_64_08_l(f2[j2]<<32)>>5) | (word_cling_64_08_l(f3[j2]<<32)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<32)>>3) | (word_cling_64_08_l(f5[j2]<<32)>>2) | (word_cling_64_08_l(f6[j2]<<32)>>1);
        tmp    = (word_cling_64_08_l(f0[j2]<<24)>>7) | (word_cling_64_08_l(f1[j2]<<24)>>6) | (word_cling_64_08_l(f2[j2]<<24)>>5) | (word_cling_64_08_l(f3[j2]<<24)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<24)>>3) | (word_cling_64_08_l(f5[j2]<<24)>>2) | (word_cling_64_08_l(f6[j2]<<24)>>1);
        t[j+4] = (t[j+4] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 4:
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3) | (word_cling_64_08_l(f5[j2]<<56)>>2) | (word_cling_64_08_l(f6[j2]<<56)>>1);
        t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3) | (word_cling_64_08_l(f5[j2]<<48)>>2) | (word_cling_64_08_l(f6[j2]<<48)>>1);
        t[j+2] = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<40)>>3) | (word_cling_64_08_l(f5[j2]<<40)>>2) | (word_cling_64_08_l(f6[j2]<<40)>>1);
        tmp    = (word_cling_64_08_l(f0[j2]<<32)>>7) | (word_cling_64_08_l(f1[j2]<<32)>>6) | (word_cling_64_08_l(f2[j2]<<32)>>5) | (word_cling_64_08_l(f3[j2]<<32)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<32)>>3) | (word_cling_64_08_l(f5[j2]<<32)>>2) | (word_cling_64_08_l(f6[j2]<<32)>>1);
        t[j+3] = (t[j+3] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 3:
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3) | (word_cling_64_08_l(f5[j2]<<56)>>2) | (word_cling_64_08_l(f6[j2]<<56)>>1);
        t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3) | (word_cling_64_08_l(f5[j2]<<48)>>2) | (word_cling_64_08_l(f6[j2]<<48)>>1);
        tmp    = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<40)>>3) | (word_cling_64_08_l(f5[j2]<<40)>>2) | (word_cling_64_08_l(f6[j2]<<40)>>1);
        t[j+2] = (t[j+2] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 2:
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3) | (word_cling_64_08_l(f5[j2]<<56)>>2) | (word_cling_64_08_l(f6[j2]<<56)>>1);
        tmp    = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3) | (word_cling_64_08_l(f5[j2]<<48)>>2) | (word_cling_64_08_l(f6[j2]<<48)>>1);
        t[j+1] = (t[j+1] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 1:
        tmp    = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3) | (word_cling_64_08_l(f5[j2]<<56)>>2) | (word_cling_64_08_l(f6[j2]<<56)>>1);
        t[j+0] = (t[j+0] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      default:
        m4ri_die("impossible");
      }
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

      register word tmp=0;
      switch(T->x->width - j) {
      case 8:
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
        tmp    = (word_cling_64_08_l(f0[j2]<< 0)>>7) | (word_cling_64_08_l(f1[j2]<< 0)>>6) | (word_cling_64_08_l(f2[j2]<< 0)>>5) | (word_cling_64_08_l(f3[j2]<< 0)>>4) \
          |      (word_cling_64_08_l(f4[j2]<< 0)>>3) | (word_cling_64_08_l(f5[j2]<< 0)>>2);
        t[j+7] = (t[j+7] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 7:
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
        tmp    = (word_cling_64_08_l(f0[j2]<< 8)>>7) | (word_cling_64_08_l(f1[j2]<< 8)>>6) | (word_cling_64_08_l(f2[j2]<< 8)>>5) | (word_cling_64_08_l(f3[j2]<< 8)>>4) \
          |      (word_cling_64_08_l(f4[j2]<< 8)>>3) | (word_cling_64_08_l(f5[j2]<< 8)>>2);
        t[j+6] = (t[j+6] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 6:
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
        tmp    = (word_cling_64_08_l(f0[j2]<<16)>>7) | (word_cling_64_08_l(f1[j2]<<16)>>6) | (word_cling_64_08_l(f2[j2]<<16)>>5) | (word_cling_64_08_l(f3[j2]<<16)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<16)>>3) | (word_cling_64_08_l(f5[j2]<<16)>>2);
        t[j+5] = (t[j+5] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 5:
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3) | (word_cling_64_08_l(f5[j2]<<56)>>2);
        t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3) | (word_cling_64_08_l(f5[j2]<<48)>>2);
        t[j+2] = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<40)>>3) | (word_cling_64_08_l(f5[j2]<<40)>>2);
        t[j+3] = (word_cling_64_08_l(f0[j2]<<32)>>7) | (word_cling_64_08_l(f1[j2]<<32)>>6) | (word_cling_64_08_l(f2[j2]<<32)>>5) | (word_cling_64_08_l(f3[j2]<<32)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<32)>>3) | (word_cling_64_08_l(f5[j2]<<32)>>2);
        tmp    = (word_cling_64_08_l(f0[j2]<<24)>>7) | (word_cling_64_08_l(f1[j2]<<24)>>6) | (word_cling_64_08_l(f2[j2]<<24)>>5) | (word_cling_64_08_l(f3[j2]<<24)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<24)>>3) | (word_cling_64_08_l(f5[j2]<<24)>>2);
        t[j+4] = (t[j+4] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 4:
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3) | (word_cling_64_08_l(f5[j2]<<56)>>2);
        t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3) | (word_cling_64_08_l(f5[j2]<<48)>>2);
        t[j+2] = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<40)>>3) | (word_cling_64_08_l(f5[j2]<<40)>>2);
        tmp    = (word_cling_64_08_l(f0[j2]<<32)>>7) | (word_cling_64_08_l(f1[j2]<<32)>>6) | (word_cling_64_08_l(f2[j2]<<32)>>5) | (word_cling_64_08_l(f3[j2]<<32)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<32)>>3) | (word_cling_64_08_l(f5[j2]<<32)>>2);
        t[j+3] = (t[j+3] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 3:
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3) | (word_cling_64_08_l(f5[j2]<<56)>>2);
        t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3) | (word_cling_64_08_l(f5[j2]<<48)>>2);
        tmp    = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<40)>>3) | (word_cling_64_08_l(f5[j2]<<40)>>2);
        t[j+2] = (t[j+2] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 2:
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3) | (word_cling_64_08_l(f5[j2]<<56)>>2);
        tmp    = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3) | (word_cling_64_08_l(f5[j2]<<48)>>2);
        t[j+1] = (t[j+1] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 1:
        tmp    = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3) | (word_cling_64_08_l(f5[j2]<<56)>>2);
        t[j+0] = (t[j+0] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      default:
        m4ri_die("impossible");
      }
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
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3);
        t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3);
        t[j+2] = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<40)>>3);
        t[j+3] = (word_cling_64_08_l(f0[j2]<<32)>>7) | (word_cling_64_08_l(f1[j2]<<32)>>6) | (word_cling_64_08_l(f2[j2]<<32)>>5) | (word_cling_64_08_l(f3[j2]<<32)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<32)>>3);
        t[j+4] = (word_cling_64_08_l(f0[j2]<<24)>>7) | (word_cling_64_08_l(f1[j2]<<24)>>6) | (word_cling_64_08_l(f2[j2]<<24)>>5) | (word_cling_64_08_l(f3[j2]<<24)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<24)>>3);
        t[j+5] = (word_cling_64_08_l(f0[j2]<<16)>>7) | (word_cling_64_08_l(f1[j2]<<16)>>6) | (word_cling_64_08_l(f2[j2]<<16)>>5) | (word_cling_64_08_l(f3[j2]<<16)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<16)>>3);
        t[j+6] = (word_cling_64_08_l(f0[j2]<< 8)>>7) | (word_cling_64_08_l(f1[j2]<< 8)>>6) | (word_cling_64_08_l(f2[j2]<< 8)>>5) | (word_cling_64_08_l(f3[j2]<< 8)>>4) \
          |      (word_cling_64_08_l(f4[j2]<< 8)>>3);
        t[j+7] = (word_cling_64_08_l(f0[j2]<< 0)>>7) | (word_cling_64_08_l(f1[j2]<< 0)>>6) | (word_cling_64_08_l(f2[j2]<< 0)>>5) | (word_cling_64_08_l(f3[j2]<< 0)>>4) \
          |      (word_cling_64_08_l(f4[j2]<< 0)>>3);
      }

      register word tmp=0;
      switch(T->x->width - j) {
      case 8:
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3);
        t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3);
        t[j+2] = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<40)>>3);
        t[j+3] = (word_cling_64_08_l(f0[j2]<<32)>>7) | (word_cling_64_08_l(f1[j2]<<32)>>6) | (word_cling_64_08_l(f2[j2]<<32)>>5) | (word_cling_64_08_l(f3[j2]<<32)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<32)>>3);
        t[j+4] = (word_cling_64_08_l(f0[j2]<<24)>>7) | (word_cling_64_08_l(f1[j2]<<24)>>6) | (word_cling_64_08_l(f2[j2]<<24)>>5) | (word_cling_64_08_l(f3[j2]<<24)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<24)>>3);
        t[j+5] = (word_cling_64_08_l(f0[j2]<<16)>>7) | (word_cling_64_08_l(f1[j2]<<16)>>6) | (word_cling_64_08_l(f2[j2]<<16)>>5) | (word_cling_64_08_l(f3[j2]<<16)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<16)>>3);
        t[j+6] = (word_cling_64_08_l(f0[j2]<< 8)>>7) | (word_cling_64_08_l(f1[j2]<< 8)>>6) | (word_cling_64_08_l(f2[j2]<< 8)>>5) | (word_cling_64_08_l(f3[j2]<< 8)>>4) \
          |      (word_cling_64_08_l(f4[j2]<< 8)>>3);
        tmp    = (word_cling_64_08_l(f0[j2]<< 0)>>7) | (word_cling_64_08_l(f1[j2]<< 0)>>6) | (word_cling_64_08_l(f2[j2]<< 0)>>5) | (word_cling_64_08_l(f3[j2]<< 0)>>4) \
          |      (word_cling_64_08_l(f4[j2]<< 0)>>3);
        t[j+7] = (t[j+7] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 7:
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3);
        t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3);
        t[j+2] = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<40)>>3);
        t[j+3] = (word_cling_64_08_l(f0[j2]<<32)>>7) | (word_cling_64_08_l(f1[j2]<<32)>>6) | (word_cling_64_08_l(f2[j2]<<32)>>5) | (word_cling_64_08_l(f3[j2]<<32)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<32)>>3);
        t[j+4] = (word_cling_64_08_l(f0[j2]<<24)>>7) | (word_cling_64_08_l(f1[j2]<<24)>>6) | (word_cling_64_08_l(f2[j2]<<24)>>5) | (word_cling_64_08_l(f3[j2]<<24)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<24)>>3);
        t[j+5] = (word_cling_64_08_l(f0[j2]<<16)>>7) | (word_cling_64_08_l(f1[j2]<<16)>>6) | (word_cling_64_08_l(f2[j2]<<16)>>5) | (word_cling_64_08_l(f3[j2]<<16)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<16)>>3);
        tmp    = (word_cling_64_08_l(f0[j2]<< 8)>>7) | (word_cling_64_08_l(f1[j2]<< 8)>>6) | (word_cling_64_08_l(f2[j2]<< 8)>>5) | (word_cling_64_08_l(f3[j2]<< 8)>>4) \
          |      (word_cling_64_08_l(f4[j2]<< 8)>>3);
        t[j+6] = (t[j+6] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 6:
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3);
        t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3);
        t[j+2] = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<40)>>3);
        t[j+3] = (word_cling_64_08_l(f0[j2]<<32)>>7) | (word_cling_64_08_l(f1[j2]<<32)>>6) | (word_cling_64_08_l(f2[j2]<<32)>>5) | (word_cling_64_08_l(f3[j2]<<32)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<32)>>3);
        t[j+4] = (word_cling_64_08_l(f0[j2]<<24)>>7) | (word_cling_64_08_l(f1[j2]<<24)>>6) | (word_cling_64_08_l(f2[j2]<<24)>>5) | (word_cling_64_08_l(f3[j2]<<24)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<24)>>3);
        tmp    = (word_cling_64_08_l(f0[j2]<<16)>>7) | (word_cling_64_08_l(f1[j2]<<16)>>6) | (word_cling_64_08_l(f2[j2]<<16)>>5) | (word_cling_64_08_l(f3[j2]<<16)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<16)>>3);
        t[j+5] = (t[j+5] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 5:
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3);
        t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3);
        t[j+2] = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<40)>>3);
        t[j+3] = (word_cling_64_08_l(f0[j2]<<32)>>7) | (word_cling_64_08_l(f1[j2]<<32)>>6) | (word_cling_64_08_l(f2[j2]<<32)>>5) | (word_cling_64_08_l(f3[j2]<<32)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<32)>>3);
        tmp    = (word_cling_64_08_l(f0[j2]<<24)>>7) | (word_cling_64_08_l(f1[j2]<<24)>>6) | (word_cling_64_08_l(f2[j2]<<24)>>5) | (word_cling_64_08_l(f3[j2]<<24)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<24)>>3);
        t[j+4] = (t[j+4] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 4:
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3);
        t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3);
        t[j+2] = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<40)>>3);
        tmp    = (word_cling_64_08_l(f0[j2]<<32)>>7) | (word_cling_64_08_l(f1[j2]<<32)>>6) | (word_cling_64_08_l(f2[j2]<<32)>>5) | (word_cling_64_08_l(f3[j2]<<32)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<32)>>3);
        t[j+3] = (t[j+3] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 3:
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3);
        t[j+1] = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3);
        tmp    = (word_cling_64_08_l(f0[j2]<<40)>>7) | (word_cling_64_08_l(f1[j2]<<40)>>6) | (word_cling_64_08_l(f2[j2]<<40)>>5) | (word_cling_64_08_l(f3[j2]<<40)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<40)>>3);
        t[j+2] = (t[j+2] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 2:
        t[j+0] = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3);
        tmp    = (word_cling_64_08_l(f0[j2]<<48)>>7) | (word_cling_64_08_l(f1[j2]<<48)>>6) | (word_cling_64_08_l(f2[j2]<<48)>>5) | (word_cling_64_08_l(f3[j2]<<48)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<48)>>3);
        t[j+1] = (t[j+1] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      case 1:
        tmp    = (word_cling_64_08_l(f0[j2]<<56)>>7) | (word_cling_64_08_l(f1[j2]<<56)>>6) | (word_cling_64_08_l(f2[j2]<<56)>>5) | (word_cling_64_08_l(f3[j2]<<56)>>4) \
          |      (word_cling_64_08_l(f4[j2]<<56)>>3);
        t[j+0] = (t[j+0] & ~bitmask_end) | (tmp & bitmask_end);
        break;
      default:
        m4ri_die("impossible");
      }
    }
  }
    break;
  default:
    m4ri_die("impossible");
  } 
  return T;
}

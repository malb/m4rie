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

#include "bitslice.h"

static const word x88888888 = 0x8888888888888888ULL;
static const word xaaaaaaaa = 0xaaaaaaaaaaaaaaaaULL;
static const word xcccccccc = 0xccccccccccccccccULL;
static const word xf0f0f0f0 = 0xf0f0f0f0f0f0f0f0ULL;
static const word xff00ff00 = 0xff00ff00ff00ff00ULL;
static const word xffff0000 = 0xffff0000ffff0000ULL;
static const word xffffffff = 0xffffffff00000000ULL;
static const word x__left16 = 0xffff000000000000ULL;

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

static inline word word_cling_64_02_l(word a) {
  a = (a & xffff0000) | (a & xffff0000>>16)>>16;
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
  case  5:
  case  6:
  case  7:
  case  8:
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
  case  5:
  case  6:
  case  7:
  case  8:
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
  register word r0,r1,r2,r3;

  const word bitmask_end = __M4RI_LEFT_BITMASK((T->x->offset + T->x->ncols) % m4ri_radix);

  if (mzd_slice_is_zero(F))
    return T;

  for(size_t i=0; i<T->nrows; i++) {
    word *f0 = F->x[0]->rows[i];
    word *f1 = F->x[1]->rows[i];
    word *t  = T->x->rows[i];

    for(j=0, j2=0; j+2 < T->x->width; j+=2, j2++) {
      if (!(f0[j2] | f1[j2]) )
        continue;
      r0 = word_cling_64_02_l(f0[j2]<<32 & xffffffff)>>1;
      r1 = word_cling_64_02_l(f0[j2]     & xffffffff)>>1;
      r2 = word_cling_64_02_l(f1[j2]<<32 & xffffffff)>>0;
      r3 = word_cling_64_02_l(f1[j2]     & xffffffff)>>0;

      t[j+0] = r0|r2;
      t[j+1] = r1|r3;
    }
    switch(T->x->width - j) {
    case 2:
      r0 = word_cling_64_02_l(f0[j2]<<32 & xffffffff)>>1;
      r1 = word_cling_64_02_l(f0[j2]     & xffffffff)>>1;
      r2 = word_cling_64_02_l(f1[j2]<<32 & xffffffff)>>0;
      r3 = word_cling_64_02_l(f1[j2]     & xffffffff)>>0;

      t[j+0] = r0|r2;      
      t[j+1] = (t[j+1] & ~bitmask_end) | ((r1|r3) & bitmask_end);
      break;
    case 1:
      r0 = word_cling_64_02_l(f0[j2]<<32 & xffffffff)>>1;
      r2 = word_cling_64_02_l(f1[j2]<<32 & xffffffff)>>0;

      t[j+0] = (t[j+0] & ~bitmask_end) | ((r0|r2) & bitmask_end);
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

  /* A0 */
  for(size_t i=0; i<T->nrows; i++) {
    word *t0 = T->x[0]->rows[i];
    word *t1 = T->x[1]->rows[i];
    word *t2 = T->x[2]->rows[i];
    const word const *f  = F->x->rows[i];

    /* bulk of work */
    for(j=0, j2=0; j+4 < F->x->width; j+=4,j2++) {
      r0 = word_slice_64_04_l(f[j+0]<<3 & x88888888)>>48;
      r1 = word_slice_64_04_l(f[j+1]<<3 & x88888888)>>32; 
      r2 = word_slice_64_04_l(f[j+2]<<3 & x88888888)>>16; 
      r3 = word_slice_64_04_l(f[j+3]<<3 & x88888888)>> 0; 
      t0[j2] = r0|r1|r2|r3;

      r0 = word_slice_64_04_l(f[j+0]<<2 & x88888888)>>48;
      r1 = word_slice_64_04_l(f[j+1]<<2 & x88888888)>>32; 
      r2 = word_slice_64_04_l(f[j+2]<<2 & x88888888)>>16; 
      r3 = word_slice_64_04_l(f[j+3]<<2 & x88888888)>> 0; 
      t1[j2] = r0|r1|r2|r3;

      r0 = word_slice_64_04_l(f[j+0]<<1 & x88888888)>>48;
      r1 = word_slice_64_04_l(f[j+1]<<1 & x88888888)>>32; 
      r2 = word_slice_64_04_l(f[j+2]<<1 & x88888888)>>16; 
      r3 = word_slice_64_04_l(f[j+3]<<1 & x88888888)>> 0; 
      t2[j2] = r0|r1|r2|r3;
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

  if(T->depth == 3)
    return T;

  /* A3 */
  for(size_t i=0; i<T->nrows; i++) {
    word *t3 = T->x[3]->rows[i];
    const word const *f  = F->x->rows[i];

    /* bulk of work */
    for(j=0, j2=0; j+4 < F->x->width; j+=4,j2++) {
      if ( !(f[j+0] | f[j+1] | f[j+2] | f[j+3]) )
        continue;
      r0 = word_slice_64_04_l(f[j+0]<<0 & x88888888)>>48;
      r1 = word_slice_64_04_l(f[j+1]<<0 & x88888888)>>32; 
      r2 = word_slice_64_04_l(f[j+2]<<0 & x88888888)>>16; 
      r3 = word_slice_64_04_l(f[j+3]<<0 & x88888888)>> 0; 
      t3[j2] = r0|r1|r2|r3;
    }
    r3 = 0;
    switch(F->x->width - j) {
    case 4:
      r3 |= word_slice_64_04_l(f[j+3]<<0 & x88888888)>> 0; 
    case 3:
      r3 |= word_slice_64_04_l(f[j+2]<<0 & x88888888)>>16; 
    case 2:
      r3 |= word_slice_64_04_l(f[j+1]<<0 & x88888888)>>32; 
    case 1:
      r3 |= word_slice_64_04_l(f[j+0]<<0 & x88888888)>>48;
      break;
    default:
      m4ri_die("impossible");
    }
    t3[j2] |= r3 & bitmask_end;    
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

void mzd_slice_set_ui(mzd_slice_t *A, word value) {
  for(int i=0; i<A->depth; i++) {
    mzd_set_ui(A->x[i], (value>>i)&1);
  }
}

mzd_slice_t *_mzd_slice_mul_karatsuba2(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B) {
  // two temporaries
  if (C == NULL)
    C = mzd_slice_init(A->finite_field, A->nrows, B->ncols);

  mzd_addmul(C->x[0], A->x[1], B->x[1], 0);  /* C0 += A1*B1 */

  mzd_t *T0 = mzd_addmul(NULL, A->x[0], B->x[0], 0);  /* A0B0 = A0*B0 */
  mzd_add(C->x[0], C->x[0], T0); /*C0 += A0*B0 */
  mzd_add(C->x[1], C->x[1], T0); /*C1 += A0*B0 */
  mzd_free(T0);

  T0 = mzd_add(NULL, A->x[1], A->x[0]); /*T0 = A1 + A0 */

  mzd_t *T1 = mzd_add(NULL, B->x[1], B->x[0]); /*T1 = B1 + B0 */

  mzd_addmul(C->x[1], T0, T1, 0); /* C1 += A0*B0 + T0*T1 */

  mzd_free(T0);  mzd_free(T1);

  return C;
}

mzd_slice_t *_mzd_slice_mul_karatsuba3(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B) {
  /* using three temporary matrices */
  if (C == NULL)
    C = mzd_slice_init(A->finite_field, A->nrows, B->ncols);

  C = _mzd_slice_adapt_depth(C,4);
  
  const mzd_t *a0 = A->x[0];
  const mzd_t *a1 = A->x[1];
  const mzd_t *a2 = A->x[2];

  const mzd_t *b0 = B->x[0];
  const mzd_t *b1 = B->x[1];
  const mzd_t *b2 = B->x[2];

  mzd_t *t0 = mzd_init(a0->nrows, a0->ncols);
  mzd_t *t1 = mzd_init(b0->nrows, b0->ncols);

  mzd_t **X = C->x;

  mzd_add(t0, a0, a1);
  mzd_add(t1, b0, b1);
  mzd_addmul(X[1], t0, t1, 0); /* + (a0+a1)(b0+b1)X */

  mzd_add(t0, a0, a2);
  mzd_add(t1, b0, b2);
  mzd_addmul(X[2], t0, t1, 0); /* + (a0+a2)(b0+b2)X^2 */

  mzd_add(t0, a1, a2);
  mzd_add(t1, b1, b2);
  mzd_addmul(X[3], t0, t1, 0); /* + (a1+a2)(b1+b2)X^3 */

  mzd_free(t0);
  mzd_free(t1);

  t0 = mzd_init(a0->nrows, b0->ncols);

  mzd_mul(t0, a0, b0, 0); /* + a0b0(1-X-X^2) */
  mzd_add(X[0], X[0], t0);
  mzd_add(X[1], X[1], t0);
  mzd_add(X[2], X[2], t0);

  mzd_mul(t0, a1, b1, 0); /* + a1b1(X+X^2-X^3) */
  mzd_add(X[1], X[1], t0);
  mzd_add(X[2], X[2], t0);
  mzd_add(X[3], X[3], t0);

  mzd_mul(t0, a2, b2, 0); /* + a2b2(-X^2-X^3+X^4) */
  
  /* modular reductions and final additions */

  if( (A->finite_field->minpoly & 1<<2) == 0)
    mzd_add(X[3], X[3], t0);
  else
    mzd_add(X[2], X[2], t0);
  mzd_add(X[1], X[1], t0);

  if(A->finite_field->minpoly & 1<<2) 
    mzd_add(X[2],X[2],X[3]);
  else  //if (A->finite_field->minpoly & 1<<1) {=
    mzd_add(X[1],X[1],X[3]);
  mzd_add(X[0],X[0],X[3]);

  mzd_free(t0);
  _mzd_slice_adapt_depth(C,3);

  return C;
}

static void _poly2_addmul(mzd_t **X, const mzd_t **a, const mzd_t **b) {
  mzd_t *t0 = mzd_init(a[0]->nrows, a[0]->ncols);
  mzd_t *t1 = mzd_init(b[0]->nrows, b[0]->ncols);

  mzd_add(t0, a[0], a[1]);
  mzd_add(t1, b[0], b[1]);

  mzd_addmul(X[1], t0, t1, 0); /* + (a0+a1)(b0+b1)X */

  mzd_free(t0);
  mzd_free(t1);

  t0 = mzd_init(a[0]->nrows, b[0]->ncols);

  mzd_mul(t0, a[0], b[0], 0); /* + a0b0(1-X) */
  mzd_add(X[0], X[0], t0);
  mzd_add(X[1], X[1], t0);

  mzd_mul(t0, a[1], b[1], 0); /* + a1b1(X+X^2) */
  mzd_add(X[1], X[1], t0);
  mzd_add(X[2], X[2], t0);

  mzd_free(t0);
}

static void _poly_add(mzd_t **c, const mzd_t **a, const mzd_t **b,const int length) {
  switch(length) {
  case 16: mzd_add(c[15], a[15], b[15]);
  case 15: mzd_add(c[14], a[14], b[14]);
  case 14: mzd_add(c[13], a[13], b[13]);
  case 13: mzd_add(c[12], a[12], b[12]);
  case 12: mzd_add(c[11], a[11], b[11]);
  case 11: mzd_add(c[10], a[10], b[10]);
  case 10: mzd_add(c[ 9], a[ 9], b[ 9]);
  case  9: mzd_add(c[ 8], a[ 8], b[ 8]);
  case  8: mzd_add(c[ 7], a[ 7], b[ 7]);
  case  7: mzd_add(c[ 6], a[ 6], b[ 6]);
  case  6: mzd_add(c[ 5], a[ 5], b[ 5]);
  case  5: mzd_add(c[ 4], a[ 4], b[ 4]);
  case  4: mzd_add(c[ 3], a[ 3], b[ 3]);
  case  3: mzd_add(c[ 2], a[ 2], b[ 2]);
  case  2: mzd_add(c[ 1], a[ 1], b[ 1]);
  case  1: mzd_add(c[ 0], a[ 0], b[ 0]);
    break;
  case 0:
  default:
    m4ri_die("this should never happen.");
  } 
}

mzd_slice_t *_mzd_slice_mul_karatsuba4(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B) {
  /* using five + two = 7 temporary matrices */

  if (C == NULL)
    C = mzd_slice_init(A->finite_field, A->nrows, B->ncols);

  C = _mzd_slice_adapt_depth(C,5);
  
  const mzd_t *a0[2] = {A->x[0],A->x[1]};
  const mzd_t *a1[2] = {A->x[2],A->x[3]};
  const mzd_t *b0[2] = {B->x[0],B->x[1]};
  const mzd_t *b1[2] = {B->x[2],B->x[3]};

  mzd_t *X[2][3] = { {C->x[0],C->x[1],C->x[2]}, 
                     {C->x[2],C->x[3],C->x[4]} };

  mzd_t *t0[3];
  mzd_t *t1[2];

  t0[0] = mzd_init(A->nrows, A->ncols);
  t0[1] = mzd_init(A->nrows, A->ncols);
  t1[0] = mzd_init(B->nrows, B->ncols);
  t1[1] = mzd_init(B->nrows, B->ncols);

  _poly_add(t0, a0, a1, 2);
  _poly_add(t1, b0, b1, 2);

  _poly2_addmul(X[1], (const mzd_t**)t0, (const mzd_t**)t1);

  mzd_free(t0[0]);  mzd_free(t0[1]);
  mzd_free(t1[0]);  mzd_free(t1[1]);

  t0[0] = mzd_init(A->nrows, B->ncols);
  t0[1] = mzd_init(A->nrows, B->ncols);
  t0[2] = mzd_init(A->nrows, B->ncols);

  _poly2_addmul(t0, a0, b0);
  _poly_add(X[0], (const mzd_t**)X[0], (const mzd_t**)t0, 3);
  _poly_add(X[1], (const mzd_t**)X[1], (const mzd_t**)t0, 3);

  mzd_set_ui(t0[0], 0);
  mzd_set_ui(t0[1], 0);
  mzd_set_ui(t0[2], 0);
  
  _poly2_addmul(t0, a1, b1);
  _poly_add(X[1], (const mzd_t**)X[1], (const mzd_t**)t0, 3);

  /* we would now do 
   *
   *    _poly_add(X[2], (const mzd_t**)X[2], (const mzd_t**)t0, 3);
   *
   * but we want avoid C->x[6] and C->x[5], hence we combine it with the
   * modular reduction.
   */

  mzd_add(C->x[4], C->x[4], t0[0]);
  if(A->finite_field->minpoly & 1<<3) {
    mzd_add(t0[1], t0[1], t0[2]);
    mzd_add(C->x[4], C->x[4], t0[1]);
    mzd_add(C->x[3], C->x[3], C->x[4]);
  }
  if(A->finite_field->minpoly & 1<<2) {
    mzd_add(C->x[4], C->x[4], t0[2]);
    mzd_add(C->x[3], C->x[3], t0[1]);
    mzd_add(C->x[2], C->x[2], C->x[4]);
  }
  if(A->finite_field->minpoly & 1<<1) {
    mzd_add(C->x[3], C->x[3], t0[2]);
    mzd_add(C->x[2], C->x[2], t0[1]);
    mzd_add(C->x[1], C->x[1], C->x[4]);
  }
  mzd_add(C->x[2], C->x[2], t0[2]);
  mzd_add(C->x[1], C->x[1], t0[1]);
  mzd_add(C->x[0], C->x[0], C->x[4]);

  mzd_free(t0[0]); mzd_free(t0[1]); mzd_free(t0[2]); 
  _mzd_slice_adapt_depth(C,4);

  return C;
}

void mzd_slice_print(const mzd_slice_t *A) {
  char formatstr[10];
  int width = gf2e_degree_to_w(A->finite_field)/4;
  if (gf2e_degree_to_w(A->finite_field)%4) 
    width += 1;
  sprintf(formatstr,"%%%dx",width);

  for (rci_t i=0; i < A->nrows; ++i) {
    printf("[");
    for (rci_t j=0; j < A->ncols; j++) {
      word tmp = mzd_slice_read_elem(A,i,j);
      printf(formatstr,(int)tmp);
      if(j<A->ncols-1)
        printf(" ");
    }
    printf("]\n");
  }
}

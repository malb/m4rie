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

mzd_slice_t *mzed_slice(mzd_slice_t *A, const mzed_t *Z) {
  if (A == NULL) {
    assert(Z->x->offset == 0);
    A = mzd_slice_init(Z->finite_field, Z->nrows, Z->ncols);
  } else {
    assert(Z->x->offset == (Z->w*A->x[0]->offset));
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
    assert(A->x->offset == (A->w*Z->x[0]->offset));
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

mzd_slice_t *_mzed_slice2(mzd_slice_t *A, const mzed_t *Z) {
  assert(A && (A->depth >= 2));
  size_t j, j2 = 0;
  register word t0,t1;

  const word bitmask_begin = __M4RI_RIGHT_BITMASK(m4ri_radix - A->x[0]->offset%m4ri_radix);
  const word bitmask_end = __M4RI_LEFT_BITMASK((A->x[0]->offset + A->ncols) % m4ri_radix);
  const word one = m4ri_one;

  for(size_t i=0; i<A->nrows; i++) {
    word *a0 = A->x[0]->rows[i];
    word *a1 = A->x[1]->rows[i];
    const word *z  = Z->x->rows[i];    

    const word a0_fix = a0[0]; 
    const word a1_fix = a1[0]; 

    /* bulk of work */
    for(j=0, j2=0; j+2 < Z->x->width; j+=2,j2++) {
      if ( !(z[j+0] | z[j+1]) )
        continue;
      t0 = t1 = 0;

      t0 |= (z[j+0] & one<< 0) >>  0; t0 |= (z[j+0] & one<< 2) >>  1; t0 |= (z[j+0] & one<< 4) >>  2; t0 |= (z[j+0] & one<< 6) >>  3; 
      t0 |= (z[j+0] & one<< 8) >>  4; t0 |= (z[j+0] & one<<10) >>  5; t0 |= (z[j+0] & one<<12) >>  6; t0 |= (z[j+0] & one<<14) >>  7; 
      t0 |= (z[j+0] & one<<16) >>  8; t0 |= (z[j+0] & one<<18) >>  9; t0 |= (z[j+0] & one<<20) >> 10; t0 |= (z[j+0] & one<<22) >> 11; 
      t0 |= (z[j+0] & one<<24) >> 12; t0 |= (z[j+0] & one<<26) >> 13; t0 |= (z[j+0] & one<<28) >> 14; t0 |= (z[j+0] & one<<30) >> 15; 
      t0 |= (z[j+0] & one<<32) >> 16; t0 |= (z[j+0] & one<<34) >> 17; t0 |= (z[j+0] & one<<36) >> 18; t0 |= (z[j+0] & one<<38) >> 19; 
      t0 |= (z[j+0] & one<<40) >> 20; t0 |= (z[j+0] & one<<42) >> 21; t0 |= (z[j+0] & one<<44) >> 22; t0 |= (z[j+0] & one<<46) >> 23; 
      t0 |= (z[j+0] & one<<48) >> 24; t0 |= (z[j+0] & one<<50) >> 25; t0 |= (z[j+0] & one<<52) >> 26; t0 |= (z[j+0] & one<<54) >> 27; 
      t0 |= (z[j+0] & one<<56) >> 28; t0 |= (z[j+0] & one<<58) >> 29; t0 |= (z[j+0] & one<<60) >> 30; t0 |= (z[j+0] & one<<62) >> 31; 

      t0 |= (z[j+1] & one<< 0) << 32; t0 |= (z[j+1] & one<< 2) << 31; t0 |= (z[j+1] & one<< 4) << 30; t0 |= (z[j+1] & one<< 6) << 29; 
      t0 |= (z[j+1] & one<< 8) << 28; t0 |= (z[j+1] & one<<10) << 27; t0 |= (z[j+1] & one<<12) << 26; t0 |= (z[j+1] & one<<14) << 25; 
      t0 |= (z[j+1] & one<<16) << 24; t0 |= (z[j+1] & one<<18) << 23; t0 |= (z[j+1] & one<<20) << 22; t0 |= (z[j+1] & one<<22) << 21; 
      t0 |= (z[j+1] & one<<24) << 20; t0 |= (z[j+1] & one<<26) << 19; t0 |= (z[j+1] & one<<28) << 18; t0 |= (z[j+1] & one<<30) << 17; 
      t0 |= (z[j+1] & one<<32) << 16; t0 |= (z[j+1] & one<<34) << 15; t0 |= (z[j+1] & one<<36) << 14; t0 |= (z[j+1] & one<<38) << 13; 
      t0 |= (z[j+1] & one<<40) << 12; t0 |= (z[j+1] & one<<42) << 11; t0 |= (z[j+1] & one<<44) << 10; t0 |= (z[j+1] & one<<46) <<  9; 
      t0 |= (z[j+1] & one<<48) <<  8; t0 |= (z[j+1] & one<<50) <<  7; t0 |= (z[j+1] & one<<52) <<  6; t0 |= (z[j+1] & one<<54) <<  5; 
      t0 |= (z[j+1] & one<<56) <<  4; t0 |= (z[j+1] & one<<58) <<  3; t0 |= (z[j+1] & one<<60) <<  2; t0 |= (z[j+1] & one<<62) <<  1;  

      t1 |= (z[j+0] & one<<( 0 + 1)) >>  1; t1 |= (z[j+0] & one<<( 2 + 1)) >>  2; t1 |= (z[j+0] & one<<( 4 + 1)) >>  3; t1 |= (z[j+0] & one<<( 6 + 1)) >>  4; 
      t1 |= (z[j+0] & one<<( 8 + 1)) >>  5; t1 |= (z[j+0] & one<<(10 + 1)) >>  6; t1 |= (z[j+0] & one<<(12 + 1)) >>  7; t1 |= (z[j+0] & one<<(14 + 1)) >>  8; 
      t1 |= (z[j+0] & one<<(16 + 1)) >>  9; t1 |= (z[j+0] & one<<(18 + 1)) >> 10; t1 |= (z[j+0] & one<<(20 + 1)) >> 11; t1 |= (z[j+0] & one<<(22 + 1)) >> 12; 
      t1 |= (z[j+0] & one<<(24 + 1)) >> 13; t1 |= (z[j+0] & one<<(26 + 1)) >> 14; t1 |= (z[j+0] & one<<(28 + 1)) >> 15; t1 |= (z[j+0] & one<<(30 + 1)) >> 16; 
      t1 |= (z[j+0] & one<<(32 + 1)) >> 17; t1 |= (z[j+0] & one<<(34 + 1)) >> 18; t1 |= (z[j+0] & one<<(36 + 1)) >> 19; t1 |= (z[j+0] & one<<(38 + 1)) >> 20; 
      t1 |= (z[j+0] & one<<(40 + 1)) >> 21; t1 |= (z[j+0] & one<<(42 + 1)) >> 22; t1 |= (z[j+0] & one<<(44 + 1)) >> 23; t1 |= (z[j+0] & one<<(46 + 1)) >> 24; 
      t1 |= (z[j+0] & one<<(48 + 1)) >> 25; t1 |= (z[j+0] & one<<(50 + 1)) >> 26; t1 |= (z[j+0] & one<<(52 + 1)) >> 27; t1 |= (z[j+0] & one<<(54 + 1)) >> 28; 
      t1 |= (z[j+0] & one<<(56 + 1)) >> 29; t1 |= (z[j+0] & one<<(58 + 1)) >> 30; t1 |= (z[j+0] & one<<(60 + 1)) >> 31; t1 |= (z[j+0] & one<<(62 + 1)) >> 32; 

      t1 |= (z[j+1] & one<<( 0 + 1)) << 31; t1 |= (z[j+1] & one<<( 2 + 1)) << 30; t1 |= (z[j+1] & one<<( 4 + 1)) << 29; t1 |= (z[j+1] & one<<( 6 + 1)) << 28; 
      t1 |= (z[j+1] & one<<( 8 + 1)) << 27; t1 |= (z[j+1] & one<<(10 + 1)) << 26; t1 |= (z[j+1] & one<<(12 + 1)) << 25; t1 |= (z[j+1] & one<<(14 + 1)) << 24; 
      t1 |= (z[j+1] & one<<(16 + 1)) << 23; t1 |= (z[j+1] & one<<(18 + 1)) << 22; t1 |= (z[j+1] & one<<(20 + 1)) << 21; t1 |= (z[j+1] & one<<(22 + 1)) << 20; 
      t1 |= (z[j+1] & one<<(24 + 1)) << 19; t1 |= (z[j+1] & one<<(26 + 1)) << 18; t1 |= (z[j+1] & one<<(28 + 1)) << 17; t1 |= (z[j+1] & one<<(30 + 1)) << 16; 
      t1 |= (z[j+1] & one<<(32 + 1)) << 15; t1 |= (z[j+1] & one<<(34 + 1)) << 14; t1 |= (z[j+1] & one<<(36 + 1)) << 13; t1 |= (z[j+1] & one<<(38 + 1)) << 12; 
      t1 |= (z[j+1] & one<<(40 + 1)) << 11; t1 |= (z[j+1] & one<<(42 + 1)) << 10; t1 |= (z[j+1] & one<<(44 + 1)) <<  9; t1 |= (z[j+1] & one<<(46 + 1)) <<  8; 
      t1 |= (z[j+1] & one<<(48 + 1)) <<  7; t1 |= (z[j+1] & one<<(50 + 1)) <<  6; t1 |= (z[j+1] & one<<(52 + 1)) <<  5; t1 |= (z[j+1] & one<<(54 + 1)) <<  4; 
      t1 |= (z[j+1] & one<<(56 + 1)) <<  3; t1 |= (z[j+1] & one<<(58 + 1)) <<  2; t1 |= (z[j+1] & one<<(60 + 1)) <<  1; t1 |= (z[j+1] & one<<(62 + 1)) <<  0; 

      a0[j2] = t0;
      a1[j2] = t1;
    }

    t0 = 0;
    t1 = 0;
    switch(Z->x->width - j) {
    case 2:
      t0 |= (z[j+1] & one<< 0) << 32; t0 |= (z[j+1] & one<< 2) << 31; t0 |= (z[j+1] & one<< 4) << 30; t0 |= (z[j+1] & one<< 6) << 29; 
      t0 |= (z[j+1] & one<< 8) << 28; t0 |= (z[j+1] & one<<10) << 27; t0 |= (z[j+1] & one<<12) << 26; t0 |= (z[j+1] & one<<14) << 25; 
      t0 |= (z[j+1] & one<<16) << 24; t0 |= (z[j+1] & one<<18) << 23; t0 |= (z[j+1] & one<<20) << 22; t0 |= (z[j+1] & one<<22) << 21; 
      t0 |= (z[j+1] & one<<24) << 20; t0 |= (z[j+1] & one<<26) << 19; t0 |= (z[j+1] & one<<28) << 18; t0 |= (z[j+1] & one<<30) << 17; 
      t0 |= (z[j+1] & one<<32) << 16; t0 |= (z[j+1] & one<<34) << 15; t0 |= (z[j+1] & one<<36) << 14; t0 |= (z[j+1] & one<<38) << 13; 
      t0 |= (z[j+1] & one<<40) << 12; t0 |= (z[j+1] & one<<42) << 11; t0 |= (z[j+1] & one<<44) << 10; t0 |= (z[j+1] & one<<46) <<  9; 
      t0 |= (z[j+1] & one<<48) <<  8; t0 |= (z[j+1] & one<<50) <<  7; t0 |= (z[j+1] & one<<52) <<  6; t0 |= (z[j+1] & one<<54) <<  5; 
      t0 |= (z[j+1] & one<<56) <<  4; t0 |= (z[j+1] & one<<58) <<  3; t0 |= (z[j+1] & one<<60) <<  2; t0 |= (z[j+1] & one<<62) <<  1;  

      t1 |= (z[j+1] & one<<( 0 + 1)) << 31; t1 |= (z[j+1] & one<<( 2 + 1)) << 30; t1 |= (z[j+1] & one<<( 4 + 1)) << 29; t1 |= (z[j+1] & one<<( 6 + 1)) << 28; 
      t1 |= (z[j+1] & one<<( 8 + 1)) << 27; t1 |= (z[j+1] & one<<(10 + 1)) << 26; t1 |= (z[j+1] & one<<(12 + 1)) << 25; t1 |= (z[j+1] & one<<(14 + 1)) << 24; 
      t1 |= (z[j+1] & one<<(16 + 1)) << 23; t1 |= (z[j+1] & one<<(18 + 1)) << 22; t1 |= (z[j+1] & one<<(20 + 1)) << 21; t1 |= (z[j+1] & one<<(22 + 1)) << 20; 
      t1 |= (z[j+1] & one<<(24 + 1)) << 19; t1 |= (z[j+1] & one<<(26 + 1)) << 18; t1 |= (z[j+1] & one<<(28 + 1)) << 17; t1 |= (z[j+1] & one<<(30 + 1)) << 16; 
      t1 |= (z[j+1] & one<<(32 + 1)) << 15; t1 |= (z[j+1] & one<<(34 + 1)) << 14; t1 |= (z[j+1] & one<<(36 + 1)) << 13; t1 |= (z[j+1] & one<<(38 + 1)) << 12; 
      t1 |= (z[j+1] & one<<(40 + 1)) << 11; t1 |= (z[j+1] & one<<(42 + 1)) << 10; t1 |= (z[j+1] & one<<(44 + 1)) <<  9; t1 |= (z[j+1] & one<<(46 + 1)) <<  8; 
      t1 |= (z[j+1] & one<<(48 + 1)) <<  7; t1 |= (z[j+1] & one<<(50 + 1)) <<  6; t1 |= (z[j+1] & one<<(52 + 1)) <<  5; t1 |= (z[j+1] & one<<(54 + 1)) <<  4; 
      t1 |= (z[j+1] & one<<(56 + 1)) <<  3; t1 |= (z[j+1] & one<<(58 + 1)) <<  2; t1 |= (z[j+1] & one<<(60 + 1)) <<  1; t1 |= (z[j+1] & one<<(62 + 1)) <<  0; 
    case 1:
      t0 |= (z[j+0] & one<< 0) >>  0; t0 |= (z[j+0] & one<< 2) >>  1; t0 |= (z[j+0] & one<< 4) >>  2; t0 |= (z[j+0] & one<< 6) >>  3; 
      t0 |= (z[j+0] & one<< 8) >>  4; t0 |= (z[j+0] & one<<10) >>  5; t0 |= (z[j+0] & one<<12) >>  6; t0 |= (z[j+0] & one<<14) >>  7; 
      t0 |= (z[j+0] & one<<16) >>  8; t0 |= (z[j+0] & one<<18) >>  9; t0 |= (z[j+0] & one<<20) >> 10; t0 |= (z[j+0] & one<<22) >> 11; 
      t0 |= (z[j+0] & one<<24) >> 12; t0 |= (z[j+0] & one<<26) >> 13; t0 |= (z[j+0] & one<<28) >> 14; t0 |= (z[j+0] & one<<30) >> 15; 
      t0 |= (z[j+0] & one<<32) >> 16; t0 |= (z[j+0] & one<<34) >> 17; t0 |= (z[j+0] & one<<36) >> 18; t0 |= (z[j+0] & one<<38) >> 19; 
      t0 |= (z[j+0] & one<<40) >> 20; t0 |= (z[j+0] & one<<42) >> 21; t0 |= (z[j+0] & one<<44) >> 22; t0 |= (z[j+0] & one<<46) >> 23; 
      t0 |= (z[j+0] & one<<48) >> 24; t0 |= (z[j+0] & one<<50) >> 25; t0 |= (z[j+0] & one<<52) >> 26; t0 |= (z[j+0] & one<<54) >> 27; 
      t0 |= (z[j+0] & one<<56) >> 28; t0 |= (z[j+0] & one<<58) >> 29; t0 |= (z[j+0] & one<<60) >> 30; t0 |= (z[j+0] & one<<62) >> 31; 
      
      t1 |= (z[j+0] & one<<( 0 + 1)) >>  1; t1 |= (z[j+0] & one<<( 2 + 1)) >>  2; t1 |= (z[j+0] & one<<( 4 + 1)) >>  3; t1 |= (z[j+0] & one<<( 6 + 1)) >>  4; 
      t1 |= (z[j+0] & one<<( 8 + 1)) >>  5; t1 |= (z[j+0] & one<<(10 + 1)) >>  6; t1 |= (z[j+0] & one<<(12 + 1)) >>  7; t1 |= (z[j+0] & one<<(14 + 1)) >>  8; 
      t1 |= (z[j+0] & one<<(16 + 1)) >>  9; t1 |= (z[j+0] & one<<(18 + 1)) >> 10; t1 |= (z[j+0] & one<<(20 + 1)) >> 11; t1 |= (z[j+0] & one<<(22 + 1)) >> 12; 
      t1 |= (z[j+0] & one<<(24 + 1)) >> 13; t1 |= (z[j+0] & one<<(26 + 1)) >> 14; t1 |= (z[j+0] & one<<(28 + 1)) >> 15; t1 |= (z[j+0] & one<<(30 + 1)) >> 16; 
      t1 |= (z[j+0] & one<<(32 + 1)) >> 17; t1 |= (z[j+0] & one<<(34 + 1)) >> 18; t1 |= (z[j+0] & one<<(36 + 1)) >> 19; t1 |= (z[j+0] & one<<(38 + 1)) >> 20; 
      t1 |= (z[j+0] & one<<(40 + 1)) >> 21; t1 |= (z[j+0] & one<<(42 + 1)) >> 22; t1 |= (z[j+0] & one<<(44 + 1)) >> 23; t1 |= (z[j+0] & one<<(46 + 1)) >> 24; 
      t1 |= (z[j+0] & one<<(48 + 1)) >> 25; t1 |= (z[j+0] & one<<(50 + 1)) >> 26; t1 |= (z[j+0] & one<<(52 + 1)) >> 27; t1 |= (z[j+0] & one<<(54 + 1)) >> 28; 
      t1 |= (z[j+0] & one<<(56 + 1)) >> 29; t1 |= (z[j+0] & one<<(58 + 1)) >> 30; t1 |= (z[j+0] & one<<(60 + 1)) >> 31; t1 |= (z[j+0] & one<<(62 + 1)) >> 32; 
      break;
    default:
      m4ri_die("impossible");
    }
    a0[j2] &= ~bitmask_end;
    a0[j2] |= t0 & bitmask_end;
    a1[j2] &= ~bitmask_end;
    a1[j2] |= t1 & bitmask_end;

    /* fix bits before offset */
    a0[0] = (a0[0] & bitmask_begin) | (a0_fix & ~bitmask_begin);
    a1[0] = (a1[0] & bitmask_begin) | (a1_fix & ~bitmask_begin);
  }
  
  return A;
}

mzed_t *_mzed_cling2(mzed_t *A, const mzd_slice_t *Z) {
  size_t j,j2 = 0;
  register word aw0;
  register word aw1;

  const word bitmask_begin = __M4RI_RIGHT_BITMASK(m4ri_radix - A->x->offset%m4ri_radix);
  const word bitmask_end = __M4RI_LEFT_BITMASK((A->x->offset + A->x->ncols) % m4ri_radix);

  /** A0 **/
  for(size_t i=0; i<A->nrows; i++) {
    word *z0 = Z->x[0]->rows[i];
    word *a  = A->x->rows[i];
    const word a_fix = a[0];

    for(j=0, j2=0; j+2 < A->x->width; j+=2, j2++) {
      if (!z0[j2] )
        continue;
      aw0 = aw1 = 0;
      aw0 |= (z0[j2] & m4ri_one<< 0) <<  0;    aw1 |= (z0[j2] & m4ri_one<<32) >> 32;
      aw0 |= (z0[j2] & m4ri_one<< 1) <<  1;    aw1 |= (z0[j2] & m4ri_one<<33) >> 31;
      aw0 |= (z0[j2] & m4ri_one<< 2) <<  2;    aw1 |= (z0[j2] & m4ri_one<<34) >> 30;
      aw0 |= (z0[j2] & m4ri_one<< 3) <<  3;    aw1 |= (z0[j2] & m4ri_one<<35) >> 29;
      aw0 |= (z0[j2] & m4ri_one<< 4) <<  4;    aw1 |= (z0[j2] & m4ri_one<<36) >> 28;
      aw0 |= (z0[j2] & m4ri_one<< 5) <<  5;    aw1 |= (z0[j2] & m4ri_one<<37) >> 27;
      aw0 |= (z0[j2] & m4ri_one<< 6) <<  6;    aw1 |= (z0[j2] & m4ri_one<<38) >> 26;
      aw0 |= (z0[j2] & m4ri_one<< 7) <<  7;    aw1 |= (z0[j2] & m4ri_one<<39) >> 25;
      aw0 |= (z0[j2] & m4ri_one<< 8) <<  8;    aw1 |= (z0[j2] & m4ri_one<<40) >> 24;
      aw0 |= (z0[j2] & m4ri_one<< 9) <<  9;    aw1 |= (z0[j2] & m4ri_one<<41) >> 23;
      aw0 |= (z0[j2] & m4ri_one<<10) << 10;    aw1 |= (z0[j2] & m4ri_one<<42) >> 22;
      aw0 |= (z0[j2] & m4ri_one<<11) << 11;    aw1 |= (z0[j2] & m4ri_one<<43) >> 21;
      aw0 |= (z0[j2] & m4ri_one<<12) << 12;    aw1 |= (z0[j2] & m4ri_one<<44) >> 20;
      aw0 |= (z0[j2] & m4ri_one<<13) << 13;    aw1 |= (z0[j2] & m4ri_one<<45) >> 19;
      aw0 |= (z0[j2] & m4ri_one<<14) << 14;    aw1 |= (z0[j2] & m4ri_one<<46) >> 18;
      aw0 |= (z0[j2] & m4ri_one<<15) << 15;    aw1 |= (z0[j2] & m4ri_one<<47) >> 17;
      aw0 |= (z0[j2] & m4ri_one<<16) << 16;    aw1 |= (z0[j2] & m4ri_one<<48) >> 16;
      aw0 |= (z0[j2] & m4ri_one<<17) << 17;    aw1 |= (z0[j2] & m4ri_one<<49) >> 15;
      aw0 |= (z0[j2] & m4ri_one<<18) << 18;    aw1 |= (z0[j2] & m4ri_one<<50) >> 14;
      aw0 |= (z0[j2] & m4ri_one<<19) << 19;    aw1 |= (z0[j2] & m4ri_one<<51) >> 13;
      aw0 |= (z0[j2] & m4ri_one<<20) << 20;    aw1 |= (z0[j2] & m4ri_one<<52) >> 12;
      aw0 |= (z0[j2] & m4ri_one<<21) << 21;    aw1 |= (z0[j2] & m4ri_one<<53) >> 11;
      aw0 |= (z0[j2] & m4ri_one<<22) << 22;    aw1 |= (z0[j2] & m4ri_one<<54) >> 10;
      aw0 |= (z0[j2] & m4ri_one<<23) << 23;    aw1 |= (z0[j2] & m4ri_one<<55) >>  9;
      aw0 |= (z0[j2] & m4ri_one<<24) << 24;    aw1 |= (z0[j2] & m4ri_one<<56) >>  8;
      aw0 |= (z0[j2] & m4ri_one<<25) << 25;    aw1 |= (z0[j2] & m4ri_one<<57) >>  7;
      aw0 |= (z0[j2] & m4ri_one<<26) << 26;    aw1 |= (z0[j2] & m4ri_one<<58) >>  6;
      aw0 |= (z0[j2] & m4ri_one<<27) << 27;    aw1 |= (z0[j2] & m4ri_one<<59) >>  5;
      aw0 |= (z0[j2] & m4ri_one<<28) << 28;    aw1 |= (z0[j2] & m4ri_one<<60) >>  4;
      aw0 |= (z0[j2] & m4ri_one<<29) << 29;    aw1 |= (z0[j2] & m4ri_one<<61) >>  3;
      aw0 |= (z0[j2] & m4ri_one<<30) << 30;    aw1 |= (z0[j2] & m4ri_one<<62) >>  2;
      aw0 |= (z0[j2] & m4ri_one<<31) << 31;    aw1 |= (z0[j2] & m4ri_one<<63) >>  1;
      a[j+0] = aw0;
      a[j+1] = aw1;
    }

    if(j+2 == A->x->width) {  /* we have to deal with two words */
      a[j+0] = 0;
      a[j+0] |= (z0[j2] & m4ri_one<<  0 ) <<  0;
      a[j+0] |= (z0[j2] & m4ri_one<<  1 ) <<  1;
      a[j+0] |= (z0[j2] & m4ri_one<<  2 ) <<  2;
      a[j+0] |= (z0[j2] & m4ri_one<<  3 ) <<  3;
      a[j+0] |= (z0[j2] & m4ri_one<<  4 ) <<  4;
      a[j+0] |= (z0[j2] & m4ri_one<<  5 ) <<  5;
      a[j+0] |= (z0[j2] & m4ri_one<<  6 ) <<  6;
      a[j+0] |= (z0[j2] & m4ri_one<<  7 ) <<  7;
      a[j+0] |= (z0[j2] & m4ri_one<<  8 ) <<  8;
      a[j+0] |= (z0[j2] & m4ri_one<<  9 ) <<  9;
      a[j+0] |= (z0[j2] & m4ri_one<< 10 ) << 10;
      a[j+0] |= (z0[j2] & m4ri_one<< 11 ) << 11;
      a[j+0] |= (z0[j2] & m4ri_one<< 12 ) << 12;
      a[j+0] |= (z0[j2] & m4ri_one<< 13 ) << 13;
      a[j+0] |= (z0[j2] & m4ri_one<< 14 ) << 14;
      a[j+0] |= (z0[j2] & m4ri_one<< 15 ) << 15;
      a[j+0] |= (z0[j2] & m4ri_one<< 16 ) << 16;
      a[j+0] |= (z0[j2] & m4ri_one<< 17 ) << 17;
      a[j+0] |= (z0[j2] & m4ri_one<< 18 ) << 18;
      a[j+0] |= (z0[j2] & m4ri_one<< 19 ) << 19;
      a[j+0] |= (z0[j2] & m4ri_one<< 20 ) << 20;
      a[j+0] |= (z0[j2] & m4ri_one<< 21 ) << 21;
      a[j+0] |= (z0[j2] & m4ri_one<< 22 ) << 22;
      a[j+0] |= (z0[j2] & m4ri_one<< 23 ) << 23;
      a[j+0] |= (z0[j2] & m4ri_one<< 24 ) << 24;
      a[j+0] |= (z0[j2] & m4ri_one<< 25 ) << 25;
      a[j+0] |= (z0[j2] & m4ri_one<< 26 ) << 26;
      a[j+0] |= (z0[j2] & m4ri_one<< 27 ) << 27;
      a[j+0] |= (z0[j2] & m4ri_one<< 28 ) << 28;
      a[j+0] |= (z0[j2] & m4ri_one<< 29 ) << 29;
      a[j+0] |= (z0[j2] & m4ri_one<< 30 ) << 30;
      a[j+0] |= (z0[j2] & m4ri_one<< 31 ) << 31;

      a[j+1] &= ~bitmask_end;
      switch((A->x->offset+A->x->ncols) % m4ri_radix) {
      case  0:      a[j+1] |= (z0[j2] & m4ri_one<<63) >>  1;
      case 62:      a[j+1] |= (z0[j2] & m4ri_one<<62) >>  2;
      case 60:      a[j+1] |= (z0[j2] & m4ri_one<<61) >>  3;
      case 58:      a[j+1] |= (z0[j2] & m4ri_one<<60) >>  4;
      case 56:      a[j+1] |= (z0[j2] & m4ri_one<<59) >>  5;
      case 54:      a[j+1] |= (z0[j2] & m4ri_one<<58) >>  6;
      case 52:      a[j+1] |= (z0[j2] & m4ri_one<<57) >>  7;
      case 50:      a[j+1] |= (z0[j2] & m4ri_one<<56) >>  8;
      case 48:      a[j+1] |= (z0[j2] & m4ri_one<<55) >>  9;
      case 46:      a[j+1] |= (z0[j2] & m4ri_one<<54) >> 10;
      case 44:      a[j+1] |= (z0[j2] & m4ri_one<<53) >> 11;
      case 42:      a[j+1] |= (z0[j2] & m4ri_one<<52) >> 12;
      case 40:      a[j+1] |= (z0[j2] & m4ri_one<<51) >> 13;
      case 38:      a[j+1] |= (z0[j2] & m4ri_one<<50) >> 14;
      case 36:      a[j+1] |= (z0[j2] & m4ri_one<<49) >> 15;
      case 34:      a[j+1] |= (z0[j2] & m4ri_one<<48) >> 16;
      case 32:      a[j+1] |= (z0[j2] & m4ri_one<<47) >> 17;
      case 30:      a[j+1] |= (z0[j2] & m4ri_one<<46) >> 18;
      case 28:      a[j+1] |= (z0[j2] & m4ri_one<<45) >> 19;
      case 26:      a[j+1] |= (z0[j2] & m4ri_one<<44) >> 20;
      case 24:      a[j+1] |= (z0[j2] & m4ri_one<<43) >> 21;
      case 22:      a[j+1] |= (z0[j2] & m4ri_one<<42) >> 22;
      case 20:      a[j+1] |= (z0[j2] & m4ri_one<<41) >> 23;
      case 18:      a[j+1] |= (z0[j2] & m4ri_one<<40) >> 24;
      case 16:      a[j+1] |= (z0[j2] & m4ri_one<<39) >> 25;
      case 14:      a[j+1] |= (z0[j2] & m4ri_one<<38) >> 26;
      case 12:      a[j+1] |= (z0[j2] & m4ri_one<<37) >> 27;
      case 10:      a[j+1] |= (z0[j2] & m4ri_one<<36) >> 28;
      case  8:      a[j+1] |= (z0[j2] & m4ri_one<<35) >> 29;
      case  6:      a[j+1] |= (z0[j2] & m4ri_one<<34) >> 30;
      case  4:      a[j+1] |= (z0[j2] & m4ri_one<<33) >> 31;
      case  2:      a[j+1] |= (z0[j2] & m4ri_one<<32) >> 32;
      }
    
    } else {  /* only one word */
      a[j+0] &= ~bitmask_end;
      switch((A->x->offset+A->x->ncols) % m4ri_radix) {
      case  0:      a[j+0] |= (z0[j2] & m4ri_one<<31) << 31;
      case 62:      a[j+0] |= (z0[j2] & m4ri_one<<30) << 30;
      case 60:      a[j+0] |= (z0[j2] & m4ri_one<<29) << 29;
      case 58:      a[j+0] |= (z0[j2] & m4ri_one<<28) << 28;
      case 56:      a[j+0] |= (z0[j2] & m4ri_one<<27) << 27;
      case 54:      a[j+0] |= (z0[j2] & m4ri_one<<26) << 26;
      case 52:      a[j+0] |= (z0[j2] & m4ri_one<<25) << 25;
      case 50:      a[j+0] |= (z0[j2] & m4ri_one<<24) << 24;
      case 48:      a[j+0] |= (z0[j2] & m4ri_one<<23) << 23;
      case 46:      a[j+0] |= (z0[j2] & m4ri_one<<22) << 22;
      case 44:      a[j+0] |= (z0[j2] & m4ri_one<<21) << 21;
      case 42:      a[j+0] |= (z0[j2] & m4ri_one<<20) << 20;
      case 40:      a[j+0] |= (z0[j2] & m4ri_one<<19) << 19;
      case 38:      a[j+0] |= (z0[j2] & m4ri_one<<18) << 18;
      case 36:      a[j+0] |= (z0[j2] & m4ri_one<<17) << 17;
      case 34:      a[j+0] |= (z0[j2] & m4ri_one<<16) << 16;
      case 32:      a[j+0] |= (z0[j2] & m4ri_one<<15) << 15;
      case 30:      a[j+0] |= (z0[j2] & m4ri_one<<14) << 14;
      case 28:      a[j+0] |= (z0[j2] & m4ri_one<<13) << 13;
      case 26:      a[j+0] |= (z0[j2] & m4ri_one<<12) << 12;
      case 24:      a[j+0] |= (z0[j2] & m4ri_one<<11) << 11;
      case 22:      a[j+0] |= (z0[j2] & m4ri_one<<10) << 10;
      case 20:      a[j+0] |= (z0[j2] & m4ri_one<< 9) <<  9;
      case 18:      a[j+0] |= (z0[j2] & m4ri_one<< 8) <<  8;
      case 16:      a[j+0] |= (z0[j2] & m4ri_one<< 7) <<  7;
      case 14:      a[j+0] |= (z0[j2] & m4ri_one<< 6) <<  6;
      case 12:      a[j+0] |= (z0[j2] & m4ri_one<< 5) <<  5;
      case 10:      a[j+0] |= (z0[j2] & m4ri_one<< 4) <<  4;
      case  8:      a[j+0] |= (z0[j2] & m4ri_one<< 3) <<  3;
      case  6:      a[j+0] |= (z0[j2] & m4ri_one<< 2) <<  2;
      case  4:      a[j+0] |= (z0[j2] & m4ri_one<< 1) <<  1;
      case  2:      a[j+0] |= (z0[j2] & m4ri_one<< 0) <<  0;
      }
    }
    ;
    a[0] = (a[0] & bitmask_begin) | (a_fix & ~bitmask_begin);
  }

  /** A1 **/
  for(size_t i=0; i<A->nrows; i++) {
    word *z1 = Z->x[1]->rows[i];
    word *a  = A->x->rows[i];    
    const word a_fix = a[0];

    for(j=0, j2=0; j+2 < A->x->width; j+=2, j2++) {
      if (!z1[j2] )
        continue;
      aw0 = a[j+0]; /** we wrote A0 already **/
      aw1 = a[j+1];
      aw0 |= (z1[j2] & m4ri_one<< 0) <<  1;    aw1 |= (z1[j2] & m4ri_one<<32) >> 31;
      aw0 |= (z1[j2] & m4ri_one<< 1) <<  2;    aw1 |= (z1[j2] & m4ri_one<<33) >> 30;
      aw0 |= (z1[j2] & m4ri_one<< 2) <<  3;    aw1 |= (z1[j2] & m4ri_one<<34) >> 29;
      aw0 |= (z1[j2] & m4ri_one<< 3) <<  4;    aw1 |= (z1[j2] & m4ri_one<<35) >> 28;
      aw0 |= (z1[j2] & m4ri_one<< 4) <<  5;    aw1 |= (z1[j2] & m4ri_one<<36) >> 27;
      aw0 |= (z1[j2] & m4ri_one<< 5) <<  6;    aw1 |= (z1[j2] & m4ri_one<<37) >> 26;
      aw0 |= (z1[j2] & m4ri_one<< 6) <<  7;    aw1 |= (z1[j2] & m4ri_one<<38) >> 25;
      aw0 |= (z1[j2] & m4ri_one<< 7) <<  8;    aw1 |= (z1[j2] & m4ri_one<<39) >> 24;
      aw0 |= (z1[j2] & m4ri_one<< 8) <<  9;    aw1 |= (z1[j2] & m4ri_one<<40) >> 23;
      aw0 |= (z1[j2] & m4ri_one<< 9) << 10;    aw1 |= (z1[j2] & m4ri_one<<41) >> 22;
      aw0 |= (z1[j2] & m4ri_one<<10) << 11;    aw1 |= (z1[j2] & m4ri_one<<42) >> 21;
      aw0 |= (z1[j2] & m4ri_one<<11) << 12;    aw1 |= (z1[j2] & m4ri_one<<43) >> 20;
      aw0 |= (z1[j2] & m4ri_one<<12) << 13;    aw1 |= (z1[j2] & m4ri_one<<44) >> 19;
      aw0 |= (z1[j2] & m4ri_one<<13) << 14;    aw1 |= (z1[j2] & m4ri_one<<45) >> 18;
      aw0 |= (z1[j2] & m4ri_one<<14) << 15;    aw1 |= (z1[j2] & m4ri_one<<46) >> 17;
      aw0 |= (z1[j2] & m4ri_one<<15) << 16;    aw1 |= (z1[j2] & m4ri_one<<47) >> 16;
      aw0 |= (z1[j2] & m4ri_one<<16) << 17;    aw1 |= (z1[j2] & m4ri_one<<48) >> 15;
      aw0 |= (z1[j2] & m4ri_one<<17) << 18;    aw1 |= (z1[j2] & m4ri_one<<49) >> 14;
      aw0 |= (z1[j2] & m4ri_one<<18) << 19;    aw1 |= (z1[j2] & m4ri_one<<50) >> 13;
      aw0 |= (z1[j2] & m4ri_one<<19) << 20;    aw1 |= (z1[j2] & m4ri_one<<51) >> 12;
      aw0 |= (z1[j2] & m4ri_one<<20) << 21;    aw1 |= (z1[j2] & m4ri_one<<52) >> 11;
      aw0 |= (z1[j2] & m4ri_one<<21) << 22;    aw1 |= (z1[j2] & m4ri_one<<53) >> 10;
      aw0 |= (z1[j2] & m4ri_one<<22) << 23;    aw1 |= (z1[j2] & m4ri_one<<54) >>  9;
      aw0 |= (z1[j2] & m4ri_one<<23) << 24;    aw1 |= (z1[j2] & m4ri_one<<55) >>  8;
      aw0 |= (z1[j2] & m4ri_one<<24) << 25;    aw1 |= (z1[j2] & m4ri_one<<56) >>  7;
      aw0 |= (z1[j2] & m4ri_one<<25) << 26;    aw1 |= (z1[j2] & m4ri_one<<57) >>  6;
      aw0 |= (z1[j2] & m4ri_one<<26) << 27;    aw1 |= (z1[j2] & m4ri_one<<58) >>  5;
      aw0 |= (z1[j2] & m4ri_one<<27) << 28;    aw1 |= (z1[j2] & m4ri_one<<59) >>  4;
      aw0 |= (z1[j2] & m4ri_one<<28) << 29;    aw1 |= (z1[j2] & m4ri_one<<60) >>  3;
      aw0 |= (z1[j2] & m4ri_one<<29) << 30;    aw1 |= (z1[j2] & m4ri_one<<61) >>  2;
      aw0 |= (z1[j2] & m4ri_one<<30) << 31;    aw1 |= (z1[j2] & m4ri_one<<62) >>  1;
      aw0 |= (z1[j2] & m4ri_one<<31) << 32;    aw1 |= (z1[j2] & m4ri_one<<63) >>  0;
      a[j+0] = aw0;
      a[j+1] = aw1;
    }

    if(j+2 == A->x->width) {  /* we have to deal with two words */
      a[j+0] |= (z1[j2] & m4ri_one<< 0) <<  1;
      a[j+0] |= (z1[j2] & m4ri_one<< 1) <<  2;
      a[j+0] |= (z1[j2] & m4ri_one<< 2) <<  3;
      a[j+0] |= (z1[j2] & m4ri_one<< 3) <<  4;
      a[j+0] |= (z1[j2] & m4ri_one<< 4) <<  5;
      a[j+0] |= (z1[j2] & m4ri_one<< 5) <<  6;
      a[j+0] |= (z1[j2] & m4ri_one<< 6) <<  7;
      a[j+0] |= (z1[j2] & m4ri_one<< 7) <<  8;
      a[j+0] |= (z1[j2] & m4ri_one<< 8) <<  9;
      a[j+0] |= (z1[j2] & m4ri_one<< 9) << 10;
      a[j+0] |= (z1[j2] & m4ri_one<<10) << 11;
      a[j+0] |= (z1[j2] & m4ri_one<<11) << 12;
      a[j+0] |= (z1[j2] & m4ri_one<<12) << 13;
      a[j+0] |= (z1[j2] & m4ri_one<<13) << 14;
      a[j+0] |= (z1[j2] & m4ri_one<<14) << 15;
      a[j+0] |= (z1[j2] & m4ri_one<<15) << 16;
      a[j+0] |= (z1[j2] & m4ri_one<<16) << 17;
      a[j+0] |= (z1[j2] & m4ri_one<<17) << 18;
      a[j+0] |= (z1[j2] & m4ri_one<<18) << 19;
      a[j+0] |= (z1[j2] & m4ri_one<<19) << 20;
      a[j+0] |= (z1[j2] & m4ri_one<<20) << 21;
      a[j+0] |= (z1[j2] & m4ri_one<<21) << 22;
      a[j+0] |= (z1[j2] & m4ri_one<<22) << 23;
      a[j+0] |= (z1[j2] & m4ri_one<<23) << 24;
      a[j+0] |= (z1[j2] & m4ri_one<<24) << 25;
      a[j+0] |= (z1[j2] & m4ri_one<<25) << 26;
      a[j+0] |= (z1[j2] & m4ri_one<<26) << 27;
      a[j+0] |= (z1[j2] & m4ri_one<<27) << 28;
      a[j+0] |= (z1[j2] & m4ri_one<<28) << 29;
      a[j+0] |= (z1[j2] & m4ri_one<<29) << 30;
      a[j+0] |= (z1[j2] & m4ri_one<<30) << 31;
      a[j+0] |= (z1[j2] & m4ri_one<<31) << 32;
 
      switch((A->x->offset+A->x->ncols) % m4ri_radix) {
      case  0:      a[j+1] |= (z1[j2] & m4ri_one<<63) >>  0;
      case 62:      a[j+1] |= (z1[j2] & m4ri_one<<62) >>  1;
      case 60:      a[j+1] |= (z1[j2] & m4ri_one<<61) >>  2;
      case 58:      a[j+1] |= (z1[j2] & m4ri_one<<60) >>  3;
      case 56:      a[j+1] |= (z1[j2] & m4ri_one<<59) >>  4;
      case 54:      a[j+1] |= (z1[j2] & m4ri_one<<58) >>  5;
      case 52:      a[j+1] |= (z1[j2] & m4ri_one<<57) >>  6;
      case 50:      a[j+1] |= (z1[j2] & m4ri_one<<56) >>  7;
      case 48:      a[j+1] |= (z1[j2] & m4ri_one<<55) >>  8;
      case 46:      a[j+1] |= (z1[j2] & m4ri_one<<54) >>  9;
      case 44:      a[j+1] |= (z1[j2] & m4ri_one<<53) >> 10;
      case 42:      a[j+1] |= (z1[j2] & m4ri_one<<52) >> 11;
      case 40:      a[j+1] |= (z1[j2] & m4ri_one<<51) >> 12;
      case 38:      a[j+1] |= (z1[j2] & m4ri_one<<50) >> 13;
      case 36:      a[j+1] |= (z1[j2] & m4ri_one<<49) >> 14;
      case 34:      a[j+1] |= (z1[j2] & m4ri_one<<48) >> 15;
      case 32:      a[j+1] |= (z1[j2] & m4ri_one<<47) >> 16;
      case 30:      a[j+1] |= (z1[j2] & m4ri_one<<46) >> 17;
      case 28:      a[j+1] |= (z1[j2] & m4ri_one<<45) >> 18;
      case 26:      a[j+1] |= (z1[j2] & m4ri_one<<44) >> 19;
      case 24:      a[j+1] |= (z1[j2] & m4ri_one<<43) >> 20;
      case 22:      a[j+1] |= (z1[j2] & m4ri_one<<42) >> 21;
      case 20:      a[j+1] |= (z1[j2] & m4ri_one<<41) >> 22;
      case 18:      a[j+1] |= (z1[j2] & m4ri_one<<40) >> 23;
      case 16:      a[j+1] |= (z1[j2] & m4ri_one<<39) >> 24;
      case 14:      a[j+1] |= (z1[j2] & m4ri_one<<38) >> 25;
      case 12:      a[j+1] |= (z1[j2] & m4ri_one<<37) >> 26;
      case 10:      a[j+1] |= (z1[j2] & m4ri_one<<36) >> 27;
      case  8:      a[j+1] |= (z1[j2] & m4ri_one<<35) >> 28;
      case  6:      a[j+1] |= (z1[j2] & m4ri_one<<34) >> 29;
      case  4:      a[j+1] |= (z1[j2] & m4ri_one<<33) >> 30;
      case  2:      a[j+1] |= (z1[j2] & m4ri_one<<32) >> 31;
      }

    } else { /* only one word */
      switch((A->x->offset+A->x->ncols) % m4ri_radix) {
      case  0:      a[j+0] |= (z1[j2] & m4ri_one<<31) << 32;
      case 62:      a[j+0] |= (z1[j2] & m4ri_one<<30) << 31;
      case 60:      a[j+0] |= (z1[j2] & m4ri_one<<29) << 30;
      case 58:      a[j+0] |= (z1[j2] & m4ri_one<<28) << 29;
      case 56:      a[j+0] |= (z1[j2] & m4ri_one<<27) << 28;
      case 54:      a[j+0] |= (z1[j2] & m4ri_one<<26) << 27;
      case 52:      a[j+0] |= (z1[j2] & m4ri_one<<25) << 26;
      case 50:      a[j+0] |= (z1[j2] & m4ri_one<<24) << 25;
      case 48:      a[j+0] |= (z1[j2] & m4ri_one<<23) << 24;
      case 46:      a[j+0] |= (z1[j2] & m4ri_one<<22) << 23;
      case 44:      a[j+0] |= (z1[j2] & m4ri_one<<21) << 22;
      case 42:      a[j+0] |= (z1[j2] & m4ri_one<<20) << 21;
      case 40:      a[j+0] |= (z1[j2] & m4ri_one<<19) << 20;
      case 38:      a[j+0] |= (z1[j2] & m4ri_one<<18) << 19;
      case 36:      a[j+0] |= (z1[j2] & m4ri_one<<17) << 18;
      case 34:      a[j+0] |= (z1[j2] & m4ri_one<<16) << 17;
      case 32:      a[j+0] |= (z1[j2] & m4ri_one<<15) << 16;
      case 30:      a[j+0] |= (z1[j2] & m4ri_one<<14) << 15;
      case 28:      a[j+0] |= (z1[j2] & m4ri_one<<13) << 14;
      case 26:      a[j+0] |= (z1[j2] & m4ri_one<<12) << 13;
      case 24:      a[j+0] |= (z1[j2] & m4ri_one<<11) << 12;
      case 22:      a[j+0] |= (z1[j2] & m4ri_one<<10) << 11;
      case 20:      a[j+0] |= (z1[j2] & m4ri_one<< 9) << 10;
      case 18:      a[j+0] |= (z1[j2] & m4ri_one<< 8) <<  9;
      case 16:      a[j+0] |= (z1[j2] & m4ri_one<< 7) <<  8;
      case 14:      a[j+0] |= (z1[j2] & m4ri_one<< 6) <<  7;
      case 12:      a[j+0] |= (z1[j2] & m4ri_one<< 5) <<  6;
      case 10:      a[j+0] |= (z1[j2] & m4ri_one<< 4) <<  5;
      case  8:      a[j+0] |= (z1[j2] & m4ri_one<< 3) <<  4;
      case  6:      a[j+0] |= (z1[j2] & m4ri_one<< 2) <<  3;
      case  4:      a[j+0] |= (z1[j2] & m4ri_one<< 1) <<  2;
      case  2:      a[j+0] |= (z1[j2] & m4ri_one<< 0) <<  1;
      }
    }
    a[0] = (a[0] & bitmask_begin) | (a_fix & ~bitmask_begin);
  }
  return A;
}

mzd_slice_t *_mzed_slice4(mzd_slice_t *A, const mzed_t *Z) {
  assert(A && (A->depth == 3 || A->depth == 4));
  size_t j, j2 = 0;
  register word t0,t1,t2,t3 = 0;

  const word one = m4ri_one;

  const word bitmask_begin = __M4RI_RIGHT_BITMASK(m4ri_radix - A->x[0]->offset%m4ri_radix);
  const word bitmask_end = __M4RI_LEFT_BITMASK((A->x[0]->offset + A->ncols) % m4ri_radix);

  /* A0 */
  for(size_t i=0; i<A->nrows; i++) {
    word *a0 = A->x[0]->rows[i];
    word *a1 = A->x[1]->rows[i];
    word *a2 = A->x[2]->rows[i];
    const word const *z  = Z->x->rows[i];

    const word a0_fix = a0[0]; 
    const word a1_fix = a1[0]; 
    const word a2_fix = a2[0]; 

    /* bulk of work */
    for(j=0, j2=0; j+4 < Z->x->width; j+=4,j2++) {
      t0 = t1 = t2 = 0;

      t0 |= (z[j+3] & (one<< 0)) <<48; t0 |= (z[j+3] & (one<< 4)) <<45; t0 |= (z[j+3] & (one<< 8)) <<42; t0 |= (z[j+3] & (one<<12)) <<39;
      t0 |= (z[j+3] & (one<<16)) <<36; t0 |= (z[j+3] & (one<<20)) <<33; t0 |= (z[j+3] & (one<<24)) <<30; t0 |= (z[j+3] & (one<<28)) <<27;
      t0 |= (z[j+3] & (one<<32)) <<24; t0 |= (z[j+3] & (one<<36)) <<21; t0 |= (z[j+3] & (one<<40)) <<18; t0 |= (z[j+3] & (one<<44)) <<15;
      t0 |= (z[j+3] & (one<<48)) <<12; t0 |= (z[j+3] & (one<<52)) << 9; t0 |= (z[j+3] & (one<<56)) << 6; t0 |= (z[j+3] & (one<<60)) << 3;

      t0 |= (z[j+2] & (one<< 0)) <<32; t0 |= (z[j+2] & (one<< 4)) <<29; t0 |= (z[j+2] & (one<< 8)) <<26; t0 |= (z[j+2] & (one<<12)) <<23;
      t0 |= (z[j+2] & (one<<16)) <<20; t0 |= (z[j+2] & (one<<20)) <<17; t0 |= (z[j+2] & (one<<24)) <<14; t0 |= (z[j+2] & (one<<28)) <<11;
      t0 |= (z[j+2] & (one<<32)) << 8; t0 |= (z[j+2] & (one<<36)) << 5; t0 |= (z[j+2] & (one<<40)) << 2; t0 |= (z[j+2] & (one<<44)) >> 1;
      t0 |= (z[j+2] & (one<<48)) >> 4; t0 |= (z[j+2] & (one<<52)) >> 7; t0 |= (z[j+2] & (one<<56)) >>10; t0 |= (z[j+2] & (one<<60)) >>13;

      t0 |= (z[j+1] & (one<< 0)) <<16; t0 |= (z[j+1] & (one<< 4)) <<13; t0 |= (z[j+1] & (one<< 8)) <<10; t0 |= (z[j+1] & (one<<12)) << 7;
      t0 |= (z[j+1] & (one<<16)) << 4; t0 |= (z[j+1] & (one<<20)) << 1; t0 |= (z[j+1] & (one<<24)) >> 2; t0 |= (z[j+1] & (one<<28)) >> 5;
      t0 |= (z[j+1] & (one<<32)) >> 8; t0 |= (z[j+1] & (one<<36)) >>11; t0 |= (z[j+1] & (one<<40)) >>14; t0 |= (z[j+1] & (one<<44)) >>17;
      t0 |= (z[j+1] & (one<<48)) >>20; t0 |= (z[j+1] & (one<<52)) >>23; t0 |= (z[j+1] & (one<<56)) >>26; t0 |= (z[j+1] & (one<<60)) >>29;

      t0 |= (z[j+0] & (one<< 0)) >> 0; t0 |= (z[j+0] & (one<< 4)) >> 3; t0 |= (z[j+0] & (one<< 8)) >> 6; t0 |= (z[j+0] & (one<<12)) >> 9;
      t0 |= (z[j+0] & (one<<16)) >>12; t0 |= (z[j+0] & (one<<20)) >>15; t0 |= (z[j+0] & (one<<24)) >>18; t0 |= (z[j+0] & (one<<28)) >>21;
      t0 |= (z[j+0] & (one<<32)) >>24; t0 |= (z[j+0] & (one<<36)) >>27; t0 |= (z[j+0] & (one<<40)) >>30; t0 |= (z[j+0] & (one<<44)) >>33;
      t0 |= (z[j+0] & (one<<48)) >>36; t0 |= (z[j+0] & (one<<52)) >>39; t0 |= (z[j+0] & (one<<56)) >>42; t0 |= (z[j+0] & (one<<60)) >>45;

      t1 |= (z[j+3] & (one<< 1)) <<47; t1 |= (z[j+3] & (one<< 5)) <<44; t1 |= (z[j+3] & (one<< 9)) <<41; t1 |= (z[j+3] & (one<<13)) <<38;
      t1 |= (z[j+3] & (one<<17)) <<35; t1 |= (z[j+3] & (one<<21)) <<32; t1 |= (z[j+3] & (one<<25)) <<29; t1 |= (z[j+3] & (one<<29)) <<26;
      t1 |= (z[j+3] & (one<<33)) <<23; t1 |= (z[j+3] & (one<<37)) <<20; t1 |= (z[j+3] & (one<<41)) <<17; t1 |= (z[j+3] & (one<<45)) <<14;
      t1 |= (z[j+3] & (one<<49)) <<11; t1 |= (z[j+3] & (one<<53)) << 8; t1 |= (z[j+3] & (one<<57)) << 5; t1 |= (z[j+3] & (one<<61)) << 2;

      t1 |= (z[j+2] & (one<< 1)) <<31; t1 |= (z[j+2] & (one<< 5)) <<28; t1 |= (z[j+2] & (one<< 9)) <<25; t1 |= (z[j+2] & (one<<13)) <<22;
      t1 |= (z[j+2] & (one<<17)) <<19; t1 |= (z[j+2] & (one<<21)) <<16; t1 |= (z[j+2] & (one<<25)) <<13; t1 |= (z[j+2] & (one<<29)) <<10;
      t1 |= (z[j+2] & (one<<33)) << 7; t1 |= (z[j+2] & (one<<37)) << 4; t1 |= (z[j+2] & (one<<41)) << 1; t1 |= (z[j+2] & (one<<45)) >> 2;
      t1 |= (z[j+2] & (one<<49)) >> 5; t1 |= (z[j+2] & (one<<53)) >> 8; t1 |= (z[j+2] & (one<<57)) >>11; t1 |= (z[j+2] & (one<<61)) >>14;

      t1 |= (z[j+1] & (one<< 1)) <<15; t1 |= (z[j+1] & (one<< 5)) <<12; t1 |= (z[j+1] & (one<< 9)) << 9; t1 |= (z[j+1] & (one<<13)) << 6;
      t1 |= (z[j+1] & (one<<17)) << 3; t1 |= (z[j+1] & (one<<21)) << 0; t1 |= (z[j+1] & (one<<25)) >> 3; t1 |= (z[j+1] & (one<<29)) >> 6;
      t1 |= (z[j+1] & (one<<33)) >> 9; t1 |= (z[j+1] & (one<<37)) >>12; t1 |= (z[j+1] & (one<<41)) >>15; t1 |= (z[j+1] & (one<<45)) >>18;
      t1 |= (z[j+1] & (one<<49)) >>21; t1 |= (z[j+1] & (one<<53)) >>24; t1 |= (z[j+1] & (one<<57)) >>27; t1 |= (z[j+1] & (one<<61)) >>30;

      t1 |= (z[j+0] & (one<< 1)) >> 1; t1 |= (z[j+0] & (one<< 5)) >> 4; t1 |= (z[j+0] & (one<< 9)) >> 7; t1 |= (z[j+0] & (one<<13)) >>10;
      t1 |= (z[j+0] & (one<<17)) >>13; t1 |= (z[j+0] & (one<<21)) >>16; t1 |= (z[j+0] & (one<<25)) >>19; t1 |= (z[j+0] & (one<<29)) >>22;
      t1 |= (z[j+0] & (one<<33)) >>25; t1 |= (z[j+0] & (one<<37)) >>28; t1 |= (z[j+0] & (one<<41)) >>31; t1 |= (z[j+0] & (one<<45)) >>34;
      t1 |= (z[j+0] & (one<<49)) >>37; t1 |= (z[j+0] & (one<<53)) >>40; t1 |= (z[j+0] & (one<<57)) >>43; t1 |= (z[j+0] & (one<<61)) >>46;

      t2 |= (z[j+3] & (one<< 2)) <<46; t2 |= (z[j+3] & (one<< 6)) <<43; t2 |= (z[j+3] & (one<<10)) <<40; t2 |= (z[j+3] & (one<<14)) <<37;
      t2 |= (z[j+3] & (one<<18)) <<34; t2 |= (z[j+3] & (one<<22)) <<31; t2 |= (z[j+3] & (one<<26)) <<28; t2 |= (z[j+3] & (one<<30)) <<25;
      t2 |= (z[j+3] & (one<<34)) <<22; t2 |= (z[j+3] & (one<<38)) <<19; t2 |= (z[j+3] & (one<<42)) <<16; t2 |= (z[j+3] & (one<<46)) <<13;
      t2 |= (z[j+3] & (one<<50)) <<10; t2 |= (z[j+3] & (one<<54)) << 7; t2 |= (z[j+3] & (one<<58)) << 4; t2 |= (z[j+3] & (one<<62)) << 1;
                                      
      t2 |= (z[j+2] & (one<< 2)) <<30; t2 |= (z[j+2] & (one<< 6)) <<27; t2 |= (z[j+2] & (one<<10)) <<24; t2 |= (z[j+2] & (one<<14)) <<21;
      t2 |= (z[j+2] & (one<<18)) <<18; t2 |= (z[j+2] & (one<<22)) <<15; t2 |= (z[j+2] & (one<<26)) <<12; t2 |= (z[j+2] & (one<<30)) << 9;
      t2 |= (z[j+2] & (one<<34)) << 6; t2 |= (z[j+2] & (one<<38)) << 3; t2 |= (z[j+2] & (one<<42)) << 0; t2 |= (z[j+2] & (one<<46)) >> 3;
      t2 |= (z[j+2] & (one<<50)) >> 6; t2 |= (z[j+2] & (one<<54)) >> 9; t2 |= (z[j+2] & (one<<58)) >>12; t2 |= (z[j+2] & (one<<62)) >>15;
                                      
      t2 |= (z[j+1] & (one<< 2)) <<14; t2 |= (z[j+1] & (one<< 6)) <<11; t2 |= (z[j+1] & (one<<10)) << 8; t2 |= (z[j+1] & (one<<14)) << 5;
      t2 |= (z[j+1] & (one<<18)) << 2; t2 |= (z[j+1] & (one<<22)) >> 1; t2 |= (z[j+1] & (one<<26)) >> 4; t2 |= (z[j+1] & (one<<30)) >> 7;
      t2 |= (z[j+1] & (one<<34)) >>10; t2 |= (z[j+1] & (one<<38)) >>13; t2 |= (z[j+1] & (one<<42)) >>16; t2 |= (z[j+1] & (one<<46)) >>19;
      t2 |= (z[j+1] & (one<<50)) >>22; t2 |= (z[j+1] & (one<<54)) >>25; t2 |= (z[j+1] & (one<<58)) >>28; t2 |= (z[j+1] & (one<<62)) >>31;
                                      
      t2 |= (z[j+0] & (one<< 2)) >> 2; t2 |= (z[j+0] & (one<< 6)) >> 5; t2 |= (z[j+0] & (one<<10)) >> 8; t2 |= (z[j+0] & (one<<14)) >>11;
      t2 |= (z[j+0] & (one<<18)) >>14; t2 |= (z[j+0] & (one<<22)) >>17; t2 |= (z[j+0] & (one<<26)) >>20; t2 |= (z[j+0] & (one<<30)) >>23;
      t2 |= (z[j+0] & (one<<34)) >>26; t2 |= (z[j+0] & (one<<38)) >>29; t2 |= (z[j+0] & (one<<42)) >>32; t2 |= (z[j+0] & (one<<46)) >>35;
      t2 |= (z[j+0] & (one<<50)) >>38; t2 |= (z[j+0] & (one<<54)) >>41; t2 |= (z[j+0] & (one<<58)) >>44; t2 |= (z[j+0] & (one<<62)) >>47;

      a0[j2] = t0;
      a1[j2] = t1;
      a2[j2] = t2;

    }
    t0 = t1 = t2 = 0;
    switch(Z->x->width - j) {
    case 4:
      t0 |= (z[j+3] & (one<< 0)) <<48; t0 |= (z[j+3] & (one<< 4)) <<45; t0 |= (z[j+3] & (one<< 8)) <<42; t0 |= (z[j+3] & (one<<12)) <<39;
      t0 |= (z[j+3] & (one<<16)) <<36; t0 |= (z[j+3] & (one<<20)) <<33; t0 |= (z[j+3] & (one<<24)) <<30; t0 |= (z[j+3] & (one<<28)) <<27;
      t0 |= (z[j+3] & (one<<32)) <<24; t0 |= (z[j+3] & (one<<36)) <<21; t0 |= (z[j+3] & (one<<40)) <<18; t0 |= (z[j+3] & (one<<44)) <<15;
      t0 |= (z[j+3] & (one<<48)) <<12; t0 |= (z[j+3] & (one<<52)) << 9; t0 |= (z[j+3] & (one<<56)) << 6; t0 |= (z[j+3] & (one<<60)) << 3;

      t1 |= (z[j+3] & (one<< 1)) <<47; t1 |= (z[j+3] & (one<< 5)) <<44; t1 |= (z[j+3] & (one<< 9)) <<41; t1 |= (z[j+3] & (one<<13)) <<38;
      t1 |= (z[j+3] & (one<<17)) <<35; t1 |= (z[j+3] & (one<<21)) <<32; t1 |= (z[j+3] & (one<<25)) <<29; t1 |= (z[j+3] & (one<<29)) <<26;
      t1 |= (z[j+3] & (one<<33)) <<23; t1 |= (z[j+3] & (one<<37)) <<20; t1 |= (z[j+3] & (one<<41)) <<17; t1 |= (z[j+3] & (one<<45)) <<14;
      t1 |= (z[j+3] & (one<<49)) <<11; t1 |= (z[j+3] & (one<<53)) << 8; t1 |= (z[j+3] & (one<<57)) << 5; t1 |= (z[j+3] & (one<<61)) << 2;

      t2 |= (z[j+3] & (one<< 2)) <<46; t2 |= (z[j+3] & (one<< 6)) <<43; t2 |= (z[j+3] & (one<<10)) <<40; t2 |= (z[j+3] & (one<<14)) <<37;
      t2 |= (z[j+3] & (one<<18)) <<34; t2 |= (z[j+3] & (one<<22)) <<31; t2 |= (z[j+3] & (one<<26)) <<28; t2 |= (z[j+3] & (one<<30)) <<25;
      t2 |= (z[j+3] & (one<<34)) <<22; t2 |= (z[j+3] & (one<<38)) <<19; t2 |= (z[j+3] & (one<<42)) <<16; t2 |= (z[j+3] & (one<<46)) <<13;
      t2 |= (z[j+3] & (one<<50)) <<10; t2 |= (z[j+3] & (one<<54)) << 7; t2 |= (z[j+3] & (one<<58)) << 4; t2 |= (z[j+3] & (one<<62)) << 1;
    case 3:
      t0 |= (z[j+2] & (one<< 0)) <<32; t0 |= (z[j+2] & (one<< 4)) <<29; t0 |= (z[j+2] & (one<< 8)) <<26; t0 |= (z[j+2] & (one<<12)) <<23;
      t0 |= (z[j+2] & (one<<16)) <<20; t0 |= (z[j+2] & (one<<20)) <<17; t0 |= (z[j+2] & (one<<24)) <<14; t0 |= (z[j+2] & (one<<28)) <<11;
      t0 |= (z[j+2] & (one<<32)) << 8; t0 |= (z[j+2] & (one<<36)) << 5; t0 |= (z[j+2] & (one<<40)) << 2; t0 |= (z[j+2] & (one<<44)) >> 1;
      t0 |= (z[j+2] & (one<<48)) >> 4; t0 |= (z[j+2] & (one<<52)) >> 7; t0 |= (z[j+2] & (one<<56)) >>10; t0 |= (z[j+2] & (one<<60)) >>13;

      t1 |= (z[j+2] & (one<< 1)) <<31; t1 |= (z[j+2] & (one<< 5)) <<28; t1 |= (z[j+2] & (one<< 9)) <<25; t1 |= (z[j+2] & (one<<13)) <<22;
      t1 |= (z[j+2] & (one<<17)) <<19; t1 |= (z[j+2] & (one<<21)) <<16; t1 |= (z[j+2] & (one<<25)) <<13; t1 |= (z[j+2] & (one<<29)) <<10;
      t1 |= (z[j+2] & (one<<33)) << 7; t1 |= (z[j+2] & (one<<37)) << 4; t1 |= (z[j+2] & (one<<41)) << 1; t1 |= (z[j+2] & (one<<45)) >> 2;
      t1 |= (z[j+2] & (one<<49)) >> 5; t1 |= (z[j+2] & (one<<53)) >> 8; t1 |= (z[j+2] & (one<<57)) >>11; t1 |= (z[j+2] & (one<<61)) >>14;

      t2 |= (z[j+2] & (one<< 2)) <<30; t2 |= (z[j+2] & (one<< 6)) <<27; t2 |= (z[j+2] & (one<<10)) <<24; t2 |= (z[j+2] & (one<<14)) <<21;
      t2 |= (z[j+2] & (one<<18)) <<18; t2 |= (z[j+2] & (one<<22)) <<15; t2 |= (z[j+2] & (one<<26)) <<12; t2 |= (z[j+2] & (one<<30)) << 9;
      t2 |= (z[j+2] & (one<<34)) << 6; t2 |= (z[j+2] & (one<<38)) << 3; t2 |= (z[j+2] & (one<<42)) << 0; t2 |= (z[j+2] & (one<<46)) >> 3;
      t2 |= (z[j+2] & (one<<50)) >> 6; t2 |= (z[j+2] & (one<<54)) >> 9; t2 |= (z[j+2] & (one<<58)) >>12; t2 |= (z[j+2] & (one<<62)) >>15;
    case 2:
      t0 |= (z[j+1] & (one<< 0)) <<16; t0 |= (z[j+1] & (one<< 4)) <<13; t0 |= (z[j+1] & (one<< 8)) <<10; t0 |= (z[j+1] & (one<<12)) << 7;
      t0 |= (z[j+1] & (one<<16)) << 4; t0 |= (z[j+1] & (one<<20)) << 1; t0 |= (z[j+1] & (one<<24)) >> 2; t0 |= (z[j+1] & (one<<28)) >> 5;
      t0 |= (z[j+1] & (one<<32)) >> 8; t0 |= (z[j+1] & (one<<36)) >>11; t0 |= (z[j+1] & (one<<40)) >>14; t0 |= (z[j+1] & (one<<44)) >>17;
      t0 |= (z[j+1] & (one<<48)) >>20; t0 |= (z[j+1] & (one<<52)) >>23; t0 |= (z[j+1] & (one<<56)) >>26; t0 |= (z[j+1] & (one<<60)) >>29;

      t1 |= (z[j+1] & (one<< 1)) <<15; t1 |= (z[j+1] & (one<< 5)) <<12; t1 |= (z[j+1] & (one<< 9)) << 9; t1 |= (z[j+1] & (one<<13)) << 6;
      t1 |= (z[j+1] & (one<<17)) << 3; t1 |= (z[j+1] & (one<<21)) << 0; t1 |= (z[j+1] & (one<<25)) >> 3; t1 |= (z[j+1] & (one<<29)) >> 6;
      t1 |= (z[j+1] & (one<<33)) >> 9; t1 |= (z[j+1] & (one<<37)) >>12; t1 |= (z[j+1] & (one<<41)) >>15; t1 |= (z[j+1] & (one<<45)) >>18;
      t1 |= (z[j+1] & (one<<49)) >>21; t1 |= (z[j+1] & (one<<53)) >>24; t1 |= (z[j+1] & (one<<57)) >>27; t1 |= (z[j+1] & (one<<61)) >>30;

      t2 |= (z[j+1] & (one<< 2)) <<14; t2 |= (z[j+1] & (one<< 6)) <<11; t2 |= (z[j+1] & (one<<10)) << 8; t2 |= (z[j+1] & (one<<14)) << 5;
      t2 |= (z[j+1] & (one<<18)) << 2; t2 |= (z[j+1] & (one<<22)) >> 1; t2 |= (z[j+1] & (one<<26)) >> 4; t2 |= (z[j+1] & (one<<30)) >> 7;
      t2 |= (z[j+1] & (one<<34)) >>10; t2 |= (z[j+1] & (one<<38)) >>13; t2 |= (z[j+1] & (one<<42)) >>16; t2 |= (z[j+1] & (one<<46)) >>19;
      t2 |= (z[j+1] & (one<<50)) >>22; t2 |= (z[j+1] & (one<<54)) >>25; t2 |= (z[j+1] & (one<<58)) >>28; t2 |= (z[j+1] & (one<<62)) >>31;
    case 1:
      t0 |= (z[j+0] & (one<< 0)) >> 0; t0 |= (z[j+0] & (one<< 4)) >> 3; t0 |= (z[j+0] & (one<< 8)) >> 6; t0 |= (z[j+0] & (one<<12)) >> 9;
      t0 |= (z[j+0] & (one<<16)) >>12; t0 |= (z[j+0] & (one<<20)) >>15; t0 |= (z[j+0] & (one<<24)) >>18; t0 |= (z[j+0] & (one<<28)) >>21;
      t0 |= (z[j+0] & (one<<32)) >>24; t0 |= (z[j+0] & (one<<36)) >>27; t0 |= (z[j+0] & (one<<40)) >>30; t0 |= (z[j+0] & (one<<44)) >>33;
      t0 |= (z[j+0] & (one<<48)) >>36; t0 |= (z[j+0] & (one<<52)) >>39; t0 |= (z[j+0] & (one<<56)) >>42; t0 |= (z[j+0] & (one<<60)) >>45;

      t1 |= (z[j+0] & (one<< 1)) >> 1; t1 |= (z[j+0] & (one<< 5)) >> 4; t1 |= (z[j+0] & (one<< 9)) >> 7; t1 |= (z[j+0] & (one<<13)) >>10;
      t1 |= (z[j+0] & (one<<17)) >>13; t1 |= (z[j+0] & (one<<21)) >>16; t1 |= (z[j+0] & (one<<25)) >>19; t1 |= (z[j+0] & (one<<29)) >>22;
      t1 |= (z[j+0] & (one<<33)) >>25; t1 |= (z[j+0] & (one<<37)) >>28; t1 |= (z[j+0] & (one<<41)) >>31; t1 |= (z[j+0] & (one<<45)) >>34;
      t1 |= (z[j+0] & (one<<49)) >>37; t1 |= (z[j+0] & (one<<53)) >>40; t1 |= (z[j+0] & (one<<57)) >>43; t1 |= (z[j+0] & (one<<61)) >>46;

      t2 |= (z[j+0] & (one<< 2)) >> 2; t2 |= (z[j+0] & (one<< 6)) >> 5; t2 |= (z[j+0] & (one<<10)) >> 8; t2 |= (z[j+0] & (one<<14)) >>11;
      t2 |= (z[j+0] & (one<<18)) >>14; t2 |= (z[j+0] & (one<<22)) >>17; t2 |= (z[j+0] & (one<<26)) >>20; t2 |= (z[j+0] & (one<<30)) >>23;
      t2 |= (z[j+0] & (one<<34)) >>26; t2 |= (z[j+0] & (one<<38)) >>29; t2 |= (z[j+0] & (one<<42)) >>32; t2 |= (z[j+0] & (one<<46)) >>35;
      t2 |= (z[j+0] & (one<<50)) >>38; t2 |= (z[j+0] & (one<<54)) >>41; t2 |= (z[j+0] & (one<<58)) >>44; t2 |= (z[j+0] & (one<<62)) >>47;
      break;
    default:
      m4ri_die("impossible");
    }
    a0[j2] |= t0 & bitmask_end;
    a1[j2] |= t1 & bitmask_end;
    a2[j2] |= t2 & bitmask_end;

    /* fix first bits before offset */
    a0[0] = (a0[0] & bitmask_begin) | (a0_fix & ~bitmask_begin);
    a1[0] = (a1[0] & bitmask_begin) | (a1_fix & ~bitmask_begin);
    a2[0] = (a2[0] & bitmask_begin) | (a2_fix & ~bitmask_begin);
  }

  if(A->depth == 3)
    return A;

  /* A3 */
  for(size_t i=0; i<A->nrows; i++) {
    word *a3 = A->x[3]->rows[i];
    const word const *z  = Z->x->rows[i];
    const word a3_fix = a3[0]; 

    /* bulk of work */
    for(j=0, j2=0; j+4 < Z->x->width; j+=4,j2++) {
      if ( !(z[j+0] | z[j+1] | z[j+2] | z[j+3]) )
        continue;
      t3 = 0;

      t3 |= (z[j+3] & (one<< 3)) <<45; t3 |= (z[j+3] & (one<< 7)) <<42; t3 |= (z[j+3] & (one<<11)) <<39; t3 |= (z[j+3] & (one<<15)) <<36;
      t3 |= (z[j+3] & (one<<19)) <<33; t3 |= (z[j+3] & (one<<23)) <<30; t3 |= (z[j+3] & (one<<27)) <<27; t3 |= (z[j+3] & (one<<31)) <<24;
      t3 |= (z[j+3] & (one<<35)) <<21; t3 |= (z[j+3] & (one<<39)) <<18; t3 |= (z[j+3] & (one<<43)) <<15; t3 |= (z[j+3] & (one<<47)) <<12;
      t3 |= (z[j+3] & (one<<51)) << 9; t3 |= (z[j+3] & (one<<55)) << 6; t3 |= (z[j+3] & (one<<59)) << 3; t3 |= (z[j+3] & (one<<63)) << 0;
                                      
      t3 |= (z[j+2] & (one<< 3)) <<29; t3 |= (z[j+2] & (one<< 7)) <<26; t3 |= (z[j+2] & (one<<11)) <<23; t3 |= (z[j+2] & (one<<15)) <<20;
      t3 |= (z[j+2] & (one<<19)) <<17; t3 |= (z[j+2] & (one<<23)) <<14; t3 |= (z[j+2] & (one<<27)) <<11; t3 |= (z[j+2] & (one<<31)) << 8;
      t3 |= (z[j+2] & (one<<35)) << 5; t3 |= (z[j+2] & (one<<39)) << 2; t3 |= (z[j+2] & (one<<43)) >> 1; t3 |= (z[j+2] & (one<<47)) >> 4;
      t3 |= (z[j+2] & (one<<51)) >> 7; t3 |= (z[j+2] & (one<<55)) >>10; t3 |= (z[j+2] & (one<<59)) >>13; t3 |= (z[j+2] & (one<<63)) >>16;
                                      
      t3 |= (z[j+1] & (one<< 3)) <<13; t3 |= (z[j+1] & (one<< 7)) <<10; t3 |= (z[j+1] & (one<<11)) << 7; t3 |= (z[j+1] & (one<<15)) << 4;
      t3 |= (z[j+1] & (one<<19)) << 1; t3 |= (z[j+1] & (one<<23)) >> 2; t3 |= (z[j+1] & (one<<27)) >> 5; t3 |= (z[j+1] & (one<<31)) >> 8;
      t3 |= (z[j+1] & (one<<35)) >>11; t3 |= (z[j+1] & (one<<39)) >>14; t3 |= (z[j+1] & (one<<43)) >>17; t3 |= (z[j+1] & (one<<47)) >>20;
      t3 |= (z[j+1] & (one<<51)) >>23; t3 |= (z[j+1] & (one<<55)) >>26; t3 |= (z[j+1] & (one<<59)) >>29; t3 |= (z[j+1] & (one<<63)) >>32;
                                      
      t3 |= (z[j+0] & (one<< 3)) >> 3; t3 |= (z[j+0] & (one<< 7)) >> 6; t3 |= (z[j+0] & (one<<11)) >> 9; t3 |= (z[j+0] & (one<<15)) >>12;
      t3 |= (z[j+0] & (one<<19)) >>15; t3 |= (z[j+0] & (one<<23)) >>18; t3 |= (z[j+0] & (one<<27)) >>21; t3 |= (z[j+0] & (one<<31)) >>24;
      t3 |= (z[j+0] & (one<<35)) >>27; t3 |= (z[j+0] & (one<<39)) >>30; t3 |= (z[j+0] & (one<<43)) >>33; t3 |= (z[j+0] & (one<<47)) >>36;
      t3 |= (z[j+0] & (one<<51)) >>39; t3 |= (z[j+0] & (one<<55)) >>42; t3 |= (z[j+0] & (one<<59)) >>45; t3 |= (z[j+0] & (one<<63)) >>48;

      a3[j2] = t3;
    }
    t3 = 0;
    switch(Z->x->width - j) {
    case 4:
      t3 |= (z[j+3] & (one<< 3)) <<45; t3 |= (z[j+3] & (one<< 7)) <<42; t3 |= (z[j+3] & (one<<11)) <<39; t3 |= (z[j+3] & (one<<15)) <<36;
      t3 |= (z[j+3] & (one<<19)) <<33; t3 |= (z[j+3] & (one<<23)) <<30; t3 |= (z[j+3] & (one<<27)) <<27; t3 |= (z[j+3] & (one<<31)) <<24;
      t3 |= (z[j+3] & (one<<35)) <<21; t3 |= (z[j+3] & (one<<39)) <<18; t3 |= (z[j+3] & (one<<43)) <<15; t3 |= (z[j+3] & (one<<47)) <<12;
      t3 |= (z[j+3] & (one<<51)) << 9; t3 |= (z[j+3] & (one<<55)) << 6; t3 |= (z[j+3] & (one<<59)) << 3; t3 |= (z[j+3] & (one<<63)) << 0;
    case 3:
      t3 |= (z[j+2] & (one<< 3)) <<29; t3 |= (z[j+2] & (one<< 7)) <<26; t3 |= (z[j+2] & (one<<11)) <<23; t3 |= (z[j+2] & (one<<15)) <<20;
      t3 |= (z[j+2] & (one<<19)) <<17; t3 |= (z[j+2] & (one<<23)) <<14; t3 |= (z[j+2] & (one<<27)) <<11; t3 |= (z[j+2] & (one<<31)) << 8;
      t3 |= (z[j+2] & (one<<35)) << 5; t3 |= (z[j+2] & (one<<39)) << 2; t3 |= (z[j+2] & (one<<43)) >> 1; t3 |= (z[j+2] & (one<<47)) >> 4;
      t3 |= (z[j+2] & (one<<51)) >> 7; t3 |= (z[j+2] & (one<<55)) >>10; t3 |= (z[j+2] & (one<<59)) >>13; t3 |= (z[j+2] & (one<<63)) >>16;
    case 2:
      t3 |= (z[j+1] & (one<< 3)) <<13; t3 |= (z[j+1] & (one<< 7)) <<10; t3 |= (z[j+1] & (one<<11)) << 7; t3 |= (z[j+1] & (one<<15)) << 4;
      t3 |= (z[j+1] & (one<<19)) << 1; t3 |= (z[j+1] & (one<<23)) >> 2; t3 |= (z[j+1] & (one<<27)) >> 5; t3 |= (z[j+1] & (one<<31)) >> 8;
      t3 |= (z[j+1] & (one<<35)) >>11; t3 |= (z[j+1] & (one<<39)) >>14; t3 |= (z[j+1] & (one<<43)) >>17; t3 |= (z[j+1] & (one<<47)) >>20;
      t3 |= (z[j+1] & (one<<51)) >>23; t3 |= (z[j+1] & (one<<55)) >>26; t3 |= (z[j+1] & (one<<59)) >>29; t3 |= (z[j+1] & (one<<63)) >>32;
    case 1:
      t3 |= (z[j+0] & (one<< 3)) >> 3; t3 |= (z[j+0] & (one<< 7)) >> 6; t3 |= (z[j+0] & (one<<11)) >> 9; t3 |= (z[j+0] & (one<<15)) >>12;
      t3 |= (z[j+0] & (one<<19)) >>15; t3 |= (z[j+0] & (one<<23)) >>18; t3 |= (z[j+0] & (one<<27)) >>21; t3 |= (z[j+0] & (one<<31)) >>24;
      t3 |= (z[j+0] & (one<<35)) >>27; t3 |= (z[j+0] & (one<<39)) >>30; t3 |= (z[j+0] & (one<<43)) >>33; t3 |= (z[j+0] & (one<<47)) >>36;
      t3 |= (z[j+0] & (one<<51)) >>39; t3 |= (z[j+0] & (one<<55)) >>42; t3 |= (z[j+0] & (one<<59)) >>45; t3 |= (z[j+0] & (one<<63)) >>48;
      break;
    default:
      m4ri_die("impossible");
    }
    a3[j2] |= t3 & bitmask_end;
    a3[0] = (a3[0] & bitmask_begin) | (a3_fix & ~bitmask_begin);

 }
  return A;
}

mzed_t *_mzed_cling4(mzed_t *A, const mzd_slice_t *Z) {
  size_t j,j2 = 0;
  register word t0, t1, t2, t3;
  const word one = m4ri_one;

  const word bitmask_begin = __M4RI_RIGHT_BITMASK(m4ri_radix - A->x->offset%m4ri_radix);
  const word bitmask_end = __M4RI_LEFT_BITMASK((A->x->offset + A->x->ncols) % m4ri_radix);

  for(size_t i=0; i<A->nrows; i++) {
    word *z0 = Z->x[0]->rows[i];
    word *z1 = Z->x[1]->rows[i];
    word *z2 = Z->x[2]->rows[i];
    word *a  = A->x->rows[i];
    const word a_fix = a[0];

    for(j=0, j2=0; j+4 < A->x->width; j+=4, j2++) {
      if (!z0[j2] && !z1[j2] && !z2[j2] )
        continue;
      t0 = t1 = t2 = t3 = 0;
      t0 |= (z0[j2] & one<< 0) <<  0; t0 |= (z0[j2] & one<< 1) <<  3; t0 |= (z0[j2] & one<< 2) <<  6; t0 |= (z0[j2] & one<< 3) <<  9;
      t0 |= (z0[j2] & one<< 4) << 12; t0 |= (z0[j2] & one<< 5) << 15; t0 |= (z0[j2] & one<< 6) << 18; t0 |= (z0[j2] & one<< 7) << 21;
      t0 |= (z0[j2] & one<< 8) << 24; t0 |= (z0[j2] & one<< 9) << 27; t0 |= (z0[j2] & one<<10) << 30; t0 |= (z0[j2] & one<<11) << 33;
      t0 |= (z0[j2] & one<<12) << 36; t0 |= (z0[j2] & one<<13) << 39; t0 |= (z0[j2] & one<<14) << 42; t0 |= (z0[j2] & one<<15) << 45;

      t0 |= (z1[j2] & one<< 0) <<  1; t0 |= (z1[j2] & one<< 1) <<  4; t0 |= (z1[j2] & one<< 2) <<  7; t0 |= (z1[j2] & one<< 3) << 10;
      t0 |= (z1[j2] & one<< 4) << 13; t0 |= (z1[j2] & one<< 5) << 16; t0 |= (z1[j2] & one<< 6) << 19; t0 |= (z1[j2] & one<< 7) << 22;
      t0 |= (z1[j2] & one<< 8) << 25; t0 |= (z1[j2] & one<< 9) << 28; t0 |= (z1[j2] & one<<10) << 31; t0 |= (z1[j2] & one<<11) << 34;
      t0 |= (z1[j2] & one<<12) << 37; t0 |= (z1[j2] & one<<13) << 40; t0 |= (z1[j2] & one<<14) << 43; t0 |= (z1[j2] & one<<15) << 46;

      t0 |= (z2[j2] & one<< 0) <<  2; t0 |= (z2[j2] & one<< 1) <<  5; t0 |= (z2[j2] & one<< 2) <<  8; t0 |= (z2[j2] & one<< 3) << 11;
      t0 |= (z2[j2] & one<< 4) << 14; t0 |= (z2[j2] & one<< 5) << 17; t0 |= (z2[j2] & one<< 6) << 20; t0 |= (z2[j2] & one<< 7) << 23;
      t0 |= (z2[j2] & one<< 8) << 26; t0 |= (z2[j2] & one<< 9) << 29; t0 |= (z2[j2] & one<<10) << 32; t0 |= (z2[j2] & one<<11) << 35;
      t0 |= (z2[j2] & one<<12) << 38; t0 |= (z2[j2] & one<<13) << 41; t0 |= (z2[j2] & one<<14) << 44; t0 |= (z2[j2] & one<<15) << 47;

      t1 |= (z0[j2] & one<<16) >> 16; t1 |= (z0[j2] & one<<17) >> 13; t1 |= (z0[j2] & one<<18) >> 10; t1 |= (z0[j2] & one<<19) >>  7;
      t1 |= (z0[j2] & one<<20) >>  4; t1 |= (z0[j2] & one<<21) >>  1; t1 |= (z0[j2] & one<<22) <<  2; t1 |= (z0[j2] & one<<23) <<  5;
      t1 |= (z0[j2] & one<<24) <<  8; t1 |= (z0[j2] & one<<25) << 11; t1 |= (z0[j2] & one<<26) << 14; t1 |= (z0[j2] & one<<27) << 17;
      t1 |= (z0[j2] & one<<28) << 20; t1 |= (z0[j2] & one<<29) << 23; t1 |= (z0[j2] & one<<30) << 26; t1 |= (z0[j2] & one<<31) << 29;

      t1 |= (z1[j2] & one<<16) >> 15; t1 |= (z1[j2] & one<<17) >> 12; t1 |= (z1[j2] & one<<18) >>  9; t1 |= (z1[j2] & one<<19) >>  6;
      t1 |= (z1[j2] & one<<20) >>  3; t1 |= (z1[j2] & one<<21) <<  0; t1 |= (z1[j2] & one<<22) <<  3; t1 |= (z1[j2] & one<<23) <<  6;
      t1 |= (z1[j2] & one<<24) <<  9; t1 |= (z1[j2] & one<<25) << 12; t1 |= (z1[j2] & one<<26) << 15; t1 |= (z1[j2] & one<<27) << 18;
      t1 |= (z1[j2] & one<<28) << 21; t1 |= (z1[j2] & one<<29) << 24; t1 |= (z1[j2] & one<<30) << 27; t1 |= (z1[j2] & one<<31) << 30;

      t1 |= (z2[j2] & one<<16) >> 14; t1 |= (z2[j2] & one<<17) >> 11; t1 |= (z2[j2] & one<<18) >>  8; t1 |= (z2[j2] & one<<19) >>  5;
      t1 |= (z2[j2] & one<<20) >>  2; t1 |= (z2[j2] & one<<21) <<  1; t1 |= (z2[j2] & one<<22) <<  4; t1 |= (z2[j2] & one<<23) <<  7;
      t1 |= (z2[j2] & one<<24) << 10; t1 |= (z2[j2] & one<<25) << 13; t1 |= (z2[j2] & one<<26) << 16; t1 |= (z2[j2] & one<<27) << 19;
      t1 |= (z2[j2] & one<<28) << 22; t1 |= (z2[j2] & one<<29) << 25; t1 |= (z2[j2] & one<<30) << 28; t1 |= (z2[j2] & one<<31) << 31;

      t2 |= (z0[j2] & one<<32) >> 32; t2 |= (z0[j2] & one<<33) >> 29; t2 |= (z0[j2] & one<<34) >> 26; t2 |= (z0[j2] & one<<35) >> 23;
      t2 |= (z0[j2] & one<<36) >> 20; t2 |= (z0[j2] & one<<37) >> 17; t2 |= (z0[j2] & one<<38) >> 14; t2 |= (z0[j2] & one<<39) >> 11;
      t2 |= (z0[j2] & one<<40) >>  8; t2 |= (z0[j2] & one<<41) >>  5; t2 |= (z0[j2] & one<<42) >>  2; t2 |= (z0[j2] & one<<43) <<  1;
      t2 |= (z0[j2] & one<<44) <<  4; t2 |= (z0[j2] & one<<45) <<  7; t2 |= (z0[j2] & one<<46) << 10; t2 |= (z0[j2] & one<<47) << 13;

      t2 |= (z1[j2] & one<<32) >> 31; t2 |= (z1[j2] & one<<33) >> 28; t2 |= (z1[j2] & one<<34) >> 25; t2 |= (z1[j2] & one<<35) >> 22;
      t2 |= (z1[j2] & one<<36) >> 19; t2 |= (z1[j2] & one<<37) >> 16; t2 |= (z1[j2] & one<<38) >> 13; t2 |= (z1[j2] & one<<39) >> 10;
      t2 |= (z1[j2] & one<<40) >>  7; t2 |= (z1[j2] & one<<41) >>  4; t2 |= (z1[j2] & one<<42) >>  1; t2 |= (z1[j2] & one<<43) <<  2;
      t2 |= (z1[j2] & one<<44) <<  5; t2 |= (z1[j2] & one<<45) <<  8; t2 |= (z1[j2] & one<<46) << 11; t2 |= (z1[j2] & one<<47) << 14;

      t2 |= (z2[j2] & one<<32) >> 30; t2 |= (z2[j2] & one<<33) >> 27; t2 |= (z2[j2] & one<<34) >> 24; t2 |= (z2[j2] & one<<35) >> 21;
      t2 |= (z2[j2] & one<<36) >> 18; t2 |= (z2[j2] & one<<37) >> 15; t2 |= (z2[j2] & one<<38) >> 12; t2 |= (z2[j2] & one<<39) >>  9;
      t2 |= (z2[j2] & one<<40) >>  6; t2 |= (z2[j2] & one<<41) >>  3; t2 |= (z2[j2] & one<<42) <<  0; t2 |= (z2[j2] & one<<43) <<  3;
      t2 |= (z2[j2] & one<<44) <<  6; t2 |= (z2[j2] & one<<45) <<  9; t2 |= (z2[j2] & one<<46) << 12; t2 |= (z2[j2] & one<<47) << 15;

      t3 |= (z0[j2] & one<<48) >> 48; t3 |= (z0[j2] & one<<49) >> 45; t3 |= (z0[j2] & one<<50) >> 42; t3 |= (z0[j2] & one<<51) >> 39;
      t3 |= (z0[j2] & one<<52) >> 36; t3 |= (z0[j2] & one<<53) >> 33; t3 |= (z0[j2] & one<<54) >> 30; t3 |= (z0[j2] & one<<55) >> 27;
      t3 |= (z0[j2] & one<<56) >> 24; t3 |= (z0[j2] & one<<57) >> 21; t3 |= (z0[j2] & one<<58) >> 18; t3 |= (z0[j2] & one<<59) >> 15;
      t3 |= (z0[j2] & one<<60) >> 12; t3 |= (z0[j2] & one<<61) >>  9; t3 |= (z0[j2] & one<<62) >>  6; t3 |= (z0[j2] & one<<63) >>  3;

      t3 |= (z1[j2] & one<<48) >> 47; t3 |= (z1[j2] & one<<49) >> 44; t3 |= (z1[j2] & one<<50) >> 41; t3 |= (z1[j2] & one<<51) >> 38;
      t3 |= (z1[j2] & one<<52) >> 35; t3 |= (z1[j2] & one<<53) >> 32; t3 |= (z1[j2] & one<<54) >> 29; t3 |= (z1[j2] & one<<55) >> 26;
      t3 |= (z1[j2] & one<<56) >> 23; t3 |= (z1[j2] & one<<57) >> 20; t3 |= (z1[j2] & one<<58) >> 17; t3 |= (z1[j2] & one<<59) >> 14;
      t3 |= (z1[j2] & one<<60) >> 11; t3 |= (z1[j2] & one<<61) >>  8; t3 |= (z1[j2] & one<<62) >>  5; t3 |= (z1[j2] & one<<63) >>  2;

      t3 |= (z2[j2] & one<<48) >> 46; t3 |= (z2[j2] & one<<49) >> 43; t3 |= (z2[j2] & one<<50) >> 40; t3 |= (z2[j2] & one<<51) >> 37;
      t3 |= (z2[j2] & one<<52) >> 34; t3 |= (z2[j2] & one<<53) >> 31; t3 |= (z2[j2] & one<<54) >> 28; t3 |= (z2[j2] & one<<55) >> 25;
      t3 |= (z2[j2] & one<<56) >> 22; t3 |= (z2[j2] & one<<57) >> 19; t3 |= (z2[j2] & one<<58) >> 16; t3 |= (z2[j2] & one<<59) >> 13;
      t3 |= (z2[j2] & one<<60) >> 10; t3 |= (z2[j2] & one<<61) >>  7; t3 |= (z2[j2] & one<<62) >>  4; t3 |= (z2[j2] & one<<63) >>  1;

      a[j+0] = t0;
      a[j+1] = t1;
      a[j+2] = t2;
      a[j+3] = t3;
    }

    /*
     * This is not efficient for very small matrices since we are
     * doing a lot of useless computations.
     */

    t0 = t1 = t2 = t3 = 0;
    switch(A->x->width - j) {
    case 4:
      t3 |= (z0[j2] & one<<48) >> 48; t3 |= (z0[j2] & one<<49) >> 45; t3 |= (z0[j2] & one<<50) >> 42; t3 |= (z0[j2] & one<<51) >> 39;
      t3 |= (z0[j2] & one<<52) >> 36; t3 |= (z0[j2] & one<<53) >> 33; t3 |= (z0[j2] & one<<54) >> 30; t3 |= (z0[j2] & one<<55) >> 27;
      t3 |= (z0[j2] & one<<56) >> 24; t3 |= (z0[j2] & one<<57) >> 21; t3 |= (z0[j2] & one<<58) >> 18; t3 |= (z0[j2] & one<<59) >> 15;
      t3 |= (z0[j2] & one<<60) >> 12; t3 |= (z0[j2] & one<<61) >>  9; t3 |= (z0[j2] & one<<62) >>  6; t3 |= (z0[j2] & one<<63) >>  3;

      t3 |= (z1[j2] & one<<48) >> 47; t3 |= (z1[j2] & one<<49) >> 44; t3 |= (z1[j2] & one<<50) >> 41; t3 |= (z1[j2] & one<<51) >> 38;
      t3 |= (z1[j2] & one<<52) >> 35; t3 |= (z1[j2] & one<<53) >> 32; t3 |= (z1[j2] & one<<54) >> 29; t3 |= (z1[j2] & one<<55) >> 26;
      t3 |= (z1[j2] & one<<56) >> 23; t3 |= (z1[j2] & one<<57) >> 20; t3 |= (z1[j2] & one<<58) >> 17; t3 |= (z1[j2] & one<<59) >> 14;
      t3 |= (z1[j2] & one<<60) >> 11; t3 |= (z1[j2] & one<<61) >>  8; t3 |= (z1[j2] & one<<62) >>  5; t3 |= (z1[j2] & one<<63) >>  2;

      t3 |= (z2[j2] & one<<48) >> 46; t3 |= (z2[j2] & one<<49) >> 43; t3 |= (z2[j2] & one<<50) >> 40; t3 |= (z2[j2] & one<<51) >> 37;
      t3 |= (z2[j2] & one<<52) >> 34; t3 |= (z2[j2] & one<<53) >> 31; t3 |= (z2[j2] & one<<54) >> 28; t3 |= (z2[j2] & one<<55) >> 25;
      t3 |= (z2[j2] & one<<56) >> 22; t3 |= (z2[j2] & one<<57) >> 19; t3 |= (z2[j2] & one<<58) >> 16; t3 |= (z2[j2] & one<<59) >> 13;
      t3 |= (z2[j2] & one<<60) >> 10; t3 |= (z2[j2] & one<<61) >>  7; t3 |= (z2[j2] & one<<62) >>  4; t3 |= (z2[j2] & one<<63) >>  1;
    case 3:
      t2 |= (z0[j2] & one<<32) >> 32; t2 |= (z0[j2] & one<<33) >> 29; t2 |= (z0[j2] & one<<34) >> 26; t2 |= (z0[j2] & one<<35) >> 23;
      t2 |= (z0[j2] & one<<36) >> 20; t2 |= (z0[j2] & one<<37) >> 17; t2 |= (z0[j2] & one<<38) >> 14; t2 |= (z0[j2] & one<<39) >> 11;
      t2 |= (z0[j2] & one<<40) >>  8; t2 |= (z0[j2] & one<<41) >>  5; t2 |= (z0[j2] & one<<42) >>  2; t2 |= (z0[j2] & one<<43) <<  1;
      t2 |= (z0[j2] & one<<44) <<  4; t2 |= (z0[j2] & one<<45) <<  7; t2 |= (z0[j2] & one<<46) << 10; t2 |= (z0[j2] & one<<47) << 13;

      t2 |= (z1[j2] & one<<32) >> 31; t2 |= (z1[j2] & one<<33) >> 28; t2 |= (z1[j2] & one<<34) >> 25; t2 |= (z1[j2] & one<<35) >> 22;
      t2 |= (z1[j2] & one<<36) >> 19; t2 |= (z1[j2] & one<<37) >> 16; t2 |= (z1[j2] & one<<38) >> 13; t2 |= (z1[j2] & one<<39) >> 10;
      t2 |= (z1[j2] & one<<40) >>  7; t2 |= (z1[j2] & one<<41) >>  4; t2 |= (z1[j2] & one<<42) >>  1; t2 |= (z1[j2] & one<<43) <<  2;
      t2 |= (z1[j2] & one<<44) <<  5; t2 |= (z1[j2] & one<<45) <<  8; t2 |= (z1[j2] & one<<46) << 11; t2 |= (z1[j2] & one<<47) << 14;

      t2 |= (z2[j2] & one<<32) >> 30; t2 |= (z2[j2] & one<<33) >> 27; t2 |= (z2[j2] & one<<34) >> 24; t2 |= (z2[j2] & one<<35) >> 21;
      t2 |= (z2[j2] & one<<36) >> 18; t2 |= (z2[j2] & one<<37) >> 15; t2 |= (z2[j2] & one<<38) >> 12; t2 |= (z2[j2] & one<<39) >>  9;
      t2 |= (z2[j2] & one<<40) >>  6; t2 |= (z2[j2] & one<<41) >>  3; t2 |= (z2[j2] & one<<42) <<  0; t2 |= (z2[j2] & one<<43) <<  3;
      t2 |= (z2[j2] & one<<44) <<  6; t2 |= (z2[j2] & one<<45) <<  9; t2 |= (z2[j2] & one<<46) << 12; t2 |= (z2[j2] & one<<47) << 15;
    case 2:
      t1 |= (z0[j2] & one<<16) >> 16; t1 |= (z0[j2] & one<<17) >> 13; t1 |= (z0[j2] & one<<18) >> 10; t1 |= (z0[j2] & one<<19) >>  7;
      t1 |= (z0[j2] & one<<20) >>  4; t1 |= (z0[j2] & one<<21) >>  1; t1 |= (z0[j2] & one<<22) <<  2; t1 |= (z0[j2] & one<<23) <<  5;
      t1 |= (z0[j2] & one<<24) <<  8; t1 |= (z0[j2] & one<<25) << 11; t1 |= (z0[j2] & one<<26) << 14; t1 |= (z0[j2] & one<<27) << 17;
      t1 |= (z0[j2] & one<<28) << 20; t1 |= (z0[j2] & one<<29) << 23; t1 |= (z0[j2] & one<<30) << 26; t1 |= (z0[j2] & one<<31) << 29;

      t1 |= (z1[j2] & one<<16) >> 15; t1 |= (z1[j2] & one<<17) >> 12; t1 |= (z1[j2] & one<<18) >>  9; t1 |= (z1[j2] & one<<19) >>  6;
      t1 |= (z1[j2] & one<<20) >>  3; t1 |= (z1[j2] & one<<21) <<  0; t1 |= (z1[j2] & one<<22) <<  3; t1 |= (z1[j2] & one<<23) <<  6;
      t1 |= (z1[j2] & one<<24) <<  9; t1 |= (z1[j2] & one<<25) << 12; t1 |= (z1[j2] & one<<26) << 15; t1 |= (z1[j2] & one<<27) << 18;
      t1 |= (z1[j2] & one<<28) << 21; t1 |= (z1[j2] & one<<29) << 24; t1 |= (z1[j2] & one<<30) << 27; t1 |= (z1[j2] & one<<31) << 30;

      t1 |= (z2[j2] & one<<16) >> 14; t1 |= (z2[j2] & one<<17) >> 11; t1 |= (z2[j2] & one<<18) >>  8; t1 |= (z2[j2] & one<<19) >>  5;
      t1 |= (z2[j2] & one<<20) >>  2; t1 |= (z2[j2] & one<<21) <<  1; t1 |= (z2[j2] & one<<22) <<  4; t1 |= (z2[j2] & one<<23) <<  7;
      t1 |= (z2[j2] & one<<24) << 10; t1 |= (z2[j2] & one<<25) << 13; t1 |= (z2[j2] & one<<26) << 16; t1 |= (z2[j2] & one<<27) << 19;
      t1 |= (z2[j2] & one<<28) << 22; t1 |= (z2[j2] & one<<29) << 25; t1 |= (z2[j2] & one<<30) << 28; t1 |= (z2[j2] & one<<31) << 31;
    case 1:
      t0 |= (z0[j2] & one<< 0) <<  0; t0 |= (z0[j2] & one<< 1) <<  3; t0 |= (z0[j2] & one<< 2) <<  6; t0 |= (z0[j2] & one<< 3) <<  9;
      t0 |= (z0[j2] & one<< 4) << 12; t0 |= (z0[j2] & one<< 5) << 15; t0 |= (z0[j2] & one<< 6) << 18; t0 |= (z0[j2] & one<< 7) << 21;
      t0 |= (z0[j2] & one<< 8) << 24; t0 |= (z0[j2] & one<< 9) << 27; t0 |= (z0[j2] & one<<10) << 30; t0 |= (z0[j2] & one<<11) << 33;
      t0 |= (z0[j2] & one<<12) << 36; t0 |= (z0[j2] & one<<13) << 39; t0 |= (z0[j2] & one<<14) << 42; t0 |= (z0[j2] & one<<15) << 45;

      t0 |= (z1[j2] & one<< 0) <<  1; t0 |= (z1[j2] & one<< 1) <<  4; t0 |= (z1[j2] & one<< 2) <<  7; t0 |= (z1[j2] & one<< 3) << 10;
      t0 |= (z1[j2] & one<< 4) << 13; t0 |= (z1[j2] & one<< 5) << 16; t0 |= (z1[j2] & one<< 6) << 19; t0 |= (z1[j2] & one<< 7) << 22;
      t0 |= (z1[j2] & one<< 8) << 25; t0 |= (z1[j2] & one<< 9) << 28; t0 |= (z1[j2] & one<<10) << 31; t0 |= (z1[j2] & one<<11) << 34;
      t0 |= (z1[j2] & one<<12) << 37; t0 |= (z1[j2] & one<<13) << 40; t0 |= (z1[j2] & one<<14) << 43; t0 |= (z1[j2] & one<<15) << 46;

      t0 |= (z2[j2] & one<< 0) <<  2; t0 |= (z2[j2] & one<< 1) <<  5; t0 |= (z2[j2] & one<< 2) <<  8; t0 |= (z2[j2] & one<< 3) << 11;
      t0 |= (z2[j2] & one<< 4) << 14; t0 |= (z2[j2] & one<< 5) << 17; t0 |= (z2[j2] & one<< 6) << 20; t0 |= (z2[j2] & one<< 7) << 23;
      t0 |= (z2[j2] & one<< 8) << 26; t0 |= (z2[j2] & one<< 9) << 29; t0 |= (z2[j2] & one<<10) << 32; t0 |= (z2[j2] & one<<11) << 35;
      t0 |= (z2[j2] & one<<12) << 38; t0 |= (z2[j2] & one<<13) << 41; t0 |= (z2[j2] & one<<14) << 44; t0 |= (z2[j2] & one<<15) << 47;
      break;
    default:
      m4ri_die("impossible");
    }

    /*
     * We could avoid this second switch case by extendeding the
     * previous one with a lot of copy'n'paste. But this version is
     * easier to read and to maintain.
     */
    
    switch(A->x->width - j) {
    case 4:
      a[j+0] = t0, a[j+1] = t1, a[j+2] = t2, a[j+3] |= t3 & bitmask_end;
      break;
    case 3:
      a[j+0] = t0, a[j+1] = t1, a[j+2] |= t2 & bitmask_end;
      break;
    case 2:
      a[j+0] = t0, a[j+1] |= t1 & bitmask_end;
      break;
    case 1:
      a[j+0] |= t0 & bitmask_end;
      break;
    default:
      m4ri_die("impossible");
    }
    a[0] = (a[0] & bitmask_begin) | (a_fix & ~bitmask_begin);
  }

  if(A->finite_field->degree == 3)
    return A;

  for(size_t i=0; i<A->nrows; i++) {
    word *z3 = Z->x[3]->rows[i];
    word *a  = A->x->rows[i];
    const word a_fix = a[0];

    for(j=0, j2=0; j+4 < A->x->width; j+=4, j2++) {
      if (!z3[j2] )
        continue;
      t0 = t1 = t2 = t3 = 0;

      t0 |= (z3[j2] & one<< 0) <<  3; t0 |= (z3[j2] & one<< 1) <<  6; t0 |= (z3[j2] & one<< 2) <<  9; t0 |= (z3[j2] & one<< 3) << 12;
      t0 |= (z3[j2] & one<< 4) << 15; t0 |= (z3[j2] & one<< 5) << 18; t0 |= (z3[j2] & one<< 6) << 21; t0 |= (z3[j2] & one<< 7) << 24;
      t0 |= (z3[j2] & one<< 8) << 27; t0 |= (z3[j2] & one<< 9) << 30; t0 |= (z3[j2] & one<<10) << 33; t0 |= (z3[j2] & one<<11) << 36;
      t0 |= (z3[j2] & one<<12) << 39; t0 |= (z3[j2] & one<<13) << 42; t0 |= (z3[j2] & one<<14) << 45; t0 |= (z3[j2] & one<<15) << 48;

      t1 |= (z3[j2] & one<<16) >> 13; t1 |= (z3[j2] & one<<17) >> 10; t1 |= (z3[j2] & one<<18) >>  7; t1 |= (z3[j2] & one<<19) >>  4;
      t1 |= (z3[j2] & one<<20) >>  1; t1 |= (z3[j2] & one<<21) <<  2; t1 |= (z3[j2] & one<<22) <<  5; t1 |= (z3[j2] & one<<23) <<  8;
      t1 |= (z3[j2] & one<<24) << 11; t1 |= (z3[j2] & one<<25) << 14; t1 |= (z3[j2] & one<<26) << 17; t1 |= (z3[j2] & one<<27) << 20;
      t1 |= (z3[j2] & one<<28) << 23; t1 |= (z3[j2] & one<<29) << 26; t1 |= (z3[j2] & one<<30) << 29; t1 |= (z3[j2] & one<<31) << 32;

      t2 |= (z3[j2] & one<<32) >> 29; t2 |= (z3[j2] & one<<33) >> 26; t2 |= (z3[j2] & one<<34) >> 23; t2 |= (z3[j2] & one<<35) >> 20;
      t2 |= (z3[j2] & one<<36) >> 17; t2 |= (z3[j2] & one<<37) >> 14; t2 |= (z3[j2] & one<<38) >> 11; t2 |= (z3[j2] & one<<39) >>  8;
      t2 |= (z3[j2] & one<<40) >>  5; t2 |= (z3[j2] & one<<41) >>  2; t2 |= (z3[j2] & one<<42) <<  1; t2 |= (z3[j2] & one<<43) <<  4;
      t2 |= (z3[j2] & one<<44) <<  7; t2 |= (z3[j2] & one<<45) << 10; t2 |= (z3[j2] & one<<46) << 13; t2 |= (z3[j2] & one<<47) << 16;

      t3 |= (z3[j2] & one<<48) >> 45; t3 |= (z3[j2] & one<<49) >> 42; t3 |= (z3[j2] & one<<50) >> 39; t3 |= (z3[j2] & one<<51) >> 36;
      t3 |= (z3[j2] & one<<52) >> 33; t3 |= (z3[j2] & one<<53) >> 30; t3 |= (z3[j2] & one<<54) >> 27; t3 |= (z3[j2] & one<<55) >> 24;
      t3 |= (z3[j2] & one<<56) >> 21; t3 |= (z3[j2] & one<<57) >> 18; t3 |= (z3[j2] & one<<58) >> 15; t3 |= (z3[j2] & one<<59) >> 12;
      t3 |= (z3[j2] & one<<60) >>  9; t3 |= (z3[j2] & one<<61) >>  6; t3 |= (z3[j2] & one<<62) >>  3; t3 |= (z3[j2] & one<<63) <<  0;
      
      a[j+0] |= t0;
      a[j+1] |= t1;
      a[j+2] |= t2;
      a[j+3] |= t3;
    }

    /*
     * This is not efficient for very small matrices since we are
     * doing a lot of useless computations.
     */

    t0 = t1 = t2 = t3 = 0;
    switch(A->x->width - j) {
    case 4:
      t3 |= (z3[j2] & one<<48) >> 45; t3 |= (z3[j2] & one<<49) >> 42; t3 |= (z3[j2] & one<<50) >> 39; t3 |= (z3[j2] & one<<51) >> 36;
      t3 |= (z3[j2] & one<<52) >> 33; t3 |= (z3[j2] & one<<53) >> 30; t3 |= (z3[j2] & one<<54) >> 27; t3 |= (z3[j2] & one<<55) >> 24;
      t3 |= (z3[j2] & one<<56) >> 21; t3 |= (z3[j2] & one<<57) >> 18; t3 |= (z3[j2] & one<<58) >> 15; t3 |= (z3[j2] & one<<59) >> 12;
      t3 |= (z3[j2] & one<<60) >>  9; t3 |= (z3[j2] & one<<61) >>  6; t3 |= (z3[j2] & one<<62) >>  3; t3 |= (z3[j2] & one<<63) <<  0;
    case 3:
      t2 |= (z3[j2] & one<<32) >> 29; t2 |= (z3[j2] & one<<33) >> 26; t2 |= (z3[j2] & one<<34) >> 23; t2 |= (z3[j2] & one<<35) >> 20;
      t2 |= (z3[j2] & one<<36) >> 17; t2 |= (z3[j2] & one<<37) >> 14; t2 |= (z3[j2] & one<<38) >> 11; t2 |= (z3[j2] & one<<39) >>  8;
      t2 |= (z3[j2] & one<<40) >>  5; t2 |= (z3[j2] & one<<41) >>  2; t2 |= (z3[j2] & one<<42) <<  1; t2 |= (z3[j2] & one<<43) <<  4;
      t2 |= (z3[j2] & one<<44) <<  7; t2 |= (z3[j2] & one<<45) << 10; t2 |= (z3[j2] & one<<46) << 13; t2 |= (z3[j2] & one<<47) << 16;
    case 2:
      t1 |= (z3[j2] & one<<16) >> 13; t1 |= (z3[j2] & one<<17) >> 10; t1 |= (z3[j2] & one<<18) >>  7; t1 |= (z3[j2] & one<<19) >>  4;
      t1 |= (z3[j2] & one<<20) >>  1; t1 |= (z3[j2] & one<<21) <<  2; t1 |= (z3[j2] & one<<22) <<  5; t1 |= (z3[j2] & one<<23) <<  8;
      t1 |= (z3[j2] & one<<24) << 11; t1 |= (z3[j2] & one<<25) << 14; t1 |= (z3[j2] & one<<26) << 17; t1 |= (z3[j2] & one<<27) << 20;
      t1 |= (z3[j2] & one<<28) << 23; t1 |= (z3[j2] & one<<29) << 26; t1 |= (z3[j2] & one<<30) << 29; t1 |= (z3[j2] & one<<31) << 32;
    case 1:
      t0 |= (z3[j2] & one<< 0) <<  3; t0 |= (z3[j2] & one<< 1) <<  6; t0 |= (z3[j2] & one<< 2) <<  9; t0 |= (z3[j2] & one<< 3) << 12;
      t0 |= (z3[j2] & one<< 4) << 15; t0 |= (z3[j2] & one<< 5) << 18; t0 |= (z3[j2] & one<< 6) << 21; t0 |= (z3[j2] & one<< 7) << 24;
      t0 |= (z3[j2] & one<< 8) << 27; t0 |= (z3[j2] & one<< 9) << 30; t0 |= (z3[j2] & one<<10) << 33; t0 |= (z3[j2] & one<<11) << 36;
      t0 |= (z3[j2] & one<<12) << 39; t0 |= (z3[j2] & one<<13) << 42; t0 |= (z3[j2] & one<<14) << 45; t0 |= (z3[j2] & one<<15) << 48;
      break;
    default:
      m4ri_die("impossible");
    }

    /*
     * We could avoid this second switch case by extendeding the
     * previous one with a lot of copy'n'paste. But this version is
     * easier to read and to maintain.
     */
    
    switch(A->x->width - j) {
    case 4:
      a[j+0] |= t0, a[j+1] |= t1, a[j+2] |= t2, a[j+3] |= t3 & bitmask_end;
      break;
    case 3:
      a[j+0] |= t0, a[j+1] |= t1, a[j+2] |= t2 & bitmask_end;
      break;
    case 2:
      a[j+0] |= t0, a[j+1] |= t1 & bitmask_end;
      break;
    case 1:
      a[j+0] |= t0 & bitmask_end;
      break;
    default:
      m4ri_die("impossible");
    }
    a[0] = (a[0] & bitmask_begin) | (a_fix & ~bitmask_begin);
  }
  return A;
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

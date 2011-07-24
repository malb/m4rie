#include "bitslice.h"

mzed_t *mzed_cling(mzed_t *A, const mzd_slice_t *Z) {
  switch(Z->finite_field->degree) {
  case  2: return mzed_cling2(A,Z);
  case  3: 
  case  4:
  case  5:
  case  6:
  case  7:
  case  8:
  case  9:
  case 10:
  default:
    m4ri_die("Clinging not implemented for this degree.");
  }
  return A;
}

mzd_slice_t *mzed_slice(mzd_slice_t *A, const mzed_t *Z) {
  switch(Z->finite_field->degree) {
  case  2: return mzed_slice2(A,Z);
  case  3: 
  case  4:
  case  5:
  case  6:
  case  7:
  case  8:
  case  9:
  case 10:
  default:
    m4ri_die("Slicing not implemented for this degree.");
  }
  return A;
}

mzed_t *_mzed_cling2(mzed_t *A, const mzd_slice_t *Z) {
  size_t j,j2 = 0;
  register word aw0;
  register word aw1;

  /** A0 **/
  for(size_t i=0; i<A->nrows; i++) {
    word *z0 = Z->x[0]->rows[i];
    word *a  = A->x->rows[i];
    for(j=0, j2=0; j+2 < A->x->width; j+=2, j2++) {
      if (!z0[j2] )
        continue;
      aw0 = 0;
      aw1 = 0;
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

      a[j+1] = 0;
      switch((2*A->ncols) % m4ri_radix) {
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
      a[j+0] = 0;
      switch((2*A->ncols) % m4ri_radix) {
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
  }

  /** A1 **/
  for(size_t i=0; i<A->nrows; i++) {
    word *z1 = Z->x[1]->rows[i];
    word *a  = A->x->rows[i];    
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
 
      switch((2*A->ncols) % m4ri_radix) {
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
      switch((2*A->ncols) % m4ri_radix) {
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
  }
  return A;
}

mzd_slice_t *_mzed_slice2(mzd_slice_t *A, const mzed_t *Z) {
  size_t j, j2 = 0;
  register word tmp = 0;

  /* Z */
  for(size_t i=0; i<A->nrows; i++) {
    word *a0 = A->x[0]->rows[i];
    const word *z  = Z->x->rows[i];    

    /* bulk of work */
    for(j=0, j2=0; j+2 < Z->x->width; j+=2,j2++) {
      if ( !(z[j+0] | z[j+1]) )
        continue;
      tmp = 0;
      tmp |= (z[j+0] & m4ri_one<< 0) >>  0, tmp |= (z[j+1] & m4ri_one<< 0) << 32,
      tmp |= (z[j+0] & m4ri_one<< 2) >>  1, tmp |= (z[j+1] & m4ri_one<< 2) << 31,
      tmp |= (z[j+0] & m4ri_one<< 4) >>  2, tmp |= (z[j+1] & m4ri_one<< 4) << 30,
      tmp |= (z[j+0] & m4ri_one<< 6) >>  3, tmp |= (z[j+1] & m4ri_one<< 6) << 29,
      tmp |= (z[j+0] & m4ri_one<< 8) >>  4, tmp |= (z[j+1] & m4ri_one<< 8) << 28,
      tmp |= (z[j+0] & m4ri_one<<10) >>  5, tmp |= (z[j+1] & m4ri_one<<10) << 27,
      tmp |= (z[j+0] & m4ri_one<<12) >>  6, tmp |= (z[j+1] & m4ri_one<<12) << 26,
      tmp |= (z[j+0] & m4ri_one<<14) >>  7, tmp |= (z[j+1] & m4ri_one<<14) << 25,
      tmp |= (z[j+0] & m4ri_one<<16) >>  8, tmp |= (z[j+1] & m4ri_one<<16) << 24,
      tmp |= (z[j+0] & m4ri_one<<18) >>  9, tmp |= (z[j+1] & m4ri_one<<18) << 23,
      tmp |= (z[j+0] & m4ri_one<<20) >> 10, tmp |= (z[j+1] & m4ri_one<<20) << 22,
      tmp |= (z[j+0] & m4ri_one<<22) >> 11, tmp |= (z[j+1] & m4ri_one<<22) << 21,
      tmp |= (z[j+0] & m4ri_one<<24) >> 12, tmp |= (z[j+1] & m4ri_one<<24) << 20,
      tmp |= (z[j+0] & m4ri_one<<26) >> 13, tmp |= (z[j+1] & m4ri_one<<26) << 19,
      tmp |= (z[j+0] & m4ri_one<<28) >> 14, tmp |= (z[j+1] & m4ri_one<<28) << 18,
      tmp |= (z[j+0] & m4ri_one<<30) >> 15, tmp |= (z[j+1] & m4ri_one<<30) << 17,
      tmp |= (z[j+0] & m4ri_one<<32) >> 16, tmp |= (z[j+1] & m4ri_one<<32) << 16,
      tmp |= (z[j+0] & m4ri_one<<34) >> 17, tmp |= (z[j+1] & m4ri_one<<34) << 15,
      tmp |= (z[j+0] & m4ri_one<<36) >> 18, tmp |= (z[j+1] & m4ri_one<<36) << 14,
      tmp |= (z[j+0] & m4ri_one<<38) >> 19, tmp |= (z[j+1] & m4ri_one<<38) << 13,
      tmp |= (z[j+0] & m4ri_one<<40) >> 20, tmp |= (z[j+1] & m4ri_one<<40) << 12,
      tmp |= (z[j+0] & m4ri_one<<42) >> 21, tmp |= (z[j+1] & m4ri_one<<42) << 11,
      tmp |= (z[j+0] & m4ri_one<<44) >> 22, tmp |= (z[j+1] & m4ri_one<<44) << 10,
      tmp |= (z[j+0] & m4ri_one<<46) >> 23, tmp |= (z[j+1] & m4ri_one<<46) <<  9,
      tmp |= (z[j+0] & m4ri_one<<48) >> 24, tmp |= (z[j+1] & m4ri_one<<48) <<  8,
      tmp |= (z[j+0] & m4ri_one<<50) >> 25, tmp |= (z[j+1] & m4ri_one<<50) <<  7,
      tmp |= (z[j+0] & m4ri_one<<52) >> 26, tmp |= (z[j+1] & m4ri_one<<52) <<  6,
      tmp |= (z[j+0] & m4ri_one<<54) >> 27, tmp |= (z[j+1] & m4ri_one<<54) <<  5,
      tmp |= (z[j+0] & m4ri_one<<56) >> 28, tmp |= (z[j+1] & m4ri_one<<56) <<  4,
      tmp |= (z[j+0] & m4ri_one<<58) >> 29, tmp |= (z[j+1] & m4ri_one<<58) <<  3,
      tmp |= (z[j+0] & m4ri_one<<60) >> 30, tmp |= (z[j+1] & m4ri_one<<60) <<  2,
      tmp |= (z[j+0] & m4ri_one<<62) >> 31, tmp |= (z[j+1] & m4ri_one<<62) <<  1;
      a0[j2] = tmp;
    }
    
    /* deal with the tail */
    if(j+2 == Z->x->width) { /* we have to deal with two words */
      tmp = 0;
      tmp |= (z[j] & m4ri_one<< 0) >>  0;
      tmp |= (z[j] & m4ri_one<< 2) >>  1;
      tmp |= (z[j] & m4ri_one<< 4) >>  2;
      tmp |= (z[j] & m4ri_one<< 6) >>  3;
      tmp |= (z[j] & m4ri_one<< 8) >>  4;
      tmp |= (z[j] & m4ri_one<<10) >>  5;
      tmp |= (z[j] & m4ri_one<<12) >>  6;
      tmp |= (z[j] & m4ri_one<<14) >>  7;
      tmp |= (z[j] & m4ri_one<<16) >>  8;
      tmp |= (z[j] & m4ri_one<<18) >>  9;
      tmp |= (z[j] & m4ri_one<<20) >> 10;
      tmp |= (z[j] & m4ri_one<<22) >> 11;
      tmp |= (z[j] & m4ri_one<<24) >> 12;
      tmp |= (z[j] & m4ri_one<<26) >> 13;
      tmp |= (z[j] & m4ri_one<<28) >> 14;
      tmp |= (z[j] & m4ri_one<<30) >> 15;
      tmp |= (z[j] & m4ri_one<<32) >> 16;
      tmp |= (z[j] & m4ri_one<<34) >> 17;
      tmp |= (z[j] & m4ri_one<<36) >> 18;
      tmp |= (z[j] & m4ri_one<<38) >> 19;
      tmp |= (z[j] & m4ri_one<<40) >> 20;
      tmp |= (z[j] & m4ri_one<<42) >> 21;
      tmp |= (z[j] & m4ri_one<<44) >> 22;
      tmp |= (z[j] & m4ri_one<<46) >> 23;
      tmp |= (z[j] & m4ri_one<<48) >> 24;
      tmp |= (z[j] & m4ri_one<<50) >> 25;
      tmp |= (z[j] & m4ri_one<<52) >> 26;
      tmp |= (z[j] & m4ri_one<<54) >> 27;
      tmp |= (z[j] & m4ri_one<<56) >> 28;
      tmp |= (z[j] & m4ri_one<<58) >> 29;
      tmp |= (z[j] & m4ri_one<<60) >> 30;
      tmp |= (z[j] & m4ri_one<<62) >> 31;
      a0[j2] = tmp;

      switch((2*A->ncols) % m4ri_radix) {
      case  0:      a0[j2] |= (z[j+1] & m4ri_one<< 62) <<  1;
      case 62:      a0[j2] |= (z[j+1] & m4ri_one<< 60) <<  2;
      case 60:      a0[j2] |= (z[j+1] & m4ri_one<< 58) <<  3;
      case 58:      a0[j2] |= (z[j+1] & m4ri_one<< 56) <<  4;
      case 56:      a0[j2] |= (z[j+1] & m4ri_one<< 54) <<  5;
      case 54:      a0[j2] |= (z[j+1] & m4ri_one<< 52) <<  6;
      case 52:      a0[j2] |= (z[j+1] & m4ri_one<< 50) <<  7;
      case 50:      a0[j2] |= (z[j+1] & m4ri_one<< 48) <<  8;
      case 48:      a0[j2] |= (z[j+1] & m4ri_one<< 46) <<  9;
      case 46:      a0[j2] |= (z[j+1] & m4ri_one<< 44) << 10;
      case 44:      a0[j2] |= (z[j+1] & m4ri_one<< 42) << 11;
      case 42:      a0[j2] |= (z[j+1] & m4ri_one<< 40) << 12;
      case 40:      a0[j2] |= (z[j+1] & m4ri_one<< 38) << 13;
      case 38:      a0[j2] |= (z[j+1] & m4ri_one<< 36) << 14;
      case 36:      a0[j2] |= (z[j+1] & m4ri_one<< 34) << 15;
      case 34:      a0[j2] |= (z[j+1] & m4ri_one<< 32) << 16;
      case 32:      a0[j2] |= (z[j+1] & m4ri_one<< 30) << 17;
      case 30:      a0[j2] |= (z[j+1] & m4ri_one<< 28) << 18;
      case 28:      a0[j2] |= (z[j+1] & m4ri_one<< 26) << 19;
      case 26:      a0[j2] |= (z[j+1] & m4ri_one<< 24) << 20;
      case 24:      a0[j2] |= (z[j+1] & m4ri_one<< 22) << 21;
      case 22:      a0[j2] |= (z[j+1] & m4ri_one<< 20) << 22;
      case 20:      a0[j2] |= (z[j+1] & m4ri_one<< 18) << 23;
      case 18:      a0[j2] |= (z[j+1] & m4ri_one<< 16) << 24;
      case 16:      a0[j2] |= (z[j+1] & m4ri_one<< 14) << 25;
      case 14:      a0[j2] |= (z[j+1] & m4ri_one<< 12) << 26;
      case 12:      a0[j2] |= (z[j+1] & m4ri_one<< 10) << 27;
      case 10:      a0[j2] |= (z[j+1] & m4ri_one<<  8) << 28;
      case  8:      a0[j2] |= (z[j+1] & m4ri_one<<  6) << 29;
      case  6:      a0[j2] |= (z[j+1] & m4ri_one<<  4) << 30;
      case  4:      a0[j2] |= (z[j+1] & m4ri_one<<  2) << 31;
      case  2:      a0[j2] |= (z[j+1] & m4ri_one<<  0) << 32;
      }

    } else { /* only one word */
      a0[j2] = 0;
      switch((2*A->ncols) % m4ri_radix) {
      case  0:      a0[j2] |= (z[j] & m4ri_one<<62) >> 31;
      case 62:      a0[j2] |= (z[j] & m4ri_one<<60) >> 30;
      case 60:      a0[j2] |= (z[j] & m4ri_one<<58) >> 29;
      case 58:      a0[j2] |= (z[j] & m4ri_one<<56) >> 28;
      case 56:      a0[j2] |= (z[j] & m4ri_one<<54) >> 27;
      case 54:      a0[j2] |= (z[j] & m4ri_one<<52) >> 26;
      case 52:      a0[j2] |= (z[j] & m4ri_one<<50) >> 25;
      case 50:      a0[j2] |= (z[j] & m4ri_one<<48) >> 24;
      case 48:      a0[j2] |= (z[j] & m4ri_one<<46) >> 23;
      case 46:      a0[j2] |= (z[j] & m4ri_one<<44) >> 22;
      case 44:      a0[j2] |= (z[j] & m4ri_one<<42) >> 21;
      case 42:      a0[j2] |= (z[j] & m4ri_one<<40) >> 20;
      case 40:      a0[j2] |= (z[j] & m4ri_one<<38) >> 19;
      case 38:      a0[j2] |= (z[j] & m4ri_one<<36) >> 18;
      case 36:      a0[j2] |= (z[j] & m4ri_one<<34) >> 17;
      case 34:      a0[j2] |= (z[j] & m4ri_one<<32) >> 16;
      case 32:      a0[j2] |= (z[j] & m4ri_one<<30) >> 15;
      case 30:      a0[j2] |= (z[j] & m4ri_one<<28) >> 14;
      case 28:      a0[j2] |= (z[j] & m4ri_one<<26) >> 13;
      case 26:      a0[j2] |= (z[j] & m4ri_one<<24) >> 12;
      case 24:      a0[j2] |= (z[j] & m4ri_one<<22) >> 11;
      case 22:      a0[j2] |= (z[j] & m4ri_one<<20) >> 10;
      case 20:      a0[j2] |= (z[j] & m4ri_one<<18) >>  9;
      case 18:      a0[j2] |= (z[j] & m4ri_one<<16) >>  8;
      case 16:      a0[j2] |= (z[j] & m4ri_one<<14) >>  7;
      case 14:      a0[j2] |= (z[j] & m4ri_one<<12) >>  6;
      case 12:      a0[j2] |= (z[j] & m4ri_one<<10) >>  5;
      case 10:      a0[j2] |= (z[j] & m4ri_one<< 8) >>  4;
      case  8:      a0[j2] |= (z[j] & m4ri_one<< 6) >>  3;
      case  6:      a0[j2] |= (z[j] & m4ri_one<< 4) >>  2;
      case  4:      a0[j2] |= (z[j] & m4ri_one<< 2) >>  1;
      case  2:      a0[j2] |= (z[j] & m4ri_one<< 0) >>  0;
      }
    }
  }

  /* A1 */
  for(size_t i=0; i<A->nrows; i++) {
    word *a1 = A->x[1]->rows[i];
    const word *z  = Z->x->rows[i];    

    for(j=0, j2=0; j+2 < Z->x->width; j+=2, j2++) {
      if ( !(z[j+0] | z[j+1]) )
        continue;
      tmp = 0;
      tmp |= (z[j+1] & m4ri_one<<( 0 + 1)) << 31,   tmp |= (z[j+0] & m4ri_one<<( 0 + 1)) >>  1,
      tmp |= (z[j+1] & m4ri_one<<( 2 + 1)) << 30,   tmp |= (z[j+0] & m4ri_one<<( 2 + 1)) >>  2,
      tmp |= (z[j+1] & m4ri_one<<( 4 + 1)) << 29,   tmp |= (z[j+0] & m4ri_one<<( 4 + 1)) >>  3,
      tmp |= (z[j+1] & m4ri_one<<( 6 + 1)) << 28,   tmp |= (z[j+0] & m4ri_one<<( 6 + 1)) >>  4,
      tmp |= (z[j+1] & m4ri_one<<( 8 + 1)) << 27,   tmp |= (z[j+0] & m4ri_one<<( 8 + 1)) >>  5,
      tmp |= (z[j+1] & m4ri_one<<(10 + 1)) << 26,   tmp |= (z[j+0] & m4ri_one<<(10 + 1)) >>  6,
      tmp |= (z[j+1] & m4ri_one<<(12 + 1)) << 25,   tmp |= (z[j+0] & m4ri_one<<(12 + 1)) >>  7,
      tmp |= (z[j+1] & m4ri_one<<(14 + 1)) << 24,   tmp |= (z[j+0] & m4ri_one<<(14 + 1)) >>  8,
      tmp |= (z[j+1] & m4ri_one<<(16 + 1)) << 23,   tmp |= (z[j+0] & m4ri_one<<(16 + 1)) >>  9,
      tmp |= (z[j+1] & m4ri_one<<(18 + 1)) << 22,   tmp |= (z[j+0] & m4ri_one<<(18 + 1)) >> 10,
      tmp |= (z[j+1] & m4ri_one<<(20 + 1)) << 21,   tmp |= (z[j+0] & m4ri_one<<(20 + 1)) >> 11,
      tmp |= (z[j+1] & m4ri_one<<(22 + 1)) << 20,   tmp |= (z[j+0] & m4ri_one<<(22 + 1)) >> 12,
      tmp |= (z[j+1] & m4ri_one<<(24 + 1)) << 19,   tmp |= (z[j+0] & m4ri_one<<(24 + 1)) >> 13,
      tmp |= (z[j+1] & m4ri_one<<(26 + 1)) << 18,   tmp |= (z[j+0] & m4ri_one<<(26 + 1)) >> 14,
      tmp |= (z[j+1] & m4ri_one<<(28 + 1)) << 17,   tmp |= (z[j+0] & m4ri_one<<(28 + 1)) >> 15,
      tmp |= (z[j+1] & m4ri_one<<(30 + 1)) << 16,   tmp |= (z[j+0] & m4ri_one<<(30 + 1)) >> 16,
      tmp |= (z[j+1] & m4ri_one<<(32 + 1)) << 15,   tmp |= (z[j+0] & m4ri_one<<(32 + 1)) >> 17,
      tmp |= (z[j+1] & m4ri_one<<(34 + 1)) << 14,   tmp |= (z[j+0] & m4ri_one<<(34 + 1)) >> 18,
      tmp |= (z[j+1] & m4ri_one<<(36 + 1)) << 13,   tmp |= (z[j+0] & m4ri_one<<(36 + 1)) >> 19,
      tmp |= (z[j+1] & m4ri_one<<(38 + 1)) << 12,   tmp |= (z[j+0] & m4ri_one<<(38 + 1)) >> 20,
      tmp |= (z[j+1] & m4ri_one<<(40 + 1)) << 11,   tmp |= (z[j+0] & m4ri_one<<(40 + 1)) >> 21,
      tmp |= (z[j+1] & m4ri_one<<(42 + 1)) << 10,   tmp |= (z[j+0] & m4ri_one<<(42 + 1)) >> 22,
      tmp |= (z[j+1] & m4ri_one<<(44 + 1)) <<  9,   tmp |= (z[j+0] & m4ri_one<<(44 + 1)) >> 23,
      tmp |= (z[j+1] & m4ri_one<<(46 + 1)) <<  8,   tmp |= (z[j+0] & m4ri_one<<(46 + 1)) >> 24,
      tmp |= (z[j+1] & m4ri_one<<(48 + 1)) <<  7,   tmp |= (z[j+0] & m4ri_one<<(48 + 1)) >> 25,
      tmp |= (z[j+1] & m4ri_one<<(50 + 1)) <<  6,   tmp |= (z[j+0] & m4ri_one<<(50 + 1)) >> 26,
      tmp |= (z[j+1] & m4ri_one<<(52 + 1)) <<  5,   tmp |= (z[j+0] & m4ri_one<<(52 + 1)) >> 27,
      tmp |= (z[j+1] & m4ri_one<<(54 + 1)) <<  4,   tmp |= (z[j+0] & m4ri_one<<(54 + 1)) >> 28,
      tmp |= (z[j+1] & m4ri_one<<(56 + 1)) <<  3,   tmp |= (z[j+0] & m4ri_one<<(56 + 1)) >> 29,
      tmp |= (z[j+1] & m4ri_one<<(58 + 1)) <<  2,   tmp |= (z[j+0] & m4ri_one<<(58 + 1)) >> 30,
      tmp |= (z[j+1] & m4ri_one<<(60 + 1)) <<  1,   tmp |= (z[j+0] & m4ri_one<<(60 + 1)) >> 31,
      tmp |= (z[j+1] & m4ri_one<<(62 + 1)) <<  0,   tmp |= (z[j+0] & m4ri_one<<(62 + 1)) >> 32;
      a1[j2] = tmp;
    }
    /* deal with the tail */
    if(j+2 == Z->x->width) { /* we have to deal with two words */
      tmp = 0;
      tmp |= (z[j] & m4ri_one<<( 0 + 1)) >>  1;
      tmp |= (z[j] & m4ri_one<<( 2 + 1)) >>  2;
      tmp |= (z[j] & m4ri_one<<( 4 + 1)) >>  3;
      tmp |= (z[j] & m4ri_one<<( 6 + 1)) >>  4;
      tmp |= (z[j] & m4ri_one<<( 8 + 1)) >>  5;
      tmp |= (z[j] & m4ri_one<<(10 + 1)) >>  6;
      tmp |= (z[j] & m4ri_one<<(12 + 1)) >>  7;
      tmp |= (z[j] & m4ri_one<<(14 + 1)) >>  8;
      tmp |= (z[j] & m4ri_one<<(16 + 1)) >>  9;
      tmp |= (z[j] & m4ri_one<<(18 + 1)) >> 10;
      tmp |= (z[j] & m4ri_one<<(20 + 1)) >> 11;
      tmp |= (z[j] & m4ri_one<<(22 + 1)) >> 12;
      tmp |= (z[j] & m4ri_one<<(24 + 1)) >> 13;
      tmp |= (z[j] & m4ri_one<<(26 + 1)) >> 14;
      tmp |= (z[j] & m4ri_one<<(28 + 1)) >> 15;
      tmp |= (z[j] & m4ri_one<<(30 + 1)) >> 16;
      tmp |= (z[j] & m4ri_one<<(32 + 1)) >> 17;
      tmp |= (z[j] & m4ri_one<<(34 + 1)) >> 18;
      tmp |= (z[j] & m4ri_one<<(36 + 1)) >> 19;
      tmp |= (z[j] & m4ri_one<<(38 + 1)) >> 20;
      tmp |= (z[j] & m4ri_one<<(40 + 1)) >> 21;
      tmp |= (z[j] & m4ri_one<<(42 + 1)) >> 22;
      tmp |= (z[j] & m4ri_one<<(44 + 1)) >> 23;
      tmp |= (z[j] & m4ri_one<<(46 + 1)) >> 24;
      tmp |= (z[j] & m4ri_one<<(48 + 1)) >> 25;
      tmp |= (z[j] & m4ri_one<<(50 + 1)) >> 26;
      tmp |= (z[j] & m4ri_one<<(52 + 1)) >> 27;
      tmp |= (z[j] & m4ri_one<<(54 + 1)) >> 28;
      tmp |= (z[j] & m4ri_one<<(56 + 1)) >> 29;
      tmp |= (z[j] & m4ri_one<<(58 + 1)) >> 30;
      tmp |= (z[j] & m4ri_one<<(60 + 1)) >> 31;
      tmp |= (z[j] & m4ri_one<<(62 + 1)) >> 32;
      a1[j2] = tmp;

      switch((2*A->ncols) % m4ri_radix) {
      case  0:      a1[j2] |= (z[j+1] & (m4ri_one<<(62 + 1))) <<  0;
      case 62:      a1[j2] |= (z[j+1] & (m4ri_one<<(60 + 1))) <<  1;
      case 60:      a1[j2] |= (z[j+1] & (m4ri_one<<(58 + 1))) <<  2;
      case 58:      a1[j2] |= (z[j+1] & (m4ri_one<<(56 + 1))) <<  3;
      case 56:      a1[j2] |= (z[j+1] & (m4ri_one<<(54 + 1))) <<  4;
      case 54:      a1[j2] |= (z[j+1] & (m4ri_one<<(52 + 1))) <<  5;
      case 52:      a1[j2] |= (z[j+1] & (m4ri_one<<(50 + 1))) <<  6;
      case 50:      a1[j2] |= (z[j+1] & (m4ri_one<<(48 + 1))) <<  7;
      case 48:      a1[j2] |= (z[j+1] & (m4ri_one<<(46 + 1))) <<  8;
      case 46:      a1[j2] |= (z[j+1] & (m4ri_one<<(44 + 1))) <<  9;
      case 44:      a1[j2] |= (z[j+1] & (m4ri_one<<(42 + 1))) << 10;
      case 42:      a1[j2] |= (z[j+1] & (m4ri_one<<(40 + 1))) << 11;
      case 40:      a1[j2] |= (z[j+1] & (m4ri_one<<(38 + 1))) << 12;
      case 38:      a1[j2] |= (z[j+1] & (m4ri_one<<(36 + 1))) << 13;
      case 36:      a1[j2] |= (z[j+1] & (m4ri_one<<(34 + 1))) << 14;
      case 34:      a1[j2] |= (z[j+1] & (m4ri_one<<(32 + 1))) << 15;
      case 32:      a1[j2] |= (z[j+1] & (m4ri_one<<(30 + 1))) << 16;
      case 30:      a1[j2] |= (z[j+1] & (m4ri_one<<(28 + 1))) << 17;
      case 28:      a1[j2] |= (z[j+1] & (m4ri_one<<(26 + 1))) << 18;
      case 26:      a1[j2] |= (z[j+1] & (m4ri_one<<(24 + 1))) << 19;
      case 24:      a1[j2] |= (z[j+1] & (m4ri_one<<(22 + 1))) << 20;
      case 22:      a1[j2] |= (z[j+1] & (m4ri_one<<(20 + 1))) << 21;
      case 20:      a1[j2] |= (z[j+1] & (m4ri_one<<(18 + 1))) << 22;
      case 18:      a1[j2] |= (z[j+1] & (m4ri_one<<(16 + 1))) << 23;
      case 16:      a1[j2] |= (z[j+1] & (m4ri_one<<(14 + 1))) << 24;
      case 14:      a1[j2] |= (z[j+1] & (m4ri_one<<(12 + 1))) << 25;
      case 12:      a1[j2] |= (z[j+1] & (m4ri_one<<(10 + 1))) << 26;
      case 10:      a1[j2] |= (z[j+1] & (m4ri_one<<( 8 + 1))) << 27;
      case  8:      a1[j2] |= (z[j+1] & (m4ri_one<<( 6 + 1))) << 28;
      case  6:      a1[j2] |= (z[j+1] & (m4ri_one<<( 4 + 1))) << 29;
      case  4:      a1[j2] |= (z[j+1] & (m4ri_one<<( 2 + 1))) << 30;
      case  2:      a1[j2] |= (z[j+1] & (m4ri_one<<( 0 + 1))) << 31;
      }

    } else { /* only one word */
      a1[j2] = 0;
      switch((2*A->ncols) % m4ri_radix) {
      case  0:      a1[j2] |= (z[j] & (m4ri_one<<(62 + 1))) >> 32;
      case 62:      a1[j2] |= (z[j] & (m4ri_one<<(60 + 1))) >> 31;
      case 60:      a1[j2] |= (z[j] & (m4ri_one<<(58 + 1))) >> 30;
      case 58:      a1[j2] |= (z[j] & (m4ri_one<<(56 + 1))) >> 29;
      case 56:      a1[j2] |= (z[j] & (m4ri_one<<(54 + 1))) >> 28;
      case 54:      a1[j2] |= (z[j] & (m4ri_one<<(52 + 1))) >> 27;
      case 52:      a1[j2] |= (z[j] & (m4ri_one<<(50 + 1))) >> 26;
      case 50:      a1[j2] |= (z[j] & (m4ri_one<<(48 + 1))) >> 25;
      case 48:      a1[j2] |= (z[j] & (m4ri_one<<(46 + 1))) >> 24;
      case 46:      a1[j2] |= (z[j] & (m4ri_one<<(44 + 1))) >> 23;
      case 44:      a1[j2] |= (z[j] & (m4ri_one<<(42 + 1))) >> 22;
      case 42:      a1[j2] |= (z[j] & (m4ri_one<<(40 + 1))) >> 21;
      case 40:      a1[j2] |= (z[j] & (m4ri_one<<(38 + 1))) >> 20;
      case 38:      a1[j2] |= (z[j] & (m4ri_one<<(36 + 1))) >> 19;
      case 36:      a1[j2] |= (z[j] & (m4ri_one<<(34 + 1))) >> 18;
      case 34:      a1[j2] |= (z[j] & (m4ri_one<<(32 + 1))) >> 17;
      case 32:      a1[j2] |= (z[j] & (m4ri_one<<(30 + 1))) >> 16;
      case 30:      a1[j2] |= (z[j] & (m4ri_one<<(28 + 1))) >> 15;
      case 28:      a1[j2] |= (z[j] & (m4ri_one<<(26 + 1))) >> 14;
      case 26:      a1[j2] |= (z[j] & (m4ri_one<<(24 + 1))) >> 13;
      case 24:      a1[j2] |= (z[j] & (m4ri_one<<(22 + 1))) >> 12;
      case 22:      a1[j2] |= (z[j] & (m4ri_one<<(20 + 1))) >> 11;
      case 20:      a1[j2] |= (z[j] & (m4ri_one<<(18 + 1))) >> 10;
      case 18:      a1[j2] |= (z[j] & (m4ri_one<<(16 + 1))) >>  9;
      case 16:      a1[j2] |= (z[j] & (m4ri_one<<(14 + 1))) >>  8;
      case 14:      a1[j2] |= (z[j] & (m4ri_one<<(12 + 1))) >>  7;
      case 12:      a1[j2] |= (z[j] & (m4ri_one<<(10 + 1))) >>  6;
      case 10:      a1[j2] |= (z[j] & (m4ri_one<<( 8 + 1))) >>  5;
      case  8:      a1[j2] |= (z[j] & (m4ri_one<<( 6 + 1))) >>  4;
      case  6:      a1[j2] |= (z[j] & (m4ri_one<<( 4 + 1))) >>  3;
      case  4:      a1[j2] |= (z[j] & (m4ri_one<<( 2 + 1))) >>  2;
      case  2:      a1[j2] |= (z[j] & (m4ri_one<<( 0 + 1))) >>  1;
      }                                                          
    }
  }
  return A;
}

mzed_t *mzed_mul_karatsuba(mzed_t *C, const mzed_t *A, const mzed_t *B) {
  if (A->ncols != B->nrows || A->finite_field != B->finite_field) {
    m4ri_die("mzed_mul: rows, columns and fields must match.\n");
  }
  if (C != NULL) {
    if (C->finite_field != A->finite_field || C->nrows != A->nrows || C->ncols != B->ncols) {
      m4ri_die("mzed_mul: rows and columns of returned matrix must match.\n");
    }
    mzed_set_ui(C,0);
  }
  switch(A->finite_field->degree) {
  case 2:
    C =  _mzed_mul_karatsuba2(C, A, B);
    break;
  default:
    m4ri_die("mzed_mul_karatsuba: only implemented for GF(2^2)");
  }
  return C; 
}

mzed_t *mzed_addmul_karatsuba(mzed_t *C, const mzed_t *A, const mzed_t *B) {
  assert(C != NULL);

  if (A->ncols != B->nrows || A->finite_field != B->finite_field) {
    m4ri_die("mzed_mul: rows, columns and fields must match.\n");
  }
  if (C->finite_field != A->finite_field || C->nrows != A->nrows || C->ncols != B->ncols) {
    m4ri_die("mzed_mul: rows and columns of returned matrix must match.\n");
  }
  switch(A->finite_field->degree) {
  case 2:
    C =  _mzed_mul_karatsuba2(C, A, B);
    break;
  default:
    m4ri_die("mzed_mul_karatsuba: only implemented for GF(2^2)");
  }
  return C; 
}

mzed_t *_mzed_mul_karatsuba2(mzed_t *C, const mzed_t *A, const mzed_t *B) {
  mzd_slice_t *As, *Bs, *Cs;
  if (C != NULL)
    Cs = mzed_slice2(NULL, C);
  else
    Cs = mzd_slice_init(A->finite_field, A->nrows, B->ncols);

  As = mzed_slice2(NULL, A);
  Bs = mzed_slice2(NULL, B);

  /* compute */

  mzd_t *T0 = mzd_mul(NULL, As->x[1], Bs->x[1], 0);  /* A1B1 = A1*B1 */
  mzd_t *T1 = mzd_mul(NULL, As->x[0], Bs->x[0], 0);  /* A0B0 = A0*B0 */

  mzd_add(Cs->x[0], Cs->x[0], T0);
  mzd_add(Cs->x[0], Cs->x[0], T1); /*C0 += A1*B1 + A0*B0 */

  mzd_t *T2 = mzd_add(NULL, As->x[1], As->x[0]); /*T2 = A1 + A0 */
  mzd_t *T3 = mzd_add(NULL, Bs->x[1], Bs->x[0]); /*T3 = B1 + B0 */

  mzd_add(Cs->x[1], Cs->x[1], T1);
  mzd_addmul(Cs->x[1], T2, T3, 0); /* C1 += T1 + T2*T3 */

  /* pack */
  C = mzed_cling2(C, Cs);

  /* clean */

  mzd_free(T0);  mzd_free(T1);  mzd_free(T2);  mzd_free(T3);
  mzd_slice_free(As);
  mzd_slice_free(Bs);
  mzd_slice_free(Cs);

  return C;
}

void mzd_slice_set_ui(mzd_slice_t *A, word value) {
  for(int i=0; i<A->depth; i++) {
    mzd_set_ui(A->x[i], (value>>i)&1);
  }
}

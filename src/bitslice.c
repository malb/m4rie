#include "bitslice.h"

void _mzed_cling2(mzed_t *A, const mzd_t *A1, const mzd_t *A0) {
  size_t j,j2 = 0;
  register word aw0;
  register word aw1;

  /** A0 **/
  for(size_t i=0; i<A->nrows; i++) {
    word *a0 = A0->rows[i];
    word *a  = A->x->rows[i];
    for(j=0, j2=0; j+2 < A->x->width; j+=2, j2++) {
      if (!a0[j2] )
        continue;
      aw0 = a[j+0];
      aw1 = a[j+1];
      aw0 |= (a0[j2] & ONE<< 0) <<  0;    aw1 |= (a0[j2] & ONE<<32) >> 32;
      aw0 |= (a0[j2] & ONE<< 1) <<  1;    aw1 |= (a0[j2] & ONE<<33) >> 31;
      aw0 |= (a0[j2] & ONE<< 2) <<  2;    aw1 |= (a0[j2] & ONE<<34) >> 30;
      aw0 |= (a0[j2] & ONE<< 3) <<  3;    aw1 |= (a0[j2] & ONE<<35) >> 29;
      aw0 |= (a0[j2] & ONE<< 4) <<  4;    aw1 |= (a0[j2] & ONE<<36) >> 28;
      aw0 |= (a0[j2] & ONE<< 5) <<  5;    aw1 |= (a0[j2] & ONE<<37) >> 27;
      aw0 |= (a0[j2] & ONE<< 6) <<  6;    aw1 |= (a0[j2] & ONE<<38) >> 26;
      aw0 |= (a0[j2] & ONE<< 7) <<  7;    aw1 |= (a0[j2] & ONE<<39) >> 25;
      aw0 |= (a0[j2] & ONE<< 8) <<  8;    aw1 |= (a0[j2] & ONE<<40) >> 24;
      aw0 |= (a0[j2] & ONE<< 9) <<  9;    aw1 |= (a0[j2] & ONE<<41) >> 23;
      aw0 |= (a0[j2] & ONE<<10) << 10;    aw1 |= (a0[j2] & ONE<<42) >> 22;
      aw0 |= (a0[j2] & ONE<<11) << 11;    aw1 |= (a0[j2] & ONE<<43) >> 21;
      aw0 |= (a0[j2] & ONE<<12) << 12;    aw1 |= (a0[j2] & ONE<<44) >> 20;
      aw0 |= (a0[j2] & ONE<<13) << 13;    aw1 |= (a0[j2] & ONE<<45) >> 19;
      aw0 |= (a0[j2] & ONE<<14) << 14;    aw1 |= (a0[j2] & ONE<<46) >> 18;
      aw0 |= (a0[j2] & ONE<<15) << 15;    aw1 |= (a0[j2] & ONE<<47) >> 17;
      aw0 |= (a0[j2] & ONE<<16) << 16;    aw1 |= (a0[j2] & ONE<<48) >> 16;
      aw0 |= (a0[j2] & ONE<<17) << 17;    aw1 |= (a0[j2] & ONE<<49) >> 15;
      aw0 |= (a0[j2] & ONE<<18) << 18;    aw1 |= (a0[j2] & ONE<<50) >> 14;
      aw0 |= (a0[j2] & ONE<<19) << 19;    aw1 |= (a0[j2] & ONE<<51) >> 13;
      aw0 |= (a0[j2] & ONE<<20) << 20;    aw1 |= (a0[j2] & ONE<<52) >> 12;
      aw0 |= (a0[j2] & ONE<<21) << 21;    aw1 |= (a0[j2] & ONE<<53) >> 11;
      aw0 |= (a0[j2] & ONE<<22) << 22;    aw1 |= (a0[j2] & ONE<<54) >> 10;
      aw0 |= (a0[j2] & ONE<<23) << 23;    aw1 |= (a0[j2] & ONE<<55) >>  9;
      aw0 |= (a0[j2] & ONE<<24) << 24;    aw1 |= (a0[j2] & ONE<<56) >>  8;
      aw0 |= (a0[j2] & ONE<<25) << 25;    aw1 |= (a0[j2] & ONE<<57) >>  7;
      aw0 |= (a0[j2] & ONE<<26) << 26;    aw1 |= (a0[j2] & ONE<<58) >>  6;
      aw0 |= (a0[j2] & ONE<<27) << 27;    aw1 |= (a0[j2] & ONE<<59) >>  5;
      aw0 |= (a0[j2] & ONE<<28) << 28;    aw1 |= (a0[j2] & ONE<<60) >>  4;
      aw0 |= (a0[j2] & ONE<<29) << 29;    aw1 |= (a0[j2] & ONE<<61) >>  3;
      aw0 |= (a0[j2] & ONE<<30) << 30;    aw1 |= (a0[j2] & ONE<<62) >>  2;
      aw0 |= (a0[j2] & ONE<<31) << 31;    aw1 |= (a0[j2] & ONE<<63) >>  1;
      a[j+0] = aw0;
      a[j+1] = aw1;
    }

    if(j+2 == A->x->width) {  /* we have to deal with two words */
      a[j+0] |= (a0[j2] & ONE<<  0 ) <<  0;
      a[j+0] |= (a0[j2] & ONE<<  1 ) <<  1;
      a[j+0] |= (a0[j2] & ONE<<  2 ) <<  2;
      a[j+0] |= (a0[j2] & ONE<<  3 ) <<  3;
      a[j+0] |= (a0[j2] & ONE<<  4 ) <<  4;
      a[j+0] |= (a0[j2] & ONE<<  5 ) <<  5;
      a[j+0] |= (a0[j2] & ONE<<  6 ) <<  6;
      a[j+0] |= (a0[j2] & ONE<<  7 ) <<  7;
      a[j+0] |= (a0[j2] & ONE<<  8 ) <<  8;
      a[j+0] |= (a0[j2] & ONE<<  9 ) <<  9;
      a[j+0] |= (a0[j2] & ONE<< 10 ) << 10;
      a[j+0] |= (a0[j2] & ONE<< 11 ) << 11;
      a[j+0] |= (a0[j2] & ONE<< 12 ) << 12;
      a[j+0] |= (a0[j2] & ONE<< 13 ) << 13;
      a[j+0] |= (a0[j2] & ONE<< 14 ) << 14;
      a[j+0] |= (a0[j2] & ONE<< 15 ) << 15;
      a[j+0] |= (a0[j2] & ONE<< 16 ) << 16;
      a[j+0] |= (a0[j2] & ONE<< 17 ) << 17;
      a[j+0] |= (a0[j2] & ONE<< 18 ) << 18;
      a[j+0] |= (a0[j2] & ONE<< 19 ) << 19;
      a[j+0] |= (a0[j2] & ONE<< 20 ) << 20;
      a[j+0] |= (a0[j2] & ONE<< 21 ) << 21;
      a[j+0] |= (a0[j2] & ONE<< 22 ) << 22;
      a[j+0] |= (a0[j2] & ONE<< 23 ) << 23;
      a[j+0] |= (a0[j2] & ONE<< 24 ) << 24;
      a[j+0] |= (a0[j2] & ONE<< 25 ) << 25;
      a[j+0] |= (a0[j2] & ONE<< 26 ) << 26;
      a[j+0] |= (a0[j2] & ONE<< 27 ) << 27;
      a[j+0] |= (a0[j2] & ONE<< 28 ) << 28;
      a[j+0] |= (a0[j2] & ONE<< 29 ) << 29;
      a[j+0] |= (a0[j2] & ONE<< 30 ) << 30;
      a[j+0] |= (a0[j2] & ONE<< 31 ) << 31;

      switch((2*A->ncols) % RADIX) {
      case  0:      a[j+1] |= (a0[j2] & ONE<<63) >>  1;
      case 62:      a[j+1] |= (a0[j2] & ONE<<62) >>  2;
      case 60:      a[j+1] |= (a0[j2] & ONE<<61) >>  3;
      case 58:      a[j+1] |= (a0[j2] & ONE<<60) >>  4;
      case 56:      a[j+1] |= (a0[j2] & ONE<<59) >>  5;
      case 54:      a[j+1] |= (a0[j2] & ONE<<58) >>  6;
      case 52:      a[j+1] |= (a0[j2] & ONE<<57) >>  7;
      case 50:      a[j+1] |= (a0[j2] & ONE<<56) >>  8;
      case 48:      a[j+1] |= (a0[j2] & ONE<<55) >>  9;
      case 46:      a[j+1] |= (a0[j2] & ONE<<54) >> 10;
      case 44:      a[j+1] |= (a0[j2] & ONE<<53) >> 11;
      case 42:      a[j+1] |= (a0[j2] & ONE<<52) >> 12;
      case 40:      a[j+1] |= (a0[j2] & ONE<<51) >> 13;
      case 38:      a[j+1] |= (a0[j2] & ONE<<50) >> 14;
      case 36:      a[j+1] |= (a0[j2] & ONE<<49) >> 15;
      case 34:      a[j+1] |= (a0[j2] & ONE<<48) >> 16;
      case 32:      a[j+1] |= (a0[j2] & ONE<<47) >> 17;
      case 30:      a[j+1] |= (a0[j2] & ONE<<46) >> 18;
      case 28:      a[j+1] |= (a0[j2] & ONE<<45) >> 19;
      case 26:      a[j+1] |= (a0[j2] & ONE<<44) >> 20;
      case 24:      a[j+1] |= (a0[j2] & ONE<<43) >> 21;
      case 22:      a[j+1] |= (a0[j2] & ONE<<42) >> 22;
      case 20:      a[j+1] |= (a0[j2] & ONE<<41) >> 23;
      case 18:      a[j+1] |= (a0[j2] & ONE<<40) >> 24;
      case 16:      a[j+1] |= (a0[j2] & ONE<<39) >> 25;
      case 14:      a[j+1] |= (a0[j2] & ONE<<38) >> 26;
      case 12:      a[j+1] |= (a0[j2] & ONE<<37) >> 27;
      case 10:      a[j+1] |= (a0[j2] & ONE<<36) >> 28;
      case  8:      a[j+1] |= (a0[j2] & ONE<<35) >> 29;
      case  6:      a[j+1] |= (a0[j2] & ONE<<34) >> 30;
      case  4:      a[j+1] |= (a0[j2] & ONE<<33) >> 31;
      case  2:      a[j+1] |= (a0[j2] & ONE<<32) >> 32;
      }
    
    } else {  /* only one word */
      switch((2*A->ncols) % RADIX) {
      case  0:      a[j+0] |= (a0[j2] & ONE<<31) << 31;
      case 62:      a[j+0] |= (a0[j2] & ONE<<30) << 30;
      case 60:      a[j+0] |= (a0[j2] & ONE<<29) << 29;
      case 58:      a[j+0] |= (a0[j2] & ONE<<28) << 28;
      case 56:      a[j+0] |= (a0[j2] & ONE<<27) << 27;
      case 54:      a[j+0] |= (a0[j2] & ONE<<26) << 26;
      case 52:      a[j+0] |= (a0[j2] & ONE<<25) << 25;
      case 50:      a[j+0] |= (a0[j2] & ONE<<24) << 24;
      case 48:      a[j+0] |= (a0[j2] & ONE<<23) << 23;
      case 46:      a[j+0] |= (a0[j2] & ONE<<22) << 22;
      case 44:      a[j+0] |= (a0[j2] & ONE<<21) << 21;
      case 42:      a[j+0] |= (a0[j2] & ONE<<20) << 20;
      case 40:      a[j+0] |= (a0[j2] & ONE<<19) << 19;
      case 38:      a[j+0] |= (a0[j2] & ONE<<18) << 18;
      case 36:      a[j+0] |= (a0[j2] & ONE<<17) << 17;
      case 34:      a[j+0] |= (a0[j2] & ONE<<16) << 16;
      case 32:      a[j+0] |= (a0[j2] & ONE<<15) << 15;
      case 30:      a[j+0] |= (a0[j2] & ONE<<14) << 14;
      case 28:      a[j+0] |= (a0[j2] & ONE<<13) << 13;
      case 26:      a[j+0] |= (a0[j2] & ONE<<12) << 12;
      case 24:      a[j+0] |= (a0[j2] & ONE<<11) << 11;
      case 22:      a[j+0] |= (a0[j2] & ONE<<10) << 10;
      case 20:      a[j+0] |= (a0[j2] & ONE<< 9) <<  9;
      case 18:      a[j+0] |= (a0[j2] & ONE<< 8) <<  8;
      case 16:      a[j+0] |= (a0[j2] & ONE<< 7) <<  7;
      case 14:      a[j+0] |= (a0[j2] & ONE<< 6) <<  6;
      case 12:      a[j+0] |= (a0[j2] & ONE<< 5) <<  5;
      case 10:      a[j+0] |= (a0[j2] & ONE<< 4) <<  4;
      case  8:      a[j+0] |= (a0[j2] & ONE<< 3) <<  3;
      case  6:      a[j+0] |= (a0[j2] & ONE<< 2) <<  2;
      case  4:      a[j+0] |= (a0[j2] & ONE<< 1) <<  1;
      case  2:      a[j+0] |= (a0[j2] & ONE<< 0) <<  0;
      }
    }
  }

  /** A1 **/
  for(size_t i=0; i<A->nrows; i++) {
    word *a1 = A1->rows[i];
    word *a  = A->x->rows[i];    
    for(j=0, j2=0; j+2 < A->x->width; j+=2, j2++) {
      if (!a1[j2] )
        continue;
      aw0 = a[j+0];
      aw1 = a[j+1];
      aw0 |= (a1[j2] & ONE<< 0) <<  1;    aw1 |= (a1[j2] & ONE<<32) >> 31;
      aw0 |= (a1[j2] & ONE<< 1) <<  2;    aw1 |= (a1[j2] & ONE<<33) >> 30;
      aw0 |= (a1[j2] & ONE<< 2) <<  3;    aw1 |= (a1[j2] & ONE<<34) >> 29;
      aw0 |= (a1[j2] & ONE<< 3) <<  4;    aw1 |= (a1[j2] & ONE<<35) >> 28;
      aw0 |= (a1[j2] & ONE<< 4) <<  5;    aw1 |= (a1[j2] & ONE<<36) >> 27;
      aw0 |= (a1[j2] & ONE<< 5) <<  6;    aw1 |= (a1[j2] & ONE<<37) >> 26;
      aw0 |= (a1[j2] & ONE<< 6) <<  7;    aw1 |= (a1[j2] & ONE<<38) >> 25;
      aw0 |= (a1[j2] & ONE<< 7) <<  8;    aw1 |= (a1[j2] & ONE<<39) >> 24;
      aw0 |= (a1[j2] & ONE<< 8) <<  9;    aw1 |= (a1[j2] & ONE<<40) >> 23;
      aw0 |= (a1[j2] & ONE<< 9) << 10;    aw1 |= (a1[j2] & ONE<<41) >> 22;
      aw0 |= (a1[j2] & ONE<<10) << 11;    aw1 |= (a1[j2] & ONE<<42) >> 21;
      aw0 |= (a1[j2] & ONE<<11) << 12;    aw1 |= (a1[j2] & ONE<<43) >> 20;
      aw0 |= (a1[j2] & ONE<<12) << 13;    aw1 |= (a1[j2] & ONE<<44) >> 19;
      aw0 |= (a1[j2] & ONE<<13) << 14;    aw1 |= (a1[j2] & ONE<<45) >> 18;
      aw0 |= (a1[j2] & ONE<<14) << 15;    aw1 |= (a1[j2] & ONE<<46) >> 17;
      aw0 |= (a1[j2] & ONE<<15) << 16;    aw1 |= (a1[j2] & ONE<<47) >> 16;
      aw0 |= (a1[j2] & ONE<<16) << 17;    aw1 |= (a1[j2] & ONE<<48) >> 15;
      aw0 |= (a1[j2] & ONE<<17) << 18;    aw1 |= (a1[j2] & ONE<<49) >> 14;
      aw0 |= (a1[j2] & ONE<<18) << 19;    aw1 |= (a1[j2] & ONE<<50) >> 13;
      aw0 |= (a1[j2] & ONE<<19) << 20;    aw1 |= (a1[j2] & ONE<<51) >> 12;
      aw0 |= (a1[j2] & ONE<<20) << 21;    aw1 |= (a1[j2] & ONE<<52) >> 11;
      aw0 |= (a1[j2] & ONE<<21) << 22;    aw1 |= (a1[j2] & ONE<<53) >> 10;
      aw0 |= (a1[j2] & ONE<<22) << 23;    aw1 |= (a1[j2] & ONE<<54) >>  9;
      aw0 |= (a1[j2] & ONE<<23) << 24;    aw1 |= (a1[j2] & ONE<<55) >>  8;
      aw0 |= (a1[j2] & ONE<<24) << 25;    aw1 |= (a1[j2] & ONE<<56) >>  7;
      aw0 |= (a1[j2] & ONE<<25) << 26;    aw1 |= (a1[j2] & ONE<<57) >>  6;
      aw0 |= (a1[j2] & ONE<<26) << 27;    aw1 |= (a1[j2] & ONE<<58) >>  5;
      aw0 |= (a1[j2] & ONE<<27) << 28;    aw1 |= (a1[j2] & ONE<<59) >>  4;
      aw0 |= (a1[j2] & ONE<<28) << 29;    aw1 |= (a1[j2] & ONE<<60) >>  3;
      aw0 |= (a1[j2] & ONE<<29) << 30;    aw1 |= (a1[j2] & ONE<<61) >>  2;
      aw0 |= (a1[j2] & ONE<<30) << 31;    aw1 |= (a1[j2] & ONE<<62) >>  1;
      aw0 |= (a1[j2] & ONE<<31) << 32;    aw1 |= (a1[j2] & ONE<<63) >>  0;
      a[j+0] = aw0;
      a[j+1] = aw1;
    }

    if(j+2 == A->x->width) {  /* we have to deal with two words */
      a[j+0] |= (a1[j2] & ONE<< 0) <<  1;
      a[j+0] |= (a1[j2] & ONE<< 1) <<  2;
      a[j+0] |= (a1[j2] & ONE<< 2) <<  3;
      a[j+0] |= (a1[j2] & ONE<< 3) <<  4;
      a[j+0] |= (a1[j2] & ONE<< 4) <<  5;
      a[j+0] |= (a1[j2] & ONE<< 5) <<  6;
      a[j+0] |= (a1[j2] & ONE<< 6) <<  7;
      a[j+0] |= (a1[j2] & ONE<< 7) <<  8;
      a[j+0] |= (a1[j2] & ONE<< 8) <<  9;
      a[j+0] |= (a1[j2] & ONE<< 9) << 10;
      a[j+0] |= (a1[j2] & ONE<<10) << 11;
      a[j+0] |= (a1[j2] & ONE<<11) << 12;
      a[j+0] |= (a1[j2] & ONE<<12) << 13;
      a[j+0] |= (a1[j2] & ONE<<13) << 14;
      a[j+0] |= (a1[j2] & ONE<<14) << 15;
      a[j+0] |= (a1[j2] & ONE<<15) << 16;
      a[j+0] |= (a1[j2] & ONE<<16) << 17;
      a[j+0] |= (a1[j2] & ONE<<17) << 18;
      a[j+0] |= (a1[j2] & ONE<<18) << 19;
      a[j+0] |= (a1[j2] & ONE<<19) << 20;
      a[j+0] |= (a1[j2] & ONE<<20) << 21;
      a[j+0] |= (a1[j2] & ONE<<21) << 22;
      a[j+0] |= (a1[j2] & ONE<<22) << 23;
      a[j+0] |= (a1[j2] & ONE<<23) << 24;
      a[j+0] |= (a1[j2] & ONE<<24) << 25;
      a[j+0] |= (a1[j2] & ONE<<25) << 26;
      a[j+0] |= (a1[j2] & ONE<<26) << 27;
      a[j+0] |= (a1[j2] & ONE<<27) << 28;
      a[j+0] |= (a1[j2] & ONE<<28) << 29;
      a[j+0] |= (a1[j2] & ONE<<29) << 30;
      a[j+0] |= (a1[j2] & ONE<<30) << 31;
      a[j+0] |= (a1[j2] & ONE<<31) << 32;
 
      switch((2*A->ncols) % RADIX) {
      case  0:      a[j+1] |= (a1[j2] & ONE<<63) >>  0;
      case 62:      a[j+1] |= (a1[j2] & ONE<<62) >>  1;
      case 60:      a[j+1] |= (a1[j2] & ONE<<61) >>  2;
      case 58:      a[j+1] |= (a1[j2] & ONE<<60) >>  3;
      case 56:      a[j+1] |= (a1[j2] & ONE<<59) >>  4;
      case 54:      a[j+1] |= (a1[j2] & ONE<<58) >>  5;
      case 52:      a[j+1] |= (a1[j2] & ONE<<57) >>  6;
      case 50:      a[j+1] |= (a1[j2] & ONE<<56) >>  7;
      case 48:      a[j+1] |= (a1[j2] & ONE<<55) >>  8;
      case 46:      a[j+1] |= (a1[j2] & ONE<<54) >>  9;
      case 44:      a[j+1] |= (a1[j2] & ONE<<53) >> 10;
      case 42:      a[j+1] |= (a1[j2] & ONE<<52) >> 11;
      case 40:      a[j+1] |= (a1[j2] & ONE<<51) >> 12;
      case 38:      a[j+1] |= (a1[j2] & ONE<<50) >> 13;
      case 36:      a[j+1] |= (a1[j2] & ONE<<49) >> 14;
      case 34:      a[j+1] |= (a1[j2] & ONE<<48) >> 15;
      case 32:      a[j+1] |= (a1[j2] & ONE<<47) >> 16;
      case 30:      a[j+1] |= (a1[j2] & ONE<<46) >> 17;
      case 28:      a[j+1] |= (a1[j2] & ONE<<45) >> 18;
      case 26:      a[j+1] |= (a1[j2] & ONE<<44) >> 19;
      case 24:      a[j+1] |= (a1[j2] & ONE<<43) >> 20;
      case 22:      a[j+1] |= (a1[j2] & ONE<<42) >> 21;
      case 20:      a[j+1] |= (a1[j2] & ONE<<41) >> 22;
      case 18:      a[j+1] |= (a1[j2] & ONE<<40) >> 23;
      case 16:      a[j+1] |= (a1[j2] & ONE<<39) >> 24;
      case 14:      a[j+1] |= (a1[j2] & ONE<<38) >> 25;
      case 12:      a[j+1] |= (a1[j2] & ONE<<37) >> 26;
      case 10:      a[j+1] |= (a1[j2] & ONE<<36) >> 27;
      case  8:      a[j+1] |= (a1[j2] & ONE<<35) >> 28;
      case  6:      a[j+1] |= (a1[j2] & ONE<<34) >> 29;
      case  4:      a[j+1] |= (a1[j2] & ONE<<33) >> 30;
      case  2:      a[j+1] |= (a1[j2] & ONE<<32) >> 31;
      }

    } else { /* only one word */
      switch((2*A->ncols) % RADIX) {
      case  0:      a[j+0] |= (a1[j2] & ONE<<31) << 32;
      case 62:      a[j+0] |= (a1[j2] & ONE<<30) << 31;
      case 60:      a[j+0] |= (a1[j2] & ONE<<29) << 30;
      case 58:      a[j+0] |= (a1[j2] & ONE<<28) << 29;
      case 56:      a[j+0] |= (a1[j2] & ONE<<27) << 28;
      case 54:      a[j+0] |= (a1[j2] & ONE<<26) << 27;
      case 52:      a[j+0] |= (a1[j2] & ONE<<25) << 26;
      case 50:      a[j+0] |= (a1[j2] & ONE<<24) << 25;
      case 48:      a[j+0] |= (a1[j2] & ONE<<23) << 24;
      case 46:      a[j+0] |= (a1[j2] & ONE<<22) << 23;
      case 44:      a[j+0] |= (a1[j2] & ONE<<21) << 22;
      case 42:      a[j+0] |= (a1[j2] & ONE<<20) << 21;
      case 40:      a[j+0] |= (a1[j2] & ONE<<19) << 20;
      case 38:      a[j+0] |= (a1[j2] & ONE<<18) << 19;
      case 36:      a[j+0] |= (a1[j2] & ONE<<17) << 18;
      case 34:      a[j+0] |= (a1[j2] & ONE<<16) << 17;
      case 32:      a[j+0] |= (a1[j2] & ONE<<15) << 16;
      case 30:      a[j+0] |= (a1[j2] & ONE<<14) << 15;
      case 28:      a[j+0] |= (a1[j2] & ONE<<13) << 14;
      case 26:      a[j+0] |= (a1[j2] & ONE<<12) << 13;
      case 24:      a[j+0] |= (a1[j2] & ONE<<11) << 12;
      case 22:      a[j+0] |= (a1[j2] & ONE<<10) << 11;
      case 20:      a[j+0] |= (a1[j2] & ONE<< 9) << 10;
      case 18:      a[j+0] |= (a1[j2] & ONE<< 8) <<  9;
      case 16:      a[j+0] |= (a1[j2] & ONE<< 7) <<  8;
      case 14:      a[j+0] |= (a1[j2] & ONE<< 6) <<  7;
      case 12:      a[j+0] |= (a1[j2] & ONE<< 5) <<  6;
      case 10:      a[j+0] |= (a1[j2] & ONE<< 4) <<  5;
      case  8:      a[j+0] |= (a1[j2] & ONE<< 3) <<  4;
      case  6:      a[j+0] |= (a1[j2] & ONE<< 2) <<  3;
      case  4:      a[j+0] |= (a1[j2] & ONE<< 1) <<  2;
      case  2:      a[j+0] |= (a1[j2] & ONE<< 0) <<  1;
      }
    }
  }
}

void _mzed_slice2(mzd_t *A1, mzd_t *A0, const mzed_t *A) {
  size_t j, j2 = 0;
  register word tmp = 0;

  /* A0 */
  for(size_t i=0; i<A->nrows; i++) {
    word *a0 = A0->rows[i];
    const word *a  = A->x->rows[i];    

    /* bulk of work */
    for(j=0, j2=0; j+2 < A->x->width; j+=2,j2++) {
      if ( !(a[j+0] | a[j+1]) )
        continue;
      tmp =0;
      tmp |= (a[j+0] & ONE<< 0) >>  0, tmp |= (a[j+1] & ONE<< 0) << 32,
      tmp |= (a[j+0] & ONE<< 2) >>  1, tmp |= (a[j+1] & ONE<< 2) << 31,
      tmp |= (a[j+0] & ONE<< 4) >>  2, tmp |= (a[j+1] & ONE<< 4) << 30,
      tmp |= (a[j+0] & ONE<< 6) >>  3, tmp |= (a[j+1] & ONE<< 6) << 29,
      tmp |= (a[j+0] & ONE<< 8) >>  4, tmp |= (a[j+1] & ONE<< 8) << 28,
      tmp |= (a[j+0] & ONE<<10) >>  5, tmp |= (a[j+1] & ONE<<10) << 27,
      tmp |= (a[j+0] & ONE<<12) >>  6, tmp |= (a[j+1] & ONE<<12) << 26,
      tmp |= (a[j+0] & ONE<<14) >>  7, tmp |= (a[j+1] & ONE<<14) << 25,
      tmp |= (a[j+0] & ONE<<16) >>  8, tmp |= (a[j+1] & ONE<<16) << 24,
      tmp |= (a[j+0] & ONE<<18) >>  9, tmp |= (a[j+1] & ONE<<18) << 23,
      tmp |= (a[j+0] & ONE<<20) >> 10, tmp |= (a[j+1] & ONE<<20) << 22,
      tmp |= (a[j+0] & ONE<<22) >> 11, tmp |= (a[j+1] & ONE<<22) << 21,
      tmp |= (a[j+0] & ONE<<24) >> 12, tmp |= (a[j+1] & ONE<<24) << 20,
      tmp |= (a[j+0] & ONE<<26) >> 13, tmp |= (a[j+1] & ONE<<26) << 19,
      tmp |= (a[j+0] & ONE<<28) >> 14, tmp |= (a[j+1] & ONE<<28) << 18,
      tmp |= (a[j+0] & ONE<<30) >> 15, tmp |= (a[j+1] & ONE<<30) << 17,
      tmp |= (a[j+0] & ONE<<32) >> 16, tmp |= (a[j+1] & ONE<<32) << 16,
      tmp |= (a[j+0] & ONE<<34) >> 17, tmp |= (a[j+1] & ONE<<34) << 15,
      tmp |= (a[j+0] & ONE<<36) >> 18, tmp |= (a[j+1] & ONE<<36) << 14,
      tmp |= (a[j+0] & ONE<<38) >> 19, tmp |= (a[j+1] & ONE<<38) << 13,
      tmp |= (a[j+0] & ONE<<40) >> 20, tmp |= (a[j+1] & ONE<<40) << 12,
      tmp |= (a[j+0] & ONE<<42) >> 21, tmp |= (a[j+1] & ONE<<42) << 11,
      tmp |= (a[j+0] & ONE<<44) >> 22, tmp |= (a[j+1] & ONE<<44) << 10,
      tmp |= (a[j+0] & ONE<<46) >> 23, tmp |= (a[j+1] & ONE<<46) <<  9,
      tmp |= (a[j+0] & ONE<<48) >> 24, tmp |= (a[j+1] & ONE<<48) <<  8,
      tmp |= (a[j+0] & ONE<<50) >> 25, tmp |= (a[j+1] & ONE<<50) <<  7,
      tmp |= (a[j+0] & ONE<<52) >> 26, tmp |= (a[j+1] & ONE<<52) <<  6,
      tmp |= (a[j+0] & ONE<<54) >> 27, tmp |= (a[j+1] & ONE<<54) <<  5,
      tmp |= (a[j+0] & ONE<<56) >> 28, tmp |= (a[j+1] & ONE<<56) <<  4,
      tmp |= (a[j+0] & ONE<<58) >> 29, tmp |= (a[j+1] & ONE<<58) <<  3,
      tmp |= (a[j+0] & ONE<<60) >> 30, tmp |= (a[j+1] & ONE<<60) <<  2,
      tmp |= (a[j+0] & ONE<<62) >> 31, tmp |= (a[j+1] & ONE<<62) <<  1;
      a0[j2] = tmp;
    }
    
    /* deal with the tail */
    if(j+2 == A->x->width) { /* we have to deal with two words */
      tmp = 0;
      tmp |= (a[j] & ONE<< 0) >>  0;
      tmp |= (a[j] & ONE<< 2) >>  1;
      tmp |= (a[j] & ONE<< 4) >>  2;
      tmp |= (a[j] & ONE<< 6) >>  3;
      tmp |= (a[j] & ONE<< 8) >>  4;
      tmp |= (a[j] & ONE<<10) >>  5;
      tmp |= (a[j] & ONE<<12) >>  6;
      tmp |= (a[j] & ONE<<14) >>  7;
      tmp |= (a[j] & ONE<<16) >>  8;
      tmp |= (a[j] & ONE<<18) >>  9;
      tmp |= (a[j] & ONE<<20) >> 10;
      tmp |= (a[j] & ONE<<22) >> 11;
      tmp |= (a[j] & ONE<<24) >> 12;
      tmp |= (a[j] & ONE<<26) >> 13;
      tmp |= (a[j] & ONE<<28) >> 14;
      tmp |= (a[j] & ONE<<30) >> 15;
      tmp |= (a[j] & ONE<<32) >> 16;
      tmp |= (a[j] & ONE<<34) >> 17;
      tmp |= (a[j] & ONE<<36) >> 18;
      tmp |= (a[j] & ONE<<38) >> 19;
      tmp |= (a[j] & ONE<<40) >> 20;
      tmp |= (a[j] & ONE<<42) >> 21;
      tmp |= (a[j] & ONE<<44) >> 22;
      tmp |= (a[j] & ONE<<46) >> 23;
      tmp |= (a[j] & ONE<<48) >> 24;
      tmp |= (a[j] & ONE<<50) >> 25;
      tmp |= (a[j] & ONE<<52) >> 26;
      tmp |= (a[j] & ONE<<54) >> 27;
      tmp |= (a[j] & ONE<<56) >> 28;
      tmp |= (a[j] & ONE<<58) >> 29;
      tmp |= (a[j] & ONE<<60) >> 30;
      tmp |= (a[j] & ONE<<62) >> 31;
      a0[j2] = tmp;

      switch((2*A->ncols) % RADIX) {
      case  0:      a0[j2] |= (a[j+1] & ONE<< 62) <<  1;
      case 62:      a0[j2] |= (a[j+1] & ONE<< 60) <<  2;
      case 60:      a0[j2] |= (a[j+1] & ONE<< 58) <<  3;
      case 58:      a0[j2] |= (a[j+1] & ONE<< 56) <<  4;
      case 56:      a0[j2] |= (a[j+1] & ONE<< 54) <<  5;
      case 54:      a0[j2] |= (a[j+1] & ONE<< 52) <<  6;
      case 52:      a0[j2] |= (a[j+1] & ONE<< 50) <<  7;
      case 50:      a0[j2] |= (a[j+1] & ONE<< 48) <<  8;
      case 48:      a0[j2] |= (a[j+1] & ONE<< 46) <<  9;
      case 46:      a0[j2] |= (a[j+1] & ONE<< 44) << 10;
      case 44:      a0[j2] |= (a[j+1] & ONE<< 42) << 11;
      case 42:      a0[j2] |= (a[j+1] & ONE<< 40) << 12;
      case 40:      a0[j2] |= (a[j+1] & ONE<< 38) << 13;
      case 38:      a0[j2] |= (a[j+1] & ONE<< 36) << 14;
      case 36:      a0[j2] |= (a[j+1] & ONE<< 34) << 15;
      case 34:      a0[j2] |= (a[j+1] & ONE<< 32) << 16;
      case 32:      a0[j2] |= (a[j+1] & ONE<< 30) << 17;
      case 30:      a0[j2] |= (a[j+1] & ONE<< 28) << 18;
      case 28:      a0[j2] |= (a[j+1] & ONE<< 26) << 19;
      case 26:      a0[j2] |= (a[j+1] & ONE<< 24) << 20;
      case 24:      a0[j2] |= (a[j+1] & ONE<< 22) << 21;
      case 22:      a0[j2] |= (a[j+1] & ONE<< 20) << 22;
      case 20:      a0[j2] |= (a[j+1] & ONE<< 18) << 23;
      case 18:      a0[j2] |= (a[j+1] & ONE<< 16) << 24;
      case 16:      a0[j2] |= (a[j+1] & ONE<< 14) << 25;
      case 14:      a0[j2] |= (a[j+1] & ONE<< 12) << 26;
      case 12:      a0[j2] |= (a[j+1] & ONE<< 10) << 27;
      case 10:      a0[j2] |= (a[j+1] & ONE<<  8) << 28;
      case  8:      a0[j2] |= (a[j+1] & ONE<<  6) << 29;
      case  6:      a0[j2] |= (a[j+1] & ONE<<  4) << 30;
      case  4:      a0[j2] |= (a[j+1] & ONE<<  2) << 31;
      case  2:      a0[j2] |= (a[j+1] & ONE<<  0) << 32;
      }

    } else { /* only one word */
      switch((2*A->ncols) % RADIX) {
      case  0:      a0[j2] |= (a[j] & ONE<<62) >> 31;
      case 62:      a0[j2] |= (a[j] & ONE<<60) >> 30;
      case 60:      a0[j2] |= (a[j] & ONE<<58) >> 29;
      case 58:      a0[j2] |= (a[j] & ONE<<56) >> 28;
      case 56:      a0[j2] |= (a[j] & ONE<<54) >> 27;
      case 54:      a0[j2] |= (a[j] & ONE<<52) >> 26;
      case 52:      a0[j2] |= (a[j] & ONE<<50) >> 25;
      case 50:      a0[j2] |= (a[j] & ONE<<48) >> 24;
      case 48:      a0[j2] |= (a[j] & ONE<<46) >> 23;
      case 46:      a0[j2] |= (a[j] & ONE<<44) >> 22;
      case 44:      a0[j2] |= (a[j] & ONE<<42) >> 21;
      case 42:      a0[j2] |= (a[j] & ONE<<40) >> 20;
      case 40:      a0[j2] |= (a[j] & ONE<<38) >> 19;
      case 38:      a0[j2] |= (a[j] & ONE<<36) >> 18;
      case 36:      a0[j2] |= (a[j] & ONE<<34) >> 17;
      case 34:      a0[j2] |= (a[j] & ONE<<32) >> 16;
      case 32:      a0[j2] |= (a[j] & ONE<<30) >> 15;
      case 30:      a0[j2] |= (a[j] & ONE<<28) >> 14;
      case 28:      a0[j2] |= (a[j] & ONE<<26) >> 13;
      case 26:      a0[j2] |= (a[j] & ONE<<24) >> 12;
      case 24:      a0[j2] |= (a[j] & ONE<<22) >> 11;
      case 22:      a0[j2] |= (a[j] & ONE<<20) >> 10;
      case 20:      a0[j2] |= (a[j] & ONE<<18) >>  9;
      case 18:      a0[j2] |= (a[j] & ONE<<16) >>  8;
      case 16:      a0[j2] |= (a[j] & ONE<<14) >>  7;
      case 14:      a0[j2] |= (a[j] & ONE<<12) >>  6;
      case 12:      a0[j2] |= (a[j] & ONE<<10) >>  5;
      case 10:      a0[j2] |= (a[j] & ONE<< 8) >>  4;
      case  8:      a0[j2] |= (a[j] & ONE<< 6) >>  3;
      case  6:      a0[j2] |= (a[j] & ONE<< 4) >>  2;
      case  4:      a0[j2] |= (a[j] & ONE<< 2) >>  1;
      case  2:      a0[j2] |= (a[j] & ONE<< 0) >>  0;
      }
    }
  }

  /* A1 */
  for(size_t i=0; i<A->nrows; i++) {
    word *a1 = A1->rows[i];
    const word *a  = A->x->rows[i];    

    for(j=0, j2=0; j+2 < A->x->width; j+=2, j2++) {
      if ( !(a[j+0] | a[j+1]) )
        continue;
      tmp = 0;
      tmp |= (a[j+1] & ONE<<( 0 + 1)) << 31,   tmp |= (a[j+0] & ONE<<( 0 + 1)) >>  1,
      tmp |= (a[j+1] & ONE<<( 2 + 1)) << 30,   tmp |= (a[j+0] & ONE<<( 2 + 1)) >>  2,
      tmp |= (a[j+1] & ONE<<( 4 + 1)) << 29,   tmp |= (a[j+0] & ONE<<( 4 + 1)) >>  3,
      tmp |= (a[j+1] & ONE<<( 6 + 1)) << 28,   tmp |= (a[j+0] & ONE<<( 6 + 1)) >>  4,
      tmp |= (a[j+1] & ONE<<( 8 + 1)) << 27,   tmp |= (a[j+0] & ONE<<( 8 + 1)) >>  5,
      tmp |= (a[j+1] & ONE<<(10 + 1)) << 26,   tmp |= (a[j+0] & ONE<<(10 + 1)) >>  6,
      tmp |= (a[j+1] & ONE<<(12 + 1)) << 25,   tmp |= (a[j+0] & ONE<<(12 + 1)) >>  7,
      tmp |= (a[j+1] & ONE<<(14 + 1)) << 24,   tmp |= (a[j+0] & ONE<<(14 + 1)) >>  8,
      tmp |= (a[j+1] & ONE<<(16 + 1)) << 23,   tmp |= (a[j+0] & ONE<<(16 + 1)) >>  9,
      tmp |= (a[j+1] & ONE<<(18 + 1)) << 22,   tmp |= (a[j+0] & ONE<<(18 + 1)) >> 10,
      tmp |= (a[j+1] & ONE<<(20 + 1)) << 21,   tmp |= (a[j+0] & ONE<<(20 + 1)) >> 11,
      tmp |= (a[j+1] & ONE<<(22 + 1)) << 20,   tmp |= (a[j+0] & ONE<<(22 + 1)) >> 12,
      tmp |= (a[j+1] & ONE<<(24 + 1)) << 19,   tmp |= (a[j+0] & ONE<<(24 + 1)) >> 13,
      tmp |= (a[j+1] & ONE<<(26 + 1)) << 18,   tmp |= (a[j+0] & ONE<<(26 + 1)) >> 14,
      tmp |= (a[j+1] & ONE<<(28 + 1)) << 17,   tmp |= (a[j+0] & ONE<<(28 + 1)) >> 15,
      tmp |= (a[j+1] & ONE<<(30 + 1)) << 16,   tmp |= (a[j+0] & ONE<<(30 + 1)) >> 16,
      tmp |= (a[j+1] & ONE<<(32 + 1)) << 15,   tmp |= (a[j+0] & ONE<<(32 + 1)) >> 17,
      tmp |= (a[j+1] & ONE<<(34 + 1)) << 14,   tmp |= (a[j+0] & ONE<<(34 + 1)) >> 18,
      tmp |= (a[j+1] & ONE<<(36 + 1)) << 13,   tmp |= (a[j+0] & ONE<<(36 + 1)) >> 19,
      tmp |= (a[j+1] & ONE<<(38 + 1)) << 12,   tmp |= (a[j+0] & ONE<<(38 + 1)) >> 20,
      tmp |= (a[j+1] & ONE<<(40 + 1)) << 11,   tmp |= (a[j+0] & ONE<<(40 + 1)) >> 21,
      tmp |= (a[j+1] & ONE<<(42 + 1)) << 10,   tmp |= (a[j+0] & ONE<<(42 + 1)) >> 22,
      tmp |= (a[j+1] & ONE<<(44 + 1)) <<  9,   tmp |= (a[j+0] & ONE<<(44 + 1)) >> 23,
      tmp |= (a[j+1] & ONE<<(46 + 1)) <<  8,   tmp |= (a[j+0] & ONE<<(46 + 1)) >> 24,
      tmp |= (a[j+1] & ONE<<(48 + 1)) <<  7,   tmp |= (a[j+0] & ONE<<(48 + 1)) >> 25,
      tmp |= (a[j+1] & ONE<<(50 + 1)) <<  6,   tmp |= (a[j+0] & ONE<<(50 + 1)) >> 26,
      tmp |= (a[j+1] & ONE<<(52 + 1)) <<  5,   tmp |= (a[j+0] & ONE<<(52 + 1)) >> 27,
      tmp |= (a[j+1] & ONE<<(54 + 1)) <<  4,   tmp |= (a[j+0] & ONE<<(54 + 1)) >> 28,
      tmp |= (a[j+1] & ONE<<(56 + 1)) <<  3,   tmp |= (a[j+0] & ONE<<(56 + 1)) >> 29,
      tmp |= (a[j+1] & ONE<<(58 + 1)) <<  2,   tmp |= (a[j+0] & ONE<<(58 + 1)) >> 30,
      tmp |= (a[j+1] & ONE<<(60 + 1)) <<  1,   tmp |= (a[j+0] & ONE<<(60 + 1)) >> 31,
      tmp |= (a[j+1] & ONE<<(62 + 1)) <<  0,   tmp |= (a[j+0] & ONE<<(62 + 1)) >> 32;
      a1[j2] = tmp;
    }
    /* deal with the tail */
    if(j+2 == A->x->width) { /* we have to deal with two words */
      tmp = 0;
      tmp |= (a[j] & ONE<<( 0 + 1)) >>  1;
      tmp |= (a[j] & ONE<<( 2 + 1)) >>  2;
      tmp |= (a[j] & ONE<<( 4 + 1)) >>  3;
      tmp |= (a[j] & ONE<<( 6 + 1)) >>  4;
      tmp |= (a[j] & ONE<<( 8 + 1)) >>  5;
      tmp |= (a[j] & ONE<<(10 + 1)) >>  6;
      tmp |= (a[j] & ONE<<(12 + 1)) >>  7;
      tmp |= (a[j] & ONE<<(14 + 1)) >>  8;
      tmp |= (a[j] & ONE<<(16 + 1)) >>  9;
      tmp |= (a[j] & ONE<<(18 + 1)) >> 10;
      tmp |= (a[j] & ONE<<(20 + 1)) >> 11;
      tmp |= (a[j] & ONE<<(22 + 1)) >> 12;
      tmp |= (a[j] & ONE<<(24 + 1)) >> 13;
      tmp |= (a[j] & ONE<<(26 + 1)) >> 14;
      tmp |= (a[j] & ONE<<(28 + 1)) >> 15;
      tmp |= (a[j] & ONE<<(30 + 1)) >> 16;
      tmp |= (a[j] & ONE<<(32 + 1)) >> 17;
      tmp |= (a[j] & ONE<<(34 + 1)) >> 18;
      tmp |= (a[j] & ONE<<(36 + 1)) >> 19;
      tmp |= (a[j] & ONE<<(38 + 1)) >> 20;
      tmp |= (a[j] & ONE<<(40 + 1)) >> 21;
      tmp |= (a[j] & ONE<<(42 + 1)) >> 22;
      tmp |= (a[j] & ONE<<(44 + 1)) >> 23;
      tmp |= (a[j] & ONE<<(46 + 1)) >> 24;
      tmp |= (a[j] & ONE<<(48 + 1)) >> 25;
      tmp |= (a[j] & ONE<<(50 + 1)) >> 26;
      tmp |= (a[j] & ONE<<(52 + 1)) >> 27;
      tmp |= (a[j] & ONE<<(54 + 1)) >> 28;
      tmp |= (a[j] & ONE<<(56 + 1)) >> 29;
      tmp |= (a[j] & ONE<<(58 + 1)) >> 30;
      tmp |= (a[j] & ONE<<(60 + 1)) >> 31;
      tmp |= (a[j] & ONE<<(62 + 1)) >> 32;
      a1[j2] = tmp;

      switch((2*A->ncols) % RADIX) {
      case  0:      a1[j2] |= (a[j+1] & (ONE<<(62 + 1))) <<  0;
      case 62:      a1[j2] |= (a[j+1] & (ONE<<(60 + 1))) <<  1;
      case 60:      a1[j2] |= (a[j+1] & (ONE<<(58 + 1))) <<  2;
      case 58:      a1[j2] |= (a[j+1] & (ONE<<(56 + 1))) <<  3;
      case 56:      a1[j2] |= (a[j+1] & (ONE<<(54 + 1))) <<  4;
      case 54:      a1[j2] |= (a[j+1] & (ONE<<(52 + 1))) <<  5;
      case 52:      a1[j2] |= (a[j+1] & (ONE<<(50 + 1))) <<  6;
      case 50:      a1[j2] |= (a[j+1] & (ONE<<(48 + 1))) <<  7;
      case 48:      a1[j2] |= (a[j+1] & (ONE<<(46 + 1))) <<  8;
      case 46:      a1[j2] |= (a[j+1] & (ONE<<(44 + 1))) <<  9;
      case 44:      a1[j2] |= (a[j+1] & (ONE<<(42 + 1))) << 10;
      case 42:      a1[j2] |= (a[j+1] & (ONE<<(40 + 1))) << 11;
      case 40:      a1[j2] |= (a[j+1] & (ONE<<(38 + 1))) << 12;
      case 38:      a1[j2] |= (a[j+1] & (ONE<<(36 + 1))) << 13;
      case 36:      a1[j2] |= (a[j+1] & (ONE<<(34 + 1))) << 14;
      case 34:      a1[j2] |= (a[j+1] & (ONE<<(32 + 1))) << 15;
      case 32:      a1[j2] |= (a[j+1] & (ONE<<(30 + 1))) << 16;
      case 30:      a1[j2] |= (a[j+1] & (ONE<<(28 + 1))) << 17;
      case 28:      a1[j2] |= (a[j+1] & (ONE<<(26 + 1))) << 18;
      case 26:      a1[j2] |= (a[j+1] & (ONE<<(24 + 1))) << 19;
      case 24:      a1[j2] |= (a[j+1] & (ONE<<(22 + 1))) << 20;
      case 22:      a1[j2] |= (a[j+1] & (ONE<<(20 + 1))) << 21;
      case 20:      a1[j2] |= (a[j+1] & (ONE<<(18 + 1))) << 22;
      case 18:      a1[j2] |= (a[j+1] & (ONE<<(16 + 1))) << 23;
      case 16:      a1[j2] |= (a[j+1] & (ONE<<(14 + 1))) << 24;
      case 14:      a1[j2] |= (a[j+1] & (ONE<<(12 + 1))) << 25;
      case 12:      a1[j2] |= (a[j+1] & (ONE<<(10 + 1))) << 26;
      case 10:      a1[j2] |= (a[j+1] & (ONE<<( 8 + 1))) << 27;
      case  8:      a1[j2] |= (a[j+1] & (ONE<<( 6 + 1))) << 28;
      case  6:      a1[j2] |= (a[j+1] & (ONE<<( 4 + 1))) << 29;
      case  4:      a1[j2] |= (a[j+1] & (ONE<<( 2 + 1))) << 30;
      case  2:      a1[j2] |= (a[j+1] & (ONE<<( 0 + 1))) << 31;
      }

    } else { /* only one word */
      switch((2*A->ncols) % RADIX) {
      case  0:      a1[j2] |= (a[j] & (ONE<<(62 + 1))) >> 32;
      case 62:      a1[j2] |= (a[j] & (ONE<<(60 + 1))) >> 31;
      case 60:      a1[j2] |= (a[j] & (ONE<<(58 + 1))) >> 30;
      case 58:      a1[j2] |= (a[j] & (ONE<<(56 + 1))) >> 29;
      case 56:      a1[j2] |= (a[j] & (ONE<<(54 + 1))) >> 28;
      case 54:      a1[j2] |= (a[j] & (ONE<<(52 + 1))) >> 27;
      case 52:      a1[j2] |= (a[j] & (ONE<<(50 + 1))) >> 26;
      case 50:      a1[j2] |= (a[j] & (ONE<<(48 + 1))) >> 25;
      case 48:      a1[j2] |= (a[j] & (ONE<<(46 + 1))) >> 24;
      case 46:      a1[j2] |= (a[j] & (ONE<<(44 + 1))) >> 23;
      case 44:      a1[j2] |= (a[j] & (ONE<<(42 + 1))) >> 22;
      case 42:      a1[j2] |= (a[j] & (ONE<<(40 + 1))) >> 21;
      case 40:      a1[j2] |= (a[j] & (ONE<<(38 + 1))) >> 20;
      case 38:      a1[j2] |= (a[j] & (ONE<<(36 + 1))) >> 19;
      case 36:      a1[j2] |= (a[j] & (ONE<<(34 + 1))) >> 18;
      case 34:      a1[j2] |= (a[j] & (ONE<<(32 + 1))) >> 17;
      case 32:      a1[j2] |= (a[j] & (ONE<<(30 + 1))) >> 16;
      case 30:      a1[j2] |= (a[j] & (ONE<<(28 + 1))) >> 15;
      case 28:      a1[j2] |= (a[j] & (ONE<<(26 + 1))) >> 14;
      case 26:      a1[j2] |= (a[j] & (ONE<<(24 + 1))) >> 13;
      case 24:      a1[j2] |= (a[j] & (ONE<<(22 + 1))) >> 12;
      case 22:      a1[j2] |= (a[j] & (ONE<<(20 + 1))) >> 11;
      case 20:      a1[j2] |= (a[j] & (ONE<<(18 + 1))) >> 10;
      case 18:      a1[j2] |= (a[j] & (ONE<<(16 + 1))) >>  9;
      case 16:      a1[j2] |= (a[j] & (ONE<<(14 + 1))) >>  8;
      case 14:      a1[j2] |= (a[j] & (ONE<<(12 + 1))) >>  7;
      case 12:      a1[j2] |= (a[j] & (ONE<<(10 + 1))) >>  6;
      case 10:      a1[j2] |= (a[j] & (ONE<<( 8 + 1))) >>  5;
      case  8:      a1[j2] |= (a[j] & (ONE<<( 6 + 1))) >>  4;
      case  6:      a1[j2] |= (a[j] & (ONE<<( 4 + 1))) >>  3;
      case  4:      a1[j2] |= (a[j] & (ONE<<( 2 + 1))) >>  2;
      case  2:      a1[j2] |= (a[j] & (ONE<<( 0 + 1))) >>  1;
      }                                                          
    }
  }
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

mzed_t *_mzed_mul_karatsuba2(mzed_t *C, const mzed_t *A, const mzed_t *B) {
  mzd_t *C1 = mzd_init(A->nrows, B->ncols);
  mzd_t *C0 = mzd_init(A->nrows, B->ncols);
  if (C != NULL)
    _mzed_slice2(C1, C0, C);

  mzd_t *A1 = mzd_init(A->nrows, A->ncols);
  mzd_t *A0 = mzd_init(A->nrows, A->ncols);
  _mzed_slice2(A1, A0, A);

  mzd_t *B1 = mzd_init(B->nrows, B->ncols);
  mzd_t *B0 = mzd_init(B->nrows, B->ncols);
  _mzed_slice2(B1, B0, B);

  /* compute */

  mzd_t *T0 = mzd_mul(NULL, A1, B1, 0);  /* A1B1 = A1*B1 */
  mzd_t *T1 = mzd_mul(NULL, A0, B0, 0);  /* A0B0 = A0*B0 */

  mzd_add(C0, C0, T0);
  mzd_add(C0, C0, T1); /*C0 += A1*B1 + A0*B0 */

  mzd_t *T2 = mzd_add(NULL, A1, A0); /*T2 = A1 + A0 */
  mzd_t *T3 = mzd_add(NULL, B1, B0); /*T3 = B1 + B0 */

  mzd_add(C1, C1, T1);
  mzd_addmul(C1, T2, T3, 0); /* C1 += T1 + T2*T3 */

  /* pack */
  if (C != NULL) {
    mzed_set_ui(C, 0);
  } else {
    C = mzed_init(A->finite_field, A->nrows, B->ncols);
  }
  _mzed_cling2(C, C1, C0);

  /* clean */

  mzd_free(T0);  mzd_free(T1);  mzd_free(T2);  mzd_free(T3);
  mzd_free(A1);  mzd_free(A0);  mzd_free(B1);  mzd_free(B0);
  mzd_free(C1);  mzd_free(C0);

  return C;
}


#include "bitslice.h"

void _mzed_cling2(mzed_t *A, const mzd_t *A0, const mzd_t *A1) {
  size_t j,j2 = 0;
  register word aw0;
  register word aw1;

  for(size_t i=0; i<A->nrows; i++) {
    word *a0 = A0->rows[i];
    word *a  = A->x->rows[i];
    for(j=0, j2=0; j+2 < A->x->width; j+=2, j2++) {
      if (!a0[j2] )
        continue;
      aw0 = a[j+0];
      aw1 = a[j+1];
      aw0 |= (a0[j2] & ONE<<(RADIX -  0 - 1))>> 0;    aw1 |= (a0[j2] & ONE<<(RADIX - 32 - 1))<<32;
      aw0 |= (a0[j2] & ONE<<(RADIX -  1 - 1))>> 1;    aw1 |= (a0[j2] & ONE<<(RADIX - 33 - 1))<<31;
      aw0 |= (a0[j2] & ONE<<(RADIX -  2 - 1))>> 2;    aw1 |= (a0[j2] & ONE<<(RADIX - 34 - 1))<<30;
      aw0 |= (a0[j2] & ONE<<(RADIX -  3 - 1))>> 3;    aw1 |= (a0[j2] & ONE<<(RADIX - 35 - 1))<<29;
      aw0 |= (a0[j2] & ONE<<(RADIX -  4 - 1))>> 4;    aw1 |= (a0[j2] & ONE<<(RADIX - 36 - 1))<<28;
      aw0 |= (a0[j2] & ONE<<(RADIX -  5 - 1))>> 5;    aw1 |= (a0[j2] & ONE<<(RADIX - 37 - 1))<<27;
      aw0 |= (a0[j2] & ONE<<(RADIX -  6 - 1))>> 6;    aw1 |= (a0[j2] & ONE<<(RADIX - 38 - 1))<<26;
      aw0 |= (a0[j2] & ONE<<(RADIX -  7 - 1))>> 7;    aw1 |= (a0[j2] & ONE<<(RADIX - 39 - 1))<<25;
      aw0 |= (a0[j2] & ONE<<(RADIX -  8 - 1))>> 8;    aw1 |= (a0[j2] & ONE<<(RADIX - 40 - 1))<<24;
      aw0 |= (a0[j2] & ONE<<(RADIX -  9 - 1))>> 9;    aw1 |= (a0[j2] & ONE<<(RADIX - 41 - 1))<<23;
      aw0 |= (a0[j2] & ONE<<(RADIX - 10 - 1))>>10;    aw1 |= (a0[j2] & ONE<<(RADIX - 42 - 1))<<22;
      aw0 |= (a0[j2] & ONE<<(RADIX - 11 - 1))>>11;    aw1 |= (a0[j2] & ONE<<(RADIX - 43 - 1))<<21;
      aw0 |= (a0[j2] & ONE<<(RADIX - 12 - 1))>>12;    aw1 |= (a0[j2] & ONE<<(RADIX - 44 - 1))<<20;
      aw0 |= (a0[j2] & ONE<<(RADIX - 13 - 1))>>13;    aw1 |= (a0[j2] & ONE<<(RADIX - 45 - 1))<<19;
      aw0 |= (a0[j2] & ONE<<(RADIX - 14 - 1))>>14;    aw1 |= (a0[j2] & ONE<<(RADIX - 46 - 1))<<18;
      aw0 |= (a0[j2] & ONE<<(RADIX - 15 - 1))>>15;    aw1 |= (a0[j2] & ONE<<(RADIX - 47 - 1))<<17;
      aw0 |= (a0[j2] & ONE<<(RADIX - 16 - 1))>>16;    aw1 |= (a0[j2] & ONE<<(RADIX - 48 - 1))<<16;
      aw0 |= (a0[j2] & ONE<<(RADIX - 17 - 1))>>17;    aw1 |= (a0[j2] & ONE<<(RADIX - 49 - 1))<<15;
      aw0 |= (a0[j2] & ONE<<(RADIX - 18 - 1))>>18;    aw1 |= (a0[j2] & ONE<<(RADIX - 50 - 1))<<14;
      aw0 |= (a0[j2] & ONE<<(RADIX - 19 - 1))>>19;    aw1 |= (a0[j2] & ONE<<(RADIX - 51 - 1))<<13;
      aw0 |= (a0[j2] & ONE<<(RADIX - 20 - 1))>>20;    aw1 |= (a0[j2] & ONE<<(RADIX - 52 - 1))<<12;
      aw0 |= (a0[j2] & ONE<<(RADIX - 21 - 1))>>21;    aw1 |= (a0[j2] & ONE<<(RADIX - 53 - 1))<<11;
      aw0 |= (a0[j2] & ONE<<(RADIX - 22 - 1))>>22;    aw1 |= (a0[j2] & ONE<<(RADIX - 54 - 1))<<10;
      aw0 |= (a0[j2] & ONE<<(RADIX - 23 - 1))>>23;    aw1 |= (a0[j2] & ONE<<(RADIX - 55 - 1))<< 9;
      aw0 |= (a0[j2] & ONE<<(RADIX - 24 - 1))>>24;    aw1 |= (a0[j2] & ONE<<(RADIX - 56 - 1))<< 8;
      aw0 |= (a0[j2] & ONE<<(RADIX - 25 - 1))>>25;    aw1 |= (a0[j2] & ONE<<(RADIX - 57 - 1))<< 7;
      aw0 |= (a0[j2] & ONE<<(RADIX - 26 - 1))>>26;    aw1 |= (a0[j2] & ONE<<(RADIX - 58 - 1))<< 6;
      aw0 |= (a0[j2] & ONE<<(RADIX - 27 - 1))>>27;    aw1 |= (a0[j2] & ONE<<(RADIX - 59 - 1))<< 5;
      aw0 |= (a0[j2] & ONE<<(RADIX - 28 - 1))>>28;    aw1 |= (a0[j2] & ONE<<(RADIX - 60 - 1))<< 4;
      aw0 |= (a0[j2] & ONE<<(RADIX - 29 - 1))>>29;    aw1 |= (a0[j2] & ONE<<(RADIX - 61 - 1))<< 3;
      aw0 |= (a0[j2] & ONE<<(RADIX - 30 - 1))>>30;    aw1 |= (a0[j2] & ONE<<(RADIX - 62 - 1))<< 2;
      aw0 |= (a0[j2] & ONE<<(RADIX - 31 - 1))>>31;    aw1 |= (a0[j2] & ONE<<(RADIX - 63 - 1))<< 1;
      a[j+0] = aw0;
      a[j+1] = aw1;
    }

    if(j+2 == A->x->width) {  /* we have to deal with two words */
      a[j+0] |= (a0[j2] & ONE<<(RADIX -  0 - 1))>> 0;
      a[j+0] |= (a0[j2] & ONE<<(RADIX -  1 - 1))>> 1;
      a[j+0] |= (a0[j2] & ONE<<(RADIX -  2 - 1))>> 2;
      a[j+0] |= (a0[j2] & ONE<<(RADIX -  3 - 1))>> 3;
      a[j+0] |= (a0[j2] & ONE<<(RADIX -  4 - 1))>> 4;
      a[j+0] |= (a0[j2] & ONE<<(RADIX -  5 - 1))>> 5;
      a[j+0] |= (a0[j2] & ONE<<(RADIX -  6 - 1))>> 6;
      a[j+0] |= (a0[j2] & ONE<<(RADIX -  7 - 1))>> 7;
      a[j+0] |= (a0[j2] & ONE<<(RADIX -  8 - 1))>> 8;
      a[j+0] |= (a0[j2] & ONE<<(RADIX -  9 - 1))>> 9;
      a[j+0] |= (a0[j2] & ONE<<(RADIX - 10 - 1))>>10;
      a[j+0] |= (a0[j2] & ONE<<(RADIX - 11 - 1))>>11;
      a[j+0] |= (a0[j2] & ONE<<(RADIX - 12 - 1))>>12;
      a[j+0] |= (a0[j2] & ONE<<(RADIX - 13 - 1))>>13;
      a[j+0] |= (a0[j2] & ONE<<(RADIX - 14 - 1))>>14;
      a[j+0] |= (a0[j2] & ONE<<(RADIX - 15 - 1))>>15;
      a[j+0] |= (a0[j2] & ONE<<(RADIX - 16 - 1))>>16;
      a[j+0] |= (a0[j2] & ONE<<(RADIX - 17 - 1))>>17;
      a[j+0] |= (a0[j2] & ONE<<(RADIX - 18 - 1))>>18;
      a[j+0] |= (a0[j2] & ONE<<(RADIX - 19 - 1))>>19;
      a[j+0] |= (a0[j2] & ONE<<(RADIX - 20 - 1))>>20;
      a[j+0] |= (a0[j2] & ONE<<(RADIX - 21 - 1))>>21;
      a[j+0] |= (a0[j2] & ONE<<(RADIX - 22 - 1))>>22;
      a[j+0] |= (a0[j2] & ONE<<(RADIX - 23 - 1))>>23;
      a[j+0] |= (a0[j2] & ONE<<(RADIX - 24 - 1))>>24;
      a[j+0] |= (a0[j2] & ONE<<(RADIX - 25 - 1))>>25;
      a[j+0] |= (a0[j2] & ONE<<(RADIX - 26 - 1))>>26;
      a[j+0] |= (a0[j2] & ONE<<(RADIX - 27 - 1))>>27;
      a[j+0] |= (a0[j2] & ONE<<(RADIX - 28 - 1))>>28;
      a[j+0] |= (a0[j2] & ONE<<(RADIX - 29 - 1))>>29;
      a[j+0] |= (a0[j2] & ONE<<(RADIX - 30 - 1))>>30;
      a[j+0] |= (a0[j2] & ONE<<(RADIX - 31 - 1))>>31;

      switch((2*A->ncols) % RADIX) {
      case  0:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 63 - 1))<< 1;
      case 62:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 62 - 1))<< 2;
      case 60:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 61 - 1))<< 3;
      case 58:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 60 - 1))<< 4;
      case 56:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 59 - 1))<< 5;
      case 54:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 58 - 1))<< 6;
      case 52:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 57 - 1))<< 7;
      case 50:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 56 - 1))<< 8;
      case 48:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 55 - 1))<< 9;
      case 46:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 54 - 1))<<10;
      case 44:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 53 - 1))<<11;
      case 42:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 52 - 1))<<12;
      case 40:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 51 - 1))<<13;
      case 38:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 50 - 1))<<14;
      case 36:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 49 - 1))<<15;
      case 34:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 48 - 1))<<16;
      case 32:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 47 - 1))<<17;
      case 30:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 46 - 1))<<18;
      case 28:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 45 - 1))<<19;
      case 26:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 44 - 1))<<20;
      case 24:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 43 - 1))<<21;
      case 22:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 42 - 1))<<22;
      case 20:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 41 - 1))<<23;
      case 18:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 40 - 1))<<24;
      case 16:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 39 - 1))<<25;
      case 14:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 38 - 1))<<26;
      case 12:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 37 - 1))<<27;
      case 10:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 36 - 1))<<28;
      case  8:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 35 - 1))<<29;
      case  6:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 34 - 1))<<30;
      case  4:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 33 - 1))<<31;
      case  2:      a[j+1] |= (a0[j2] & ONE<<(RADIX - 32 - 1))<<32;
      }
    
    } else {  /* only one word */
      switch((2*A->ncols) % RADIX) {
      case  0:      a[j+0] |= (a0[j2] & ONE<<(RADIX - 31 - 1))>>31;
      case 62:      a[j+0] |= (a0[j2] & ONE<<(RADIX - 30 - 1))>>30;
      case 60:      a[j+0] |= (a0[j2] & ONE<<(RADIX - 29 - 1))>>29;
      case 58:      a[j+0] |= (a0[j2] & ONE<<(RADIX - 28 - 1))>>28;
      case 56:      a[j+0] |= (a0[j2] & ONE<<(RADIX - 27 - 1))>>27;
      case 54:      a[j+0] |= (a0[j2] & ONE<<(RADIX - 26 - 1))>>26;
      case 52:      a[j+0] |= (a0[j2] & ONE<<(RADIX - 25 - 1))>>25;
      case 50:      a[j+0] |= (a0[j2] & ONE<<(RADIX - 24 - 1))>>24;
      case 48:      a[j+0] |= (a0[j2] & ONE<<(RADIX - 23 - 1))>>23;
      case 46:      a[j+0] |= (a0[j2] & ONE<<(RADIX - 22 - 1))>>22;
      case 44:      a[j+0] |= (a0[j2] & ONE<<(RADIX - 21 - 1))>>21;
      case 42:      a[j+0] |= (a0[j2] & ONE<<(RADIX - 20 - 1))>>20;
      case 40:      a[j+0] |= (a0[j2] & ONE<<(RADIX - 19 - 1))>>19;
      case 38:      a[j+0] |= (a0[j2] & ONE<<(RADIX - 18 - 1))>>18;
      case 36:      a[j+0] |= (a0[j2] & ONE<<(RADIX - 17 - 1))>>17;
      case 34:      a[j+0] |= (a0[j2] & ONE<<(RADIX - 16 - 1))>>16;
      case 32:      a[j+0] |= (a0[j2] & ONE<<(RADIX - 15 - 1))>>15;
      case 30:      a[j+0] |= (a0[j2] & ONE<<(RADIX - 14 - 1))>>14;
      case 28:      a[j+0] |= (a0[j2] & ONE<<(RADIX - 13 - 1))>>13;
      case 26:      a[j+0] |= (a0[j2] & ONE<<(RADIX - 12 - 1))>>12;
      case 24:      a[j+0] |= (a0[j2] & ONE<<(RADIX - 11 - 1))>>11;
      case 22:      a[j+0] |= (a0[j2] & ONE<<(RADIX - 10 - 1))>>10;
      case 20:      a[j+0] |= (a0[j2] & ONE<<(RADIX -  9 - 1))>> 9;
      case 18:      a[j+0] |= (a0[j2] & ONE<<(RADIX -  8 - 1))>> 8;
      case 16:      a[j+0] |= (a0[j2] & ONE<<(RADIX -  7 - 1))>> 7;
      case 14:      a[j+0] |= (a0[j2] & ONE<<(RADIX -  6 - 1))>> 6;
      case 12:      a[j+0] |= (a0[j2] & ONE<<(RADIX -  5 - 1))>> 5;
      case 10:      a[j+0] |= (a0[j2] & ONE<<(RADIX -  4 - 1))>> 4;
      case  8:      a[j+0] |= (a0[j2] & ONE<<(RADIX -  3 - 1))>> 3;
      case  6:      a[j+0] |= (a0[j2] & ONE<<(RADIX -  2 - 1))>> 2;
      case  4:      a[j+0] |= (a0[j2] & ONE<<(RADIX -  1 - 1))>> 1;
      case  2:      a[j+0] |= (a0[j2] & ONE<<(RADIX -  0 - 1))>> 0;
      }
    }
  }

  for(size_t i=0; i<A->nrows; i++) {
    word *a1 = A1->rows[i];
    word *a  = A->x->rows[i];    
    for(j=0, j2=0; j+2 < A->x->width; j+=2, j2++) {
      if (!a1[j2] )
        continue;
      aw0 = a[j+0];
      aw1 = a[j+1];
      aw0 |= (a1[j2] & ONE<<(RADIX -  0 - 1))>> 1;    aw1 |= (a1[j2] & ONE<<(RADIX - 32 - 1))<<31;
      aw0 |= (a1[j2] & ONE<<(RADIX -  1 - 1))>> 2;    aw1 |= (a1[j2] & ONE<<(RADIX - 33 - 1))<<30;
      aw0 |= (a1[j2] & ONE<<(RADIX -  2 - 1))>> 3;    aw1 |= (a1[j2] & ONE<<(RADIX - 34 - 1))<<29;
      aw0 |= (a1[j2] & ONE<<(RADIX -  3 - 1))>> 4;    aw1 |= (a1[j2] & ONE<<(RADIX - 35 - 1))<<28;
      aw0 |= (a1[j2] & ONE<<(RADIX -  4 - 1))>> 5;    aw1 |= (a1[j2] & ONE<<(RADIX - 36 - 1))<<27;
      aw0 |= (a1[j2] & ONE<<(RADIX -  5 - 1))>> 6;    aw1 |= (a1[j2] & ONE<<(RADIX - 37 - 1))<<26;
      aw0 |= (a1[j2] & ONE<<(RADIX -  6 - 1))>> 7;    aw1 |= (a1[j2] & ONE<<(RADIX - 38 - 1))<<25;
      aw0 |= (a1[j2] & ONE<<(RADIX -  7 - 1))>> 8;    aw1 |= (a1[j2] & ONE<<(RADIX - 39 - 1))<<24;
      aw0 |= (a1[j2] & ONE<<(RADIX -  8 - 1))>> 9;    aw1 |= (a1[j2] & ONE<<(RADIX - 40 - 1))<<23;
      aw0 |= (a1[j2] & ONE<<(RADIX -  9 - 1))>>10;    aw1 |= (a1[j2] & ONE<<(RADIX - 41 - 1))<<22;
      aw0 |= (a1[j2] & ONE<<(RADIX - 10 - 1))>>11;    aw1 |= (a1[j2] & ONE<<(RADIX - 42 - 1))<<21;
      aw0 |= (a1[j2] & ONE<<(RADIX - 11 - 1))>>12;    aw1 |= (a1[j2] & ONE<<(RADIX - 43 - 1))<<20;
      aw0 |= (a1[j2] & ONE<<(RADIX - 12 - 1))>>13;    aw1 |= (a1[j2] & ONE<<(RADIX - 44 - 1))<<19;
      aw0 |= (a1[j2] & ONE<<(RADIX - 13 - 1))>>14;    aw1 |= (a1[j2] & ONE<<(RADIX - 45 - 1))<<18;
      aw0 |= (a1[j2] & ONE<<(RADIX - 14 - 1))>>15;    aw1 |= (a1[j2] & ONE<<(RADIX - 46 - 1))<<17;
      aw0 |= (a1[j2] & ONE<<(RADIX - 15 - 1))>>16;    aw1 |= (a1[j2] & ONE<<(RADIX - 47 - 1))<<16;
      aw0 |= (a1[j2] & ONE<<(RADIX - 16 - 1))>>17;    aw1 |= (a1[j2] & ONE<<(RADIX - 48 - 1))<<15;
      aw0 |= (a1[j2] & ONE<<(RADIX - 17 - 1))>>18;    aw1 |= (a1[j2] & ONE<<(RADIX - 49 - 1))<<14;
      aw0 |= (a1[j2] & ONE<<(RADIX - 18 - 1))>>19;    aw1 |= (a1[j2] & ONE<<(RADIX - 50 - 1))<<13;
      aw0 |= (a1[j2] & ONE<<(RADIX - 19 - 1))>>20;    aw1 |= (a1[j2] & ONE<<(RADIX - 51 - 1))<<12;
      aw0 |= (a1[j2] & ONE<<(RADIX - 20 - 1))>>21;    aw1 |= (a1[j2] & ONE<<(RADIX - 52 - 1))<<11;
      aw0 |= (a1[j2] & ONE<<(RADIX - 21 - 1))>>22;    aw1 |= (a1[j2] & ONE<<(RADIX - 53 - 1))<<10;
      aw0 |= (a1[j2] & ONE<<(RADIX - 22 - 1))>>23;    aw1 |= (a1[j2] & ONE<<(RADIX - 54 - 1))<< 9;
      aw0 |= (a1[j2] & ONE<<(RADIX - 23 - 1))>>24;    aw1 |= (a1[j2] & ONE<<(RADIX - 55 - 1))<< 8;
      aw0 |= (a1[j2] & ONE<<(RADIX - 24 - 1))>>25;    aw1 |= (a1[j2] & ONE<<(RADIX - 56 - 1))<< 7;
      aw0 |= (a1[j2] & ONE<<(RADIX - 25 - 1))>>26;    aw1 |= (a1[j2] & ONE<<(RADIX - 57 - 1))<< 6;
      aw0 |= (a1[j2] & ONE<<(RADIX - 26 - 1))>>27;    aw1 |= (a1[j2] & ONE<<(RADIX - 58 - 1))<< 5;
      aw0 |= (a1[j2] & ONE<<(RADIX - 27 - 1))>>28;    aw1 |= (a1[j2] & ONE<<(RADIX - 59 - 1))<< 4;
      aw0 |= (a1[j2] & ONE<<(RADIX - 28 - 1))>>29;    aw1 |= (a1[j2] & ONE<<(RADIX - 60 - 1))<< 3;
      aw0 |= (a1[j2] & ONE<<(RADIX - 29 - 1))>>30;    aw1 |= (a1[j2] & ONE<<(RADIX - 61 - 1))<< 2;
      aw0 |= (a1[j2] & ONE<<(RADIX - 30 - 1))>>31;    aw1 |= (a1[j2] & ONE<<(RADIX - 62 - 1))<< 1;
      aw0 |= (a1[j2] & ONE<<(RADIX - 31 - 1))>>32;    aw1 |= (a1[j2] & ONE<<(RADIX - 63 - 1))<< 0;
      a[j+0] = aw0;
      a[j+1] = aw1;
    }

    if(j+2 == A->x->width) {  /* we have to deal with two words */
      a[j+0] |= (a1[j2] & ONE<<(RADIX -  0 - 1))>> 1;
      a[j+0] |= (a1[j2] & ONE<<(RADIX -  1 - 1))>> 2;
      a[j+0] |= (a1[j2] & ONE<<(RADIX -  2 - 1))>> 3;
      a[j+0] |= (a1[j2] & ONE<<(RADIX -  3 - 1))>> 4;
      a[j+0] |= (a1[j2] & ONE<<(RADIX -  4 - 1))>> 5;
      a[j+0] |= (a1[j2] & ONE<<(RADIX -  5 - 1))>> 6;
      a[j+0] |= (a1[j2] & ONE<<(RADIX -  6 - 1))>> 7;
      a[j+0] |= (a1[j2] & ONE<<(RADIX -  7 - 1))>> 8;
      a[j+0] |= (a1[j2] & ONE<<(RADIX -  8 - 1))>> 9;
      a[j+0] |= (a1[j2] & ONE<<(RADIX -  9 - 1))>>10;
      a[j+0] |= (a1[j2] & ONE<<(RADIX - 10 - 1))>>11;
      a[j+0] |= (a1[j2] & ONE<<(RADIX - 11 - 1))>>12;
      a[j+0] |= (a1[j2] & ONE<<(RADIX - 12 - 1))>>13;
      a[j+0] |= (a1[j2] & ONE<<(RADIX - 13 - 1))>>14;
      a[j+0] |= (a1[j2] & ONE<<(RADIX - 14 - 1))>>15;
      a[j+0] |= (a1[j2] & ONE<<(RADIX - 15 - 1))>>16;
      a[j+0] |= (a1[j2] & ONE<<(RADIX - 16 - 1))>>17;
      a[j+0] |= (a1[j2] & ONE<<(RADIX - 17 - 1))>>18;
      a[j+0] |= (a1[j2] & ONE<<(RADIX - 18 - 1))>>19;
      a[j+0] |= (a1[j2] & ONE<<(RADIX - 19 - 1))>>20;
      a[j+0] |= (a1[j2] & ONE<<(RADIX - 20 - 1))>>21;
      a[j+0] |= (a1[j2] & ONE<<(RADIX - 21 - 1))>>22;
      a[j+0] |= (a1[j2] & ONE<<(RADIX - 22 - 1))>>23;
      a[j+0] |= (a1[j2] & ONE<<(RADIX - 23 - 1))>>24;
      a[j+0] |= (a1[j2] & ONE<<(RADIX - 24 - 1))>>25;
      a[j+0] |= (a1[j2] & ONE<<(RADIX - 25 - 1))>>26;
      a[j+0] |= (a1[j2] & ONE<<(RADIX - 26 - 1))>>27;
      a[j+0] |= (a1[j2] & ONE<<(RADIX - 27 - 1))>>28;
      a[j+0] |= (a1[j2] & ONE<<(RADIX - 28 - 1))>>29;
      a[j+0] |= (a1[j2] & ONE<<(RADIX - 29 - 1))>>30;
      a[j+0] |= (a1[j2] & ONE<<(RADIX - 30 - 1))>>31;
      a[j+0] |= (a1[j2] & ONE<<(RADIX - 31 - 1))>>32;
 
      switch((2*A->ncols) % RADIX) {
      case  0:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 63 - 1))<< 0;
      case 62:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 62 - 1))<< 1;
      case 60:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 61 - 1))<< 2;
      case 58:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 60 - 1))<< 3;
      case 56:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 59 - 1))<< 4;
      case 54:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 58 - 1))<< 5;
      case 52:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 57 - 1))<< 6;
      case 50:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 56 - 1))<< 7;
      case 48:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 55 - 1))<< 8;
      case 46:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 54 - 1))<< 9;
      case 44:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 53 - 1))<<10;
      case 42:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 52 - 1))<<11;
      case 40:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 51 - 1))<<12;
      case 38:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 50 - 1))<<13;
      case 36:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 49 - 1))<<14;
      case 34:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 48 - 1))<<15;
      case 32:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 47 - 1))<<16;
      case 30:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 46 - 1))<<17;
      case 28:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 45 - 1))<<18;
      case 26:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 44 - 1))<<19;
      case 24:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 43 - 1))<<20;
      case 22:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 42 - 1))<<21;
      case 20:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 41 - 1))<<22;
      case 18:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 40 - 1))<<23;
      case 16:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 39 - 1))<<24;
      case 14:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 38 - 1))<<25;
      case 12:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 37 - 1))<<26;
      case 10:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 36 - 1))<<27;
      case  8:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 35 - 1))<<28;
      case  6:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 34 - 1))<<29;
      case  4:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 33 - 1))<<30;
      case  2:      a[j+1] |= (a1[j2] & ONE<<(RADIX - 32 - 1))<<31;
      }

    } else { /* only one word */
      switch((2*A->ncols) % RADIX) {
      case  0:      a[j+0] |= (a1[j2] & ONE<<(RADIX - 31 - 1))>>32;
      case 62:      a[j+0] |= (a1[j2] & ONE<<(RADIX - 30 - 1))>>31;
      case 60:      a[j+0] |= (a1[j2] & ONE<<(RADIX - 29 - 1))>>30;
      case 58:      a[j+0] |= (a1[j2] & ONE<<(RADIX - 28 - 1))>>29;
      case 56:      a[j+0] |= (a1[j2] & ONE<<(RADIX - 27 - 1))>>28;
      case 54:      a[j+0] |= (a1[j2] & ONE<<(RADIX - 26 - 1))>>27;
      case 52:      a[j+0] |= (a1[j2] & ONE<<(RADIX - 25 - 1))>>26;
      case 50:      a[j+0] |= (a1[j2] & ONE<<(RADIX - 24 - 1))>>25;
      case 48:      a[j+0] |= (a1[j2] & ONE<<(RADIX - 23 - 1))>>24;
      case 46:      a[j+0] |= (a1[j2] & ONE<<(RADIX - 22 - 1))>>23;
      case 44:      a[j+0] |= (a1[j2] & ONE<<(RADIX - 21 - 1))>>22;
      case 42:      a[j+0] |= (a1[j2] & ONE<<(RADIX - 20 - 1))>>21;
      case 40:      a[j+0] |= (a1[j2] & ONE<<(RADIX - 19 - 1))>>20;
      case 38:      a[j+0] |= (a1[j2] & ONE<<(RADIX - 18 - 1))>>19;
      case 36:      a[j+0] |= (a1[j2] & ONE<<(RADIX - 17 - 1))>>18;
      case 34:      a[j+0] |= (a1[j2] & ONE<<(RADIX - 16 - 1))>>17;
      case 32:      a[j+0] |= (a1[j2] & ONE<<(RADIX - 15 - 1))>>16;
      case 30:      a[j+0] |= (a1[j2] & ONE<<(RADIX - 14 - 1))>>15;
      case 28:      a[j+0] |= (a1[j2] & ONE<<(RADIX - 13 - 1))>>14;
      case 26:      a[j+0] |= (a1[j2] & ONE<<(RADIX - 12 - 1))>>13;
      case 24:      a[j+0] |= (a1[j2] & ONE<<(RADIX - 11 - 1))>>12;
      case 22:      a[j+0] |= (a1[j2] & ONE<<(RADIX - 10 - 1))>>11;
      case 20:      a[j+0] |= (a1[j2] & ONE<<(RADIX -  9 - 1))>>10;
      case 18:      a[j+0] |= (a1[j2] & ONE<<(RADIX -  8 - 1))>> 9;
      case 16:      a[j+0] |= (a1[j2] & ONE<<(RADIX -  7 - 1))>> 8;
      case 14:      a[j+0] |= (a1[j2] & ONE<<(RADIX -  6 - 1))>> 7;
      case 12:      a[j+0] |= (a1[j2] & ONE<<(RADIX -  5 - 1))>> 6;
      case 10:      a[j+0] |= (a1[j2] & ONE<<(RADIX -  4 - 1))>> 5;
      case  8:      a[j+0] |= (a1[j2] & ONE<<(RADIX -  3 - 1))>> 4;
      case  6:      a[j+0] |= (a1[j2] & ONE<<(RADIX -  2 - 1))>> 3;
      case  4:      a[j+0] |= (a1[j2] & ONE<<(RADIX -  1 - 1))>> 2;
      case  2:      a[j+0] |= (a1[j2] & ONE<<(RADIX -  0 - 1))>> 1;
      }
    }
  }
}

void _mzed_slice2(mzd_t *A0, mzd_t *A1, const mzed_t *A) {
  size_t j, j2 = 0;
  register word tmp = 0;

  for(size_t i=0; i<A->nrows; i++) {
    word *a0 = A0->rows[i];
    const word *a  = A->x->rows[i];    

    /* bulk of work */
    for(j=0, j2=0; j+2 < A->x->width; j+=2,j2++) {
      if ( !(a[j+0] | a[j+1]) )
        continue;
      tmp =0;
      tmp |= (a[j+0] & ONE<<(RADIX -  0 - 1))<< 0, tmp |= (a[j+1] & ONE<<(RADIX -  0 - 1))>>32,
      tmp |= (a[j+0] & ONE<<(RADIX -  2 - 1))<< 1, tmp |= (a[j+1] & ONE<<(RADIX -  2 - 1))>>31,
      tmp |= (a[j+0] & ONE<<(RADIX -  4 - 1))<< 2, tmp |= (a[j+1] & ONE<<(RADIX -  4 - 1))>>30,
      tmp |= (a[j+0] & ONE<<(RADIX -  6 - 1))<< 3, tmp |= (a[j+1] & ONE<<(RADIX -  6 - 1))>>29,
      tmp |= (a[j+0] & ONE<<(RADIX -  8 - 1))<< 4, tmp |= (a[j+1] & ONE<<(RADIX -  8 - 1))>>28,
      tmp |= (a[j+0] & ONE<<(RADIX - 10 - 1))<< 5, tmp |= (a[j+1] & ONE<<(RADIX - 10 - 1))>>27,
      tmp |= (a[j+0] & ONE<<(RADIX - 12 - 1))<< 6, tmp |= (a[j+1] & ONE<<(RADIX - 12 - 1))>>26,
      tmp |= (a[j+0] & ONE<<(RADIX - 14 - 1))<< 7, tmp |= (a[j+1] & ONE<<(RADIX - 14 - 1))>>25,
      tmp |= (a[j+0] & ONE<<(RADIX - 16 - 1))<< 8, tmp |= (a[j+1] & ONE<<(RADIX - 16 - 1))>>24,
      tmp |= (a[j+0] & ONE<<(RADIX - 18 - 1))<< 9, tmp |= (a[j+1] & ONE<<(RADIX - 18 - 1))>>23,
      tmp |= (a[j+0] & ONE<<(RADIX - 20 - 1))<<10, tmp |= (a[j+1] & ONE<<(RADIX - 20 - 1))>>22,
      tmp |= (a[j+0] & ONE<<(RADIX - 22 - 1))<<11, tmp |= (a[j+1] & ONE<<(RADIX - 22 - 1))>>21,
      tmp |= (a[j+0] & ONE<<(RADIX - 24 - 1))<<12, tmp |= (a[j+1] & ONE<<(RADIX - 24 - 1))>>20,
      tmp |= (a[j+0] & ONE<<(RADIX - 26 - 1))<<13, tmp |= (a[j+1] & ONE<<(RADIX - 26 - 1))>>19,
      tmp |= (a[j+0] & ONE<<(RADIX - 28 - 1))<<14, tmp |= (a[j+1] & ONE<<(RADIX - 28 - 1))>>18,
      tmp |= (a[j+0] & ONE<<(RADIX - 30 - 1))<<15, tmp |= (a[j+1] & ONE<<(RADIX - 30 - 1))>>17,
      tmp |= (a[j+0] & ONE<<(RADIX - 32 - 1))<<16, tmp |= (a[j+1] & ONE<<(RADIX - 32 - 1))>>16,
      tmp |= (a[j+0] & ONE<<(RADIX - 34 - 1))<<17, tmp |= (a[j+1] & ONE<<(RADIX - 34 - 1))>>15,
      tmp |= (a[j+0] & ONE<<(RADIX - 36 - 1))<<18, tmp |= (a[j+1] & ONE<<(RADIX - 36 - 1))>>14,
      tmp |= (a[j+0] & ONE<<(RADIX - 38 - 1))<<19, tmp |= (a[j+1] & ONE<<(RADIX - 38 - 1))>>13,
      tmp |= (a[j+0] & ONE<<(RADIX - 40 - 1))<<20, tmp |= (a[j+1] & ONE<<(RADIX - 40 - 1))>>12,
      tmp |= (a[j+0] & ONE<<(RADIX - 42 - 1))<<21, tmp |= (a[j+1] & ONE<<(RADIX - 42 - 1))>>11,
      tmp |= (a[j+0] & ONE<<(RADIX - 44 - 1))<<22, tmp |= (a[j+1] & ONE<<(RADIX - 44 - 1))>>10,
      tmp |= (a[j+0] & ONE<<(RADIX - 46 - 1))<<23, tmp |= (a[j+1] & ONE<<(RADIX - 46 - 1))>> 9,
      tmp |= (a[j+0] & ONE<<(RADIX - 48 - 1))<<24, tmp |= (a[j+1] & ONE<<(RADIX - 48 - 1))>> 8,
      tmp |= (a[j+0] & ONE<<(RADIX - 50 - 1))<<25, tmp |= (a[j+1] & ONE<<(RADIX - 50 - 1))>> 7,
      tmp |= (a[j+0] & ONE<<(RADIX - 52 - 1))<<26, tmp |= (a[j+1] & ONE<<(RADIX - 52 - 1))>> 6,
      tmp |= (a[j+0] & ONE<<(RADIX - 54 - 1))<<27, tmp |= (a[j+1] & ONE<<(RADIX - 54 - 1))>> 5,
      tmp |= (a[j+0] & ONE<<(RADIX - 56 - 1))<<28, tmp |= (a[j+1] & ONE<<(RADIX - 56 - 1))>> 4,
      tmp |= (a[j+0] & ONE<<(RADIX - 58 - 1))<<29, tmp |= (a[j+1] & ONE<<(RADIX - 58 - 1))>> 3,
      tmp |= (a[j+0] & ONE<<(RADIX - 60 - 1))<<30, tmp |= (a[j+1] & ONE<<(RADIX - 60 - 1))>> 2,
      tmp |= (a[j+0] & ONE<<(RADIX - 62 - 1))<<31, tmp |= (a[j+1] & ONE<<(RADIX - 62 - 1))>> 1;
      a0[j2] = tmp;
    }
    
    /* deal with the tail */
    if(j+2 == A->x->width) { /* we have to deal with two words */
      tmp = 0;
      tmp |= (a[j] & ONE<<(RADIX -  0 - 1))<< 0;
      tmp |= (a[j] & ONE<<(RADIX -  2 - 1))<< 1;
      tmp |= (a[j] & ONE<<(RADIX -  4 - 1))<< 2;
      tmp |= (a[j] & ONE<<(RADIX -  6 - 1))<< 3;
      tmp |= (a[j] & ONE<<(RADIX -  8 - 1))<< 4;
      tmp |= (a[j] & ONE<<(RADIX - 10 - 1))<< 5;
      tmp |= (a[j] & ONE<<(RADIX - 12 - 1))<< 6;
      tmp |= (a[j] & ONE<<(RADIX - 14 - 1))<< 7;
      tmp |= (a[j] & ONE<<(RADIX - 16 - 1))<< 8;
      tmp |= (a[j] & ONE<<(RADIX - 18 - 1))<< 9;
      tmp |= (a[j] & ONE<<(RADIX - 20 - 1))<<10;
      tmp |= (a[j] & ONE<<(RADIX - 22 - 1))<<11;
      tmp |= (a[j] & ONE<<(RADIX - 24 - 1))<<12;
      tmp |= (a[j] & ONE<<(RADIX - 26 - 1))<<13;
      tmp |= (a[j] & ONE<<(RADIX - 28 - 1))<<14;
      tmp |= (a[j] & ONE<<(RADIX - 30 - 1))<<15;
      tmp |= (a[j] & ONE<<(RADIX - 32 - 1))<<16;
      tmp |= (a[j] & ONE<<(RADIX - 34 - 1))<<17;
      tmp |= (a[j] & ONE<<(RADIX - 36 - 1))<<18;
      tmp |= (a[j] & ONE<<(RADIX - 38 - 1))<<19;
      tmp |= (a[j] & ONE<<(RADIX - 40 - 1))<<20;
      tmp |= (a[j] & ONE<<(RADIX - 42 - 1))<<21;
      tmp |= (a[j] & ONE<<(RADIX - 44 - 1))<<22;
      tmp |= (a[j] & ONE<<(RADIX - 46 - 1))<<23;
      tmp |= (a[j] & ONE<<(RADIX - 48 - 1))<<24;
      tmp |= (a[j] & ONE<<(RADIX - 50 - 1))<<25;
      tmp |= (a[j] & ONE<<(RADIX - 52 - 1))<<26;
      tmp |= (a[j] & ONE<<(RADIX - 54 - 1))<<27;
      tmp |= (a[j] & ONE<<(RADIX - 56 - 1))<<28;
      tmp |= (a[j] & ONE<<(RADIX - 58 - 1))<<29;
      tmp |= (a[j] & ONE<<(RADIX - 60 - 1))<<30;
      tmp |= (a[j] & ONE<<(RADIX - 62 - 1))<<31;
      a0[j2] = tmp;

      switch((2*A->ncols) % RADIX) {
      case  0:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 62 - 1))>> 1;
      case 62:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 60 - 1))>> 2;
      case 60:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 58 - 1))>> 3;
      case 58:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 56 - 1))>> 4;
      case 56:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 54 - 1))>> 5;
      case 54:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 52 - 1))>> 6;
      case 52:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 50 - 1))>> 7;
      case 50:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 48 - 1))>> 8;
      case 48:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 46 - 1))>> 9;
      case 46:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 44 - 1))>>10;
      case 44:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 42 - 1))>>11;
      case 42:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 40 - 1))>>12;
      case 40:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 38 - 1))>>13;
      case 38:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 36 - 1))>>14;
      case 36:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 34 - 1))>>15;
      case 34:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 32 - 1))>>16;
      case 32:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 30 - 1))>>17;
      case 30:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 28 - 1))>>18;
      case 28:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 26 - 1))>>19;
      case 26:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 24 - 1))>>20;
      case 24:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 22 - 1))>>21;
      case 22:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 20 - 1))>>22;
      case 20:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 18 - 1))>>23;
      case 18:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 16 - 1))>>24;
      case 16:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 14 - 1))>>25;
      case 14:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 12 - 1))>>26;
      case 12:      a0[j2] |= (a[j+1] & ONE<<(RADIX - 10 - 1))>>27;
      case 10:      a0[j2] |= (a[j+1] & ONE<<(RADIX -  8 - 1))>>28;
      case  8:      a0[j2] |= (a[j+1] & ONE<<(RADIX -  6 - 1))>>29;
      case  6:      a0[j2] |= (a[j+1] & ONE<<(RADIX -  4 - 1))>>30;
      case  4:      a0[j2] |= (a[j+1] & ONE<<(RADIX -  2 - 1))>>31;
      case  2:      a0[j2] |= (a[j+1] & ONE<<(RADIX -  0 - 1))>>32;
      }

    } else { /* only one word */
      switch((2*A->ncols) % RADIX) {
      case  0:      a0[j2] |= (a[j] & ONE<<(RADIX - 62 - 1))<<31;
      case 62:      a0[j2] |= (a[j] & ONE<<(RADIX - 60 - 1))<<30;
      case 60:      a0[j2] |= (a[j] & ONE<<(RADIX - 58 - 1))<<29;
      case 58:      a0[j2] |= (a[j] & ONE<<(RADIX - 56 - 1))<<28;
      case 56:      a0[j2] |= (a[j] & ONE<<(RADIX - 54 - 1))<<27;
      case 54:      a0[j2] |= (a[j] & ONE<<(RADIX - 52 - 1))<<26;
      case 52:      a0[j2] |= (a[j] & ONE<<(RADIX - 50 - 1))<<25;
      case 50:      a0[j2] |= (a[j] & ONE<<(RADIX - 48 - 1))<<24;
      case 48:      a0[j2] |= (a[j] & ONE<<(RADIX - 46 - 1))<<23;
      case 46:      a0[j2] |= (a[j] & ONE<<(RADIX - 44 - 1))<<22;
      case 44:      a0[j2] |= (a[j] & ONE<<(RADIX - 42 - 1))<<21;
      case 42:      a0[j2] |= (a[j] & ONE<<(RADIX - 40 - 1))<<20;
      case 40:      a0[j2] |= (a[j] & ONE<<(RADIX - 38 - 1))<<19;
      case 38:      a0[j2] |= (a[j] & ONE<<(RADIX - 36 - 1))<<18;
      case 36:      a0[j2] |= (a[j] & ONE<<(RADIX - 34 - 1))<<17;
      case 34:      a0[j2] |= (a[j] & ONE<<(RADIX - 32 - 1))<<16;
      case 32:      a0[j2] |= (a[j] & ONE<<(RADIX - 30 - 1))<<15;
      case 30:      a0[j2] |= (a[j] & ONE<<(RADIX - 28 - 1))<<14;
      case 28:      a0[j2] |= (a[j] & ONE<<(RADIX - 26 - 1))<<13;
      case 26:      a0[j2] |= (a[j] & ONE<<(RADIX - 24 - 1))<<12;
      case 24:      a0[j2] |= (a[j] & ONE<<(RADIX - 22 - 1))<<11;
      case 22:      a0[j2] |= (a[j] & ONE<<(RADIX - 20 - 1))<<10;
      case 20:      a0[j2] |= (a[j] & ONE<<(RADIX - 18 - 1))<< 9;
      case 18:      a0[j2] |= (a[j] & ONE<<(RADIX - 16 - 1))<< 8;
      case 16:      a0[j2] |= (a[j] & ONE<<(RADIX - 14 - 1))<< 7;
      case 14:      a0[j2] |= (a[j] & ONE<<(RADIX - 12 - 1))<< 6;
      case 12:      a0[j2] |= (a[j] & ONE<<(RADIX - 10 - 1))<< 5;
      case 10:      a0[j2] |= (a[j] & ONE<<(RADIX -  8 - 1))<< 4;
      case  8:      a0[j2] |= (a[j] & ONE<<(RADIX -  6 - 1))<< 3;
      case  6:      a0[j2] |= (a[j] & ONE<<(RADIX -  4 - 1))<< 2;
      case  4:      a0[j2] |= (a[j] & ONE<<(RADIX -  2 - 1))<< 1;
      case  2:      a0[j2] |= (a[j] & ONE<<(RADIX -  0 - 1))<< 0;
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
      tmp |= (a[j+1] & ONE<<(RADIX -  0 - 2))>>31,   tmp |= (a[j+0] & ONE<<(RADIX -  0 - 2))<< 1,
      tmp |= (a[j+1] & ONE<<(RADIX -  2 - 2))>>30,   tmp |= (a[j+0] & ONE<<(RADIX -  2 - 2))<< 2,
      tmp |= (a[j+1] & ONE<<(RADIX -  4 - 2))>>29,   tmp |= (a[j+0] & ONE<<(RADIX -  4 - 2))<< 3,
      tmp |= (a[j+1] & ONE<<(RADIX -  6 - 2))>>28,   tmp |= (a[j+0] & ONE<<(RADIX -  6 - 2))<< 4,
      tmp |= (a[j+1] & ONE<<(RADIX -  8 - 2))>>27,   tmp |= (a[j+0] & ONE<<(RADIX -  8 - 2))<< 5,
      tmp |= (a[j+1] & ONE<<(RADIX - 10 - 2))>>26,   tmp |= (a[j+0] & ONE<<(RADIX - 10 - 2))<< 6,
      tmp |= (a[j+1] & ONE<<(RADIX - 12 - 2))>>25,   tmp |= (a[j+0] & ONE<<(RADIX - 12 - 2))<< 7,
      tmp |= (a[j+1] & ONE<<(RADIX - 14 - 2))>>24,   tmp |= (a[j+0] & ONE<<(RADIX - 14 - 2))<< 8,
      tmp |= (a[j+1] & ONE<<(RADIX - 16 - 2))>>23,   tmp |= (a[j+0] & ONE<<(RADIX - 16 - 2))<< 9,
      tmp |= (a[j+1] & ONE<<(RADIX - 18 - 2))>>22,   tmp |= (a[j+0] & ONE<<(RADIX - 18 - 2))<<10,
      tmp |= (a[j+1] & ONE<<(RADIX - 20 - 2))>>21,   tmp |= (a[j+0] & ONE<<(RADIX - 20 - 2))<<11,
      tmp |= (a[j+1] & ONE<<(RADIX - 22 - 2))>>20,   tmp |= (a[j+0] & ONE<<(RADIX - 22 - 2))<<12,
      tmp |= (a[j+1] & ONE<<(RADIX - 24 - 2))>>19,   tmp |= (a[j+0] & ONE<<(RADIX - 24 - 2))<<13,
      tmp |= (a[j+1] & ONE<<(RADIX - 26 - 2))>>18,   tmp |= (a[j+0] & ONE<<(RADIX - 26 - 2))<<14,
      tmp |= (a[j+1] & ONE<<(RADIX - 28 - 2))>>17,   tmp |= (a[j+0] & ONE<<(RADIX - 28 - 2))<<15,
      tmp |= (a[j+1] & ONE<<(RADIX - 30 - 2))>>16,   tmp |= (a[j+0] & ONE<<(RADIX - 30 - 2))<<16,
      tmp |= (a[j+1] & ONE<<(RADIX - 32 - 2))>>15,   tmp |= (a[j+0] & ONE<<(RADIX - 32 - 2))<<17,
      tmp |= (a[j+1] & ONE<<(RADIX - 34 - 2))>>14,   tmp |= (a[j+0] & ONE<<(RADIX - 34 - 2))<<18,
      tmp |= (a[j+1] & ONE<<(RADIX - 36 - 2))>>13,   tmp |= (a[j+0] & ONE<<(RADIX - 36 - 2))<<19,
      tmp |= (a[j+1] & ONE<<(RADIX - 38 - 2))>>12,   tmp |= (a[j+0] & ONE<<(RADIX - 38 - 2))<<20,
      tmp |= (a[j+1] & ONE<<(RADIX - 40 - 2))>>11,   tmp |= (a[j+0] & ONE<<(RADIX - 40 - 2))<<21,
      tmp |= (a[j+1] & ONE<<(RADIX - 42 - 2))>>10,   tmp |= (a[j+0] & ONE<<(RADIX - 42 - 2))<<22,
      tmp |= (a[j+1] & ONE<<(RADIX - 44 - 2))>> 9,   tmp |= (a[j+0] & ONE<<(RADIX - 44 - 2))<<23,
      tmp |= (a[j+1] & ONE<<(RADIX - 46 - 2))>> 8,   tmp |= (a[j+0] & ONE<<(RADIX - 46 - 2))<<24,
      tmp |= (a[j+1] & ONE<<(RADIX - 48 - 2))>> 7,   tmp |= (a[j+0] & ONE<<(RADIX - 48 - 2))<<25,
      tmp |= (a[j+1] & ONE<<(RADIX - 50 - 2))>> 6,   tmp |= (a[j+0] & ONE<<(RADIX - 50 - 2))<<26,
      tmp |= (a[j+1] & ONE<<(RADIX - 52 - 2))>> 5,   tmp |= (a[j+0] & ONE<<(RADIX - 52 - 2))<<27,
      tmp |= (a[j+1] & ONE<<(RADIX - 54 - 2))>> 4,   tmp |= (a[j+0] & ONE<<(RADIX - 54 - 2))<<28,
      tmp |= (a[j+1] & ONE<<(RADIX - 56 - 2))>> 3,   tmp |= (a[j+0] & ONE<<(RADIX - 56 - 2))<<29,
      tmp |= (a[j+1] & ONE<<(RADIX - 58 - 2))>> 2,   tmp |= (a[j+0] & ONE<<(RADIX - 58 - 2))<<30,
      tmp |= (a[j+1] & ONE<<(RADIX - 60 - 2))>> 1,   tmp |= (a[j+0] & ONE<<(RADIX - 60 - 2))<<31,
      tmp |= (a[j+1] & ONE<<(RADIX - 62 - 2))>> 0,   tmp |= (a[j+0] & ONE<<(RADIX - 62 - 2))<<32;
      a1[j2] = tmp;
    }
    /* deal with the tail */
    if(j+2 == A->x->width) { /* we have to deal with two words */
      tmp = 0;
      tmp |= (a[j] & ONE<<(RADIX -  0 - 2))<< 1;
      tmp |= (a[j] & ONE<<(RADIX -  2 - 2))<< 2;
      tmp |= (a[j] & ONE<<(RADIX -  4 - 2))<< 3;
      tmp |= (a[j] & ONE<<(RADIX -  6 - 2))<< 4;
      tmp |= (a[j] & ONE<<(RADIX -  8 - 2))<< 5;
      tmp |= (a[j] & ONE<<(RADIX - 10 - 2))<< 6;
      tmp |= (a[j] & ONE<<(RADIX - 12 - 2))<< 7;
      tmp |= (a[j] & ONE<<(RADIX - 14 - 2))<< 8;
      tmp |= (a[j] & ONE<<(RADIX - 16 - 2))<< 9;
      tmp |= (a[j] & ONE<<(RADIX - 18 - 2))<<10;
      tmp |= (a[j] & ONE<<(RADIX - 20 - 2))<<11;
      tmp |= (a[j] & ONE<<(RADIX - 22 - 2))<<12;
      tmp |= (a[j] & ONE<<(RADIX - 24 - 2))<<13;
      tmp |= (a[j] & ONE<<(RADIX - 26 - 2))<<14;
      tmp |= (a[j] & ONE<<(RADIX - 28 - 2))<<15;
      tmp |= (a[j] & ONE<<(RADIX - 30 - 2))<<16;
      tmp |= (a[j] & ONE<<(RADIX - 32 - 2))<<17;
      tmp |= (a[j] & ONE<<(RADIX - 34 - 2))<<18;
      tmp |= (a[j] & ONE<<(RADIX - 36 - 2))<<19;
      tmp |= (a[j] & ONE<<(RADIX - 38 - 2))<<20;
      tmp |= (a[j] & ONE<<(RADIX - 40 - 2))<<21;
      tmp |= (a[j] & ONE<<(RADIX - 42 - 2))<<22;
      tmp |= (a[j] & ONE<<(RADIX - 44 - 2))<<23;
      tmp |= (a[j] & ONE<<(RADIX - 46 - 2))<<24;
      tmp |= (a[j] & ONE<<(RADIX - 48 - 2))<<25;
      tmp |= (a[j] & ONE<<(RADIX - 50 - 2))<<26;
      tmp |= (a[j] & ONE<<(RADIX - 52 - 2))<<27;
      tmp |= (a[j] & ONE<<(RADIX - 54 - 2))<<28;
      tmp |= (a[j] & ONE<<(RADIX - 56 - 2))<<29;
      tmp |= (a[j] & ONE<<(RADIX - 58 - 2))<<30;
      tmp |= (a[j] & ONE<<(RADIX - 60 - 2))<<31;
      tmp |= (a[j] & ONE<<(RADIX - 62 - 2))<<32;
      a1[j2] = tmp;

      switch((2*A->ncols) % RADIX) {
      case  0:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 62 - 2)))>> 0;
      case 62:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 60 - 2)))>> 1;
      case 60:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 58 - 2)))>> 2;
      case 58:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 56 - 2)))>> 3;
      case 56:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 54 - 2)))>> 4;
      case 54:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 52 - 2)))>> 5;
      case 52:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 50 - 2)))>> 6;
      case 50:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 48 - 2)))>> 7;
      case 48:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 46 - 2)))>> 8;
      case 46:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 44 - 2)))>> 9;
      case 44:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 42 - 2)))>>10;
      case 42:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 40 - 2)))>>11;
      case 40:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 38 - 2)))>>12;
      case 38:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 36 - 2)))>>13;
      case 36:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 34 - 2)))>>14;
      case 34:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 32 - 2)))>>15;
      case 32:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 30 - 2)))>>16;
      case 30:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 28 - 2)))>>17;
      case 28:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 26 - 2)))>>18;
      case 26:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 24 - 2)))>>19;
      case 24:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 22 - 2)))>>20;
      case 22:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 20 - 2)))>>21;
      case 20:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 18 - 2)))>>22;
      case 18:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 16 - 2)))>>23;
      case 16:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 14 - 2)))>>24;
      case 14:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 12 - 2)))>>25;
      case 12:      a1[j2] |= (a[j+1] & (ONE<<(RADIX - 10 - 2)))>>26;
      case 10:      a1[j2] |= (a[j+1] & (ONE<<(RADIX -  8 - 2)))>>27;
      case  8:      a1[j2] |= (a[j+1] & (ONE<<(RADIX -  6 - 2)))>>28;
      case  6:      a1[j2] |= (a[j+1] & (ONE<<(RADIX -  4 - 2)))>>29;
      case  4:      a1[j2] |= (a[j+1] & (ONE<<(RADIX -  2 - 2)))>>30;
      case  2:      a1[j2] |= (a[j+1] & (ONE<<(RADIX -  0 - 2)))>>31;
      }

    } else { /* only one word */
      switch((2*A->ncols) % RADIX) {
      case  0:      a1[j2] |= (a[j] & (ONE<<(RADIX - 62 - 2)))<<32;
      case 62:      a1[j2] |= (a[j] & (ONE<<(RADIX - 60 - 2)))<<31;
      case 60:      a1[j2] |= (a[j] & (ONE<<(RADIX - 58 - 2)))<<30;
      case 58:      a1[j2] |= (a[j] & (ONE<<(RADIX - 56 - 2)))<<29;
      case 56:      a1[j2] |= (a[j] & (ONE<<(RADIX - 54 - 2)))<<28;
      case 54:      a1[j2] |= (a[j] & (ONE<<(RADIX - 52 - 2)))<<27;
      case 52:      a1[j2] |= (a[j] & (ONE<<(RADIX - 50 - 2)))<<26;
      case 50:      a1[j2] |= (a[j] & (ONE<<(RADIX - 48 - 2)))<<25;
      case 48:      a1[j2] |= (a[j] & (ONE<<(RADIX - 46 - 2)))<<24;
      case 46:      a1[j2] |= (a[j] & (ONE<<(RADIX - 44 - 2)))<<23;
      case 44:      a1[j2] |= (a[j] & (ONE<<(RADIX - 42 - 2)))<<22;
      case 42:      a1[j2] |= (a[j] & (ONE<<(RADIX - 40 - 2)))<<21;
      case 40:      a1[j2] |= (a[j] & (ONE<<(RADIX - 38 - 2)))<<20;
      case 38:      a1[j2] |= (a[j] & (ONE<<(RADIX - 36 - 2)))<<19;
      case 36:      a1[j2] |= (a[j] & (ONE<<(RADIX - 34 - 2)))<<18;
      case 34:      a1[j2] |= (a[j] & (ONE<<(RADIX - 32 - 2)))<<17;
      case 32:      a1[j2] |= (a[j] & (ONE<<(RADIX - 30 - 2)))<<16;
      case 30:      a1[j2] |= (a[j] & (ONE<<(RADIX - 28 - 2)))<<15;
      case 28:      a1[j2] |= (a[j] & (ONE<<(RADIX - 26 - 2)))<<14;
      case 26:      a1[j2] |= (a[j] & (ONE<<(RADIX - 24 - 2)))<<13;
      case 24:      a1[j2] |= (a[j] & (ONE<<(RADIX - 22 - 2)))<<12;
      case 22:      a1[j2] |= (a[j] & (ONE<<(RADIX - 20 - 2)))<<11;
      case 20:      a1[j2] |= (a[j] & (ONE<<(RADIX - 18 - 2)))<<10;
      case 18:      a1[j2] |= (a[j] & (ONE<<(RADIX - 16 - 2)))<< 9;
      case 16:      a1[j2] |= (a[j] & (ONE<<(RADIX - 14 - 2)))<< 8;
      case 14:      a1[j2] |= (a[j] & (ONE<<(RADIX - 12 - 2)))<< 7;
      case 12:      a1[j2] |= (a[j] & (ONE<<(RADIX - 10 - 2)))<< 6;
      case 10:      a1[j2] |= (a[j] & (ONE<<(RADIX -  8 - 2)))<< 5;
      case  8:      a1[j2] |= (a[j] & (ONE<<(RADIX -  6 - 2)))<< 4;
      case  6:      a1[j2] |= (a[j] & (ONE<<(RADIX -  4 - 2)))<< 3;
      case  4:      a1[j2] |= (a[j] & (ONE<<(RADIX -  2 - 2)))<< 2;
      case  2:      a1[j2] |= (a[j] & (ONE<<(RADIX -  0 - 2)))<< 1;
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
    _mzed_slice2(C0, C1, C);

  mzd_t *A0 = mzd_init(A->nrows, A->ncols);
  mzd_t *A1 = mzd_init(A->nrows, A->ncols);
  _mzed_slice2(A0, A1, A);

  mzd_t *B0 = mzd_init(B->nrows, B->ncols);
  mzd_t *B1 = mzd_init(B->nrows, B->ncols);
  _mzed_slice2(B0, B1, B);

  /* compute */

  mzd_t *T0 = mzd_mul(NULL, A0, B0, 0);  /* A0B0 = A0*B0 */
  mzd_t *T1 = mzd_mul(NULL, A1, B1, 0);  /* A1B1 = A1*B1 */

  mzd_add(C1, C1, T0);
  mzd_add(C1, C1, T1); /*C1 += A0*B0 + A1*B1 */

  mzd_t *T2 = mzd_add(NULL, A0, A1); /*T2 = A0 + A1 */
  mzd_t *T3 = mzd_add(NULL, B0, B1); /*T3 = B0 + B1 */

  mzd_add(C0, C0, T1);
  mzd_addmul(C0, T2, T3, 0); /* C0 += T1 + T2*T3 */

  /* pack */
  if (C != NULL) {
    mzed_set_ui(C, 0);
  } else {
    C = mzed_init(A->finite_field, A->nrows, B->ncols);
  }
  _mzed_cling2(C, C0, C1);

  /* clean */

  mzd_free(T0);  mzd_free(T1);  mzd_free(T2);  mzd_free(T3);
  mzd_free(A0);  mzd_free(A1);  mzd_free(B0);  mzd_free(B1);
  mzd_free(C0);  mzd_free(C1);

  return C;
}


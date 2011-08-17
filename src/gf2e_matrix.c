/******************************************************************************
*
*            M4RIE: Linear Algebra over GF(2^e)
*
*    Copyright (C) 2010 Martin Albrecht <martinralbrecht@googlemail.com>
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

#include <stdlib.h>

#include "config.h"
#include "gf2e_matrix.h"
#include "strassen.h"
#include "bitslice.h"

mzed_t *mzed_init(gf2e* k, rci_t m, rci_t n) {
  mzed_t *A = (mzed_t *)m4ri_mm_malloc(sizeof(mzed_t));

  A->finite_field = k;
  A->w = gf2e_degree_to_w(A->finite_field);
  A->nrows = m;
  A->ncols = n;
  A->x = mzd_init(m, A->w*n);
  return A;
}

void mzed_free(mzed_t *A) {
  mzd_free(A->x);
  m4ri_mm_free(A);
}

void mzed_randomize(mzed_t *A) {
  int bitmask = (1<<A->finite_field->degree)-1;
  for(rci_t r=0; r<A->nrows; r++) {
    for(rci_t c=0; c<A->ncols; c++) {
      mzed_write_elem(A,r,c, random()&bitmask);
    }
  }
}

mzed_t *mzed_add(mzed_t *C, const mzed_t *A, const mzed_t *B) {
  if (A->nrows != B->nrows || A->ncols != B->ncols || A->finite_field != B->finite_field) {
    m4ri_die("mzed_add: rows, columns and fields must match.\n");
  }
  if (C == NULL) {
    C = mzed_init(A->finite_field, A->nrows, A->ncols);
  } else if (C != A) {
    if (C->finite_field != A->finite_field || C->nrows != A->nrows || C->ncols != A->ncols) {
      m4ri_die("mzed_add: rows and columns of returned matrix must match.\n");
    }
  }
  mzd_add(C->x, A->x, B->x);
  return C;
}

mzed_t *_mzed_add(mzed_t *C, const mzed_t *A, const mzed_t *B) {
  mzd_add(C->x, A->x, B->x);
  return C;
}


mzed_t *_mzed_mul_init(mzed_t *C, const mzed_t *A, const mzed_t *B, int clear) {
  if (A->ncols != B->nrows || A->finite_field != B->finite_field) {
    m4ri_die("mzed_mul: rows, columns and fields must match.\n");
  }
  if (C == NULL) {
    C = mzed_init(A->finite_field, A->nrows, B->ncols);
  } else {
    if (C->finite_field != A->finite_field || C->nrows != A->nrows || C->ncols != B->ncols) {
      m4ri_die("mzed_mul: rows and columns of returned matrix must match.\n");
    }
    if (clear)
      mzed_set_ui(C,0);
  }
  return C;
}


mzed_t *mzed_mul(mzed_t *C, const mzed_t *A, const mzed_t *B) {
  C = _mzed_mul_init(C,A,B, TRUE);
  _mzed_mul(C, A, B);
  return C;
}

mzed_t *mzed_addmul(mzed_t *C, const mzed_t *A, const mzed_t *B) {
  C = _mzed_mul_init(C,A,B, FALSE);
  _mzed_addmul(C, A, B);
  return C;
}

mzed_t *_mzed_mul(mzed_t *C, const mzed_t *A, const mzed_t *B) {
  if (A->finite_field->degree <= 3)
    return _mzed_mul_karatsuba(C, A, B);

  const rci_t cutoff = _mzed_strassen_cutoff(C, A, B);
  return _mzed_mul_strassen(C, A, B, cutoff);
}

mzed_t *_mzed_addmul(mzed_t *C, const mzed_t *A, const mzed_t *B) {
  if (A->finite_field->degree <= 3)
    return _mzed_mul_karatsuba(C, A, B);

  const rci_t cutoff = _mzed_strassen_cutoff(C, A, B);
  return _mzed_addmul_strassen(C, A, B, cutoff);
}

mzed_t *mzed_mul_naive(mzed_t *C, const mzed_t *A, const mzed_t *B) {
  C = _mzed_mul_init(C,A,B, TRUE);
  return _mzed_mul_naive(C, A, B);
}

mzed_t *mzed_addmul_naive(mzed_t *C, const mzed_t *A, const mzed_t *B) {
  C = _mzed_mul_init(C,A,B, FALSE);
  return _mzed_mul_naive(C, A, B);
}

mzed_t *_mzed_mul_naive(mzed_t *C, const mzed_t *A, const mzed_t *B) {
  gf2e* ff = C->finite_field;
  for (rci_t i=0; i<C->nrows; ++i) {
    for (rci_t j=0; j<C->ncols; ++j) {
      for (rci_t k=0; k<A->ncols; ++k) {
        mzed_add_elem(C, i, j, ff->mul[mzed_read_elem(A,i, k)][mzed_read_elem(B, k, j)]);
      }
    }
  }
  return C;
}

mzed_t *mzed_mul_scalar(mzed_t *C, const word a, const mzed_t *B) {
  /**
   * The algorithm proceeds as follows:
   */
  if(C == NULL)
    C = mzed_init(B->finite_field, B->nrows, B->ncols);

  const gf2e *ff = B->finite_field;
  const word *x = ff->mul[a];

  /**
   * 0) If a direct approach would need less lookups we use that.
   */

  if(ff->degree > 8 || B->nrows*B->ncols < 1<<17) {
    mzed_copy(C, B);
    for(rci_t i=0; i<B->nrows; i++)
      mzed_rescale_row(C, i, 0, x);
    return C;
  }

  word *mul = calloc(1<<16, sizeof(word));

  const word mask_16 = (1<<16)-1;
  const word mask_w = (1<<B->w)-1;

  /**
   * 1) We generate a lookup table of 16-bit wide entries
   *
   * @todo: this is a bit of overkill, we could do better
   */
  for(word i=0; i<1<<16; i++) {
    switch(C->w) {
    case 2:
      mul[i]  = x[(i&mask_w)] | x[((i>>2)&mask_w)]<<2 | x[((i>>4)&mask_w)]<<4 | x[((i>>6)&mask_w)]<<6;
      mul[i] |= x[((i>>8)&mask_w)]<<8 | x[((i>>10)&mask_w)]<<10 | x[((i>>12)&mask_w)]<<12 | x[((i>>14)&mask_w)]<<14;
      break;
    case 4:
      mul[i]  = x[(i&mask_w)] | x[((i>>4)&mask_w)]<<4 | x[((i>>8)&mask_w)]<<9 | x[((i>>12)&mask_w)]<<12;
      break;
    case 8:
      mul[i]  = x[(i&mask_w)] | x[((i>>8)&mask_w)]<<8;
      break;
    case 16:
      mul[i]  = x[(i&mask_w)];
      break;
    };
  }

  /**
   * 2) We use that lookup table to do 4 lookups per word
   */

  for(rci_t i=0; i<C->nrows; i++) {
    word *c_row = C->x->rows[i];
    const word *b_row = B->x->rows[i];
    for(wi_t j=0; j<C->x->width-1; j++) {
      const word tmp = b_row[j];
      const word a0 = tmp & mask_16;
      const word a1 = tmp>>16 & mask_16;
      const word a2 = tmp>>32 & mask_16;
      const word a3 = tmp>>48 & mask_16;
      c_row[j] = mul[a3]<<48 | mul[a2]<<32 | mul[a1]<<16 | mul[a0];      
    }
    /* deal with rest */
    const word tmp = b_row[B->x->width-1] & B->x->high_bitmask;
    const word a0 = tmp & mask_16;
    const word a1 = tmp>>16 & mask_16;
    const word a2 = tmp>>32 & mask_16;
    const word a3 = tmp>>48 & mask_16;
    c_row[C->x->width-1] &= ~B->x->high_bitmask;
    c_row[C->x->width-1] |= mul[a3]<<48 | mul[a2]<<32 | mul[a1]<<16 | mul[a0];
  }
  free(mul);
  return C;
}


mzed_t *mzed_copy(mzed_t *A, const mzed_t *B) {
  if (A == B)
    return A;
  if (A == NULL)
    A = mzed_init(B->finite_field, B->nrows, B->ncols);
  if (A->finite_field != B->finite_field || A->nrows != B->nrows || A->ncols != B->ncols) {
    m4ri_die("mzed_copy: target matrix has wrong dimensions or base field.");
  }
  mzd_copy(A->x, B->x);
  return A;
}

rci_t mzed_echelonize_naive(mzed_t *A, int full) {
  rci_t start_row,r,c,i,elim_start;
  word x = 0;

  rci_t nr = A->nrows;
  rci_t nc = A->ncols;

  gf2e *ff = A->finite_field;

  start_row = 0;

  for(c=0; c<nc; c++) {
    for(r=start_row; r<nr; r++) {
      x = mzed_read_elem(A, r, c);
      if (x) {
        mzed_rescale_row(A, r, c, ff->mul[ff->inv[x]]);
        mzd_row_swap(A->x, r, start_row);
        if (full)
          elim_start = 0;
        else
          elim_start = start_row + 1;
        for(i=elim_start; i<nr; i++) {
          if (i==start_row) 
            continue;
          x = mzed_read_elem(A,i,c);
          if(!x) continue;
          /* clear row */
          mzed_add_multiple_of_row(A, i, A, start_row, ff->mul[x], c);
        }
        start_row++;
        break;
      }
    }
  }
  return start_row;
}

void mzed_set_ui(mzed_t *A, word value) {
  mzd_set_ui(A->x, 0);
  if(!value)
    return;
  for(rci_t i=0; i< MIN(A->ncols,A->nrows); i++) {
    mzed_write_elem(A, i, i, value);
  }
}

void mzed_print(const mzed_t *A) {
  for (rci_t i=0; i < A->nrows; ++i) {
    printf("[");
    for (rci_t j=0; j < A->ncols; j++) {
      word tmp = mzed_read_elem(A,i,j);
      printf("[");
      switch(A->finite_field->degree) {
      case 16:  (tmp&(1ULL<<15)) ? printf("1") : printf("0");
      case 15:  (tmp&(1ULL<<14)) ? printf("1") : printf("0");
      case 14:  (tmp&(1ULL<<13)) ? printf("1") : printf("0");
      case 13:  (tmp&(1ULL<<12)) ? printf("1") : printf("0");
      case 12:  (tmp&(1ULL<<11)) ? printf("1") : printf("0");
      case 11:  (tmp&(1ULL<<10)) ? printf("1") : printf("0");
      case 10:  (tmp&(1ULL<< 9)) ? printf("1") : printf("0");
      case  9:  (tmp&(1ULL<< 8)) ? printf("1") : printf("0");
      case  8:  (tmp&(1ULL<< 7)) ? printf("1") : printf("0");
      case  7:  (tmp&(1ULL<< 6)) ? printf("1") : printf("0");
      case  6:  (tmp&(1ULL<< 5)) ? printf("1") : printf("0");
      case  5:  (tmp&(1ULL<< 4)) ? printf("1") : printf("0");
      case  4:  (tmp&(1ULL<< 3)) ? printf("1") : printf("0");
      case  3:  (tmp&(1ULL<< 2)) ? printf("1") : printf("0");
      case  2:  (tmp&(1ULL<< 1)) ? printf("1") : printf("0");
      case  1:  (tmp&(1ULL<< 0)) ? printf("1") : printf("0");
        break;
      default: m4ri_die("degree %lz too big\n", A->finite_field->degree);
      }
      printf("]");
      if(j<A->ncols-1)
        printf(" ");
    }
    printf("]\n");
  }
}


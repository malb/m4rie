#ifndef GF2E_MATRIX_H
#define GF2E_MATRIX_H

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

#include <m4ri/m4ri.h>
#include "finite_field.h"
#include "m4ri_functions.h"
 
typedef struct {

  /**
   * The internal representation of our matrices
   */
  mzd_t *x;

  /**
   * The finite field GF(2^n)
   */

  gf2e *finite_field;

  /**
   * Number of rows.
   */

  size_t nrows;

  /**
   * Number of columns.
   */

  size_t ncols;

  /**
   * The internal width of elements (must divide 64)
   */

  size_t width;

} mzed_t;


mzed_t *mzed_init(gf2e *ff, const size_t m, const size_t n);
void mzed_free(mzed_t *a);

void mzed_randomize(mzed_t *a);

mzed_t *mzed_add(mzed_t *c, const mzed_t *a, const mzed_t *b);
mzed_t *mzed_copy(mzed_t *o, const mzed_t *i);

static inline int mzed_cmp(mzed_t *l, mzed_t *r) {
  return mzd_cmp(l->x,r->x);
}

static inline int mzed_read_elem(const mzed_t *a, const size_t row, const size_t col) {
  return (int)__mzd_read_bits(a->x, row, a->width*col, a->width);
}

static inline void mzed_add_elem(mzed_t *a, const size_t row, const size_t col, const int elem) {
  __mzd_xor_bits(a->x, row, a->width*col, a->width, elem);
}

static inline void mzed_write_elem(mzed_t *a, const size_t row, const size_t col, const int elem) {
  __mzd_clear_bits(a->x, row, a->width*col, a->width);
  __mzd_xor_bits(a->x, row, a->width*col, a->width, elem);
}

static inline void mzed_add_multiple_of_row(mzed_t *to, size_t to_row, mzed_t *from, size_t from_row, word *X, size_t start_col) {
  assert(to->ncols == from->ncols && to->finite_field == from->finite_field);
  assert(to->x->offset == 0 && from->x->offset == 0);
  if(to->width == 4) {
    size_t startblock = start_col/RADIX;
    mzd_t *from_x = from->x;
    mzd_t *to_x = to->x;
    word *_f = from_x->rows[from_row];
    word *_t = to_x->rows[to_row];
    size_t j;
    for(j=startblock; j<to_x->width -1; j++) {
      _t[j] ^= ((word)X[(int)((_f[j] & 0xF000000000000000ULL)>>60)])<<60;
      _t[j] ^= ((word)X[(int)((_f[j] & 0x0F00000000000000ULL)>>56)])<<56;
      _t[j] ^= ((word)X[(int)((_f[j] & 0x00F0000000000000ULL)>>52)])<<52;
      _t[j] ^= ((word)X[(int)((_f[j] & 0x000F000000000000ULL)>>48)])<<48;
      _t[j] ^= ((word)X[(int)((_f[j] & 0x0000F00000000000ULL)>>44)])<<44;
      _t[j] ^= ((word)X[(int)((_f[j] & 0x00000F0000000000ULL)>>40)])<<40;
      _t[j] ^= ((word)X[(int)((_f[j] & 0x000000F000000000ULL)>>36)])<<36;
      _t[j] ^= ((word)X[(int)((_f[j] & 0x0000000F00000000ULL)>>32)])<<32;
      _t[j] ^= ((word)X[(int)((_f[j] & 0x00000000F0000000ULL)>>28)])<<28;
      _t[j] ^= ((word)X[(int)((_f[j] & 0x000000000F000000ULL)>>24)])<<24;
      _t[j] ^= ((word)X[(int)((_f[j] & 0x0000000000F00000ULL)>>20)])<<20;
      _t[j] ^= ((word)X[(int)((_f[j] & 0x00000000000F0000ULL)>>16)])<<16;
      _t[j] ^= ((word)X[(int)((_f[j] & 0x000000000000F000ULL)>>12)])<<12;
      _t[j] ^= ((word)X[(int)((_f[j] & 0x0000000000000F00ULL)>> 8)])<< 8;
      _t[j] ^= ((word)X[(int)((_f[j] & 0x00000000000000F0ULL)>> 4)])<< 4;
      _t[j] ^= ((word)X[(int)((_f[j] & 0x000000000000000FULL)>> 0)])<< 0;
    }
    switch(to_x->ncols % RADIX) {
    case  0: _t[j] ^= ((word)X[(int)((_f[j] & 0x000000000000000FULL)>> 0)])<< 0;
    case 60: _t[j] ^= ((word)X[(int)((_f[j] & 0x00000000000000F0ULL)>> 4)])<< 4;
    case 56: _t[j] ^= ((word)X[(int)((_f[j] & 0x0000000000000F00ULL)>> 8)])<< 8;
    case 52: _t[j] ^= ((word)X[(int)((_f[j] & 0x000000000000F000ULL)>>12)])<<12;
    case 48: _t[j] ^= ((word)X[(int)((_f[j] & 0x00000000000F0000ULL)>>16)])<<16;
    case 44: _t[j] ^= ((word)X[(int)((_f[j] & 0x0000000000F00000ULL)>>20)])<<20;
    case 40: _t[j] ^= ((word)X[(int)((_f[j] & 0x000000000F000000ULL)>>24)])<<24;
    case 36: _t[j] ^= ((word)X[(int)((_f[j] & 0x00000000F0000000ULL)>>28)])<<28;
    case 32: _t[j] ^= ((word)X[(int)((_f[j] & 0x0000000F00000000ULL)>>32)])<<32;
    case 28: _t[j] ^= ((word)X[(int)((_f[j] & 0x000000F000000000ULL)>>36)])<<36;
    case 24: _t[j] ^= ((word)X[(int)((_f[j] & 0x00000F0000000000ULL)>>40)])<<40;
    case 20: _t[j] ^= ((word)X[(int)((_f[j] & 0x0000F00000000000ULL)>>44)])<<44;
    case 16: _t[j] ^= ((word)X[(int)((_f[j] & 0x000F000000000000ULL)>>48)])<<48;
    case 12: _t[j] ^= ((word)X[(int)((_f[j] & 0x00F0000000000000ULL)>>52)])<<52;
    case  8: _t[j] ^= ((word)X[(int)((_f[j] & 0x0F00000000000000ULL)>>56)])<<56;
    case  4: _t[j] ^= ((word)X[(int)((_f[j] & 0xF000000000000000ULL)>>60)])<<60;
    };

  } else if (to->width == 8) {
    size_t startblock = start_col/RADIX;
    mzd_t *from_x = from->x;
    mzd_t *to_x = to->x;
    word *_f = from_x->rows[from_row];
    word *_t = to_x->rows[to_row];
    size_t j;
    register word __t, __f;
    /* for(j=startblock; j<to_x->width -1; j++) { */
    /*   tmp = _f[j]; */
    /*   _t[j] ^= ((word)X[(int)((tmp & 0xFF00000000000000ULL)>>56)])<<56; */
    /*   _t[j] ^= ((word)X[(int)((tmp & 0x00FF000000000000ULL)>>48)])<<48; */
    /*   _t[j] ^= ((word)X[(int)((tmp & 0x0000FF0000000000ULL)>>40)])<<40; */
    /*   _t[j] ^= ((word)X[(int)((tmp & 0x000000FF00000000ULL)>>32)])<<32; */
    /*   _t[j] ^= ((word)X[(int)((tmp & 0x00000000FF000000ULL)>>24)])<<24; */
    /*   _t[j] ^= ((word)X[(int)((tmp & 0x0000000000FF0000ULL)>>16)])<<16; */
    /*   _t[j] ^= ((word)X[(int)((tmp & 0x000000000000FF00ULL)>> 8)])<< 8; */
    /*   _t[j] ^= ((word)X[(int)((tmp & 0x00000000000000FFULL)>> 0)])<< 0; */
    /* } */

    for(j=startblock; j<to_x->width -1; j++) {
      __f = _f[j], __t = _t[j];
      __t ^= (X[((__f)& 0x00000000000000FFULL)])<<0;
      __f >>= 8;
      __t ^= (X[((__f)& 0x00000000000000FFULL)])<<8;
      __f >>= 8;
      __t ^= (X[((__f)& 0x00000000000000FFULL)])<<16;
      __f >>= 8;
      __t ^= (X[((__f)& 0x00000000000000FFULL)])<<24;
      __f >>= 8;
      __t ^= (X[((__f)& 0x00000000000000FFULL)])<<32;
      __f >>= 8;
      __t ^= (X[((__f)& 0x00000000000000FFULL)])<<40;
      __f >>= 8;
      __t ^= (X[((__f)& 0x00000000000000FFULL)])<<48;
      __f >>= 8;
      __t ^= (X[((__f)& 0x00000000000000FFULL)])<<56;
      __f >>= 8;
      _t[j] = __t;
    }
    
    switch(to_x->ncols % RADIX) {
    case  0: _t[j] ^= ((word)X[(int)((_f[j] & 0x00000000000000FFULL)>> 0)])<< 0;
    case 56: _t[j] ^= ((word)X[(int)((_f[j] & 0x000000000000FF00ULL)>> 8)])<< 8;
    case 48: _t[j] ^= ((word)X[(int)((_f[j] & 0x0000000000FF0000ULL)>>16)])<<16;
    case 40: _t[j] ^= ((word)X[(int)((_f[j] & 0x00000000FF000000ULL)>>24)])<<24;
    case 32: _t[j] ^= ((word)X[(int)((_f[j] & 0x000000FF00000000ULL)>>32)])<<32;
    case 24: _t[j] ^= ((word)X[(int)((_f[j] & 0x0000FF0000000000ULL)>>40)])<<40;
    case 16: _t[j] ^= ((word)X[(int)((_f[j] & 0x00FF000000000000ULL)>>48)])<<48;
    case  8: _t[j] ^= ((word)X[(int)((_f[j] & 0xFF00000000000000ULL)>>56)])<<56;
    };

  }  else {
    for(size_t j=start_col; j<from->ncols; j++) {
      mzed_add_elem(to, to_row, j, X[mzed_read_elem(from, from_row, j)]);
    }
  }
}

size_t mzed_echelonize_naive(mzed_t *a, int full);
size_t mzed_echelonize_table(mzed_t *a, int full);

mzed_t *mzed_mul_table(mzed_t *C, mzed_t *A, mzed_t *B);

#endif //MATRIX_H

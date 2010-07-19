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

#include <givaro/givgfq.h>
#include <m4rie/m4rie.h>

namespace M4RIE {

  class FiniteField: public GFqDom<int> {
  public: 
    unsigned int  log2pol(int x) { return _log2pol[x]; };
    unsigned int  pol2log(int x) { return _pol2log[x]; };
  }; 
};

static inline gf2e *gf2e_init_givgfq(M4RIE::FiniteField *givgfq) {
  gf2e *ff = (gf2e*)m4ri_mm_malloc(sizeof(gf2e));
  ff->degree = givgfq->exponent();

  ff->mul = (word **)m4ri_mm_calloc(TWOPOW(givgfq->exponent()), sizeof(word *));
  for(size_t i = 0; i<TWOPOW(givgfq->exponent()); i++) {
    ff->mul[i] = (word *)m4ri_mm_calloc(TWOPOW(givgfq->exponent()),sizeof(word));
    for(size_t j=0; j<TWOPOW(givgfq->exponent()); j++) {
      int prod = givgfq->mul(prod, givgfq->pol2log(i) , givgfq->pol2log(j));
      ff->mul[i][j] = givgfq->log2pol(prod);
    }
  }
  ff->inv = (word*)m4ri_mm_calloc(TWOPOW(givgfq->exponent()), sizeof(word));
  for(size_t i = 0; i<TWOPOW(givgfq->exponent()); i++) {
    int tmp = givgfq->inv(tmp, givgfq->pol2log(i));
    ff->inv[i] = givgfq->log2pol(tmp);
  }
  return ff;
}
 
static inline int mzed_read_elem_log(const mzed_t *a, const size_t row, const size_t col, M4RIE::FiniteField *ff) {
  return ff->pol2log((int)__mzd_read_bits(a->x, row, a->w*col, a->w));
};

static inline void mzed_write_elem_log(mzed_t *a, const size_t row, const size_t col, const int elem, M4RIE::FiniteField *ff) {
  __mzd_clear_bits(a->x, row, a->w*col, a->w);
  __mzd_xor_bits(a->x, row, a->w*col, a->w, ff->log2pol(elem));
};

static inline void mzed_add_elem_log(mzed_t *a, const size_t row, const size_t col, const int elem, M4RIE::FiniteField *ff) {
  __mzd_xor_bits(a->x, row, a->w*col, a->w, ff->log2pol(elem));
};

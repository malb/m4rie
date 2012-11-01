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

#include <givaro/givconfig.h>
#include <givaro/givgfq.h>
#include <m4rie.h>

namespace M4RIE {

#if GIVARO_VERSION  <  30400 || GIVARO_VERSION >= 196608 // old Givaro versions used 0x03xxyy
  class FiniteField: public GFqDom<int> {
  public:
  FiniteField(const unsigned int e) : GFqDom<int>(2, e){};
#else
  class FiniteField: public Givaro::GFqDom<int> {
  public:
  FiniteField(const unsigned int e) : Givaro::GFqDom<int>(2, e){};
#endif
      unsigned int  log2pol(int x) { return _log2pol[x]; };
      unsigned int  pol2log(int x) { return _pol2log[x]; };
  };
};

static inline gf2e *gf2e_init_givgfq(M4RIE::FiniteField *givgfq) {
  word minpoly = givgfq->pol2log(1);
  unsigned int degree = givgfq->exponent();
  for(unsigned int i = 0; i<degree; i++) {
    minpoly = givgfq->mul((int&)minpoly, (int)givgfq->pol2log(2) , (int)minpoly);
  }
  minpoly = givgfq->log2pol(minpoly);
  minpoly = minpoly ^ (1<<degree);
  return gf2e_init(minpoly);
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

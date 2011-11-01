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
#include <stdarg.h>

static inline mzd_t *mzd_add_accumulate(mzd_t *A, const int n, ...) {
  /**
   * adds n elements to A
   */
  assert(n>1);
  va_list b_list;
  va_start( b_list, n );

  mzd_add(A, va_arg(b_list, mzd_t *), va_arg(b_list, mzd_t *));

  for( int i = 0 ; i < n-2; i++ ) {
    mzd_t *B = va_arg(b_list, mzd_t *);
    mzd_add(A, A, B);
  }

  va_end( b_list );
  return A;
}

static inline mzd_t *mzd_add_to_all(mzd_t *A, const int n, ...) {
  /**
   * adds A to n elements
   */
  assert(n>1);
  va_list b_list;
  va_start( b_list, n );

  for( int i = 0 ; i < n; i++ ) {
    mzd_t *B = va_arg(b_list, mzd_t *);
    mzd_add(B, B, A);
  }

  va_end( b_list );
  return A;
}

void mzd_slice_set_ui(mzd_slice_t *A, word value) {
  for(int i=0; i<A->depth; i++) {
    mzd_set_ui(A->x[i], (value>>i)&1);
  }
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

mzd_slice_t *_mzd_slice_mul_naive(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B) {
  if (C == NULL)
    C = mzd_slice_init(A->finite_field, A->nrows, B->ncols);

  const unsigned int e = A->finite_field->degree;
  const word minpoly = A->finite_field->minpoly;

  C = _mzd_slice_adapt_depth(C,2*e-1);

  for(unsigned int i=0; i<e; i++) {
    for(unsigned int j=0; j<e; j++) {
      mzd_addmul(C->x[i+j], A->x[i], B->x[j], 0);
    }
  }

  for(unsigned int i=2*e-2; i>= e; i--)
    for(unsigned int j=0; j<e; j++)
      if (minpoly & 1<<j)
        mzd_add(C->x[i-e+j], C->x[i-e+j], C->x[i]);
  _mzd_slice_adapt_depth(C,e);
  return C;
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

static void _poly_addmul2(mzd_t **X, const mzd_t **a, const mzd_t **b) {
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

static void _poly_addmul4(mzd_t **c, const mzd_t **a, const mzd_t **b) {
  const mzd_t *a0[2] = {a[0],a[1]};
  const mzd_t *a1[2] = {a[2],a[3]};
  const mzd_t *b0[2] = {b[0],b[1]};
  const mzd_t *b1[2] = {b[2],b[3]};

  mzd_t *X[3][3] = { {c[0],c[1],c[2]}, 
                     {c[2],c[3],c[4]},
                     {c[4],c[5],c[6]} };

  mzd_t *t0[3];
  mzd_t *t1[2];

  t0[0] = mzd_init(a[0]->nrows, a[0]->ncols);
  t0[1] = mzd_init(a[0]->nrows, a[0]->ncols);
  t1[0] = mzd_init(b[0]->nrows, b[0]->ncols);
  t1[1] = mzd_init(b[0]->nrows, b[0]->ncols);

  _poly_add(t0, a0, a1, 2);
  _poly_add(t1, b0, b1, 2);

  _poly_addmul2(X[1], (const mzd_t**)t0, (const mzd_t**)t1);

  mzd_free(t0[0]);  mzd_free(t0[1]);
  mzd_free(t1[0]);  mzd_free(t1[1]);

  t0[0] = mzd_init(a[0]->nrows, b[0]->ncols);
  t0[1] = mzd_init(a[0]->nrows, b[0]->ncols);
  t0[2] = mzd_init(a[0]->nrows, b[0]->ncols);

  _poly_addmul2(t0, a0, b0);
  _poly_add(X[0], (const mzd_t**)X[0], (const mzd_t**)t0, 3);
  _poly_add(X[1], (const mzd_t**)X[1], (const mzd_t**)t0, 3);

  mzd_set_ui(t0[0], 0);
  mzd_set_ui(t0[1], 0);
  mzd_set_ui(t0[2], 0);
  
  _poly_addmul2(t0, a1, b1);
  _poly_add(X[1], (const mzd_t**)X[1], (const mzd_t**)t0, 3);
  _poly_add(X[2], (const mzd_t**)X[2], (const mzd_t**)t0, 3);

  mzd_free(t0[0]); mzd_free(t0[1]); mzd_free(t0[2]); 
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

  _poly_addmul2(X[1], (const mzd_t**)t0, (const mzd_t**)t1);

  mzd_free(t0[0]);  mzd_free(t0[1]);
  mzd_free(t1[0]);  mzd_free(t1[1]);

  t0[0] = mzd_init(A->nrows, B->ncols);
  t0[1] = mzd_init(A->nrows, B->ncols);
  t0[2] = mzd_init(A->nrows, B->ncols);

  _poly_addmul2(t0, a0, b0);
  _poly_add(X[0], (const mzd_t**)X[0], (const mzd_t**)t0, 3);
  _poly_add(X[1], (const mzd_t**)X[1], (const mzd_t**)t0, 3);

  mzd_set_ui(t0[0], 0);
  mzd_set_ui(t0[1], 0);
  mzd_set_ui(t0[2], 0);
  
  _poly_addmul2(t0, a1, b1);
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

static inline void mzd_add_modred5(mzd_t **X, const int t, mzd_t *A, const word minpoly) {
  /**
   * @warning the reduction is not complete, i.e. not applied
   * recursively
   */
  if (mzd_is_zero(A))
    return;

  mzd_add(X[t-5+0], X[t-5+0], A);

  switch(minpoly) {
  case 37: //x^5 + x^2 + 1
    mzd_add(X[t-5+2], X[t-5+2], A);
    break;
  case 41: //x^5 + x^3 + 1
    mzd_add(X[t-5+3], X[t-5+3], A);
    break;
  case 47: //x^5 + x^3 + x^2 + x + 1
    mzd_add(X[t-5+1], X[t-5+1], A);
    mzd_add(X[t-5+2], X[t-5+2], A);
    mzd_add(X[t-5+3], X[t-5+3], A);
    break;
  case 55: //x^5 + x^4 + x^2 + x + 1
    mzd_add(X[t-5+1], X[t-5+1], A);
    mzd_add(X[t-5+2], X[t-5+2], A);
    mzd_add(X[t-5+4], X[t-5+4], A);
    break;
  case 59: //x^5 + x^4 + x^3 + x + 1
    mzd_add(X[t-5+1], X[t-5+1], A);
    mzd_add(X[t-5+3], X[t-5+3], A);
    mzd_add(X[t-5+4], X[t-5+4], A);
    break;
  case 61: //x^5 + x^4 + x^3 + x^2 + 1
    mzd_add(X[t-5+2], X[t-5+2], A);
    mzd_add(X[t-5+3], X[t-5+3], A);
    mzd_add(X[t-5+4], X[t-5+4], A);
    break;
  default:
    m4ri_die("minpoly 0x%03x is not an irreducible polynomial of degree 5 over GF(2)",(unsigned int)minpoly);
  }
}

mzd_slice_t *_mzd_slice_mul_karatsuba5(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B) {
  /* using three + three temporary matrices */
  if (C == NULL)
    C = mzd_slice_init(A->finite_field, A->nrows, B->ncols);

  const word minpoly = A->finite_field->minpoly;
  C = _mzd_slice_adapt_depth(C, 8);

  const mzd_t *a0 = A->x[0];
  const mzd_t *a1 = A->x[1];
  const mzd_t *a2 = A->x[2];
  const mzd_t *a3 = A->x[3];
  const mzd_t *a4 = A->x[4];

  const mzd_t *b0 = B->x[0];
  const mzd_t *b1 = B->x[1];
  const mzd_t *b2 = B->x[2];
  const mzd_t *b3 = B->x[3];
  const mzd_t *b4 = B->x[4];

  mzd_t **X = C->x;

  mzd_t *t0 = mzd_init(a0->nrows, b0->ncols);

  /* a0b0(X^6 - 3X^5 + 3X^4 - 2X^3 + X^2 - X + 1) => a0b0(X^6 + X^5 + X^4 + X^2 + X + 1)*/
  mzd_mul(t0, a0, b0, 0);
  mzd_add_modred5(X, 6, t0, minpoly);
  mzd_add_modred5(X, 5, t0, minpoly);
  mzd_add(X[4], X[4], t0);
  mzd_add(X[2], X[2], t0);
  mzd_add(X[1], X[1], t0);
  mzd_add(X[0], X[0], t0);

  /* a1b1(X^4 - 2X^3 + 2X^2 - X) => a1b1(X^4 + X)*/
  mzd_mul(t0, a1, b1, 0);
  mzd_add(X[4], X[4], t0);
  mzd_add(X[1], X[1], t0);

  /* a3b3(-X^7 +2X^6 - 2X^5 + X^4) => a3b3(X^7 + X^4) */
  mzd_mul(t0, a3, b3, 0);
  mzd_add_modred5(X, 7, t0, minpoly);
  mzd_add(X[4], X[4], t0);

  /* (a4b4)(X^8 - X^7 + X^6 - 2X^5 + 3X^4 - 3X^3 + X^2) = > (a4b4)(X^8 + X^7 + X^6 + X^4 + X^3 + X^2) */
  mzd_mul(t0, a4, b4, 0);
  mzd_add_modred5(X, 8, t0, minpoly);
  mzd_add_modred5(X, 7, t0, minpoly);
  mzd_add_modred5(X, 6, t0, minpoly);
  mzd_add(X[4], X[4], t0);
  mzd_add(X[3], X[3], t0);
  mzd_add(X[2], X[2], t0);

  mzd_t *t1 = mzd_init(a0->nrows, a1->ncols);
  mzd_t *t2 = mzd_init(b0->nrows, b1->ncols);

  /* (a0-a4)(b0-b4)(-X^6 + 3X^5 - 3X^3 + X^2) => (a0+a4)(b0+b4)(X^6 + X^5 + X^3 + X^2) */
  mzd_add(t1, a0, a4);
  mzd_add(t2, b0, b4);
  mzd_mul(t0, t1, t2, 0);
  mzd_add_modred5(X, 6, t0, minpoly);
  mzd_add_modred5(X, 5, t0, minpoly);
  mzd_add(X[3], X[3], t0);
  mzd_add(X[2], X[2], t0);

  /* (a0+a1)(b0+b1)(-X^5 + X^4 - X^2 + X) => (a0+a1)(b0+b1)(X^5 + X^4 + X^2 + X) */
  mzd_add(t1, a0, a1);
  mzd_add(t2, b0, b1);
  mzd_mul(t0, t1, t2, 0);
  mzd_add_modred5(X, 5, t0, minpoly);
  mzd_add(X[4], X[4], t0);
  mzd_add(X[2], X[2], t0);
  mzd_add(X[1], X[1], t0);

  /* (a3+a4)(b3+b4)(X^7 -X^6 + X^4 - X^3) => (a3+a4)(b3+b4)(X^7 + X^6 + X^4 + X^3) */
  mzd_add(t1, a3, a4);
  mzd_add(t2, b3, b4);
  mzd_mul(t0, t1, t2, 0);
  mzd_add_modred5(X, 7, t0, minpoly);
  mzd_add_modred5(X, 6, t0, minpoly);
  mzd_add(X[4], X[4], t0);
  mzd_add(X[3], X[3], t0);

  /* (a1+a2-a4)(a1+a2-a4)(-X^4 + 2X^3 - X^2) => (a1+a2+a4)(a1+a2+a4)(X^4 + X^2) */
  mzd_add(t1, a1, a2);
  mzd_add(t1, t1, a4);
  mzd_add(t2, b1, b2);
  mzd_add(t2, t2, b4);
  mzd_mul(t0, t1, t2, 0);
  mzd_add(X[4], X[4], t0);
  mzd_add(X[2], X[2], t0);

  /* (a0-a2-a3)(b0-b2-b3)(-X^6 + 2X^5 - X^4) => (a0+a2+a3)(b0+b2+b3)(X^6 + X^4) */
  mzd_add(t1, a0, a2);
  mzd_add(t1, t1, a3);
  mzd_add(t2, b0, b2);
  mzd_add(t2, t2, b3);
  mzd_mul(t0, t1, t2, 0);
  mzd_add_modred5(X, 6, t0, minpoly);
  mzd_add(X[4], X[4], t0);

  /* (a0+a1-a3-a4)(b0+b1-b3-b4)(X^5 - 2X^4 + X^3) => (a0+a1+a3+a4)(b0+b1+b3+b4)(X^5 + X^3) */
  mzd_add(t1, a0, a1);
  mzd_add(t1, t1, a3);
  mzd_add(t1, t1, a4);
  mzd_add(t2, b0, b1);
  mzd_add(t2, t2, b3);
  mzd_add(t2, t2, b4);
  mzd_mul(t0, t1, t2, 0);
  mzd_add_modred5(X, 5, t0, minpoly);
  mzd_add(X[3], X[3], t0);

  /* (a0+a1+a2-a4)(b0+b1+b2-b4)(-X^5 + 2X^4 - 2X^3 + X^2) => (a0+a1+a2+a4)(b0+b1+b2+b4)(X^5 + X^2) */
  mzd_add(t1, a0, a1);
  mzd_add(t1, t1, a2);
  mzd_add(t1, t1, a4);
  mzd_add(t2, b0, b1);
  mzd_add(t2, t2, b2);
  mzd_add(t2, t2, b4);
  mzd_mul(t0, t1, t2, 0);
  mzd_add_modred5(X, 5, t0, minpoly);
  mzd_add(X[2], X[2], t0);

  /* (a0-a2-a3-a4)(b0-b2-b3-b4)(X^6 - 2X^5 + 2X^4 - X^3) => (a0+a2+a3+a4)(b0+b2+b3+b4)(X^6 + X^3)*/
  mzd_add(t1, a0, a2);
  mzd_add(t1, t1, a3);
  mzd_add(t1, t1, a4);
  mzd_add(t2, b0, b2);
  mzd_add(t2, t2, b3);
  mzd_add(t2, t2, b4);
  mzd_mul(t0, t1, t2, 0);
  mzd_add_modred5(X, 6, t0, minpoly);
  mzd_add(X[3], X[3], t0);

  /* (a0+a1+a2+a3+a4)(b0+b1+b2+b3+b4)(X^5 - X^4 + X^3) = > (a0+a1+a2+a3+a4)(b0+b1+b2+b3+b4)(X^5 + X^4 + X^3)*/
  mzd_add(t1, a0, a1);
  mzd_add(t1, t1, a2);
  mzd_add(t1, t1, a3);
  mzd_add(t1, t1, a4);
  mzd_add(t2, b0, b1);
  mzd_add(t2, t2, b2);
  mzd_add(t2, t2, b3);
  mzd_add(t2, t2, b4);
  mzd_mul(t0, t1, t2, 0);
  mzd_add_modred5(X, 5, t0, minpoly);
  mzd_add(X[4], X[4], t0);
  mzd_add(X[3], X[3], t0);

  mzd_free(t0);
  mzd_free(t1);
  mzd_free(t2);

  mzd_add_modred5(X, 7, X[7], minpoly);
  mzd_add_modred5(X, 6, X[6], minpoly);
  mzd_add_modred5(X, 5, X[5], minpoly);

  _mzd_slice_adapt_depth(C,5);
  return C;
}

static inline void mzd_add_modred6(mzd_t **X, const int t, mzd_t *A, const word minpoly) {
  /**
   * @warning the reduction is not complete, i.e. not applied
   * recursively
   */
  if (mzd_is_zero(A))
    return;

  mzd_add(X[t-6+0], X[t-6+0], A);

  switch(minpoly) {
  case  67: // x^6 + x + 1
    mzd_add(X[t-6+1], X[t-6+1], A);
    break;
  case  73: // x^6 + x^3 + 1
    mzd_add(X[t-6+3], X[t-6+3], A);
    break;
  case  87: // x^6 + x^4 + x^2 + x + 1
    mzd_add_to_all(A, 3, X[t-6+4], X[t-6+2], X[t-6+1]);
    break;
  case  91: // x^6 + x^4 + x^3 + x + 1
    mzd_add_to_all(A, 3, X[t-6+4], X[t-6+3], X[t-6+1]);
    break;
  case  97: // x^6 + x^5 + 1
    mzd_add(X[t-6+5], X[t-6+5], A);
    break;
  case 103: // x^6 + x^5 + x^2 + x + 1
    mzd_add_to_all(A, 3, X[t-6+5], X[t-6+2], X[t-6+1]);
    break;
  case 109: // x^6 + x^5 + x^3 + x^2 + 1
    mzd_add_to_all(A, 3, X[t-6+5], X[t-6+3], X[t-6+2]);
    break;
  case 115: // x^6 + x^5 + x^4 + x + 1
    mzd_add_to_all(A, 3, X[t-6+5], X[t-6+4], X[t-6+1]);
    break;
  case 117: // x^6 + x^5 + x^4 + x^2 + 1
    mzd_add_to_all(A, 3, X[t-6+5], X[t-6+4], X[t-6+2]);
    break;
  default:
    m4ri_die("minpoly 0x%03x is not an irreducible polynomial of degree 5 over GF(2)",(unsigned int)minpoly);
  }
}


mzd_slice_t *_mzd_slice_mul_karatsuba6(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B) {
  /* using three + four temporary matrices */
  if (C == NULL)
    C = mzd_slice_init(A->finite_field, A->nrows, B->ncols);

  const word minpoly = A->finite_field->minpoly;
  C = _mzd_slice_adapt_depth(C, 10);

  const mzd_t *a0 = A->x[0];
  const mzd_t *a1 = A->x[1];
  const mzd_t *a2 = A->x[2];
  const mzd_t *a3 = A->x[3];
  const mzd_t *a4 = A->x[4];
  const mzd_t *a5 = A->x[5];

  const mzd_t *b0 = B->x[0];
  const mzd_t *b1 = B->x[1];
  const mzd_t *b2 = B->x[2];
  const mzd_t *b3 = B->x[3];
  const mzd_t *b4 = B->x[4];
  const mzd_t *b5 = B->x[5];

  mzd_t **X = C->x;

  mzd_t *t0 = mzd_init(a0->nrows, b0->ncols);
  mzd_t *t1 = mzd_init(a0->nrows, a1->ncols);
  mzd_t *t2 = mzd_init(b0->nrows, b1->ncols);

  /* a5b5 (X^10 + X^9+ X^6 + X^5) */
  mzd_mul(t0, a5, b5, 0);
  mzd_add_modred6(X, 10, t0, minpoly);
  mzd_add_to_all(t0, 3,   X[ 9],X[ 6],X[ 5]);

  /* a4b4 (X^9 + X^7 + X^5 + X^3) */
  mzd_mul(t0, a4, b4, 0);
  mzd_add_to_all(t0, 4,   X[ 9],X[ 7],X[ 5],X[ 3]);

  /* a1b1 (X^7 + X^6 + X^5 + X^4 + X^3 + X) */
  mzd_mul(t0, a1, b1, 0);
  mzd_add_to_all(t0, 6,   X[ 7],X[ 6],X[ 5],X[ 4],X[ 3],X[ 1]);

  /* a0b0 (X^6 + X^5 + X + 1) */
  mzd_mul(t0, a0, b0, 0);
  mzd_add_to_all(t0, 4,   X[ 6],X[ 5],X[ 1],X[ 0]);

  /* (a4 + a5)(b4 + b5) (X^9 + X^8 + X^4+ X^3) */
  mzd_add_accumulate(t1, 2, a4, a5);
  mzd_add_accumulate(t2, 2, b4, b5);
  mzd_mul(t0, t1, t2, 0);
  mzd_add_to_all(t0, 4,    X[ 9], X[ 8], X[ 4], X[ 3]);

  /* (a0 + a1)(b0 + b1) (X^7 + X^4 + X^2 + X) */
  mzd_add_accumulate(t1, 2, a0, a1);
  mzd_add_accumulate(t2, 2, b0, b1);
  mzd_mul(t0, t1, t2, 0);
  mzd_add_to_all(t0, 4,   X[ 7], X[ 4], X[ 2], X[ 1]);

  /* (a3 + a4)(b3 + b4)(X^8 + X^7 + X^6 + X^3) */
  mzd_add_accumulate(t1, 2, a3, a4);
  mzd_add_accumulate(t2, 2, b3, b4);
  mzd_mul(t0, t1, t2, 0);
  mzd_add_to_all(t0, 4,   X[ 8], X[ 7], X[ 6], X[ 3]);

  /* (a1 + a2)(b1 + b2) (X^7 + X^6 + X^3 + X^2) */
  mzd_add_accumulate(t1, 2, a1, a2);
  mzd_add_accumulate(t2, 2, b1, b2);
  mzd_mul(t0, t1, t2, 0);
  mzd_add_to_all(t0, 4,   X[ 7], X[ 6], X[ 3], X[ 2]);

  /* (a1 + a4)(b1 + b4) (X^4 + X^5) */
  mzd_add_accumulate(t1, 2, a1, a4);
  mzd_add_accumulate(t2, 2, b1, b4);
  mzd_mul(t0, t1, t2, 0);
  mzd_add_to_all(t0, 2,   X[ 4], X[ 5]);

  /* (a2 + a3)(b2 + b3) (X^7 + X^6 + X^4 + X^3) */
  mzd_add_accumulate(t1, 2, a2, a3);
  mzd_add_accumulate(t2, 2, b2, b3);
  mzd_mul(t0, t1, t2, 0);
  mzd_add_to_all(t0, 4,   X[ 7], X[ 6], X[ 4], X[ 3]);

  /* (a3 + a4 + a5)(b3 + b4 + b5) (X^8 + X^6 + X^4 + X^3) */
  mzd_add_accumulate(t1, 3, a3, a4, a5);
  mzd_add_accumulate(t2, 3, b3, b4, b5);
  mzd_mul(t0, t1, t2, 0);
  mzd_add_to_all(t0, 4,   X[ 8], X[ 6], X[ 4], X[ 3]);

  /* (a0 + a1 + a2)(b0 + b1 + b2) (X^7 + X^2) */
  mzd_add_accumulate(t1, 3, a0, a1, a2);
  mzd_add_accumulate(t2, 3, b0, b1, b2);
  mzd_mul(t0, t1, t2, 0);
  mzd_add_to_all(t0, 2,   X[ 7], X[ 2]);

  /* (a0 + a3 + a5)(b0 + b3 + b5) (X^7 + X^5) */
  mzd_add_accumulate(t1, 3, a0, a3, a5);
  mzd_add_accumulate(t2, 3, b0, b3, b5);
  mzd_mul(t0, t1, t2, 0);
  mzd_add_to_all(t0, 2,   X[ 7], X[ 5]);

  /* (a0 + a2 + a5)(b0 + b2 + b5) (X^6 + X^5 + X^4 + X^3) */
  mzd_add_accumulate(t1, 3, a0, a2, a5);
  mzd_add_accumulate(t2, 3, b0, b2, b5);
  mzd_mul(t0, t1, t2, 0);
  mzd_add_to_all(t0, 4,   X[ 6], X[ 5], X[ 4], X[ 3]);

  /* (a0 + a2 + a3 + a5)(b0 + b2 + b3 + b5) (X^7 + X^5 + X^4 + X^3) */
  mzd_add_accumulate(t1, 4, a0, a2, a3, a5); // TODO: can save adds here
  mzd_add_accumulate(t2, 4, b0, b2, b3, b5); // TODO: can save adds here
  mzd_mul(t0, t1, t2, 0);
  mzd_add_to_all(t0, 4,   X[ 7], X[ 5], X[ 4], X[ 3]);

  /* (a0 + a1 + a3 + a4)(b0 + b1 + b3 + b4) (X^6 + X^4) */
  mzd_add_accumulate(t1, 4, a0, a1, a3, a4);
  mzd_add_accumulate(t2, 4, b0, b1, b3, b4);
  mzd_mul(t0, t1, t2, 0);
  mzd_add_to_all(t0, 2,   X[ 6], X[ 4]);

  /* (a0 + a1 + a2 + a3 + a4 + a5)(b0 + b1 + b2 + b3 + b4 + b5) X^6 */
  mzd_add_accumulate(t1, 6, a0, a1, a2, a3, a4, a5); // TODO: can save adds here
  mzd_add_accumulate(t2, 6, b0, b1, b2, b3, b4, b5); // TODO: can save adds here
  mzd_mul(t0, t1, t2, 0);
  mzd_add(X[ 6], X[ 6], t0);

  mzd_add_modred6(X, 9, X[ 9], minpoly);
  mzd_add_modred6(X, 8, X[ 8], minpoly);
  mzd_add_modred6(X, 7, X[ 7], minpoly);
  mzd_add_modred6(X, 6, X[ 6], minpoly);

  mzd_free(t0);
  mzd_free(t1);
  mzd_free(t2);

  _mzd_slice_adapt_depth(C,6);
  return C;
}


mzd_slice_t *_mzd_slice_mul_karatsuba8(mzd_slice_t *C, const mzd_slice_t *A, const mzd_slice_t *B) {
  /** 8 + 7 temporaries **/
  if (C == NULL)
    C = mzd_slice_init(A->finite_field, A->nrows, B->ncols);

  const word minpoly = A->finite_field->minpoly;
  C = _mzd_slice_adapt_depth(C,15);

  const mzd_t *a0[4] = {A->x[0],A->x[1],A->x[2],A->x[3]};
  const mzd_t *a1[4] = {A->x[4],A->x[5],A->x[6],A->x[7]};
  const mzd_t *b0[4] = {B->x[0],B->x[1],B->x[2],B->x[3]};
  const mzd_t *b1[4] = {B->x[4],B->x[5],B->x[6],B->x[7]};

  mzd_t *X[3][7] = { {C->x[ 0],C->x[ 1],C->x[ 2],C->x[ 3],C->x[ 4],C->x[ 5],C->x[ 6]},
                     {C->x[ 4],C->x[ 5],C->x[ 6],C->x[ 7],C->x[ 8],C->x[ 9],C->x[10]},
                     {C->x[ 8],C->x[ 9],C->x[10],C->x[11],C->x[12],C->x[13],C->x[14]} };

  mzd_t *t0[7];
  mzd_t *t1[4];

  t0[0] = mzd_init(A->nrows, A->ncols);
  t0[1] = mzd_init(A->nrows, A->ncols);
  t0[2] = mzd_init(A->nrows, A->ncols);
  t0[3] = mzd_init(A->nrows, A->ncols);

  t1[0] = mzd_init(B->nrows, B->ncols);
  t1[1] = mzd_init(B->nrows, B->ncols);
  t1[2] = mzd_init(B->nrows, B->ncols);
  t1[3] = mzd_init(B->nrows, B->ncols);

  _poly_add(t0, a0, a1, 4);
  _poly_add(t1, b0, b1, 4);

  _poly_addmul4(X[1], (const mzd_t**)t0, (const mzd_t**)t1);

  mzd_free(t0[0]);  mzd_free(t0[1]);  mzd_free(t0[2]);  mzd_free(t0[3]);
  mzd_free(t1[0]);  mzd_free(t1[1]);  mzd_free(t1[2]);  mzd_free(t1[3]);

  t0[0] = mzd_init(C->x[0]->nrows, B->x[0]->ncols);
  t0[1] = mzd_init(C->x[0]->nrows, B->x[0]->ncols);
  t0[2] = mzd_init(C->x[0]->nrows, B->x[0]->ncols);
  t0[3] = mzd_init(C->x[0]->nrows, B->x[0]->ncols);
  t0[4] = mzd_init(C->x[0]->nrows, B->x[0]->ncols);
  t0[5] = mzd_init(C->x[0]->nrows, B->x[0]->ncols);
  t0[6] = mzd_init(C->x[0]->nrows, B->x[0]->ncols);

  _poly_addmul4(t0, a0, b0);
  _poly_add(X[0], (const mzd_t**)X[0], (const mzd_t**)t0, 7);
  _poly_add(X[1], (const mzd_t**)X[1], (const mzd_t**)t0, 7);

  mzd_set_ui(t0[0], 0);
  mzd_set_ui(t0[1], 0);
  mzd_set_ui(t0[2], 0);
  mzd_set_ui(t0[3], 0);
  mzd_set_ui(t0[4], 0);
  mzd_set_ui(t0[5], 0);
  mzd_set_ui(t0[6], 0);

  _poly_addmul4(t0, a1, b1);
  _poly_add(X[1], (const mzd_t**)X[1], (const mzd_t**)t0, 7);
  _poly_add(X[2], (const mzd_t**)X[2], (const mzd_t**)t0, 7);

  mzd_free(t0[0]); mzd_free(t0[1]); mzd_free(t0[2]); mzd_free(t0[3]); 
  mzd_free(t0[4]); mzd_free(t0[5]); mzd_free(t0[6]);

  for(unsigned int i=2*8-2; i>= 8; i--)
    for(unsigned int j=0; j<8; j++)
      if (minpoly & 1<<j)
        mzd_add(C->x[i-8+j], C->x[i-8+j], C->x[i]);
  _mzd_slice_adapt_depth(C,8);
  return C;
}

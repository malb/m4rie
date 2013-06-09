#include "mzd_ptr.h"
#include "m4ri_functions.h"

void _mzd_ptr_addmul_karatsuba2(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B) {
  const mzd_t *a0 = A[0];
  const mzd_t *a1 = A[1];

  const mzd_t *b0 = B[0];
  const mzd_t *b1 = B[1];

  mzd_t *t0 = mzd_init(a0->nrows, a0->ncols);
  mzd_t *t1 = mzd_init(b0->nrows, b0->ncols);

  mzd_addmul(X[1], mzd_add(t0, a0, a1), mzd_add(t1, b0, b1), 0); /* + (a0+a1)(b0+b1)X */

  mzd_free(t0);
  mzd_free(t1);

  t0 = mzd_init(a0->nrows, b0->ncols);
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a0, b0, 0), X, 2,  0, 1); /* + a0b0(1-X) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a1, b1, 0), X, 2,  1, 2); /* + a1b1(X+X^2) */
  mzd_free(t0);
}

void _mzd_ptr_addmul_karatsuba3(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B) {
  const mzd_t *a0 = A[0];
  const mzd_t *a1 = A[1];
  const mzd_t *a2 = A[2];

  const mzd_t *b0 = B[0];
  const mzd_t *b1 = B[1];
  const mzd_t *b2 = B[2];

  mzd_t *t0 = mzd_init(a0->nrows, b0->ncols);
  mzd_t *t1 = mzd_init(a0->nrows, a0->ncols);
  mzd_t *t2 = mzd_init(b0->nrows, b0->ncols);

  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_add(t1, a0, a1), mzd_add(t2, b0, b1), 0), X, 1,   1); /* + (a0+a1)(b0+b1)X */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_add(t1, a0, a2), mzd_add(t2, b0, b2), 0), X, 1,   2); /* + (a0+a2)(b0+b2)X^2 */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_add(t1, a1, a2), mzd_add(t2, b1, b2), 0), X, 1,   3); /* + (a1+a2)(b1+b2)X^3 */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a0, b0, 0), X, 3,   0, 1, 2); /* + a0b0(1-X-X^2) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a1, b1, 0), X, 3,   1, 2, 3); /* + a1b1(X+X^2-X^3) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a2, b2, 0), X, 3,   2, 3, 4); /* + a2b2(-X^2-X^3+X^4) */

  mzd_free(t0);
  mzd_free(t1);
  mzd_free(t2);
}

void _mzd_ptr_addmul_karatsuba4(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B) {
  const mzd_t *a0 = A[0];
  const mzd_t *a1 = A[1];
  const mzd_t *a2 = A[2];
  const mzd_t *a3 = A[3];

  const mzd_t *b0 = B[0];
  const mzd_t *b1 = B[1];
  const mzd_t *b2 = B[2];
  const mzd_t *b3 = B[3];

  mzd_t *t0 = mzd_init(a0->nrows, b0->ncols);
  mzd_t *t1 = mzd_init(a0->nrows, a0->ncols);
  mzd_t *t2 = mzd_init(b0->nrows, b0->ncols);

  /* (a0 + a1 + a2 + a3)*(b0 + b1 + b2 + b3)*X^3 */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 4, a0, a1, a2, a3),
                                      mzd_sum(t2, 4, b0, b1, b2, b3), 0), X, 1,   3);
  /* (a0 + a1)*(b0 + b1)*(X   + X^3) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a0, a1),
                                      mzd_sum(t2, 2, b0, b1), 0), X, 2,   1, 3);
  /* (a0 + a2)*(b0 + b2)*(X^2 + X^3) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a0, a2),
                                      mzd_sum(t2, 2, b0, b2), 0), X, 2,   2, 3);
  /* (a1 + a3)*(b1 + b3)*(X^3 + X^4) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a1, a3),
                                      mzd_sum(t2, 2, b1, b3), 0), X, 2,   3, 4);
  /* (a2 + a3)*(b2 + b3)*(X^3 + X^5) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a2, a3),
                                      mzd_sum(t2, 2, b2, b3), 0), X, 2,   3, 5);
  /* (a0*b0)*(1 + X   + X^2 + X^3) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a0, b0, 0), X, 4,    0, 1, 2, 3);
  /* (a1*b1)*(X + X^2 + X^3 + X^4) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a1, b1, 0), X, 4,    1, 2, 3, 4);
  /* (a2*b2)*(X^2 + X^3 + X^4 + X^5) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a2, b2, 0), X, 4,    2, 3, 4, 5);
  /* (a3*b3)*(X^3 + X^4 + X^5 + X^6) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a3, b3, 0), X, 4,    3, 4, 5, 6);

  mzd_free(t0);
  mzd_free(t1);
  mzd_free(t2);
}

void _mzd_ptr_addmul_karatsuba5(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B) {
  const mzd_t *a0 = A[0];
  const mzd_t *a1 = A[1];
  const mzd_t *a2 = A[2];
  const mzd_t *a3 = A[3];
  const mzd_t *a4 = A[4];

  const mzd_t *b0 = B[0];
  const mzd_t *b1 = B[1];
  const mzd_t *b2 = B[2];
  const mzd_t *b3 = B[3];
  const mzd_t *b4 = B[4];

  mzd_t *t0 = mzd_init(a0->nrows, b0->ncols);

  /* a0b0(X^6 + X^5 + X^4 + X^2 + X + 1)*/
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a0, b0, 0), X, 6,   6, 5, 4, 2, 1, 0);
  /* a1b1(X^4 + X)*/
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a1, b1, 0), X, 2,   4, 1);
  /* a3b3(X^7 + X^4) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a3, b3, 0), X, 2,   7, 4);
  /* (a4b4)(X^8 + X^7 + X^6 + X^4 + X^3 + X^2) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a4, b4, 0), X, 6,   8, 7, 6, 4, 3, 2);

  mzd_t *t1 = mzd_init(a0->nrows, a1->ncols);
  mzd_t *t2 = mzd_init(b0->nrows, b1->ncols);

  /* (a0+a4)(b0+b4)(X^6 + X^5 + X^3 + X^2) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a0, a4),
                                      mzd_sum(t2, 2, b0, b4), 0), X, 4,   6, 5, 3, 2);

  /* (a0+a1)(b0+b1)(X^5 + X^4 + X^2 + X) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a0, a1),
                                      mzd_sum(t2, 2, b0, b1), 0), X, 4,   5, 4, 2, 1);

  /* (a3+a4)(b3+b4)(X^7 + X^6 + X^4 + X^3) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a3, a4),
                                      mzd_sum(t2, 2, b3, b4), 0), X, 4,   7, 6, 4, 3);

  /* (a1+a2+a4)(a1+a2+a4)(X^4 + X^2) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 3, a1, a2, a4),
                                      mzd_sum(t2, 3, b1, b2, b4), 0), X, 2,   4, 2);

  /* (a0+a2+a3)(b0+b2+b3)(X^6 + X^4) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 3, a0, a2, a3),
                                      mzd_sum(t2, 3, b0, b2, b3), 0), X, 2,   6, 4);

  /* (a0+a1+a3+a4)(b0+b1+b3+b4)(X^5 + X^3) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 4, a0, a1, a3, a4),
                                      mzd_sum(t2, 4, b0, b1, b3, b4), 0), X, 2,   5, 3);

  /* (a0+a1+a2+a4)(b0+b1+b2+b4)(X^5 + X^2) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 4, a0, a1, a2, a4),
                                      mzd_sum(t2, 4, b0, b1, b2, b4), 0), X, 2,   5, 2);

  /* (a0+a2+a3+a4)(b0+b2+b3+b4)(X^6 + X^3)*/
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 4, a0, a2, a3, a4),
                                      mzd_sum(t2, 4, b0, b2, b3, b4), 0), X, 2,   6, 3);

  /* (a0+a1+a2+a3+a4)(b0+b1+b2+b3+b4)(X^5 + X^4 + X^3)*/
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 5, a0, a1, a2, a3, a4),
                                      mzd_sum(t2, 5, b0, b1, b2, b3, b4), 0), X, 3,   5, 4, 3);

  mzd_free(t0);
  mzd_free(t1);
  mzd_free(t2);
}

void _mzd_ptr_addmul_karatsuba6(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B) {
  /* using 22 multiplications and three temporaries */

  const mzd_t *a0 = A[0];
  const mzd_t *a1 = A[1];
  const mzd_t *a2 = A[2];
  const mzd_t *a3 = A[3];
  const mzd_t *a4 = A[4];
  const mzd_t *a5 = A[5];

  const mzd_t *b0 = B[0];
  const mzd_t *b1 = B[1];
  const mzd_t *b2 = B[2];
  const mzd_t *b3 = B[3];
  const mzd_t *b4 = B[4];
  const mzd_t *b5 = B[5];

  mzd_t *t0 = mzd_init(a0->nrows, b0->ncols);
  mzd_t *t1 = mzd_init(a0->nrows, a1->ncols);
  mzd_t *t2 = mzd_init(b0->nrows, b1->ncols);

  /* a0b0 (X^6 + X^5 + X + 1) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a0, b0, 0), X, 4,   6, 5, 1, 0);
  /* a1b1 (X^7 + X^6 + X^5 + X^4 + X^3 + X) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a1, b1, 0), X, 6,   7, 6, 5, 4, 3, 1);
  /* a4b4 (X^9 + X^7 + X^5 + X^3) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a4, b4, 0), X, 4,   9, 7, 5, 3);
  /* a5b5 (X^10 + X^9 + X^6 + X^5) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a5, b5, 0), X, 4,   10, 9, 6, 5);

  /* (a4 + a5)(b4 + b5) (X^9 + X^8 + X^4+ X^3) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a4, a5),
                                      mzd_sum(t2, 2, b4, b5), 0), X, 4,   9, 8, 4, 3);
  /* (a0 + a1)(b0 + b1) (X^7 + X^4 + X^2 + X) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a0, a1),
                                      mzd_sum(t2, 2, b0, b1), 0), X, 4,   7, 4, 2, 1);
  /* (a3 + a4)(b3 + b4)(X^8 + X^7 + X^6 + X^3) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a3, a4),
                                      mzd_sum(t2, 2, b3, b4), 0), X, 4,   8, 7, 6, 3);
  /* (a1 + a2)(b1 + b2) (X^7 + X^6 + X^3 + X^2) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a1, a2),
                                      mzd_sum(t2, 2, b1, b2), 0), X, 4,   7, 6, 3, 2);
  /* (a1 + a4)(b1 + b4) (X^4 + X^5) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a1, a4),
                                      mzd_sum(t2, 2, b1, b4), 0), X, 2,   4, 5);
  /* (a2 + a3)(b2 + b3) (X^7 + X^6 + X^4 + X^3) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a2, a3),
                                      mzd_sum(t2, 2, b2, b3), 0), X, 4,   7, 6, 4, 3);
  /* (a3 + a4 + a5)(b3 + b4 + b5) (X^8 + X^6 + X^4 + X^3) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 3, a3, a4, a5),
                                      mzd_sum(t2, 3, b3, b4, b5), 0), X, 4,   8, 6, 4, 3);
  /* (a0 + a1 + a2)(b0 + b1 + b2) (X^7 + X^2) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 3, a0, a1, a2),
                                      mzd_sum(t2, 3, b0, b1, b2), 0), X, 2,   7, 2);
  /* (a0 + a3 + a5)(b0 + b3 + b5) (X^7 + X^5) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 3, a0, a3, a5),
                                      mzd_sum(t2, 3, b0, b3, b5), 0), X, 2,   7, 5);
  /* (a0 + a2 + a5)(b0 + b2 + b5) (X^6 + X^5 + X^4 + X^3) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 3, a0, a2, a5),
                                      mzd_sum(t2, 3, b0, b2, b5), 0), X, 4,   6, 5, 4, 3);
  /* (a0 + a2 + a3 + a5)(b0 + b2 + b3 + b5) (X^7 + X^5 + X^4 + X^3) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 4, a0, a2, a3, a5),
                                      mzd_sum(t2, 4, b0, b2, b3, b5), 0), X, 4,   7, 5, 4, 3);
  /* (a0 + a1 + a3 + a4)(b0 + b1 + b3 + b4) (X^6 + X^4) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 4, a0, a1, a3, a4),
                                      mzd_sum(t2, 4, b0, b1, b3, b4), 0), X, 2,   6, 4);
  /* (a0 + a1 + a2 + a3 + a4 + a5)(b0 + b1 + b2 + b3 + b4 + b5) X^6 */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 6, a0, a1, a2, a3, a4, a5),
                                      mzd_sum(t2, 6, b0, b1, b2, b3, b4, b5), 0), X, 1,   6);

  mzd_free(t0);
  mzd_free(t1);
  mzd_free(t2);
}

void _mzd_ptr_addmul_karatsuba7(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B) {
  /* using three temporaries */
  const mzd_t *a0 = A[0];
  const mzd_t *a1 = A[1];
  const mzd_t *a2 = A[2];
  const mzd_t *a3 = A[3];
  const mzd_t *a4 = A[4];
  const mzd_t *a5 = A[5];
  const mzd_t *a6 = A[6];

  const mzd_t *b0 = B[0];
  const mzd_t *b1 = B[1];
  const mzd_t *b2 = B[2];
  const mzd_t *b3 = B[3];
  const mzd_t *b4 = B[4];
  const mzd_t *b5 = B[5];
  const mzd_t *b6 = B[6];

  mzd_t *t0 = mzd_init(a0->nrows, b0->ncols);
  mzd_t *t1 = mzd_init(a0->nrows, a1->ncols);
  mzd_t *t2 = mzd_init(b0->nrows, b1->ncols);

  /* (a0 + a1 + a2 + a3 + a4 + a5 + a6)(b0 + b1 + b2 + b3 + b4 + b5 + b6)(X^7 + X^6 + X^5) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 7, a0, a1, a2, a3, a4, a5, a6),
                                      mzd_sum(t2, 7, b0, b1, b2, b3, b4, b5, b6), 0),
                        X, 3,   7, 6, 5);

  /* (a1 + a2 + a3 + a5 + a6)(b1 + b2 + b3 + b5 + b6)(X^9 + X^6) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 5, a1, a2, a3, a5, a6),
                                      mzd_sum(t2, 5, b1, b2, b3, b5, b6), 0),
                        X, 2,   9, 6);

  /* (a0 + a1 + a3 + a4 + a5)(b0 + b1 + b3 + b4 + b5)(X^6 + X^3) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 5, a0, a1, a3, a4, a5),
                                      mzd_sum(t2, 5, b0, b1, b3, b4, b5), 0),
                        X, 2,   6, 3);

  /* (a0 + a2 + a3 + a4 + a6)(b0 + b2 + b3 + b4 + b6)(X^9 + X^3) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 5, a0, a2, a3, a4, a6),
                                      mzd_sum(t2, 5, b0, b2, b3, b4, b6), 0),
                        X, 2,   9, 3);

  /* (a0 + a2 + a3 + a5 + a6)(b0 + b2 + b3 + b5 + b6)(X^7 + X^3) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 5, a0, a2, a3, a5, a6),
                                      mzd_sum(t2, 5, b0, b2, b3, b5, b6), 0),
                        X, 2,   7, 3);

  /* (a0 + a1 + a3 + a4 + a6)(b0 + b1 + b3 + b4 + b6)(X^9 + X^5) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 5, a0, a1, a3, a4, a6),
                                      mzd_sum(t2, 5, b0, b1, b3, b4, b6), 0),
                        X, 2,   9, 5);

  /* (a1 + a2 + a4 + a5)(b1 + b2 + b4 + b5)(X^9 + X^7 + X^5 + X^3) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 4, a1, a2, a4, a5),
                                      mzd_sum(t2, 4, b1, b2, b4, b5), 0),
                        X, 4,   9, 7, 5, 3);

  /* (a0 + a1)(b0 + b1)(X^9 + X^7 + X^3 + X) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a0, a1),
                                      mzd_sum(t2, 2, b0, b1), 0),
                        X, 4,   9, 7, 3, 1);

  /* (a0 + a2)(b0 + b2)(X^9 + X^6 + X^5 + X^2) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a0, a2),
                                      mzd_sum(t2, 2, b0, b2), 0),
                        X, 4,   9, 6, 5, 2);

  /* (a0 + a4)(b0 + b4)(X^7 + X^4) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a0, a4),
                                      mzd_sum(t2, 2, b0, b4), 0),
                        X, 2,   7, 4);

  /* (a1 + a3)(b1 + b3)(X^7 + X^6 + X^4 + X^3) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a1, a3),
                                      mzd_sum(t2, 2, b1, b3), 0),
                        X, 4,   7, 6, 4, 3);

  /* (a2 + a6)(b2 + b6)(X^8 + X^5) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a2, a6),
                                      mzd_sum(t2, 2, b2, b6), 0),
                        X, 2,   8, 5);

  /* (a3 + a5)(b3 + b5)(X^9 + X^8 + X^6 + X^5) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a3, a5),
                                      mzd_sum(t2, 2, b3, b5), 0),
                        X, 4,   9, 8, 6, 5);

  /* (a4 + a6)(b4 + b6)(X^10 + X^7 + X^6 + X^3) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a4, a6),
                                      mzd_sum(t2, 2, b4, b6), 0),
                        X, 4,  10, 7, 6, 3);

  /* (a5 + a6)(b5 + b6)(X^11 + X^9 + X^5 + X^3) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, mzd_sum(t1, 2, a5, a6),
                                      mzd_sum(t2, 2, b5, b6), 0),
                        X, 4,  11, 9, 5, 3);

  /* a0b0(X^6 + X^5 + X^4 + X^2 + X + 1) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a0, b0, 0), X, 6,   6, 5, 4, 2, 1, 0);
  /* a1b1(X^5 + X^4 + X^2 + X) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a1, b1, 0), X, 4,   5, 4, 2, 1);
  /* a2b2(X^8 + X^7 + X^6 + X^4 + X^3 + X^2) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a2, b2, 0), X, 6,   8, 7, 6, 4, 3, 2);
  /* a3b3(X^8 + X^7 + X^5 + X^4) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a3, b3, 0), X, 4,   8, 7, 5, 4);
  /* a4b4(X^10 + X^9 + X^8 + X^6 + X^5 + X^4)*/
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a4, b4, 0), X, 6,  10, 9, 8, 6, 5, 4);;
  /* a5b5(X^11 + X^10 + X^8 + X^7) */
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a5, b5, 0), X, 4,  11, 10, 8, 7);
  /* a6b6(X^12 + X^11 + X^10 + X^8 + X^7 + X^6)*/
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a6, b6, 0), X, 6,  12, 11, 10, 8, 7, 6);

  mzd_free(t0);
  mzd_free(t1);
  mzd_free(t2);

}

void _mzd_ptr_addmul_karatsuba8(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B) {
  /* using three temporaries */
  const mzd_t *a0 = A[0];
  const mzd_t *a1 = A[1];
  const mzd_t *a2 = A[2];
  const mzd_t *a3 = A[3];
  const mzd_t *a4 = A[4];
  const mzd_t *a5 = A[5];
  const mzd_t *a6 = A[6];
  const mzd_t *a7 = A[7];

  const mzd_t *b0 = B[0];
  const mzd_t *b1 = B[1];
  const mzd_t *b2 = B[2];
  const mzd_t *b3 = B[3];
  const mzd_t *b4 = B[4];
  const mzd_t *b5 = B[5];
  const mzd_t *b6 = B[6];
  const mzd_t *b7 = B[7];

  mzd_t *t0 = mzd_init(a0->nrows, b0->ncols);
  mzd_t *t1 = mzd_init(a0->nrows, a0->ncols);
  mzd_t *t2 = mzd_init(b0->nrows, b0->ncols);

  // a0*b0 * (X^0 + X^1 + X^2 + X^3 + X^4 + X^5 + X^6 + X^7) 
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a0, b0, 0), X, 8,  0, 1, 2, 3, 4, 5, 6, 7);
  // a1*b1 * (X^1 + X^2 + X^3 + X^4 + X^5 + X^6 + X^7 + X^8) 
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a1, b1, 0), X, 8,  1, 2, 3, 4, 5, 6, 7, 8);
  // a2*b2 * (X^2 + X^3 + X^4 + X^5 + X^6 + X^7 + X^8 + X^9) 
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a2, b2, 0), X, 8,  2, 3, 4, 5, 6, 7, 8, 9);
  // a3*b3 * (X^3 + X^4 + X^5 + X^6 + X^7 + X^8 + X^9 + X^10) 
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a3, b3, 0), X, 8,  3, 4, 5, 6, 7, 8, 9, 10);
  // a4*b4 * (X^4 + X^5 + X^6 + X^7 + X^8 + X^9 + X^10 + X^11) 
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a4, b4, 0), X, 8,  4, 5, 6, 7, 8, 9, 10, 11);
  // a5*b5 * (X^5 + X^6 + X^7 + X^8 + X^9 + X^10 + X^11 + X^12) 
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a5, b5, 0), X, 8,  5, 6, 7, 8, 9, 10, 11, 12);
  // a6*b6 * (X^6 + X^7 + X^8 + X^9 + X^10 + X^11 + X^12 + X^13) 
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a6, b6, 0), X, 8,  6, 7, 8, 9, 10, 11, 12, 13);
  // a7*b7 * (X^7 + X^8 + X^9 + X^10 + X^11 + X^12 + X^13 + X^14) 
  _mzd_ptr_add_to_all(ff, mzd_mul(t0, a7, b7, 0), X, 8,  7, 8, 9, 10, 11, 12, 13, 14);

  // (b0 + b2)*(a0 + a2) * (X^2 + X^3 + X^6 + X^7)
  _mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a0, a2),
                                       mzd_sum(t2, 2, b0, b2), 0), X, 4,  2, 3, 6, 7);
  // (b2 + b3)*(a2 + a3) * (X^3 + X^5 + X^7 + X^9)
  _mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a2, a3),
                                       mzd_sum(t2, 2, b2, b3), 0), X, 4,  3, 5, 7, 9);
  // (b2 + b6)*(a2 + a6) * (X^6 + X^7 + X^8 + X^9)
  _mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a2, a6),
                                       mzd_sum(t2, 2, b2, b6), 0), X, 4,  6, 7, 8, 9);
  // (b4 + b5)*(a4 + a5) * (X^5 + X^7 + X^9 + X^11)
  _mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a4, a5),
                                       mzd_sum(t2, 2, b4, b5), 0), X, 4,  5, 7, 9, 11);
  // (b6 + b7)*(a6 + a7) * (X^7 + X^9 + X^11 + X^13)
  _mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a6, a7),
                                       mzd_sum(t2, 2, b6, b7), 0), X, 4,  7, 9, 11, 13);
  // (b0 + b1 + b2 + b3)*(a0 + a1 + a2 + a3) * (X^3 + X^7)
  _mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a0, a1, a2, a3),
                                       mzd_sum(t2, 4, b0, b1, b2, b3), 0), X, 2,  3, 7);
  // (b5 + b7)*(a5 + a7) * (X^7 + X^8 + X^11 + X^12)
  _mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a5, a7),
                                       mzd_sum(t2, 2, b5, b7), 0), X, 4,  7, 8, 11, 12);
  // (b0 + b4)*(a0 + a4) * (X^4 + X^5 + X^6 + X^7)
  _mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a0, a4),
                                       mzd_sum(t2, 2, b0, b4), 0), X, 4,  4, 5, 6, 7);
  // (b1 + b3)*(a1 + a3) * (X^3 + X^4 + X^7 + X^8)
  _mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a1, a3),
                                       mzd_sum(t2, 2, b1, b3), 0), X, 4,  3, 4, 7, 8);
  // (b1 + b3 + b5 + b7)*(a1 + a3 + a5 + a7) * (X^7 + X^8)
  _mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a1, a3, a5, a7),
                                       mzd_sum(t2, 4, b1, b3, b5, b7), 0), X, 2,  7, 8);
  // (b0 + b1)*(a0 + a1) * (X^1 + X^3 + X^5 + X^7)
  _mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a0, a1),
                                       mzd_sum(t2, 2, b0, b1), 0), X, 4,  1, 3, 5, 7);
  // (b4 + b6)*(a4 + a6) * (X^6 + X^7 + X^10 + X^11)
  _mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a4, a6),
                                       mzd_sum(t2, 2, b4, b6), 0), X, 4,  6, 7, 10, 11);
  // (b1 + b5)*(a1 + a5) * (X^5 + X^6 + X^7 + X^8)
  _mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a1, a5),
                                       mzd_sum(t2, 2, b1, b5), 0), X, 4,  5, 6, 7, 8);
  // (b0 + b1 + b2 + b3 + b4 + b5 + b6 + b7)*(a0 + a1 + a2 + a3 + a4 + a5 + a6 + a7) * (X^7)
  _mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 8, a0, a1, a2, a3, a4, a5, a6, a7),
                                       mzd_sum(t2, 8, b0, b1, b2, b3, b4, b5, b6, b7), 0), X, 1,  7);
  // (b0 + b2 + b4 + b6)*(a0 + a2 + a4 + a6) * (X^6 + X^7)
  _mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a0, a2, a4, a6),
                                       mzd_sum(t2, 4, b0, b2, b4, b6), 0), X, 2,  6, 7);
  // (b0 + b1 + b4 + b5)*(a0 + a1 + a4 + a5) * (X^5 + X^7)
  _mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a0, a1, a4, a5),
                                       mzd_sum(t2, 4, b0, b1, b4, b5), 0), X, 2,  5, 7);
  // (b2 + b3 + b6 + b7)*(a2 + a3 + a6 + a7) * (X^7 + X^9)
  _mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a2, a3, a6, a7),
                                       mzd_sum(t2, 4, b2, b3, b6, b7), 0), X, 2,  7, 9);
  // (b3 + b7)*(a3 + a7) * (X^7 + X^8 + X^9 + X^10)
  _mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 2, a3, a7),
                                       mzd_sum(t2, 2, b3, b7), 0), X, 4,  7, 8, 9, 10);
  // (b4 + b5 + b6 + b7)*(a4 + a5 + a6 + a7) * (X^7 + X^11)
  _mzd_ptr_add_to_all(ff,  mzd_mul(t0, mzd_sum(t1, 4, a4, a5, a6, a7),
                                       mzd_sum(t2, 4, b4, b5, b6, b7), 0), X, 2,  7, 11);

  mzd_free(t0);
  mzd_free(t1);
  mzd_free(t2);
}


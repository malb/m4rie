#include "mzd_poly.h"

void _poly_addmul2(mzd_t **X, const mzd_t **a, const mzd_t **b) {
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

void _poly_addmul4(mzd_t **c, const mzd_t **a, const mzd_t **b) {
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

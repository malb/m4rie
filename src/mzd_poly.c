#include "mzd_poly.h"
#include "mzd_slice.h"
#include "newton_john.h"

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

static inline mzd_slice_t * mzd_slice_addmul_mzd(mzd_slice_t *C, const word a, const mzd_t *A) {
  for(int i=0; i<C->depth; i++)
    if(a & 1ULL<<i)
      mzd_add(C->x[i], C->x[i], A);
  return C;
}

mzd_poly_t *_mzd_poly_addmul1(mzd_poly_t *C, mzd_poly_t *A, mzd_poly_t *B) {
  const int d = A->depth + B->depth-1;
  const int log2d = (int)ceil(log2( (double)d ));
  assert(log2d <= 16);

  gf2e *ff = gf2e_init(irreducible_polynomials[log2d][1]);

  mzd_slice_t  *a = mzd_slice_init(ff, A->nrows, A->ncols);
  mzd_slice_t  *b = mzd_slice_init(ff, B->nrows, B->ncols);
  mzd_slice_t **c = (mzd_slice_t**)calloc(sizeof(mzd_slice_t*),d);

  mzed_t *Phi = mzed_init(ff, d, d);
  mzed_t *Rho = mzed_init(ff, d, d);

  /* evaluation at zero */
  mzed_write_elem(Phi, 0, 0, 1);
  c[0] = mzd_slice_init(ff, A->nrows, B->ncols);
  mzd_mul(c[0]->x[0], A->x[0], B->x[0], 0);

  /* evaluation at one */
  for(int i=0; i<d; i++)
    mzed_write_elem(Phi, 1, i, 1);
  for (int i = 0; i < A->depth; i++)
    mzd_add(a->x[0], a->x[0], A->x[i]);
  for (int i = 0; i < B->depth; i++)
    mzd_add(b->x[0], b->x[0], B->x[i]);

  c[1] = mzd_slice_init(ff, A->nrows, B->ncols);
  mzd_mul(c[1]->x[0], a->x[0], b->x[0], 0);

  /* evaluation at 2 ... d-1 */
  for(int i=2; i<d; i++) {
    word alpha = i;
    word acc = 1;
    mzd_slice_set_ui(a, 0);
    mzd_slice_set_ui(b, 0);
    for(int j=0; j<d; j++) {
      mzed_write_elem(Phi, i, j, acc);
      if (j < A->depth) mzd_slice_addmul_mzd(a, acc, A->x[j]);
      if (j < B->depth) mzd_slice_addmul_mzd(b, acc, B->x[j]);
      acc = ff->mul(ff, alpha, acc);
    }
    c[i] = mzd_slice_mul(NULL, a, b);
  }
  mzd_slice_free(a);
  mzd_slice_free(b);

  mzed_invert_newton_john(Rho, Phi);

  mzd_slice_t *tmp = mzd_slice_init(ff, A->nrows, B->ncols);
  for(int i=0; i<d; i++) {
    for(int j=0; j<d; j++) {
      mzd_slice_mul_scalar(tmp, mzed_read_elem(Rho, i, j), c[j]);
      mzd_add(C->x[i], C->x[i], tmp->x[0]);
    }
  }
  mzd_slice_free(tmp);

  for(int i=0; i<d; i++) {
    mzd_slice_free(c[i]);
  }
  free(c);

  mzed_free(Phi);
  mzed_free(Rho);
  gf2e_free(ff);
  return C;
}

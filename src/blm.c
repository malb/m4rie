#include <m4ri/djb.h>
#include <m4rie/blm.h>
#include <m4rie/mzd_ptr.h>


void _mzd_ptr_apply_blm_mzd(mzd_t **X, const mzd_t **A, const mzd_t **B, const blm_t *f) {
  assert((f->H!=NULL) & (f->F!=NULL) & (f->G!=NULL) &
         (f->H->ncols == f->F->nrows) &   
         (f->F->nrows == f->G->nrows));

  mzd_t *t0 = mzd_init(A[0]->nrows, B[0]->ncols);
  mzd_t *t1 = mzd_init(A[0]->nrows, A[0]->ncols);
  mzd_t *t2 = mzd_init(B[0]->nrows, B[0]->ncols);

  for(rci_t i=0; i < f->F->nrows; i++) {
    mzd_set_ui(t1, 0);
    for(rci_t j=0; j < f->F->ncols; j++) {
      if(mzd_read_bit(f->F, i, j)) {
        mzd_add(t1, t1, A[j]);
      }
    }

    mzd_set_ui(t2, 0);
    for(rci_t j=0; j < f->G->ncols; j++) {
      if(mzd_read_bit(f->G, i, j)) {
        mzd_add(t2, t2, B[j]);
      }
    }


    mzd_mul(t0, t1, t2, 0);

    for(rci_t j=0; j < f->H->nrows; j++)
      if(mzd_read_bit(f->H, j, i))
        _mzd_ptr_add_modred(NULL, t0, X, j);
  }

  mzd_free(t0);
  mzd_free(t1);
  mzd_free(t2);
}

blm_t *_blm_djb_compile(blm_t *f) {
  assert((f->f == NULL) && (f->g == NULL) && (f->h == NULL));
  assert((f->F != NULL) && (f->G != NULL) && (f->H != NULL));


  mzd_t *F = mzd_copy(NULL, f->F);
  f->f = djb_compile(F);
  mzd_free(F);
  if(mzd_equal(f->F, f->G))
    f->g = f->f;
  else {
    mzd_t *G = mzd_copy(NULL, f->G);
    f->g = djb_compile(G);
    mzd_free(G);
  }
  mzd_t *H = mzd_copy(NULL, f->H);
  f->h = djb_compile(H);
  mzd_free(H);

  return f;
}

void _mzd_ptr_apply_blm_djb(mzd_t **X, const mzd_t **A, const mzd_t **B, const blm_t *f) {
  assert((f->H!=NULL) & (f->F!=NULL) & (f->G!=NULL) &   \
         (f->H->ncols == f->F->nrows) &                 \
         (f->F->nrows == f->G->nrows));

 
  mzd_t **t0 = (mzd_t**)m4ri_mm_malloc(sizeof(mzd_t*)*f->F->nrows);
  mzd_t **t1 = (mzd_t**)m4ri_mm_malloc(sizeof(mzd_t*)*f->F->nrows);
  mzd_t **t2 = (mzd_t**)m4ri_mm_malloc(sizeof(mzd_t*)*f->F->nrows);

  for(rci_t i=0; i<f->F->nrows; i++) {
    t1[i] = mzd_init(A[0]->nrows, A[0]->ncols);
    t2[i] = mzd_init(B[0]->nrows, B[0]->ncols);
  }

  djb_apply_mzd_ptr(f->f, t1, A);
  djb_apply_mzd_ptr(f->g, t2, B);

  for(rci_t i=0; i<f->F->nrows; i++) {
    t0[i] = mzd_init(A[0]->nrows, B[0]->ncols);
    mzd_mul(t0[i], t1[i], t2[i], 0);
    mzd_free(t1[i]);
    mzd_free(t2[i]);
  }

  

  djb_apply_mzd_ptr(f->h, X, (const mzd_t**)t0);

  for(rci_t i=0; i<f->F->nrows; i++) {
    mzd_free(t0[i]);
  }

  m4ri_mm_free(t0);
  m4ri_mm_free(t1);
  m4ri_mm_free(t2);
}



int *crt_init(const deg_t f_len, const deg_t g_len) {
  int *p_best = (int*)m4ri_mm_calloc(M4RIE_CRT_LEN, sizeof(int));
  int  c_best = f_len * g_len;

  int *p = (int*)m4ri_mm_calloc(M4RIE_CRT_LEN, sizeof(int));

  for(deg_t omega=0; omega<8; omega++) {

    deg_t deg_need = f_len+g_len-1-omega;
    deg_t deg_have = 0;
    deg_t deg_poly = 1;

    p[0] = omega;
    for(deg_t d=1; d<M4RIE_CRT_LEN; d++)
      p[d] = 0;

    while (deg_have < deg_need) {
      p[deg_poly] = irreducible_polynomials[deg_poly][0];
      if (deg_have + deg_poly * p[deg_poly] < deg_need) {
        deg_have += deg_poly * p[deg_poly];
      } else {
        deg_t deg_diff = deg_need - deg_have;
        p[deg_poly] = ceil(deg_diff/(double)deg_poly);
        deg_have += deg_poly * p[deg_poly];
      }
      deg_poly ++;
    }

    deg_t deg_diff = deg_have - deg_need;
    if (deg_diff && p[deg_diff] > 0) {
      p[deg_diff]--;
      deg_have -= deg_diff;
    }

    int c = costs[p[0]];
    for(deg_t d=1; d<M4RIE_CRT_LEN; d++)
      c += costs[d] * p[d];

    if (c < c_best) {
      for(deg_t d=0; d<M4RIE_CRT_LEN; d++)
        p_best[d] = p[d];
      c_best = c;
    }
  }  
  m4ri_mm_free(p);
  return p_best;
}

mzd_t *_small_multiplication_map(const deg_t degree) {
  mzd_t *A;
  switch(degree) {
  case 1:
    A = mzd_init(1, 1);
    mzd_row(A, 0)[0] = 0x1;
    return A;
  case 2:
    A = mzd_init(3, 2);
    mzd_row(A, 0)[0] = 0x1; 
    mzd_row(A, 1)[0] = 0x3; 
    mzd_row(A, 2)[0] = 0x2;
    return A;
  case 3:
    A = mzd_init(6, 3);
    mzd_row(A, 0)[0] = 0x1; 
    mzd_row(A, 1)[0] = 0x2; 
    mzd_row(A, 2)[0] = 0x4;
    mzd_row(A, 3)[0] = 0x3; 
    mzd_row(A, 4)[0] = 0x5; 
    mzd_row(A, 5)[0] = 0x6;
    return A;
  case 4:
    A = mzd_init(9, 4);
    mzd_row(A, 0)[0] = 0x1;
    mzd_row(A, 1)[0] = 0x2;
    mzd_row(A, 2)[0] = 0x4;
    mzd_row(A, 3)[0] = 0x8; 
    mzd_row(A, 4)[0] = 0xf;
    mzd_row(A, 5)[0] = 0x3;
    mzd_row(A, 6)[0] = 0x5;
    mzd_row(A, 7)[0] = 0xa;
    mzd_row(A, 8)[0] = 0xc;
    return A;
  case 5:
    A = mzd_init(13, 5);
    mzd_row(A,  0)[0] = 0x01;
    mzd_row(A,  1)[0] = 0x02; 
    mzd_row(A,  2)[0] = 0x08; 
    mzd_row(A,  3)[0] = 0x10; 
    mzd_row(A,  4)[0] = 0x11;
    mzd_row(A,  5)[0] = 0x03; 
    mzd_row(A,  6)[0] = 0x18; 
    mzd_row(A,  7)[0] = 0x16;
    mzd_row(A,  8)[0] = 0x0d;
    mzd_row(A,  9)[0] = 0x1b; 
    mzd_row(A, 10)[0] = 0x17; 
    mzd_row(A, 11)[0] = 0x1d;
    mzd_row(A, 12)[0] = 0x1f;
    return A;
  case 6:
    A = mzd_init(17,  6);
    mzd_row(A,  0)[0] = 0x01;
    mzd_row(A,  1)[0] = 0x02;
    mzd_row(A,  2)[0] = 0x10;
    mzd_row(A,  3)[0] = 0x20;
    mzd_row(A,  4)[0] = 0x30;
    mzd_row(A,  5)[0] = 0x03;
    mzd_row(A,  6)[0] = 0x18;
    mzd_row(A,  7)[0] = 0x06;
    mzd_row(A,  8)[0] = 0x12;
    mzd_row(A,  9)[0] = 0x0c;
    mzd_row(A, 10)[0] = 0x38;
    mzd_row(A, 11)[0] = 0x07;
    mzd_row(A, 12)[0] = 0x29;
    mzd_row(A, 13)[0] = 0x25;
    mzd_row(A, 14)[0] = 0x2d;
    mzd_row(A, 15)[0] = 0x1b;
    mzd_row(A, 16)[0] = 0x3f;
    return A;
  case 7:
    A = mzd_init(22,  7);
    mzd_row(A,  0)[0] = 0x7f;
    mzd_row(A,  1)[0] = 0x6e;
    mzd_row(A,  2)[0] = 0x3b;
    mzd_row(A,  3)[0] = 0x5d;
    mzd_row(A,  4)[0] = 0x6d;
    mzd_row(A,  5)[0] = 0x5b;
    mzd_row(A,  6)[0] = 0x36;
    mzd_row(A,  7)[0] = 0x03;
    mzd_row(A,  8)[0] = 0x05;
    mzd_row(A,  9)[0] = 0x11;
    mzd_row(A, 10)[0] = 0x0a;
    mzd_row(A, 11)[0] = 0x44;
    mzd_row(A, 12)[0] = 0x28;
    mzd_row(A, 13)[0] = 0x50;
    mzd_row(A, 14)[0] = 0x60;
    mzd_row(A, 15)[0] = 0x01;
    mzd_row(A, 16)[0] = 0x02;
    mzd_row(A, 17)[0] = 0x04;
    mzd_row(A, 18)[0] = 0x08;
    mzd_row(A, 19)[0] = 0x10;
    mzd_row(A, 20)[0] = 0x20;
    mzd_row(A, 21)[0] = 0x40;
    return A;
  case 8:
    A = mzd_init(27,  8);
    mzd_row(A,  0)[0] = 0x01;
    mzd_row(A,  1)[0] = 0x02;
    mzd_row(A,  2)[0] = 0x04;     
    mzd_row(A,  3)[0] = 0x08;
    mzd_row(A,  4)[0] = 0x10;
    mzd_row(A,  5)[0] = 0x20;
    mzd_row(A,  6)[0] = 0x40;     
    mzd_row(A,  7)[0] = 0x80;
    mzd_row(A,  8)[0] = 0x05;
    mzd_row(A,  9)[0] = 0x0c;
    mzd_row(A, 10)[0] = 0x44;     
    mzd_row(A, 11)[0] = 0x30;
    mzd_row(A, 12)[0] = 0xc0;
    mzd_row(A, 13)[0] = 0x0f;
    mzd_row(A, 14)[0] = 0xa0;     
    mzd_row(A, 15)[0] = 0x11;
    mzd_row(A, 16)[0] = 0x0a;
    mzd_row(A, 17)[0] = 0xaa;
    mzd_row(A, 18)[0] = 0x03;     
    mzd_row(A, 19)[0] = 0x50;
    mzd_row(A, 20)[0] = 0x22;
    mzd_row(A, 21)[0] = 0xff;
    mzd_row(A, 22)[0] = 0x55;     
    mzd_row(A, 23)[0] = 0x33;
    mzd_row(A, 24)[0] = 0xcc;
    mzd_row(A, 25)[0] = 0x88;
    mzd_row(A, 26)[0] = 0xf0;
    return A;
  case 9:
    A = mzd_init(31,  9);
    mzd_row(A,  0)[0] = 0x100;
    mzd_row(A,  1)[0] = 0x001;
    mzd_row(A,  2)[0] = 0x002;
    mzd_row(A,  3)[0] = 0x003;
    mzd_row(A,  4)[0] = 0x155;
    mzd_row(A,  5)[0] = 0x0aa;
    mzd_row(A,  6)[0] = 0x1ff;
    mzd_row(A,  7)[0] = 0x16d;
    mzd_row(A,  8)[0] = 0x1b6;
    mzd_row(A,  9)[0] = 0x0db;
    mzd_row(A, 10)[0] = 0x0e9;
    mzd_row(A, 11)[0] = 0x13a;
    mzd_row(A, 12)[0] = 0x074;
    mzd_row(A, 13)[0] = 0x1d3;
    mzd_row(A, 14)[0] = 0x09d;
    mzd_row(A, 15)[0] = 0x14e;
    mzd_row(A, 16)[0] = 0x0b9;
    mzd_row(A, 17)[0] = 0x172;
    mzd_row(A, 18)[0] = 0x05c;
    mzd_row(A, 19)[0] = 0x1cb;
    mzd_row(A, 20)[0] = 0x0e5;
    mzd_row(A, 21)[0] = 0x12e;
    mzd_row(A, 22)[0] = 0x191;
    mzd_row(A, 23)[0] = 0x0b2;
    mzd_row(A, 24)[0] = 0x164;
    mzd_row(A, 25)[0] = 0x0c8;
    mzd_row(A, 26)[0] = 0x08f;
    mzd_row(A, 27)[0] = 0x123;
    mzd_row(A, 28)[0] = 0x0f5;
    mzd_row(A, 29)[0] = 0x07a;
    mzd_row(A, 30)[0] = 0x1ac;
    return A;
  case 10:
    A = mzd_init(36, 10);
    mzd_row(A,  0)[0] = 0x200;
    mzd_row(A,  1)[0] = 0x001;
    mzd_row(A,  2)[0] = 0x3ff;
    mzd_row(A,  3)[0] = 0x36d;
    mzd_row(A,  4)[0] = 0x1b6;
    mzd_row(A,  5)[0] = 0x2db;
    mzd_row(A,  6)[0] = 0x0e9;
    mzd_row(A,  7)[0] = 0x13a;
    mzd_row(A,  8)[0] = 0x274;
    mzd_row(A,  9)[0] = 0x1d3;
    mzd_row(A, 10)[0] = 0x29d;
    mzd_row(A, 11)[0] = 0x34e;
    mzd_row(A, 12)[0] = 0x0b9;
    mzd_row(A, 13)[0] = 0x172;
    mzd_row(A, 14)[0] = 0x25c;
    mzd_row(A, 15)[0] = 0x1cb;
    mzd_row(A, 16)[0] = 0x2e5;
    mzd_row(A, 17)[0] = 0x32e;
    mzd_row(A, 18)[0] = 0x191;
    mzd_row(A, 19)[0] = 0x2b2;
    mzd_row(A, 20)[0] = 0x164;
    mzd_row(A, 21)[0] = 0x2c8;
    mzd_row(A, 22)[0] = 0x08f;
    mzd_row(A, 23)[0] = 0x323;
    mzd_row(A, 24)[0] = 0x0f5;
    mzd_row(A, 25)[0] = 0x07a;
    mzd_row(A, 26)[0] = 0x3ac;
    mzd_row(A, 27)[0] = 0x2f1;
    mzd_row(A, 28)[0] = 0x1e2;
    mzd_row(A, 29)[0] = 0x3c4;
    mzd_row(A, 30)[0] = 0x178;
    mzd_row(A, 31)[0] = 0x1af;
    mzd_row(A, 32)[0] = 0x313;
    mzd_row(A, 33)[0] = 0x135;
    mzd_row(A, 34)[0] = 0x09a;
    mzd_row(A, 35)[0] = 0x2bc;
    return A;
  case 11:
    A = mzd_init(40, 11);
    mzd_row(A,  0)[0] = 0x001;
    mzd_row(A,  1)[0] = 0x36d;
    mzd_row(A,  2)[0] = 0x5b6;
    mzd_row(A,  3)[0] = 0x6db;
    mzd_row(A,  4)[0] = 0x555;
    mzd_row(A,  5)[0] = 0x2aa;
    mzd_row(A,  6)[0] = 0x7ff;
    mzd_row(A,  7)[0] = 0x400;
    mzd_row(A,  8)[0] = 0x200;
    mzd_row(A,  9)[0] = 0x600;
    mzd_row(A, 10)[0] = 0x4e9;
    mzd_row(A, 11)[0] = 0x53a;
    mzd_row(A, 12)[0] = 0x274;
    mzd_row(A, 13)[0] = 0x1d3; 
    mzd_row(A, 14)[0] = 0x69d; 
    mzd_row(A, 15)[0] = 0x74e;
    mzd_row(A, 16)[0] = 0x4b9;
    mzd_row(A, 17)[0] = 0x172; 
    mzd_row(A, 18)[0] = 0x65c; 
    mzd_row(A, 19)[0] = 0x5cb;
    mzd_row(A, 20)[0] = 0x2e5;
    mzd_row(A, 21)[0] = 0x72e; 
    mzd_row(A, 22)[0] = 0x591; 
    mzd_row(A, 23)[0] = 0x6b2;
    mzd_row(A, 24)[0] = 0x564;
    mzd_row(A, 25)[0] = 0x2c8; 
    mzd_row(A, 26)[0] = 0x48f; 
    mzd_row(A, 27)[0] = 0x323;
    mzd_row(A, 28)[0] = 0x0f5;
    mzd_row(A, 29)[0] = 0x47a; 
    mzd_row(A, 30)[0] = 0x7ac; 
    mzd_row(A, 31)[0] = 0x2f1;
    mzd_row(A, 32)[0] = 0x5e2;
    mzd_row(A, 33)[0] = 0x3c4; 
    mzd_row(A, 34)[0] = 0x578; 
    mzd_row(A, 35)[0] = 0x1af;
    mzd_row(A, 36)[0] = 0x713;
    mzd_row(A, 37)[0] = 0x135; 
    mzd_row(A, 38)[0] = 0x09a; 
    mzd_row(A, 39)[0] = 0x6bc;
    return A;
  case 12:
    A = mzd_init(45, 12);
    mzd_row(A,  0)[0] = 0xb6d;
    mzd_row(A,  1)[0] = 0xdb6;
    mzd_row(A,  2)[0] = 0x6db;
    mzd_row(A,  3)[0] = 0x001;
    mzd_row(A,  4)[0] = 0x002;
    mzd_row(A,  5)[0] = 0x003;
    mzd_row(A,  6)[0] = 0x555;
    mzd_row(A,  7)[0] = 0xaaa;
    mzd_row(A,  8)[0] = 0xfff;
    mzd_row(A,  9)[0] = 0x800;
    mzd_row(A, 10)[0] = 0x400; 
    mzd_row(A, 11)[0] = 0xc00;
    mzd_row(A, 12)[0] = 0x4e9;
    mzd_row(A, 13)[0] = 0xd3a; 
    mzd_row(A, 14)[0] = 0xa74; 
    mzd_row(A, 15)[0] = 0x9d3;
    mzd_row(A, 16)[0] = 0xe9d;
    mzd_row(A, 17)[0] = 0x74e; 
    mzd_row(A, 18)[0] = 0x591; 
    mzd_row(A, 19)[0] = 0xeb2;
    mzd_row(A, 20)[0] = 0xd64;
    mzd_row(A, 21)[0] = 0xac8; 
    mzd_row(A, 22)[0] = 0xc8f; 
    mzd_row(A, 23)[0] = 0xb23;
    mzd_row(A, 24)[0] = 0x8f5;
    mzd_row(A, 25)[0] = 0x47a; 
    mzd_row(A, 26)[0] = 0x7ac; 
    mzd_row(A, 27)[0] = 0xaf1;
    mzd_row(A, 28)[0] = 0x5e2;
    mzd_row(A, 29)[0] = 0xbc4; 
    mzd_row(A, 30)[0] = 0xd78; 
    mzd_row(A, 31)[0] = 0x9af;
    mzd_row(A, 32)[0] = 0xf13;
    mzd_row(A, 33)[0] = 0x135; 
    mzd_row(A, 34)[0] = 0x89a; 
    mzd_row(A, 35)[0] = 0x6bc;
    mzd_row(A, 36)[0] = 0x631;
    mzd_row(A, 37)[0] = 0xa52; 
    mzd_row(A, 38)[0] = 0x294; 
    mzd_row(A, 39)[0] = 0x318;
    mzd_row(A, 40)[0] = 0xdef;
    mzd_row(A, 41)[0] = 0xc63; 
    mzd_row(A, 42)[0] = 0x4a5; 
    mzd_row(A, 43)[0] = 0x94a;
    mzd_row(A, 44)[0] = 0x18c;
    return A;
  case 13:
    A = mzd_init(49, 13);
    mzd_row(A,  0)[0] = 0x0001;
    mzd_row(A,  1)[0] = 0x1b6d;
    mzd_row(A,  2)[0] = 0x0db6;
    mzd_row(A,  3)[0] = 0x16db;
    mzd_row(A,  4)[0] = 0x1555;
    mzd_row(A,  5)[0] = 0x0aaa;
    mzd_row(A,  6)[0] = 0x1fff;
    mzd_row(A,  7)[0] = 0x1000;
    mzd_row(A,  8)[0] = 0x0800;
    mzd_row(A,  9)[0] = 0x1800;
    mzd_row(A, 10)[0] = 0x14e9;
    mzd_row(A, 11)[0] = 0x1d3a;
    mzd_row(A, 12)[0] = 0x1a74;
    mzd_row(A, 13)[0] = 0x09d3;
    mzd_row(A, 14)[0] = 0x0e9d;
    mzd_row(A, 15)[0] = 0x074e;
    mzd_row(A, 16)[0] = 0x1cb9;
    mzd_row(A, 17)[0] = 0x1972;
    mzd_row(A, 18)[0] = 0x0e5c;
    mzd_row(A, 19)[0] = 0x05cb;
    mzd_row(A, 20)[0] = 0x12e5;
    mzd_row(A, 21)[0] = 0x172e;
    mzd_row(A, 22)[0] = 0x1591;
    mzd_row(A, 23)[0] = 0x1eb2;
    mzd_row(A, 24)[0] = 0x1d64;
    mzd_row(A, 25)[0] = 0x1ac8;
    mzd_row(A, 26)[0] = 0x0c8f;
    mzd_row(A, 27)[0] = 0x0b23;
    mzd_row(A, 28)[0] = 0x08f5;
    mzd_row(A, 29)[0] = 0x047a;
    mzd_row(A, 30)[0] = 0x07ac;
    mzd_row(A, 31)[0] = 0x1af1;
    mzd_row(A, 32)[0] = 0x15e2;
    mzd_row(A, 33)[0] = 0x0bc4;
    mzd_row(A, 34)[0] = 0x0d78;
    mzd_row(A, 35)[0] = 0x09af;
    mzd_row(A, 36)[0] = 0x0f13;
    mzd_row(A, 37)[0] = 0x1135;
    mzd_row(A, 38)[0] = 0x189a;
    mzd_row(A, 39)[0] = 0x06bc;
    mzd_row(A, 40)[0] = 0x0631;
    mzd_row(A, 41)[0] = 0x0a52;
    mzd_row(A, 42)[0] = 0x1294;
    mzd_row(A, 43)[0] = 0x0318;
    mzd_row(A, 44)[0] = 0x1def;
    mzd_row(A, 45)[0] = 0x0c63;
    mzd_row(A, 46)[0] = 0x14a5;
    mzd_row(A, 47)[0] = 0x094a;
    mzd_row(A, 48)[0] = 0x118c;
    return A;
  case 14:
    A = mzd_init(55, 14);
    mzd_row(A,  0)[0] = 0x0001;
    mzd_row(A,  1)[0] = 0x1b6d;
    mzd_row(A,  2)[0] = 0x2db6;
    mzd_row(A,  3)[0] = 0x36db;
    mzd_row(A,  4)[0] = 0x1555;
    mzd_row(A,  5)[0] = 0x2aaa;
    mzd_row(A,  6)[0] = 0x3fff;
    mzd_row(A,  7)[0] = 0x34e9;
    mzd_row(A,  8)[0] = 0x1d3a;
    mzd_row(A,  9)[0] = 0x3a74;
    mzd_row(A, 10)[0] = 0x29d3;
    mzd_row(A, 11)[0] = 0x0e9d;
    mzd_row(A, 12)[0] = 0x274e;
    mzd_row(A, 13)[0] = 0x1cb9;
    mzd_row(A, 14)[0] = 0x3972;
    mzd_row(A, 15)[0] = 0x2e5c;
    mzd_row(A, 16)[0] = 0x25cb;
    mzd_row(A, 17)[0] = 0x32e5;
    mzd_row(A, 18)[0] = 0x172e;
    mzd_row(A, 19)[0] = 0x3591;
    mzd_row(A, 20)[0] = 0x1eb2;
    mzd_row(A, 21)[0] = 0x3d64;
    mzd_row(A, 22)[0] = 0x3ac8;
    mzd_row(A, 23)[0] = 0x2c8f;
    mzd_row(A, 24)[0] = 0x2b23;
    mzd_row(A, 25)[0] = 0x08f5;
    mzd_row(A, 26)[0] = 0x247a;
    mzd_row(A, 27)[0] = 0x07ac;
    mzd_row(A, 28)[0] = 0x1af1;
    mzd_row(A, 29)[0] = 0x35e2;
    mzd_row(A, 30)[0] = 0x2bc4;
    mzd_row(A, 31)[0] = 0x0d78;
    mzd_row(A, 32)[0] = 0x09af;
    mzd_row(A, 33)[0] = 0x2f13;
    mzd_row(A, 34)[0] = 0x3135;
    mzd_row(A, 35)[0] = 0x389a;
    mzd_row(A, 36)[0] = 0x26bc;
    mzd_row(A, 37)[0] = 0x0631;
    mzd_row(A, 38)[0] = 0x0a52;
    mzd_row(A, 39)[0] = 0x1294;
    mzd_row(A, 40)[0] = 0x2318;
    mzd_row(A, 41)[0] = 0x3def;
    mzd_row(A, 42)[0] = 0x0c63;
    mzd_row(A, 43)[0] = 0x14a5;
    mzd_row(A, 44)[0] = 0x294a;
    mzd_row(A, 45)[0] = 0x318c;
    mzd_row(A, 46)[0] = 0x2000;
    mzd_row(A, 47)[0] = 0x1000;
    mzd_row(A, 48)[0] = 0x0800;
    mzd_row(A, 49)[0] = 0x0400;
    mzd_row(A, 50)[0] = 0x3c00;
    mzd_row(A, 51)[0] = 0x3000;
    mzd_row(A, 52)[0] = 0x2800;
    mzd_row(A, 53)[0] = 0x1400;
    mzd_row(A, 54)[0] = 0x0c00;
    return A;
  case 15:
    A = mzd_init(60, 15);
    mzd_row(A,  0)[0] = 0x0001;
    mzd_row(A,  1)[0] = 0x7fff;
    mzd_row(A,  2)[0] = 0x5b6d;
    mzd_row(A,  3)[0] = 0x6db6;
    mzd_row(A,  4)[0] = 0x36db;
    mzd_row(A,  5)[0] = 0x4000;
    mzd_row(A,  6)[0] = 0x2000;
    mzd_row(A,  7)[0] = 0x6000;
    mzd_row(A,  8)[0] = 0x74e9;
    mzd_row(A,  9)[0] = 0x1d3a;
    mzd_row(A, 10)[0] = 0x3a74;
    mzd_row(A, 11)[0] = 0x69d3;
    mzd_row(A, 12)[0] = 0x4e9d;
    mzd_row(A, 13)[0] = 0x274e;
    mzd_row(A, 14)[0] = 0x5cb9;
    mzd_row(A, 15)[0] = 0x3972;
    mzd_row(A, 16)[0] = 0x2e5c;
    mzd_row(A, 17)[0] = 0x65cb;
    mzd_row(A, 18)[0] = 0x72e5;
    mzd_row(A, 19)[0] = 0x172e;
    mzd_row(A, 20)[0] = 0x7591;
    mzd_row(A, 21)[0] = 0x1eb2;
    mzd_row(A, 22)[0] = 0x3d64;
    mzd_row(A, 23)[0] = 0x7ac8;
    mzd_row(A, 24)[0] = 0x2c8f;
    mzd_row(A, 25)[0] = 0x6b23;
    mzd_row(A, 26)[0] = 0x48f5;
    mzd_row(A, 27)[0] = 0x647a;
    mzd_row(A, 28)[0] = 0x47ac;
    mzd_row(A, 29)[0] = 0x1af1;
    mzd_row(A, 30)[0] = 0x35e2;
    mzd_row(A, 31)[0] = 0x6bc4;
    mzd_row(A, 32)[0] = 0x4d78;
    mzd_row(A, 33)[0] = 0x09af;
    mzd_row(A, 34)[0] = 0x2f13;
    mzd_row(A, 35)[0] = 0x7135;
    mzd_row(A, 36)[0] = 0x789a;
    mzd_row(A, 37)[0] = 0x26bc;
    mzd_row(A, 38)[0] = 0x4631;
    mzd_row(A, 39)[0] = 0x4a52;
    mzd_row(A, 40)[0] = 0x5294;
    mzd_row(A, 41)[0] = 0x6318;
    mzd_row(A, 42)[0] = 0x3def;
    mzd_row(A, 43)[0] = 0x0c63;
    mzd_row(A, 44)[0] = 0x14a5;
    mzd_row(A, 45)[0] = 0x294a;
    mzd_row(A, 46)[0] = 0x318c;
    mzd_row(A, 47)[0] = 0x4d21;
    mzd_row(A, 48)[0] = 0x1a42;
    mzd_row(A, 49)[0] = 0x7348;
    mzd_row(A, 50)[0] = 0x6690;
    mzd_row(A, 51)[0] = 0x2bb1;
    mzd_row(A, 52)[0] = 0x5763;
    mzd_row(A, 53)[0] = 0x15d8;
    mzd_row(A, 54)[0] = 0x0576;
    mzd_row(A, 55)[0] = 0x47cd;
    mzd_row(A, 56)[0] = 0x42bb;
    mzd_row(A, 57)[0] = 0x4857;
    mzd_row(A, 58)[0] = 0x215d;
    mzd_row(A, 59)[0] = 0x3b1f;
    return A;
  case 16:
    A = mzd_init(64, 16);
    mzd_row(A,  0)[0] = 0xdb6d;
    mzd_row(A,  1)[0] = 0x6db6;
    mzd_row(A,  2)[0] = 0xb6db;
    mzd_row(A,  3)[0] = 0x0001;
    mzd_row(A,  4)[0] = 0x0002;
    mzd_row(A,  5)[0] = 0x0003;
    mzd_row(A,  6)[0] = 0x5555;
    mzd_row(A,  7)[0] = 0xaaaa;
    mzd_row(A,  8)[0] = 0xffff;
    mzd_row(A,  9)[0] = 0x8000;
    mzd_row(A, 10)[0] = 0x4000;
    mzd_row(A, 11)[0] = 0xc000;
    mzd_row(A, 12)[0] = 0x74e9;
    mzd_row(A, 13)[0] = 0x9d3a;
    mzd_row(A, 14)[0] = 0x3a74;
    mzd_row(A, 15)[0] = 0xe9d3;
    mzd_row(A, 16)[0] = 0x4e9d;
    mzd_row(A, 17)[0] = 0xa74e;
    mzd_row(A, 18)[0] = 0x5cb9;
    mzd_row(A, 19)[0] = 0xb972;
    mzd_row(A, 20)[0] = 0x2e5c;
    mzd_row(A, 21)[0] = 0xe5cb;
    mzd_row(A, 22)[0] = 0x72e5;
    mzd_row(A, 23)[0] = 0x972e;
    mzd_row(A, 24)[0] = 0xf591;
    mzd_row(A, 25)[0] = 0x1eb2;
    mzd_row(A, 26)[0] = 0x3d64;
    mzd_row(A, 27)[0] = 0x7ac8;
    mzd_row(A, 28)[0] = 0xac8f;
    mzd_row(A, 29)[0] = 0xeb23;
    mzd_row(A, 30)[0] = 0xc8f5;
    mzd_row(A, 31)[0] = 0x647a;
    mzd_row(A, 32)[0] = 0x47ac;
    mzd_row(A, 33)[0] = 0x9af1;
    mzd_row(A, 34)[0] = 0x35e2;
    mzd_row(A, 35)[0] = 0x6bc4;
    mzd_row(A, 36)[0] = 0x4d78;
    mzd_row(A, 37)[0] = 0x89af;
    mzd_row(A, 38)[0] = 0xaf13;
    mzd_row(A, 39)[0] = 0xf135;
    mzd_row(A, 40)[0] = 0x789a;
    mzd_row(A, 41)[0] = 0x26bc;
    mzd_row(A, 42)[0] = 0xc631;
    mzd_row(A, 43)[0] = 0x4a52;
    mzd_row(A, 44)[0] = 0x5294;
    mzd_row(A, 45)[0] = 0x6318;
    mzd_row(A, 46)[0] = 0xbdef;
    mzd_row(A, 47)[0] = 0x8c63;
    mzd_row(A, 48)[0] = 0x94a5;
    mzd_row(A, 49)[0] = 0x294a;
    mzd_row(A, 50)[0] = 0x318c;
    mzd_row(A, 51)[0] = 0xcd21;
    mzd_row(A, 52)[0] = 0x9a42;
    mzd_row(A, 53)[0] = 0xf348;
    mzd_row(A, 54)[0] = 0xe690;
    mzd_row(A, 55)[0] = 0x2bb1;
    mzd_row(A, 56)[0] = 0x5763;
    mzd_row(A, 57)[0] = 0x15d8;
    mzd_row(A, 58)[0] = 0x8576;
    mzd_row(A, 59)[0] = 0xc7cd;
    mzd_row(A, 60)[0] = 0x42bb;
    mzd_row(A, 61)[0] = 0x4857;
    mzd_row(A, 62)[0] = 0x215d;
    mzd_row(A, 63)[0] = 0xbb1f;
    return A;

  default:
    m4ri_die("only degrees up to 16 are implemented but got degree %d\n", degree);
  }
  return NULL;
}

/**
 * \param length The length of the polynomial we want to reduce
 * \param poly A polynomial
 * \param d The degree of poly
 */

mzd_t *_crt_modred_mat(const deg_t length, const word poly, const deg_t d) {
  mzd_t *A = mzd_init(d, length);

  /* (x-infinity)^d */
  if (poly == 0) {
    for(deg_t i=0; i<d; i++) 
      mzd_write_bit(A, i, length-i-1, 1);
    return A;
  }

  mzd_t *f = mzd_init(1, length);
  mzd_t *t = mzd_init(1, length);

  for(deg_t i=0; i<length; i++) {
    mzd_set_ui(f, 0);
    mzd_row(f, 0)[i/m4ri_radix] = __M4RI_TWOPOW(i%m4ri_radix);
    word ii = i;
    while(ii >= (word)d) {
      /* f ^= gf2x_mul((1ULL<<(ii-d)), poly, length); */
      mzd_set_ui(t, 0);
      mzd_xor_bits(t, 0, ii-d, d+1, poly);

      mzd_add(f, f, t);

      /* ii = gf2x_deg(f); */
      ii = 0;
      for(wi_t j=f->width-1; j>=0; j--) {
        if (mzd_row(f, 0)[j]) {
          ii = gf2x_deg(mzd_row(f, 0)[j]) + m4ri_radix*j;
          break;
        }
      }
    }
    for(deg_t j=0; j<= (deg_t)ii; j++)
      mzd_write_bit(A, j, i, (mzd_row(f, 0)[j/m4ri_radix]>>(j%m4ri_radix)) & 0x1);
  }
  return A;
}

blm_t *_blm_finish_polymult(const gf2e *ff, blm_t *f) {
  assert( (f != NULL) & (f->F != NULL) & (f->G != NULL) );

  const rci_t m = f->F->nrows;
  const rci_t c_nrows = f->F->ncols + f->G->ncols - 1;

  mzd_t *H = mzd_init(c_nrows, m);

  mzd_t *F_T = mzd_transpose(NULL, f->F);
  mzd_t *G_T = mzd_transpose(NULL, f->G);

  mzd_t *C = mzd_init(m, m);

  /* we find a full rank matrix */

  word v = 0;
  word w = 0;
  rci_t r = 0;
  rci_t rank = 0;

  mzd_t *pivots = mzd_init(m, m4ri_radix*2);

  mzp_t *P = mzp_init(m);
  mzp_t *Q = mzp_init(m);

  while(rank < m) {
    /* x^v = (0, 0, ... , 1, ..., 0, 0) * F_T -> select row v */
    for(wi_t j=0; j< C->width; j++)
      mzd_row(C, r)[j] = mzd_row(F_T, v)[j] & mzd_row(G_T, w)[j];

    mzd_row(pivots, r)[0] = v;
    mzd_row(pivots, r)[1] = w;

    w++;
    if (w == (word)f->G->ncols) {
      v++;
      if (v == (word)f->F->ncols)
        v = 0;
      w = v;
    }
    r++;
    if (r == C->nrows) {
      mzd_t *D = mzd_copy(NULL, C);
      rank = mzd_ple(D, P, Q, 0);
      mzd_apply_p_left(pivots, P);
      mzd_apply_p_left(C, P);
      mzd_free(D);

      if (rank < m)
        r = rank;
    }
  }
  mzp_free(P);
  mzp_free(Q);

  for(r=0; r<m; r++) {
    /* x^v = (0, 0, ... , 1, ..., 0, 0) * F_T -> select row v */
    v = mzd_row(pivots, r)[0];
    w = mzd_row(pivots, r)[1];
    for(wi_t j=0; j< C->width; j++)
      mzd_row(C, r)[j] = mzd_row(F_T, v)[j] & mzd_row(G_T, w)[j];
  }
  mzd_free(F_T);
  mzd_free(G_T);

  // This should be replaced by TRSM calls
  mzd_t *D = mzd_inv_m4ri(NULL, C, 0); 
  mzd_free(C);
  mzd_t *DT = mzd_transpose(NULL, D);
  mzd_free(D);

  mzd_t *a = mzd_init(1, m);
  mzd_t *b = mzd_init(1, H->ncols);

  for(rci_t i=0; i<H->nrows; i++) {
    mzd_set_ui(a, 0);
    for(rci_t j=0; j<m; j++) {
      v = mzd_row(pivots, j)[0];
      w = mzd_row(pivots, j)[1];
      if ((v+w) == (word)i)
        mzd_write_bit(a, 0, j, 1);
    }
    mzd_mul(b, a, DT, 0);

    /* copy result to H */
    for(rci_t j=0; j<H->ncols; j++)
      mzd_write_bit(H, i, j, mzd_read_bit(b, 0, j ));
  }
  mzd_free(a);
  mzd_free(b);
  mzd_free(pivots);

  if (ff == NULL) {
    f->H = H;
  } else { 
    mzd_t *N = _crt_modred_mat(H->nrows, ff->minpoly, ff->degree);
    f->H = mzd_mul(NULL, N, H, 0);
    mzd_free(N);
    mzd_free(H);
  }
  return f;
}

const int costs[17] = {0, 1, 3, 6, 9, 13, 17, 22, 27, 31, 36, 40, 45, 49, 55, 60, 64};
//const int costs[17] = {0, 1, 3, 6, 9, 13, 15, 22, 24, 30, 33, 39, 42, 48, 51, 54, 60}; /* best possible */

blm_t *blm_init_crt(const gf2e *ff, const deg_t f_ncols, const deg_t g_ncols, const int *p, int djb) {
  blm_t *f = m4ri_mm_malloc(sizeof(blm_t));

  // iterator over irreducible polynomials
  int *p_it = (int*)m4ri_mm_calloc(sizeof(int), M4RIE_CRT_LEN); 

  mzd_t *M, *T;

  word poly = 0;

  rci_t m = costs[p[0]];
  for(int d=1; d<M4RIE_CRT_LEN; d++)
    m += costs[d] * p[d];

  f->F = mzd_init(m, f_ncols);
  f->f = NULL;
  f->G = mzd_init(m, g_ncols);
  f->g = NULL;

  rci_t r = 0;


  /**
   * 1) We construct maps F,G which combine modular reduction to degree d and the linear map
   *    required for multiplying modulo a polynomial of degree d.
   */

  /**
   * 1.1) We deal with (x+infinity)^omega first
   */

  if(p[0] != 0) {
    deg_t d = p[0];
    mzd_t *N  = _small_multiplication_map(d);
    M = _crt_modred_mat(f_ncols, poly, d);
    T = mzd_init_window(f->F, r, 0, r + costs[d], f_ncols); 
    mzd_mul(T, N, M, 0);
    mzd_free(T);
    mzd_free(M);

    M = _crt_modred_mat(g_ncols, poly, d);
    T = mzd_init_window(f->G, r, 0, r + costs[d], g_ncols); 
    mzd_mul(T, N, M, 0);
    mzd_free(T);
    mzd_free(M);
    
    mzd_free(N);

    r += costs[d];
  }

  /**
   * 1.2) We deal with regular polynomial which are co-prime
   */

  for(deg_t d=1; d<M4RIE_CRT_LEN; d++) {
    if (p[d] == 0)
      continue;

    mzd_t *N  = _small_multiplication_map(d);

    for(int i=0; i<p[d]; i++) {
      if ((word)p_it[d] < irreducible_polynomials[d][0]) {
        poly = irreducible_polynomials[d][ 1 + p_it[d] ];
        p_it[d]++;
      } else if (d/2 && (word)p_it[d/2] < irreducible_polynomials[d/2][0]) {
        /** the minimal polynomial is a square */
        poly = irreducible_polynomials[d/2][ 1 + p_it[d/2]];
        p_it[d/2]++;
        poly = gf2x_mul(poly, poly, d/2+1);
      } else if (d/4 && (word)p_it[d/4] < irreducible_polynomials[d/4][0]) {
        /** the minimal polynomial is a fourth power */
        poly = irreducible_polynomials[d/4][1 + p_it[d/4]];
        p_it[d/4]++;
        poly = gf2x_mul(poly, poly, d/4+1);
        poly = gf2x_mul(poly, poly, d/2+1);
      } else if (d/8 && (word)p_it[d/8] < irreducible_polynomials[d/8][0]) {
        /** the minimal polynomial is an eigth power */
        poly = irreducible_polynomials[d/8][p_it[d/8]+ 1];
        p_it[d/8]++;
        poly = gf2x_mul(poly, poly, d/8+1);
        poly = gf2x_mul(poly, poly, d/4+1);
        poly = gf2x_mul(poly, poly, d/2+1);
      } else {
        m4ri_die("Degree %d is not implemented\n", d);
      }
      M = _crt_modred_mat(f_ncols, poly, d);
      T = mzd_init_window(f->F, r, 0, r + costs[d], f_ncols); 
      mzd_mul(T, N, M, 0);
      mzd_free(T);
      mzd_free(M);

      M = _crt_modred_mat(g_ncols, poly, d);
      T = mzd_init_window(f->G, r, 0, r + costs[d], g_ncols); 
      mzd_mul(T, N, M, 0);
      mzd_free(T);
      mzd_free(M);

      r += costs[d];
    }

    mzd_free(N);
  }

  m4ri_mm_free(p_it);

  /**
   * 2) We solve for H as we know poly(c) and (F*vec(a) x G*vec(b)). We pick points poly(a) = x^v,
   * poly(b) = x^w (hence: poly(c) = x^(v+w)).
   */
  _blm_finish_polymult(ff, f);
  f->h = NULL;

  /**
   * 3) We compile DJB maps if asked for.
   */

  if (djb)
    _blm_djb_compile(f);

  return f;
}


void blm_free(blm_t *f) {
  mzd_free(f->F);
  mzd_free(f->G);
  mzd_free(f->H);
  if (f->f != f->g)
    djb_free(f->g);
  djb_free(f->f);
  djb_free(f->h);

  m4ri_mm_free(f);
}

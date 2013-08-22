#include "m4rie/blm.h"
#include "mzd_ptr.h"

void _mzd_ptr_apply_blm(const gf2e *ff, mzd_t **X, const mzd_t **A, const mzd_t **B, const blm_t *f) {
  assert((f->H!=NULL) & (f->F!=NULL) & (f->G!=NULL) &   \
         (f->H->ncols == f->F->nrows) &                 \
         (f->F->nrows == f->G->nrows));

  mzd_t *t0 = mzd_init(A[0]->nrows, B[0]->ncols);
  mzd_t *t1 = mzd_init(A[0]->nrows, A[0]->ncols);
  mzd_t *t2 = mzd_init(B[0]->nrows, B[0]->ncols);

  for(rci_t i=0; i<f->F->nrows; i++) {
    for(rci_t j=0; j<f->F->ncols; j++)
      if(mzd_read_bit(f->F, i, j)) 
        mzd_add(t1, t0, A[j]);

    for(rci_t j=0; j<f->G->ncols; j++)
      if(mzd_read_bit(f->G, i, j)) 
        mzd_add(t2, t0, B[j]);

    mzd_mul(t0, t1, t2, 0);

    for(rci_t j=0; j<-f->H->nrows; j++)
      if(mzd_read_bit(f->H, j, i))
        _mzd_ptr_add_modred(ff, t0, X, j);
  }

  mzd_free(t0);
  mzd_free(t1);
  mzd_free(t2);
}

mzd_t *_blm_smallmul_F(const deg_t degree) {
  mzd_t *A;
  switch(degree) {
  case 1:
    A = mzd_init(1, 1);
    A->rows[0][0] = 0x1;
    return A;
  case 2:
    A = mzd_init(3, 2);
    A->rows[0][0] = 0x1; A->rows[1][0] = 0x3; A->rows[2][0] = 0x2;
    return A;
  case 3:
    A = mzd_init(6, 3);
    A->rows[0][0] = 0x1; A->rows[1][0] = 0x2; A->rows[2][0] = 0x4;
    A->rows[3][0] = 0x3; A->rows[4][0] = 0x5; A->rows[5][0] = 0x6;
    return A;
  case 4:
    A = mzd_init(9, 4);
    A->rows[0][0] = 0x1; A->rows[1][0] = 0x2; A->rows[2][0] = 0x4; A->rows[3][0] = 0x8; 
    A->rows[4][0] = 0xf; A->rows[5][0] = 0x3; A->rows[6][0] = 0x5; A->rows[7][0] = 0xa; A->rows[8][0] = 0xc;
    return A;
  case 5:
    A = mzd_init(13, 5);
    A->rows[ 0][0] = 0x01; A->rows[ 1][0] = 0x02; A->rows[ 2][0] = 0x08; A->rows[ 3][0] = 0x10; 
    A->rows[ 4][0] = 0x11; A->rows[ 5][0] = 0x03; A->rows[ 6][0] = 0x18; A->rows[ 7][0] = 0x16;
    A->rows[ 8][0] = 0x0d; A->rows[ 9][0] = 0x1b; A->rows[10][0] = 0x17; A->rows[11][0] = 0x1d;
    A->rows[12][0] = 0x1f;
    return A;
  case 6:
    A = mzd_init(17,  6);
    A->rows[ 0][0] = 0x01;     A->rows[ 1][0] = 0x02;     A->rows[ 2][0] = 0x10;      A->rows[ 3][0] = 0x20;
    A->rows[ 4][0] = 0x30;     A->rows[ 5][0] = 0x03;     A->rows[ 6][0] = 0x18;      A->rows[ 7][0] = 0x06;
    A->rows[ 8][0] = 0x12;     A->rows[ 9][0] = 0x0c;     A->rows[10][0] = 0x38;      A->rows[11][0] = 0x07;
    A->rows[12][0] = 0x29;     A->rows[13][0] = 0x25;     A->rows[14][0] = 0x2d;      A->rows[15][0] = 0x1b;
    A->rows[16][0] = 0x3f;
    return A;
  case 7:
    A = mzd_init(22,  7);
    A->rows[ 0][0] = 0x7f;     A->rows[ 1][0] = 0x6e;     A->rows[ 2][0] = 0x3b;     A->rows[ 3][0] = 0x5d;
    A->rows[ 4][0] = 0x6d;     A->rows[ 5][0] = 0x5b;     A->rows[ 6][0] = 0x36;     A->rows[ 7][0] = 0x03;
    A->rows[ 8][0] = 0x05;     A->rows[ 9][0] = 0x11;     A->rows[10][0] = 0x0a;     A->rows[11][0] = 0x44;
    A->rows[12][0] = 0x28;     A->rows[13][0] = 0x50;     A->rows[14][0] = 0x60;     A->rows[15][0] = 0x01;
    A->rows[16][0] = 0x02;     A->rows[17][0] = 0x04;     A->rows[18][0] = 0x08;     A->rows[19][0] = 0x10;
    A->rows[20][0] = 0x20;     A->rows[21][0] = 0x40;     
    return A;
  case 8:
    A = mzd_init(27,  8);
    A->rows[ 0][0] = 0x01;     A->rows[ 1][0] = 0x02;     A->rows[ 2][0] = 0x04;     A->rows[ 3][0] = 0x08;
    A->rows[ 4][0] = 0x10;     A->rows[ 5][0] = 0x20;     A->rows[ 6][0] = 0x40;     A->rows[ 7][0] = 0x80;
    A->rows[ 8][0] = 0x05;     A->rows[ 9][0] = 0x0c;     A->rows[10][0] = 0x44;     A->rows[11][0] = 0x30;
    A->rows[12][0] = 0xc0;     A->rows[13][0] = 0x0f;     A->rows[14][0] = 0xa0;     A->rows[15][0] = 0x11;
    A->rows[16][0] = 0x0a;     A->rows[17][0] = 0xaa;     A->rows[18][0] = 0x03;     A->rows[19][0] = 0x50;
    A->rows[20][0] = 0x22;     A->rows[21][0] = 0xff;     A->rows[22][0] = 0x55;     A->rows[23][0] = 0x33;
    A->rows[24][0] = 0xcc;     A->rows[25][0] = 0x88;     A->rows[26][0] = 0xf0;
    return A;
  case 9:
    A = mzd_init(31,  9);
    A->rows[ 0][0] = 0x100;     A->rows[ 1][0] = 0x001;     A->rows[ 2][0] = 0x002;     A->rows[ 3][0] = 0x003;
    A->rows[ 4][0] = 0x155;     A->rows[ 5][0] = 0x0aa;     A->rows[ 6][0] = 0x1ff;     A->rows[ 7][0] = 0x16d;
    A->rows[ 8][0] = 0x1b6;     A->rows[ 9][0] = 0x0db;     A->rows[10][0] = 0x0e9;     A->rows[11][0] = 0x13a;
    A->rows[12][0] = 0x074;     A->rows[13][0] = 0x1d3;     A->rows[14][0] = 0x09d;     A->rows[15][0] = 0x14e;
    A->rows[16][0] = 0x0b9;     A->rows[17][0] = 0x172;     A->rows[18][0] = 0x05c;     A->rows[19][0] = 0x1cb;
    A->rows[20][0] = 0x0e5;     A->rows[21][0] = 0x12e;     A->rows[22][0] = 0x191;     A->rows[23][0] = 0x0b2;
    A->rows[24][0] = 0x164;     A->rows[25][0] = 0x0c8;     A->rows[26][0] = 0x08f;     A->rows[27][0] = 0x123;
    A->rows[28][0] = 0x0f5;     A->rows[29][0] = 0x07a;     A->rows[30][0] = 0x1ac;
    return A;
  case 10:
    A = mzd_init(36, 10);
    A->rows[ 0][0] = 0x200;     A->rows[ 1][0] = 0x001;     A->rows[ 2][0] = 0x3ff;     A->rows[ 3][0] = 0x36d;
    A->rows[ 4][0] = 0x1b6;     A->rows[ 5][0] = 0x2db;     A->rows[ 6][0] = 0x0e9;     A->rows[ 7][0] = 0x13a;
    A->rows[ 8][0] = 0x274;     A->rows[ 9][0] = 0x1d3;     A->rows[10][0] = 0x29d;     A->rows[11][0] = 0x34e;
    A->rows[12][0] = 0x0b9;     A->rows[13][0] = 0x172;     A->rows[14][0] = 0x25c;     A->rows[15][0] = 0x1cb;
    A->rows[16][0] = 0x2e5;     A->rows[17][0] = 0x32e;     A->rows[18][0] = 0x191;     A->rows[19][0] = 0x2b2;
    A->rows[20][0] = 0x164;     A->rows[21][0] = 0x2c8;     A->rows[22][0] = 0x08f;     A->rows[23][0] = 0x323;
    A->rows[24][0] = 0x0f5;     A->rows[25][0] = 0x07a;     A->rows[26][0] = 0x3ac;     A->rows[27][0] = 0x2f1;
    A->rows[28][0] = 0x1e2;     A->rows[29][0] = 0x3c4;     A->rows[30][0] = 0x178;     A->rows[31][0] = 0x1af;
    A->rows[32][0] = 0x313;     A->rows[33][0] = 0x135;     A->rows[34][0] = 0x09a;     A->rows[35][0] = 0x2bc;
    return A;
  case 11:
  case 12:
  case 13:
  case 14:
  case 15:
  case 16:
  default:
    m4ri_die("not implemented\n");
  }
  return NULL;
}

/**
 * \param length The length of the polynomial we want to reduce
 * \param mp A polynomial
 * \param d The degree of mp
 */

mzd_t *_blm_modred_mat(const deg_t length, const word mp, const deg_t d) {
  mzd_t *A = mzd_init(d, length);

  /* infinity */
  if (mp == 0) {
    for(deg_t i=0; i<d; i++) 
      mzd_write_bit(A, i, length-i-1, 1);
    return A;
  }

  mzd_t *f = mzd_init(1, length);
  mzd_t *t = mzd_init(1, length);

  for(deg_t i=0; i<length; i++) {
    mzd_set_ui(f, 0);
    f->rows[0][i/m4ri_radix] = __M4RI_TWOPOW(i%m4ri_radix);
    word ii = i;
    while(ii >= d) {
      /* f ^= gf2x_mul((1ULL<<(ii-d)), mp, length); */
      mzd_set_ui(t, 0);
      mzd_xor_bits(t, 0, ii-d, d+1, mp);

      mzd_add(f, f, t);

      /* ii = gf2x_deg(f); */
      ii = 0;
      for(wi_t j=f->width-1; j>=0; j--) {
        if (f->rows[0][j]) {
          ii = gf2x_deg(f->rows[0][j]) + m4ri_radix*j;
          break;
        }
      }
    }
    for(deg_t j=0; j<= ii; j++) 
      mzd_write_bit(A, j, i, (f->rows[0][j/m4ri_radix]>>(j%m4ri_radix)) & 0x1);
  }
  return A;
}

blm_t *_blm_finish_polymult(blm_t *f) {
  assert( (f != NULL) & (f->F != NULL) & (f->G != NULL) );

  const rci_t m = f->F->nrows;
  const rci_t c_nrows = f->F->ncols + f->G->ncols - 1;

  f->H = mzd_init(c_nrows, m);

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
      C->rows[r][j] = F_T->rows[v][j] & G_T->rows[w][j];

    pivots->rows[r][0] = v;
    pivots->rows[r][1] = w;

    w++;
    if (w == f->G->ncols) {
      v++;
      if (v == f->F->ncols)
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
  mzd_free(F_T);
  mzd_free(G_T);
  mzp_free(P);
  mzp_free(Q);

  for(r=0; r<m; r++) {
    /* x^v = (0, 0, ... , 1, ..., 0, 0) * F_T -> select row v */
    v = pivots->rows[r][0];
    w = pivots->rows[r][1];
    for(wi_t j=0; j< C->width; j++)
      C->rows[r][j] = F_T->rows[v][j] & G_T->rows[w][j];
  }

  mzd_t *D = mzd_inv_m4ri(NULL, C, 0); // This should be replaced by TRSM calls
  mzd_free(C);
  mzd_t *DT = mzd_transpose(NULL, D);
  mzd_free(D);

  mzd_t *a = mzd_init(1, m);
  mzd_t *b = mzd_init(1, f->H->ncols);

  for(rci_t i=0; i<f->H->nrows; i++) {
    mzd_set_ui(a, 0);
    for(rci_t j=0; j<m; j++) {
      v = pivots->rows[j][0];
      w = pivots->rows[j][1];
      if ((v+w) == i)
        mzd_write_bit(a, 0, j, 1);
    }
    mzd_mul(b, a, DT, 0);

    /* copy result to H */
    for(rci_t j=0; j<f->H->ncols; j++)
      mzd_write_bit(f->H, i, j, mzd_read_bit(b, 0, j ));
  }
  mzd_free(a);
  mzd_free(b);
  mzd_free(pivots);

  return f;
}

const int costs[17] = {0, 1, 3, 6, 9, 13, 17, 22, 27, 31, 36, 40, 45, 49, 55, 60, 64};

blm_t *blm_init_multimod(const deg_t f_ncols, const deg_t g_ncols, const deg_t deg, const int *primes) {
  blm_t *f = m4ri_mm_malloc(sizeof(blm_t));

  // iterator over co-primes
  int *primes_it = (int*)m4ri_mm_calloc(sizeof(int), deg); 

  mzd_t *M, *T;

  word prime = 0;
  int infinity_done = 0;

  rci_t m = 0;
  for(int d=0; d<deg; d++)
    m += costs[d] * primes[d];

  f->F = mzd_init(m, f_ncols);
  f->G = mzd_init(m, g_ncols);

  rci_t r = 0;

  assert(primes[0] == 0);

  /**
   * 1) We construct matrices F,G which combine modular reduction and the linear map required for
   * multiplying modulo a polynomial of degree d.
   */

  for(deg_t d=1; d<deg; d++) {
    mzd_t *N  = _blm_smallmul_F(d);

    for(int i=0; i<primes[d]; i++) {
      if (primes_it[d] < irreducible_polynomials[d][0]) {
        prime = irreducible_polynomials[d][ primes_it[d]+ 1 ];
        primes_it[d]++;
      } else if (d/2 && primes_it[d/2] < irreducible_polynomials[d/2][0]) {
        /** the minimal polynomial is a square */
        prime = irreducible_polynomials[d/2][primes_it[d/2]+ 1];
        primes_it[d/2]++;
        prime = gf2x_mul(prime, prime, d/2+1);
      } else if (d/4 && primes_it[d/4] < irreducible_polynomials[d/4][0]) {
        /** the minimal polynomial is a fourth power */
        prime = irreducible_polynomials[d/4][primes_it[d/4]+ 1];
        primes_it[d/4]++;
        prime = gf2x_mul(prime, prime, d/4+1);
        prime = gf2x_mul(prime, prime, d/2+1);
      } else if (!infinity_done){
        /** we assume we want to evaluate at infinity **/
        if ((d == 1) | (d == 2) | (d == 4) ) 
          prime = 0;
         else 
          m4ri_die("not imlemented\n");
        infinity_done = 1;
      } else {
        m4ri_die("not imlemented\n");
      }
      M = _blm_modred_mat(f_ncols, prime, d);
      T = mzd_init_window(f->F, r, 0, r + costs[d], f_ncols); 
      mzd_mul(T, N, M, 0);
      mzd_free(T);
      mzd_free(M);

      M = _blm_modred_mat(g_ncols, prime, d);
      T = mzd_init_window(f->G, r, 0, r + costs[d], g_ncols); 
      mzd_mul(T, N, M, 0);
      mzd_free(T);
      mzd_free(M);

      r += costs[d];
    }
    mzd_free(N);
  }

  m4ri_mm_free(primes_it);

  /**
   * 2) We solve for H as we know poly(c) and (F*vec(a) x G*vec(b)). We pick points poly(a) = x^v,
   * poly(b) = x^w (hence: poly(c) = x^(v+w)).
   */
  return _blm_finish_polymult(f);
}


void blm_free(blm_t *f) {
  mzd_free(f->F);
  mzd_free(f->G);
  mzd_free(f->H);
  m4ri_mm_free(f);
}

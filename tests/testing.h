#include <unistd.h>
#include <m4rie.h>

#define m4rie_check(expr)						\
  if (!expr) {								\
    fail_ret += 1;                                                      \
    printf("\n%s in %s:%d failed\n", #expr, __FILE__, __LINE__);  \
  } 

const word m4rie_canary = (word)0x63616E61727921ULL;

static inline word field_mask(const gf2e *ff) {
  const word mask_fld = ((1<<ff->degree)-1);
  word mask = 0;
  switch(gf2e_degree_to_w(ff)) {
  case 2:
    break;
  case 4:
    mask |= mask_fld<< 0 | mask_fld<< 4 | mask_fld<< 8 | mask_fld<<12;
    mask |= mask_fld<<16 | mask_fld<<20 | mask_fld<<24 | mask_fld<<28;
    mask |= mask_fld<<32 | mask_fld<<36 | mask_fld<<40 | mask_fld<<44;
    mask |= mask_fld<<48 | mask_fld<<52 | mask_fld<<56 | mask_fld<<60;
    break;
  case 8:
    mask |= mask_fld<< 0 | mask_fld<< 8 | mask_fld<<16 | mask_fld<<24;
    mask |= mask_fld<<32 | mask_fld<<40 | mask_fld<<48 | mask_fld<<56;
    break;
  case 16:
    mask |= mask_fld<< 0 | mask_fld<<16 | mask_fld<<32 | mask_fld<<48;
    break;
  }
  return mask;
}

static inline void mzed_set_canary(mzed_t *A) {
  const word mask_end   = __M4RI_LEFT_BITMASK(A->x->ncols % m4ri_radix);
  const word mask_field = field_mask(A->finite_field);
  const rci_t n = A->x->width-1;

  for(rci_t i=0; i<A->nrows; i++) {
    A->x->rows[i][n] = (A->x->rows[i][n] & mask_end)   | (m4rie_canary & mask_field & ~mask_end);
  }
}

static inline void mzed_clear_canary(mzed_t *A) {
  const word mask_end   = __M4RI_LEFT_BITMASK(A->x->ncols % m4ri_radix);
  const rci_t n = A->x->width-1;

  for(rci_t i=0; i<A->nrows; i++) {
    A->x->rows[i][n] &= mask_end;
  }
}

static inline int mzed_canary_is_alive(mzed_t *A) {
  const word mask_end   = __M4RI_LEFT_BITMASK(A->x->ncols % m4ri_radix);
  const word mask_field = field_mask(A->finite_field);
  const rci_t n = A->x->width-1;

  for(rci_t i=0; i<A->nrows; i++) {
    if ((A->x->rows[i][n] & ~mask_end)   !=  (m4rie_canary & mask_field & ~mask_end)) {
      return 0;
    }
  }
  return 1;
}

static inline int mzed_interval_clean(mzed_t *A) {
  const word mask_end  = __M4RI_LEFT_BITMASK(A->x->ncols % m4ri_radix);
  const word mask_field = field_mask(A->finite_field);
  for(rci_t i=0; i<A->nrows; i++) {
    for(wi_t j=0; j<A->x->width-1; j++) {
      if (A->x->rows[i][j] & mask_field)
        return 0;
    }
    if (A->x->rows[i][A->x->width-1] & mask_field & mask_end)
      return 0;
  }
  return 1;
}

static inline void mzd_slice_set_canary(mzd_slice_t *A) {
  const word mask_end   = __M4RI_LEFT_BITMASK(A->ncols % m4ri_radix);
  const rci_t n = A->x[0]->width-1;

  for(unsigned int e=0; e<A->finite_field->degree; e++) {
    for(rci_t i=0; i<A->nrows; i++) {
      A->x[e]->rows[i][n] = (A->x[e]->rows[i][n] & mask_end)   | (m4rie_canary & ~mask_end);
    }
  }
}

static inline void mzd_slice_clear_canary(mzd_slice_t *A) {
  const word mask_end   = __M4RI_LEFT_BITMASK(A->ncols % m4ri_radix);
  const rci_t n = A->x[0]->width-1;

  for(int e=0; e<A->finite_field->degree; e++) {
    for(rci_t i=0; i<A->nrows; i++) {
      A->x[e]->rows[i][n] &=mask_end;
    }
  }
}

static inline int mzd_slice_canary_is_alive(mzd_slice_t *A) {
  const word mask_end   = __M4RI_LEFT_BITMASK(A->ncols % m4ri_radix);
  const rci_t n = A->x[0]->width-1;

  for(unsigned int e=0; e<A->finite_field->degree; e++) {
    for(rci_t i=0; i<A->nrows; i++) {
      if ((A->x[e]->rows[i][n] & ~mask_end)   != (m4rie_canary & ~mask_end)) {
        return 0;
      }
    }
  }
  return 1;
}

static inline mzed_t *random_mzed_t(gf2e *ff, int m, int n) {
  mzed_t *A  = mzed_init(ff,m,n);
  mzed_randomize(A);
  mzed_set_canary(A);
  return A;
}

static inline mzd_slice_t *random_mzd_slice_t(gf2e *ff, int m, int n) {
  mzd_slice_t *A = mzd_slice_init(ff,m,n);
  mzd_slice_randomize(A);
  mzd_slice_set_canary(A);
  return A;
}

static inline mzed_t *random_mzed_t_rank(gf2e *ff, const rci_t m, const rci_t n, const rci_t r) {
  mzed_t *U = mzed_init(ff, m, n);
  mzed_t *Ur = mzed_init_window(U, 0, 0, r, U->ncols);
  mzed_t *L = mzed_init(ff, m, m);

  mzed_randomize(L);
  mzed_randomize(Ur);

  for(rci_t i=0; i<r; i++) {
    for(rci_t j=i+1; j<L->ncols; j++) {
      mzed_write_elem(L, i, j, 0);
    }
    mzed_write_elem(L, i, i, 1); 
  }
  for(rci_t i=r; i<L->nrows; i++) {
    for(rci_t j=r+1; j < L->ncols; j++) {
      mzed_write_elem(L, i, j, 0);
    }
  }

  for(rci_t i=0; i<r; i++) {
    mzed_write_elem(U, i, i, 1);
    for(rci_t j=0; j<i; j++) {
      mzed_write_elem(U, i, j, 0);
    }
  }
  mzed_t *A = mzed_mul(NULL, L, U);
  mzed_free(L);
  mzed_free_window(Ur);
  mzed_free(U);

  for(rci_t i=0; i<A->nrows; i++) {
    const rci_t ii = random() % A->nrows;
    mzed_row_swap(A, i, ii);
  };
  for(rci_t i=0; i<A->ncols; i++) {
    const rci_t ii = random() % A->ncols;
    mzed_col_swap(A, i, ii);
  };
  mzed_set_canary(A);
  return A;
}


static inline int parse_parameters(int argc, char **argv) {
  int runlong = 0;
  int c;
  while ((c = getopt(argc, argv, "l")) != -1) {
    switch (c) {
    case 'l':
      runlong = 1;
      break;
    case '?':
      printf(" -l   run long tests.\n");
      abort();
    default:
      abort();
    }
  }
  return runlong;
}

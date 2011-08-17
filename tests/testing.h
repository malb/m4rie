#include <m4rie.h>

# define m4rie_check(expr)						\
  if (!expr) {								\
    fail_ret += 1;                                                      \
    printf("%s in %s:%d failed\n",__STRING(expr), __FILE__, __LINE__);  \
  } 


static inline mzed_t *random_mzed_t(gf2e *ff, int m, int n) {
  mzed_t *A  = mzed_init(ff,m,n);
  mzed_randomize(A);
  return A;
}

static inline mzd_slice_t *random_mzd_slice_t(gf2e *ff, int m, int n) {
  mzd_slice_t *A  = mzd_slice_init(ff,m,n);
  mzd_slice_randomize(A);
  return A;
}

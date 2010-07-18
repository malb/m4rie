static inline int __mzd_read_bits(const mzd_t *a, const size_t x, const size_t y, const size_t n) {
  word temp =  a->rows[x][(y+a->offset) / RADIX]; /* get the value */
  temp <<= (y+a->offset)%RADIX; /* clear upper bits */
  temp >>= RADIX - n; /* clear lower bits and move to correct position.*/
  return temp;
}

static inline void __mzd_xor_bits(const mzd_t *a, const size_t x, const size_t y, const int n, word values) {
  word *temp =  a->rows[x] + (y+a->offset) / RADIX;
  *temp ^= values<<(RADIX-((y+a->offset)%RADIX)-n);
}

static inline void __mzd_clear_bits(const mzd_t *M, const size_t x, const size_t y, const int n) {
  word temp =  M->rows[x][(y+M->offset) / RADIX];
  temp <<= (y+M->offset)%RADIX; /* clear upper bits */
  temp >>= RADIX-n; /* clear lower bits and move to correct position.*/
  temp <<= RADIX-n - (y+M->offset)%RADIX;
  M->rows[x][(y+M->offset) / RADIX] ^= temp;
}

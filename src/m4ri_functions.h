static inline word __mzd_read_bits(const mzd_t *M, const size_t x, const size_t y, const size_t n) {
  int const spot = (y + M->offset) % RADIX;
  wi_t const block = (y + M->offset) / RADIX;
  int const spill = spot + n - RADIX;
  word temp = M->rows[x][block] << -spill;
  return temp >> (RADIX - n);
}

static inline void __mzd_xor_bits(const mzd_t *M, const size_t x, const size_t y, const size_t n, word values) {
  int const spot = (y + M->offset) % RADIX;
  wi_t const block = (y + M->offset) / RADIX;
  M->rows[x][block] ^= values << spot;
}

static inline void __mzd_clear_bits(const mzd_t *M, const size_t x, const size_t y, const size_t n) {
  word values = FFFF >> (RADIX - n);
  int const spot = (y + M->offset) % RADIX;
  wi_t const block = (y + M->offset) / RADIX;
  M->rows[x][block] &= ~(values << spot);
}

static inline word __mzd_read_bits(const mzd_t *M, const size_t x, const size_t y, const size_t n) {
  int const spot = (y + M->offset) % m4ri_radix;
  wi_t const block = (y + M->offset) / m4ri_radix;
  int const spill = spot + n - m4ri_radix;
  word temp = M->rows[x][block] << -spill;
  return temp >> (m4ri_radix - n);
}

static inline void __mzd_xor_bits(const mzd_t *M, const size_t x, const size_t y, const size_t n, word values) {
  int const spot = (y + M->offset) % m4ri_radix;
  wi_t const block = (y + M->offset) / m4ri_radix;
  M->rows[x][block] ^= values << spot;
}

static inline void __mzd_clear_bits(const mzd_t *M, const size_t x, const size_t y, const size_t n) {
  word values = m4ri_ffff >> (m4ri_radix - n);
  int const spot = (y + M->offset) % m4ri_radix;
  wi_t const block = (y + M->offset) / m4ri_radix;
  M->rows[x][block] &= ~(values << spot);
}

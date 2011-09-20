#ifndef FINITE_FIELD_H
#define FINITE_FIELD_H

#include <m4ri/m4ri.h>

typedef struct {
  word **mul;
  word *inv;
  size_t degree;
  word minpoly;
} gf2e;

void gf2e_free(gf2e *ff);


static inline size_t gf2e_degree_to_w(gf2e *ff) {
  switch(ff->degree) {
  case 2: 
    return 2; 
  case 3: 
  case 4: 
    return 4;
  case 5: 
  case 6: 
  case 7: 
  case 8: 
    return 8;
  case 9: 
  case 10: 
  case 11:
  case 12:
  case 13:
  case 14:
  case 15:
  case 16:
    return 16;
  default:
    m4ri_die("degree %d not supported.\n",ff->degree);
  }
  return 0;
}

static inline word *gf2e_t16_init(gf2e *ff, const word a) {
  word *mul = (word*)m4ri_mm_calloc(1<<16, sizeof(word));

  const int w = gf2e_degree_to_w(ff);
  const word mask_w = (1<<w)-1;
  const word *x = ff->mul[a];

  /**
   * @todo: this is a bit of overkill, we could do better
   */
  for(word i=0; i<1<<16; i++) {
    switch(w) {
    case 2:
      mul[i]  = x[(i&mask_w)] | x[((i>>2)&mask_w)]<<2 | x[((i>>4)&mask_w)]<<4 | x[((i>>6)&mask_w)]<<6;
      mul[i] |= x[((i>>8)&mask_w)]<<8 | x[((i>>10)&mask_w)]<<10 | x[((i>>12)&mask_w)]<<12 | x[((i>>14)&mask_w)]<<14;
      break;
    case 4:
      mul[i]  = x[(i&mask_w)] | x[((i>>4)&mask_w)]<<4 | x[((i>>8)&mask_w)]<<8 | x[((i>>12)&mask_w)]<<12;
      break;
    case 8:
      mul[i]  = x[(i&mask_w)] | x[((i>>8)&mask_w)]<<8;
      break;
    case 16:
      mul[i]  = x[(i&mask_w)];
      break;
    };
  }
  return mul;
}

static inline void gf2e_t16_free(word *mul) {
  m4ri_mm_free(mul);
}

#endif //FINITE_FIELD_H

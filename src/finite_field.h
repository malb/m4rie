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

#endif //FINITE_FIELD_H

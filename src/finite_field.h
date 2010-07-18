#ifndef FINITE_FIELD_H
#define FINITE_FIELD_H

#include <m4ri/m4ri.h>

typedef struct {
  word **mul;
  word *inv;
  size_t degree;
} gf2e;

void gf2e_free(gf2e *ff);


#endif //FINITE_FIELD_H

#include <m4ri/m4ri.h>
#include "finite_field.h"

void gf2e_free(gf2e *ff) {
  
  for(size_t i=0; i<TWOPOW(ff->degree); i++) {
    m4ri_mm_free(ff->mul[i]);
  }
  m4ri_mm_free(ff->mul);
  m4ri_mm_free(ff->inv);
}


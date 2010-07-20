
/******************************************************************************
*
*            M4RIE: Linear Algebra over GF(2^e)
*
*    Copyright (C) 2010 Martin Albrecht <martinralbrecht@googlemail.com>
*
*  Distributed under the terms of the GNU General Public License (GEL)
*  version 2 or higher.
*
*    This code is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*    General Public License for more details.
*
*  The full text of the GPL is available at:
*
*                  http://www.gnu.org/licenses/
******************************************************************************/

#include "config.h"

#include <m4ri/misc.h>
#include <m4ri/packedmatrix.h>
#include <m4ri/brilliantrussian.h>
#include <m4ri/xor.h>

#include "travolta.h"

/**
 * Compute C[rc,i] = C[rc,i] + T0[r0,i] + ... + T3[r3,i] for 0 <= i < ncols
 *
 * \param C Matrix
 * \apram rc Row index
 * \param T0 Matrix
 * \param r0 Row index
 * \param T1 Matrix
 * \param r1 Row index
 * \param T2 Matrix
 * \param r2 Row index
 * \param T3 Matrix
 * \param r3 Row index
 *
 * \wordoffset
 */

static inline void mzed_combine4(mzed_t *C, size_t rc, 
                                 mzed_t *T0, size_t r0, mzed_t *T1, size_t r1, mzed_t *T2, size_t r2, mzed_t *T3, size_t r3) {
#ifndef HAVE_SSE2
  size_t ii=0;
#endif
  _mzd_combine4(C->x->rows[rc], 
                T0->x->rows[r0], T1->x->rows[r1], T2->x->rows[r2], T3->x->rows[r3], 
                C->x->width);
}

/**
 * Compute C[rc,i] = C[rc,i] + T0[r0,i] + ... + T7[r7,i] for 0 <= i < ncols
 *
 * \param C Matrix
 * \apram rc Row index
 * \param T0 Matrix
 * \param r0 Row index
 * \param T1 Matrix
 * \param r1 Row index
 * \param T2 Matrix
 * \param r2 Row index
 * \param T3 Matrix
 * \param r3 Row index
 * \param T4 Matrix
 * \param r4 Row index
 * \param T5 Matrix
 * \param r5 Row index
 * \param T6 Matrix
 * \param r6 Row index
 * \param T7 Matrix
 * \param r7 Row index
 *
 * \wordoffset
 */

static inline void mzed_combine8(mzed_t *C, size_t rc, 
                                 mzed_t *T0, size_t r0, mzed_t *T1, size_t r1, mzed_t *T2, size_t r2, mzed_t *T3, size_t r3,
                                 mzed_t *T4, size_t r4, mzed_t *T5, size_t r5, mzed_t *T6, size_t r6, mzed_t *T7, size_t r7) {
#ifndef HAVE_SSE2
  size_t ii=0;
#endif
  _mzd_combine8(C->x->rows[rc], 
                T0->x->rows[r0], T1->x->rows[r1], T2->x->rows[r2], T3->x->rows[r3], 
                T4->x->rows[r4], T5->x->rows[r5], T6->x->rows[r6], T7->x->rows[r7], 
                C->x->width);
}


/**
 * \brief Perform Gaussian reduction to reduced row echelon form on a
 * submatrix.
 * 
 * The submatrix has dimension at most k starting at r x c of A. Checks
 * for pivot rows up to row endrow (exclusive). Terminates as soon as
 * finding a pivot column fails.
 *
 * \param A Matrix.
 * \param r First row.
 * \param c First column.
 * \param k Maximal dimension of identity matrix to produce.
 * \param end_row Maximal row index (exclusive) for rows to consider
 * for inclusion.
 */

size_t _mzed_gauss_submatrix_full(mzed_t *A, size_t r, size_t c, size_t end_row, int k) {
  size_t i,j,l;
  size_t start_row = r;
  int found;
  word tmp;

  gf2e *ff = A->finite_field;

  for (j=c; j<c+k; j++) {
    found = 0;
    for (i=start_row; i< end_row; i++) {
      /* first we need to clear the first columns */
      for (l=0; l<j-c; l++) {
        tmp = mzed_read_elem(A, i, c+l);
        if (tmp) mzed_add_multiple_of_row(A, i, A, r+l, ff->mul[tmp], c+l);
      }
      /* pivot? */
      const word x = mzed_read_elem(A, i, j);
      if (x) {
        mzed_rescale_row(A, i, j, ff->mul[ff->inv[x]]);
        mzd_row_swap(A->x, i, start_row);

        /* clear above */
        for (l=r; l<start_row; l++) {
          tmp = mzed_read_elem(A, l, j);
          if (tmp) mzed_add_multiple_of_row(A, l, A, start_row, ff->mul[tmp], j);
        }
        start_row++;
        found = 1;
        break;
      }
    }
    if (found==0) {
      return j - c;
    }
  }
  return j - c;
}


void mzed_make_table(const mzed_t *A, size_t r, size_t c, mzed_t *T,  size_t *L, gf2e *ff) {
  mzd_set_ui(T->x,0);

  for(size_t i=0; i< TWOPOW(ff->degree); i++) {
    word *X = ff->mul[i];
    L[i] = i;
    mzed_add_multiple_of_row(T, i, A, r, X, c);
  }
}

size_t mzed_echelonize_travolta(mzed_t *A, int full) {
  gf2e* ff = A->finite_field;

  size_t r,c;

  size_t k = ff->degree;

  /** cf. mzd_echelonize_m4ri **/
  size_t kk = (size_t)m4ri_opt_k(A->x->nrows, A->x->ncols, 0);
  if (kk>=7) 
    kk = 7;
  if ( (6*(1<<kk)*A->ncols / 8.0) > CPU_L2_CACHE / 2.0 )
    kk -= 1;
  kk = (6*kk)/k;

  /** enforcing bounds **/
  if (kk == 0)
    kk = 1;
  else if (kk > 6)
    kk = 6;

  size_t kbar = 0;

  mzed_t *T0 = mzed_init(ff, TWOPOW(k), A->ncols);
  mzed_t *T1 = mzed_init(ff, TWOPOW(k), A->ncols);
  mzed_t *T2 = mzed_init(ff, TWOPOW(k), A->ncols);
  mzed_t *T3 = mzed_init(ff, TWOPOW(k), A->ncols);
  mzed_t *T4 = mzed_init(ff, TWOPOW(k), A->ncols);
  mzed_t *T5 = mzed_init(ff, TWOPOW(k), A->ncols);
  
  /* this is dummy, we keep it for compatibility with the M4RI functions */
  size_t *L0 = (size_t *)m4ri_mm_calloc(TWOPOW(k), sizeof(size_t));
  size_t *L1 = (size_t *)m4ri_mm_calloc(TWOPOW(k), sizeof(size_t));
  size_t *L2 = (size_t *)m4ri_mm_calloc(TWOPOW(k), sizeof(size_t));
  size_t *L3 = (size_t *)m4ri_mm_calloc(TWOPOW(k), sizeof(size_t));
  size_t *L4 = (size_t *)m4ri_mm_calloc(TWOPOW(k), sizeof(size_t));
  size_t *L5 = (size_t *)m4ri_mm_calloc(TWOPOW(k), sizeof(size_t));


  r = 0;
  c = 0;
  while(c < A->ncols) {
    /**
     * \todo: If full == False we should switch over to naive once the
     *        remain matrix is small.
     */
    if(c+kk > A->ncols) kk = A->ncols - c;

    /**
     * \todo: we don't really compute the upper triangular form yet,
     *        we need to implement _mzed_gauss_submatrix and a better
     *        table creation for that.
     */ 
    kbar = _mzed_gauss_submatrix_full(A, r, c, A->nrows, kk);

    if (kbar == 6)  {
      mzed_make_table( A, r,   c,   T0, L0, ff);
      mzed_make_table( A, r+1, c+1, T1, L1, ff);
      mzed_make_table( A, r+2, c+2, T2, L2, ff);
      mzed_make_table( A, r+3, c+3, T3, L3, ff);
      mzed_make_table( A, r+4, c+4, T4, L4, ff);
      mzed_make_table( A, r+5, c+5, T5, L5, ff);
      if(kbar == kk)
        mzed_process_rows6( A, r+6, A->nrows, c, T0, L0, T1, L1, T2, L2, T3, L3, T4, L4, T5, L5);
      if(full)
        mzed_process_rows6( A,   0,        r, c, T0, L0, T1, L1, T2, L2, T3, L3, T4, L4, T5, L5);
    } else if(kbar == 5) {
      mzed_make_table( A, r,     c, T0, L0, ff);
      mzed_make_table( A, r+1, c+1, T1, L1, ff);
      mzed_make_table( A, r+2, c+2, T2, L2, ff);
      mzed_make_table( A, r+3, c+3, T3, L3, ff);
      mzed_make_table( A, r+4, c+4, T4, L4, ff);
      if(kbar == kk)
        mzed_process_rows5( A, r+5, A->nrows, c, T0, L0, T1, L1, T2, L2, T3, L3, T4, L4);
      if(full)
        mzed_process_rows5( A,   0,        r, c, T0, L0, T1, L1, T2, L2, T3, L3, T4, L4);

    } else if(kbar == 4) {
      mzed_make_table( A, r,   c,   T0, L0, ff);
      mzed_make_table( A, r+1, c+1, T1, L1, ff);
      mzed_make_table( A, r+2, c+2, T2, L2, ff);
      mzed_make_table( A, r+3, c+3, T3, L3, ff);
      if(kbar == kk)
        mzed_process_rows4( A, r+4, A->nrows, c, T0, L0, T1, L1, T2, L2, T3, L3);
      if(full)
        mzed_process_rows4( A,   0,        r, c, T0, L0, T1, L1, T2, L2, T3, L3);

    } else if(kbar == 3) {
      mzed_make_table( A, r,   c,   T0, L0, ff);
      mzed_make_table( A, r+1, c+1, T1, L1, ff);
      mzed_make_table( A, r+2, c+2, T2, L2, ff);
      if(kbar == kk)
        mzed_process_rows3( A, r+3, A->nrows, c, T0, L0, T1, L1, T2, L2);
      if(full)
        mzed_process_rows3( A,   0,        r, c, T0, L0, T1, L1, T2, L2);

    } else if(kbar == 2) {
      mzed_make_table( A, r,   c,   T0, L0, ff);
      mzed_make_table( A, r+1, c+1, T1, L1, ff);
      if(kbar == kk)
        mzed_process_rows2( A, r+2, A->nrows, c, T0, L0, T1, L1);
      if(full)
        mzed_process_rows2( A,   0,        r, c, T0, L0, T1, L1);

    } else if (kbar == 1) {
      mzed_make_table(A, r, c, T0, L0, ff);
      if(kbar == kk)
        mzed_process_rows( A, r+1, A->nrows, c, T0, L0);
      if(full)
        mzed_process_rows( A,   0,        r, c, T0, L0);

    } else {
      c++;
    }
    r += kbar;
    c += kbar;
  }

  m4ri_mm_free(L0); m4ri_mm_free(L1); m4ri_mm_free(L2);
  m4ri_mm_free(L3); m4ri_mm_free(L4); m4ri_mm_free(L5);
  mzed_free(T0); mzed_free(T1); mzed_free(T2);
  mzed_free(T3); mzed_free(T4); mzed_free(T5);
  return r;
}

mzed_t *_mzed_mul_travolta0(mzed_t *C, const mzed_t *A, const mzed_t *B) {
  mzed_t *T0 = mzed_init(C->finite_field, TWOPOW(A->finite_field->degree), B->ncols);
  size_t *L0 = (size_t*)m4ri_mm_calloc(TWOPOW(A->finite_field->degree), sizeof(size_t));

  for(size_t i=0; i < A->ncols; i++) {
    mzed_make_table(B, i, 0, T0, L0, A->finite_field);
    for(size_t j=0; j<A->nrows; j++)
      mzd_combine(C->x, j, 0, C->x, j, 0, T0->x, mzed_read_elem(A, j, i), 0);
  }
  mzed_free(T0);
  m4ri_mm_free(L0);
  return C;
}

mzed_t *_mzed_mul_travolta1(mzed_t *C, const mzed_t *A, const mzed_t *B) {
  mzed_t *T0 = mzed_init(C->finite_field, TWOPOW(A->finite_field->degree), B->ncols);
  mzed_t *T1 = mzed_init(C->finite_field, TWOPOW(A->finite_field->degree), B->ncols);
  mzed_t *T2 = mzed_init(C->finite_field, TWOPOW(A->finite_field->degree), B->ncols);
  mzed_t *T3 = mzed_init(C->finite_field, TWOPOW(A->finite_field->degree), B->ncols);
  mzed_t *T4 = mzed_init(C->finite_field, TWOPOW(A->finite_field->degree), B->ncols);
  mzed_t *T5 = mzed_init(C->finite_field, TWOPOW(A->finite_field->degree), B->ncols);
  mzed_t *T6 = mzed_init(C->finite_field, TWOPOW(A->finite_field->degree), B->ncols);
  mzed_t *T7 = mzed_init(C->finite_field, TWOPOW(A->finite_field->degree), B->ncols);

  size_t *L0 = (size_t*)m4ri_mm_calloc(TWOPOW(A->finite_field->degree), sizeof(size_t));
  size_t *L1 = (size_t*)m4ri_mm_calloc(TWOPOW(A->finite_field->degree), sizeof(size_t));
  size_t *L2 = (size_t*)m4ri_mm_calloc(TWOPOW(A->finite_field->degree), sizeof(size_t));
  size_t *L3 = (size_t*)m4ri_mm_calloc(TWOPOW(A->finite_field->degree), sizeof(size_t));
  size_t *L4 = (size_t*)m4ri_mm_calloc(TWOPOW(A->finite_field->degree), sizeof(size_t));
  size_t *L5 = (size_t*)m4ri_mm_calloc(TWOPOW(A->finite_field->degree), sizeof(size_t));
  size_t *L6 = (size_t*)m4ri_mm_calloc(TWOPOW(A->finite_field->degree), sizeof(size_t));
  size_t *L7 = (size_t*)m4ri_mm_calloc(TWOPOW(A->finite_field->degree), sizeof(size_t));
  
  const size_t kk = 8;
  const size_t end = A->ncols/kk;

  /* we are effectively disabling the cache optimisation since table
     creation is so expensive in this context. 
  */
  size_t blocksize = 1ULL<<31;

  /*
   * it seems to give some advantage for larger matrices over GF(2^2)
   * though
   */
   if (A->w == 2 && A->nrows >= 2*MZD_MUL_BLOCKSIZE) 
     blocksize = MZD_MUL_BLOCKSIZE/A->w; 

  size_t giantstep, babystep;

  for (giantstep=0; giantstep + blocksize <= A->nrows; giantstep += blocksize) {
    for(size_t i=0; i < end; i++) {
      mzed_make_table( B, kk*i  , 0, T0, L0, A->finite_field);
      mzed_make_table( B, kk*i+1, 0, T1, L1, A->finite_field);
      mzed_make_table( B, kk*i+2, 0, T2, L2, A->finite_field);
      mzed_make_table( B, kk*i+3, 0, T3, L3, A->finite_field);
      mzed_make_table( B, kk*i+4, 0, T4, L4, A->finite_field);
      mzed_make_table( B, kk*i+5, 0, T5, L5, A->finite_field);
      mzed_make_table( B, kk*i+6, 0, T6, L6, A->finite_field);
      mzed_make_table( B, kk*i+7, 0, T7, L7, A->finite_field);
      for(babystep = 0; babystep < blocksize; babystep++) {
        const size_t j = giantstep + babystep;
        const size_t x0 = mzed_read_elem(A, j, kk*  i);
        const size_t x1 = mzed_read_elem(A, j, kk*i+1);
        const size_t x2 = mzed_read_elem(A, j, kk*i+2);
        const size_t x3 = mzed_read_elem(A, j, kk*i+3);
        const size_t x4 = mzed_read_elem(A, j, kk*i+4);
        const size_t x5 = mzed_read_elem(A, j, kk*i+5);
        const size_t x6 = mzed_read_elem(A, j, kk*i+6);
        const size_t x7 = mzed_read_elem(A, j, kk*i+7);
        mzed_combine8(C, j, T0, x0, T1, x1, T2, x2, T3, x3, T4, x4, T5, x5, T6, x6, T7, x7);
      }
    }
  }

  /* last giant step */
  for(size_t i=0; i < end; i++) {
    mzed_make_table( B, kk*i  , 0, T0, L0, A->finite_field);
    mzed_make_table( B, kk*i+1, 0, T1, L1, A->finite_field);
    mzed_make_table( B, kk*i+2, 0, T2, L2, A->finite_field);
    mzed_make_table( B, kk*i+3, 0, T3, L3, A->finite_field);
    mzed_make_table( B, kk*i+4, 0, T4, L4, A->finite_field);
    mzed_make_table( B, kk*i+5, 0, T5, L5, A->finite_field);
    mzed_make_table( B, kk*i+6, 0, T6, L6, A->finite_field);
    mzed_make_table( B, kk*i+7, 0, T7, L7, A->finite_field);
    for(babystep = 0; babystep < A->nrows - giantstep; babystep++) {
      const size_t j = giantstep + babystep;
      const size_t x0 = mzed_read_elem(A, j, kk*  i);
      const size_t x1 = mzed_read_elem(A, j, kk*i+1);
      const size_t x2 = mzed_read_elem(A, j, kk*i+2);
      const size_t x3 = mzed_read_elem(A, j, kk*i+3);
      const size_t x4 = mzed_read_elem(A, j, kk*i+4);
      const size_t x5 = mzed_read_elem(A, j, kk*i+5);
      const size_t x6 = mzed_read_elem(A, j, kk*i+6);
      const size_t x7 = mzed_read_elem(A, j, kk*i+7);
      mzed_combine8(C, j, T0, x0, T1, x1, T2, x2, T3, x3, T4, x4, T5, x5, T6, x6, T7, x7);
    }
  }
  
  if (A->ncols%kk) {
    for(size_t i=kk*end; i < A->ncols; i++) {
      mzed_make_table(B, i, 0, T0, L0, A->finite_field);
      for(size_t j=0; j<A->nrows; j++)
        mzd_combine(C->x, j, 0, C->x, j, 0, T0->x, mzed_read_elem(A, j, i), 0);
    }
  }

  mzed_free(T0); mzed_free(T1);  mzed_free(T2); mzed_free(T3);
  mzed_free(T4); mzed_free(T5);  mzed_free(T6); mzed_free(T7);
  m4ri_mm_free(L0); m4ri_mm_free(L1); m4ri_mm_free(L2); m4ri_mm_free(L3);
  m4ri_mm_free(L4); m4ri_mm_free(L5); m4ri_mm_free(L6); m4ri_mm_free(L7);
  return C;
}

mzed_t *mzed_addmul_travolta(mzed_t *C, const mzed_t *A, const mzed_t *B) {
  if (C->nrows != A->nrows || C->ncols != B->ncols || C->finite_field != A->finite_field) {
    m4ri_die("mzd_mul_naive: Provided return matrix has wrong dimensions or wrong base field.\n");
  }
  return _mzed_mul_travolta1(C, A, B);
}

mzed_t *mzed_mul_travolta(mzed_t *C, const mzed_t *A, const mzed_t *B) {
  if (C==NULL) {
    C = mzed_init(A->finite_field, A->nrows, B->ncols);
  } else {
    if (C->nrows != A->nrows || C->ncols != B->ncols || C->finite_field != A->finite_field) {
      m4ri_die("mzd_mul_naive: Provided return matrix has wrong dimensions or wrong base field.\n");
    }
    mzd_set_ui(C->x, 0);
  }
  return _mzed_mul_travolta1(C, A, B);
}

mzed_t *mzed_invert_travolta(mzed_t *B, const mzed_t *A) {
  assert(A->nrows == A->ncols);
  mzed_t *I = mzed_init(A->finite_field, A->nrows, A->ncols);
  mzed_set_ui(I, 1);
  mzed_t *T = mzed_concat(NULL, A, I);
  mzed_free(I);

  size_t r = mzed_echelonize_travolta(T, 1);
  if (r != A->nrows) 
    m4ri_die("mzed_invert_travolta: input matrix does not have full rank.");
  B = mzed_submatrix(B, T, 0, A->ncols, A->nrows, T->ncols);
  mzed_free(T);
  return B;
}

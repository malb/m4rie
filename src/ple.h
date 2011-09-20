/**
 * \file ple.h
 * \brief L*E = P*A decomposition.
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 */

#ifndef M4RIE_PLE_H
#define M4RIE_PLE_H

/******************************************************************************
*
*            M4RIE: Linear Algebra over GF(2^e)
*
*    Copyright (C) 2011 Martin Albrecht <martinralbrecht@googlemail.com>
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

#include <m4ri/m4ri.h>
#include "gf2e_matrix.h"
#include "bitslice.h"

/**
 * \brief PLE decomposition: L*E = P*A 
 *
 * Modifies A in place to store lower triangular L below (and on) the
 * main diagonal and E -- a echelon form of A -- above the main
 * diagonal (pivots are stored in Q). P and Q are updated with row and
 * column permutations respectively.
 *
 * This function uses naive cubic PLE decomposition depending on the
 * size of the underlying field.
 *
 * \param A Matrix
 * \param P Permutation vector of length A->nrows
 * \param Q Permutation vector of length A->ncols
 *
 * \ingroup PLE
 * 
 * \sa mzed_ple_naive mzed_ple_travolta _mzed_ple
 * 
 * \example tests/test_ple.cc bench/bench_ple.cc
 */

rci_t mzed_ple_naive(mzed_t *A, mzp_t *P, mzp_t *Q);

/**
 * \brief PLE decomposition: L*E = P*A 
 *
 * Modifies A in place to store lower triangular L below (and on) the
 * main diagonal and E -- a echelon form of A -- above the main
 * diagonal (pivots are stored in Q). P and Q are updated with row and
 * column permutations respectively.
 *
 * This function uses either asymptotically fast PLE decomposition by
 * reducing it to matrix multiplication or naive cubic PLE
 * decomposition depending on the size of the underlying field. If
 * asymptotically fast PLE decomposition is used, then the algorithm
 * switches to mzed_ple_travolta if e * ncols * nrows is <= cutoff
 * where e is the exponent of the finite field.
 *
 * \param A Matrix
 * \param P Permutation vector of length A->nrows
 * \param Q Permutation vector of length A->ncols
 * \param cutoff Integer
 *
 * \ingroup PLE
 * 
 * \sa mzed_ple_naive mzed_ple_travolta _mzed_ple
 * 
 * \example tests/test_ple.cc bench/bench_ple.cc
 */

rci_t _mzd_slice_ple(mzd_slice_t *A, mzp_t *P, mzp_t *Q, rci_t cutoff);

/**
 * \brief PLE decomposition: L*E = P*A 
 *
 * Modifies A in place to store lower triangular L below (and on) the
 * main diagonal and E -- a echelon form of A -- above the main
 * diagonal (pivots are stored in Q). P and Q are updated with row and
 * column permutations respectively.
 *
 * This function implements asymptotically fast PLE decomposition by
 * reducing it to matrix multiplication.
 *
 * \param A Matrix
 * \param P Permutation vector of length A->nrows
 * \param Q Permutation vector of length A->ncols
 *
 * \ingroup PLE
 * 
 * \sa mzed_ple_naive mzed_ple_travolta _mzd_slice_ple
 * 
 * \example tests/test_ple.cc bench/bench_ple.cc
 */

static inline rci_t mzd_slice_ple(mzd_slice_t *A, mzp_t *P, mzp_t *Q) {
  assert(P->length == A->nrows);
  assert(Q->length == A->ncols);
  return _mzd_slice_ple(A, P, Q, 0);
}

rci_t _mzd_slice_pluq(mzd_slice_t *A, mzp_t *P, mzp_t *Q, rci_t cutoff);

static inline rci_t mzd_slice_pluq(mzd_slice_t *A, mzp_t *P, mzp_t *Q) {
  assert(P->length == A->nrows);
  assert(Q->length == A->ncols);
  return _mzd_slice_pluq(A, P, Q, 0);
}


/**
 * \brief PLE decomposition: L*E = P*A 
 *
 * Modifies A in place to store lower triangular L below (and on) the
 * main diagonal and E -- a echelon form of A -- above the main
 * diagonal (pivots are stored in Q). P and Q are updated with row and
 * column permutations respectively.
 *
 * This function uses either asymptotically fast PLE decomposition by
 * reducing it to matrix multiplication or naive cubic PLE
 * decomposition depending on the size of the underlying field. If
 * asymptotically fast PLE decomposition is used, then the algorithm
 * switches to mzed_ple_travolta if e * ncols * nrows is <= cutoff
 * where e is the exponent of the finite field.
 *
 * \param A Matrix
 * \param P Permutation vector of length A->nrows
 * \param Q Permutation vector of length A->ncols
 * \param cutoff Integer >= 0
 *
 * \ingroup PLE
 * 
 * \sa mzed_ple_naive mzed_ple_travolta _mzed_ple
 * 
 * \example tests/test_ple.cc bench/bench_ple.cc
 */

rci_t _mzed_ple(mzed_t *A, mzp_t *P, mzp_t *Q, rci_t cutoff);

/**
 * \brief PLE decomposition: L*E = P*A 
 *
 * Modifies A in place to store lower triangular L below (and on) the
 * main diagonal and E -- a echelon form of A -- above the main
 * diagonal (pivots are stored in Q). P and Q are updated with row and
 * column permutations respectively.
 *
 * This function uses either asymptotically fast PLE decomposition by
 * reducing it to matrix multiplication or naive cubic PLE
 * decomposition depending on the size of the underlying field.
 *
 * \param A Matrix
 * \param P Permutation vector of length A->nrows
 * \param Q Permutation vector of length A->ncols
 *
 * \ingroup PLE
 * 
 * \sa mzed_ple_naive mzed_ple_travolta _mzed_ple
 * 
 * \example tests/test_ple.cc bench/bench_ple.cc
 */

static inline rci_t mzed_ple(mzed_t *A, mzp_t *P, mzp_t *Q) {
  return _mzed_ple(A, P, Q, 0);
}

/**
 * Default crossover to PLE base case.
 */

#define __M4RIE_PLE_CUTOFF (__M4RI_CPU_L2_CACHE<<3) //(2048*2048)

#endif //M4RIE_PLE_H

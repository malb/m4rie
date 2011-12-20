/**
 * \file m4rie.h
 * \brief Main include file for the M4RIE library.
 *
 * \author Martin Albrecht <martinralbrecht@googlemail.com>
 *
 * \mainpage
 *
 * M4RIE is a library to do fast compations with dense matrices over
 * \GF2E for small \e. M4RIE is available under the GPLv2+.
 *
 * The two fundamental data types of this library are mzed_t and
 * mzd_slice_t. For big matrices, i.e., those which do not fit into L2
 * cache, it is recommended to use mzd_slice_t and for smaller
 * matrices mzed_t will be slightly faster and use less memory.
 *
 * Function names follow the pattern
\verbatim
[_]_[type]_[what]_[algorithm]
\endverbatim
 *
 * Function names beginning with an underscore perform less
 * consistency checks (matching dimensions, matching fields) than
 * those without, e.g., _mzed_ple() is called by mzed_ple() after some
 * checks were performed.
 *
 * For both data types almost all functions are the same, e.g., there
 * is a function mzd_slice_add() and there also should be a function
 * mzed_add() with the same signature except for the matrix type.
 *
 * Functions which do not specify an algorithm choose the best
 * available algorithm (based on some heuristic), e.g., mzed_ple()
 * might call mzed_ple_newton_john().
 *
 * \defgroup Definitions        Type definitions
 * \defgroup Constructions      Constructions
 * \defgroup Assignment         Assignment and basic manipulation
 * \defgroup RowOperations      Operations on rows
 * \defgroup StringConversions  String conversions and I/O
 * \defgroup Addition           Addition and subtraction
 * \defgroup Multiplication     Multiplication
 * \defgroup PLE                PLE and PLUQ decomposition
 * \defgroup Echelon            Echelon forms
 * \defgroup Triangular         Triangular matrices
 *
 * \example tests/test_elimination.cc tests/test_multiplication.cc
 */

#ifndef M4RIE_H
#define M4RIE_H

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

#ifdef __cplusplus
extern "C" {
#endif //__cplusplus

#include "gf2e.h"
#include "mzed.h"
#include "newton_john.h"
#include "echelonform.h"
#include "strassen.h"
#include "mzd_slice.h"
#include "trsm.h"
#include "ple.h"
#include "conversion.h"
#include "permutation.h"

#ifdef __cplusplus
}
#endif //__cplusplus


#endif //M4RIE_H

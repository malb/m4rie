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
 * \defgroup Definitions        Type definitions
 * \defgroup Construction       Constructions
 * \defgroup Assignment         Assignment and basic manipulation
 * \defgroup Addition           Addition and subtraction
 * \defgroup Multiplication     Multiplication
 * \defgroup Echelon            Echelon forms
 * \defgroup StringConversions  String conversions and I/O
 * \defgroup RowOperations      Operations on rows
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

#include "finite_field.h"
#include "gf2e_matrix.h"
#include "travolta.h"
#include "echelonform.h"
#include "strassen.h"
#include "bitslice.h"
#include "trsm.h"
#include "ple.h"
#include "conversion.h"
#include "permutation.h"

#ifdef __cplusplus
}
#endif //__cplusplus


#endif //M4RIE_H

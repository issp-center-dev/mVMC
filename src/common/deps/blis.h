/*
 * blis.h
   This header is cut out from BLIS v0.8.0,
    containing necessary type definitions & macros.
   User can use this file to avert BLIS dependency
    if their BLAS library vendors a GEMMT implementation.
   Here begins 3-clause BSD-style license of BLIS:

   BLIS
   An object-based framework for developing high-performance BLAS-like
   libraries.

   Copyright (C) 2014, The University of Texas at Austin
   Copyright (C) 2016, Hewlett Packard Enterprise Development LP
   Copyright (C) 2020, Advanced Micro Devices, Inc.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:
    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    - Neither the name(s) of the copyright holder(s) nor the names of its
      contributors may be used to endorse or promote products derived
      from this software without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef BLIS_H
#define BLIS_H


// Allow C++ users to include this header file in their source code. However,
// we make the extern "C" conditional on whether we're using a C++ compiler,
// since regular C compilers don't understand the extern "C" construct.
#ifdef __cplusplus
extern "C" {
#endif

// NOTE: PLEASE DON'T CHANGE THE ORDER IN WHICH HEADERS ARE INCLUDED UNLESS
// YOU ARE SURE THAT IT DOESN'T BREAK INTER-HEADER MACRO DEPENDENCIES.

#include <stdio.h> // skipped
#include <stdlib.h> // skipped
#include <math.h> // skipped
#include <string.h> // skipped
#include <stdarg.h> // skipped
#include <float.h> // skipped
#include <errno.h> // skipped
#include <ctype.h> // skipped
#include <assert.h> // skipped

// -- STATIC INLINE FUNCTIONS --------------------------------------------------

// C and C++ have different semantics for defining "inline" functions. In C,
// the keyword phrase "static inline" accomplishes this, though the "inline"
// is optional. In C++, the "inline" keyword is required and obviates "static"
// altogether. Why does this matter? While BLIS is compiled in C99, blis.h may
// be #included by a source file that is compiled with C++.
#ifdef __cplusplus
  #define BLIS_INLINE inline
#else
  //#define BLIS_INLINE static inline
  #define BLIS_INLINE static
#endif

// -- Common BLIS definitions --

// begin bli_type_defs.h


#ifndef BLIS_TYPE_DEFS_H
#define BLIS_TYPE_DEFS_H


//
// -- BLIS basic types ---------------------------------------------------------
//

#ifdef __cplusplus
  // For C++, include stdint.h.
#include <stdint.h> // skipped
#elif __STDC_VERSION__ >= 199901L
  // For C99 (or later), include stdint.h.
#include <stdint.h> // skipped
#include <stdbool.h> // skipped
#else
  // When stdint.h is not available, manually typedef the types we will use.
  #ifdef _WIN32
    typedef          __int32  int32_t;
    typedef unsigned __int32 uint32_t;
    typedef          __int64  int64_t;
    typedef unsigned __int64 uint64_t;
  #else
    #error "Attempting to compile on pre-C99 system without stdint.h."
  #endif
#endif

// -- General-purpose integers --

// If BLAS integers are 64 bits, mandate that BLIS integers also be 64 bits.
// NOTE: This cpp guard will only meaningfully change BLIS's behavior on
// systems where the BLIS integer size would have been automatically selected
// to be 32 bits, since explicit selection of 32 bits is prohibited at
// configure-time (and explicit or automatic selection of 64 bits is fine
// and would have had the same result).
#if BLIS_BLAS_INT_SIZE == 64
  #undef  BLIS_INT_TYPE_SIZE
  #define BLIS_INT_TYPE_SIZE 64
#endif

// Define integer types depending on what size integer was requested.
#if   BLIS_INT_TYPE_SIZE == 32
typedef           int32_t  gint_t;
typedef          uint32_t guint_t;
#elif BLIS_INT_TYPE_SIZE == 64
typedef           int64_t  gint_t;
typedef          uint64_t guint_t;
#else
typedef   signed long int  gint_t;
typedef unsigned long int guint_t;
#endif

// -- Boolean type --

// NOTE: bool_t is no longer used and has been replaced with C99's bool type.
//typedef bool bool_t;

// BLIS uses TRUE and FALSE macro constants as possible boolean values, but we
// define these macros in terms of true and false, respectively, which are
// defined by C99 in stdbool.h.
#ifndef TRUE
  #define TRUE  true
#endif

#ifndef FALSE
  #define FALSE false
#endif

// -- Special-purpose integers --

// This cpp guard provides a temporary hack to allow libflame
// interoperability with BLIS.
#ifndef _DEFINED_DIM_T
#define _DEFINED_DIM_T
typedef   gint_t dim_t;      // dimension type
#endif
typedef   gint_t inc_t;      // increment/stride type
typedef   gint_t doff_t;     // diagonal offset type
typedef  guint_t siz_t;      // byte size type
typedef uint32_t objbits_t;  // object information bit field

// -- Real types --

// Define the number of floating-point types supported, and the size of the
// largest type.
#define BLIS_NUM_FP_TYPES   4
#define BLIS_MAX_TYPE_SIZE  sizeof(dcomplex)

// There are some places where we need to use sizeof() inside of a C
// preprocessor #if conditional, and so here we define the various sizes
// for those purposes.
#define BLIS_SIZEOF_S      4  // sizeof(float)
#define BLIS_SIZEOF_D      8  // sizeof(double)
#define BLIS_SIZEOF_C      8  // sizeof(scomplex)
#define BLIS_SIZEOF_Z      16 // sizeof(dcomplex)

// -- Complex types --

#ifdef BLIS_ENABLE_C99_COMPLEX

	#if __STDC_VERSION__ >= 199901L
#include <complex.h> // skipped

		// Typedef official complex types to BLIS complex type names.
		typedef  float complex scomplex;
		typedef double complex dcomplex;
	#else
		#error "Configuration requested C99 complex types, but C99 does not appear to be supported."
	#endif

#else // ifndef BLIS_ENABLE_C99_COMPLEX

	// This cpp guard provides a temporary hack to allow libflame
	// interoperability with BLIS.
	#ifndef _DEFINED_SCOMPLEX
	#define _DEFINED_SCOMPLEX
	typedef struct
	{
		float  real;
		float  imag;
	} scomplex;
	#endif

	// This cpp guard provides a temporary hack to allow libflame
	// interoperability with BLIS.
	#ifndef _DEFINED_DCOMPLEX
	#define _DEFINED_DCOMPLEX
	typedef struct
	{
		double real;
		double imag;
	} dcomplex;
	#endif

#endif // BLIS_ENABLE_C99_COMPLEX

// -- Atom type --

// Note: atom types are used to hold "bufferless" scalar object values. Note
// that it needs to be as large as the largest possible scalar value we might
// want to hold. Thus, for now, it is a dcomplex.
typedef dcomplex atom_t;

// -- Fortran-77 types --

// Note: These types are typically only used by BLAS compatibility layer, but
// we must define them even when the compatibility layer isn't being built
// because they also occur in bli_slamch() and bli_dlamch().

// Define f77_int depending on what size of integer was requested.
#if   BLIS_BLAS_INT_TYPE_SIZE == 32
typedef int32_t   f77_int;
#elif BLIS_BLAS_INT_TYPE_SIZE == 64
typedef int64_t   f77_int;
#else
typedef long int  f77_int;
#endif

typedef char      f77_char;
typedef float     f77_float;
typedef double    f77_double;
typedef scomplex  f77_scomplex;
typedef dcomplex  f77_dcomplex;

// -- Void function pointer types --

// Note: This type should be used in any situation where the address of a
// *function* will be conveyed or stored prior to it being typecast back
// to the correct function type. It does not need to be used when conveying
// or storing the address of *data* (such as an array of float or double).

//typedef void (*void_fp)( void );
typedef void* void_fp;


//
// -- BLIS info bit field offsets ----------------------------------------------
//



// info
#define BLIS_DATATYPE_SHIFT                0
#define   BLIS_DOMAIN_SHIFT                0
#define   BLIS_PRECISION_SHIFT             1
#define BLIS_CONJTRANS_SHIFT               3
#define   BLIS_TRANS_SHIFT                 3
#define   BLIS_CONJ_SHIFT                  4
#define BLIS_UPLO_SHIFT                    5
#define   BLIS_UPPER_SHIFT                 5
#define   BLIS_DIAG_SHIFT                  6
#define   BLIS_LOWER_SHIFT                 7
#define BLIS_UNIT_DIAG_SHIFT               8
#define BLIS_INVERT_DIAG_SHIFT             9
#define BLIS_TARGET_DT_SHIFT               10
#define   BLIS_TARGET_DOMAIN_SHIFT         10
#define   BLIS_TARGET_PREC_SHIFT           11
#define BLIS_EXEC_DT_SHIFT                 13
#define   BLIS_EXEC_DOMAIN_SHIFT           13
#define   BLIS_EXEC_PREC_SHIFT             14
#define BLIS_PACK_SCHEMA_SHIFT             16
#define   BLIS_PACK_RC_SHIFT               16
#define   BLIS_PACK_PANEL_SHIFT            17
#define   BLIS_PACK_FORMAT_SHIFT           18
#define   BLIS_PACK_SHIFT                  22
#define BLIS_PACK_REV_IF_UPPER_SHIFT       23
#define BLIS_PACK_REV_IF_LOWER_SHIFT       24
#define BLIS_PACK_BUFFER_SHIFT             25
#define BLIS_STRUC_SHIFT                   27
#define BLIS_COMP_DT_SHIFT                 29
#define   BLIS_COMP_DOMAIN_SHIFT           29
#define   BLIS_COMP_PREC_SHIFT             30

// info2
#define BLIS_SCALAR_DT_SHIFT                0
#define   BLIS_SCALAR_DOMAIN_SHIFT          0
#define   BLIS_SCALAR_PREC_SHIFT            1

//
// -- BLIS info bit field masks ------------------------------------------------
//

// info
#define BLIS_DATATYPE_BITS                 ( 0x7  << BLIS_DATATYPE_SHIFT )
#define   BLIS_DOMAIN_BIT                  ( 0x1  << BLIS_DOMAIN_SHIFT )
#define   BLIS_PRECISION_BIT               ( 0x1  << BLIS_PRECISION_SHIFT )
#define BLIS_CONJTRANS_BITS                ( 0x3  << BLIS_CONJTRANS_SHIFT )
#define   BLIS_TRANS_BIT                   ( 0x1  << BLIS_TRANS_SHIFT )
#define   BLIS_CONJ_BIT                    ( 0x1  << BLIS_CONJ_SHIFT )
#define BLIS_UPLO_BITS                     ( 0x7  << BLIS_UPLO_SHIFT )
#define   BLIS_UPPER_BIT                   ( 0x1  << BLIS_UPPER_SHIFT )
#define   BLIS_DIAG_BIT                    ( 0x1  << BLIS_DIAG_SHIFT )
#define   BLIS_LOWER_BIT                   ( 0x1  << BLIS_LOWER_SHIFT )
#define BLIS_UNIT_DIAG_BIT                 ( 0x1  << BLIS_UNIT_DIAG_SHIFT )
#define BLIS_INVERT_DIAG_BIT               ( 0x1  << BLIS_INVERT_DIAG_SHIFT )
#define BLIS_TARGET_DT_BITS                ( 0x7  << BLIS_TARGET_DT_SHIFT )
#define   BLIS_TARGET_DOMAIN_BIT           ( 0x1  << BLIS_TARGET_DOMAIN_SHIFT )
#define   BLIS_TARGET_PREC_BIT             ( 0x1  << BLIS_TARGET_PREC_SHIFT )
#define BLIS_EXEC_DT_BITS                  ( 0x7  << BLIS_EXEC_DT_SHIFT )
#define   BLIS_EXEC_DOMAIN_BIT             ( 0x1  << BLIS_EXEC_DOMAIN_SHIFT )
#define   BLIS_EXEC_PREC_BIT               ( 0x1  << BLIS_EXEC_PREC_SHIFT )
#define BLIS_PACK_SCHEMA_BITS              ( 0x7F << BLIS_PACK_SCHEMA_SHIFT )
#define   BLIS_PACK_RC_BIT                 ( 0x1  << BLIS_PACK_RC_SHIFT )
#define   BLIS_PACK_PANEL_BIT              ( 0x1  << BLIS_PACK_PANEL_SHIFT )
#define   BLIS_PACK_FORMAT_BITS            ( 0xF  << BLIS_PACK_FORMAT_SHIFT )
#define   BLIS_PACK_BIT                    ( 0x1  << BLIS_PACK_SHIFT )
#define BLIS_PACK_REV_IF_UPPER_BIT         ( 0x1  << BLIS_PACK_REV_IF_UPPER_SHIFT )
#define BLIS_PACK_REV_IF_LOWER_BIT         ( 0x1  << BLIS_PACK_REV_IF_LOWER_SHIFT )
#define BLIS_PACK_BUFFER_BITS              ( 0x3  << BLIS_PACK_BUFFER_SHIFT )
#define BLIS_STRUC_BITS                    ( 0x3  << BLIS_STRUC_SHIFT )
#define BLIS_COMP_DT_BITS                  ( 0x7  << BLIS_COMP_DT_SHIFT )
#define   BLIS_COMP_DOMAIN_BIT             ( 0x1  << BLIS_COMP_DOMAIN_SHIFT )
#define   BLIS_COMP_PREC_BIT               ( 0x1  << BLIS_COMP_PREC_SHIFT )

// info2
#define BLIS_SCALAR_DT_BITS                ( 0x7  << BLIS_SCALAR_DT_SHIFT )
#define   BLIS_SCALAR_DOMAIN_BIT           ( 0x1  << BLIS_SCALAR_DOMAIN_SHIFT )
#define   BLIS_SCALAR_PREC_BIT             ( 0x1  << BLIS_SCALAR_PREC_SHIFT )


//
// -- BLIS enumerated type value definitions -----------------------------------
//

#define BLIS_BITVAL_REAL                      0x0
#define BLIS_BITVAL_COMPLEX                   BLIS_DOMAIN_BIT
#define BLIS_BITVAL_SINGLE_PREC               0x0
#define BLIS_BITVAL_DOUBLE_PREC               BLIS_PRECISION_BIT
#define   BLIS_BITVAL_FLOAT_TYPE              0x0
#define   BLIS_BITVAL_SCOMPLEX_TYPE           BLIS_DOMAIN_BIT  
#define   BLIS_BITVAL_DOUBLE_TYPE             BLIS_PRECISION_BIT
#define   BLIS_BITVAL_DCOMPLEX_TYPE         ( BLIS_DOMAIN_BIT | BLIS_PRECISION_BIT )
#define   BLIS_BITVAL_INT_TYPE                0x04
#define   BLIS_BITVAL_CONST_TYPE              0x05
#define BLIS_BITVAL_NO_TRANS                  0x0
#define BLIS_BITVAL_TRANS                     BLIS_TRANS_BIT
#define BLIS_BITVAL_NO_CONJ                   0x0
#define BLIS_BITVAL_CONJ                      BLIS_CONJ_BIT
#define BLIS_BITVAL_CONJ_TRANS              ( BLIS_CONJ_BIT | BLIS_TRANS_BIT )
#define BLIS_BITVAL_ZEROS                     0x0 
#define BLIS_BITVAL_UPPER                   ( BLIS_UPPER_BIT | BLIS_DIAG_BIT )
#define BLIS_BITVAL_LOWER                   ( BLIS_LOWER_BIT | BLIS_DIAG_BIT )
#define BLIS_BITVAL_DENSE                     BLIS_UPLO_BITS  
#define BLIS_BITVAL_NONUNIT_DIAG              0x0
#define BLIS_BITVAL_UNIT_DIAG                 BLIS_UNIT_DIAG_BIT
#define BLIS_BITVAL_INVERT_DIAG               BLIS_INVERT_DIAG_BIT
#define BLIS_BITVAL_NOT_PACKED                0x0
#define   BLIS_BITVAL_4MI                   ( 0x1  << BLIS_PACK_FORMAT_SHIFT )
#define   BLIS_BITVAL_3MI                   ( 0x2  << BLIS_PACK_FORMAT_SHIFT )
#define   BLIS_BITVAL_4MS                   ( 0x3  << BLIS_PACK_FORMAT_SHIFT )
#define   BLIS_BITVAL_3MS                   ( 0x4  << BLIS_PACK_FORMAT_SHIFT )
#define   BLIS_BITVAL_RO                    ( 0x5  << BLIS_PACK_FORMAT_SHIFT )
#define   BLIS_BITVAL_IO                    ( 0x6  << BLIS_PACK_FORMAT_SHIFT )
#define   BLIS_BITVAL_RPI                   ( 0x7  << BLIS_PACK_FORMAT_SHIFT )
#define   BLIS_BITVAL_1E                    ( 0x8  << BLIS_PACK_FORMAT_SHIFT )
#define   BLIS_BITVAL_1R                    ( 0x9  << BLIS_PACK_FORMAT_SHIFT )
#define   BLIS_BITVAL_PACKED_UNSPEC         ( BLIS_PACK_BIT                                                            )
#define   BLIS_BITVAL_PACKED_ROWS           ( BLIS_PACK_BIT                                                            )
#define   BLIS_BITVAL_PACKED_COLUMNS        ( BLIS_PACK_BIT                                         | BLIS_PACK_RC_BIT )
#define   BLIS_BITVAL_PACKED_ROW_PANELS     ( BLIS_PACK_BIT                   | BLIS_PACK_PANEL_BIT                    )
#define   BLIS_BITVAL_PACKED_COL_PANELS     ( BLIS_PACK_BIT                   | BLIS_PACK_PANEL_BIT | BLIS_PACK_RC_BIT )
#define   BLIS_BITVAL_PACKED_ROW_PANELS_4MI ( BLIS_PACK_BIT | BLIS_BITVAL_4MI | BLIS_PACK_PANEL_BIT                    )
#define   BLIS_BITVAL_PACKED_COL_PANELS_4MI ( BLIS_PACK_BIT | BLIS_BITVAL_4MI | BLIS_PACK_PANEL_BIT | BLIS_PACK_RC_BIT )
#define   BLIS_BITVAL_PACKED_ROW_PANELS_3MI ( BLIS_PACK_BIT | BLIS_BITVAL_3MI | BLIS_PACK_PANEL_BIT                    )
#define   BLIS_BITVAL_PACKED_COL_PANELS_3MI ( BLIS_PACK_BIT | BLIS_BITVAL_3MI | BLIS_PACK_PANEL_BIT | BLIS_PACK_RC_BIT )
#define   BLIS_BITVAL_PACKED_ROW_PANELS_4MS ( BLIS_PACK_BIT | BLIS_BITVAL_4MS | BLIS_PACK_PANEL_BIT                    )
#define   BLIS_BITVAL_PACKED_COL_PANELS_4MS ( BLIS_PACK_BIT | BLIS_BITVAL_4MS | BLIS_PACK_PANEL_BIT | BLIS_PACK_RC_BIT )
#define   BLIS_BITVAL_PACKED_ROW_PANELS_3MS ( BLIS_PACK_BIT | BLIS_BITVAL_3MS | BLIS_PACK_PANEL_BIT                    )
#define   BLIS_BITVAL_PACKED_COL_PANELS_3MS ( BLIS_PACK_BIT | BLIS_BITVAL_3MS | BLIS_PACK_PANEL_BIT | BLIS_PACK_RC_BIT )
#define   BLIS_BITVAL_PACKED_ROW_PANELS_RO  ( BLIS_PACK_BIT | BLIS_BITVAL_RO  | BLIS_PACK_PANEL_BIT                    )
#define   BLIS_BITVAL_PACKED_COL_PANELS_RO  ( BLIS_PACK_BIT | BLIS_BITVAL_RO  | BLIS_PACK_PANEL_BIT | BLIS_PACK_RC_BIT )
#define   BLIS_BITVAL_PACKED_ROW_PANELS_IO  ( BLIS_PACK_BIT | BLIS_BITVAL_IO  | BLIS_PACK_PANEL_BIT                    )
#define   BLIS_BITVAL_PACKED_COL_PANELS_IO  ( BLIS_PACK_BIT | BLIS_BITVAL_IO  | BLIS_PACK_PANEL_BIT | BLIS_PACK_RC_BIT )
#define   BLIS_BITVAL_PACKED_ROW_PANELS_RPI ( BLIS_PACK_BIT | BLIS_BITVAL_RPI | BLIS_PACK_PANEL_BIT                    )
#define   BLIS_BITVAL_PACKED_COL_PANELS_RPI ( BLIS_PACK_BIT | BLIS_BITVAL_RPI | BLIS_PACK_PANEL_BIT | BLIS_PACK_RC_BIT )
#define   BLIS_BITVAL_PACKED_ROW_PANELS_1E  ( BLIS_PACK_BIT | BLIS_BITVAL_1E  | BLIS_PACK_PANEL_BIT                    )
#define   BLIS_BITVAL_PACKED_COL_PANELS_1E  ( BLIS_PACK_BIT | BLIS_BITVAL_1E  | BLIS_PACK_PANEL_BIT | BLIS_PACK_RC_BIT )
#define   BLIS_BITVAL_PACKED_ROW_PANELS_1R  ( BLIS_PACK_BIT | BLIS_BITVAL_1R  | BLIS_PACK_PANEL_BIT                    )
#define   BLIS_BITVAL_PACKED_COL_PANELS_1R  ( BLIS_PACK_BIT | BLIS_BITVAL_1R  | BLIS_PACK_PANEL_BIT | BLIS_PACK_RC_BIT )
#define BLIS_BITVAL_PACK_FWD_IF_UPPER         0x0
#define BLIS_BITVAL_PACK_REV_IF_UPPER         BLIS_PACK_REV_IF_UPPER_BIT
#define BLIS_BITVAL_PACK_FWD_IF_LOWER         0x0
#define BLIS_BITVAL_PACK_REV_IF_LOWER         BLIS_PACK_REV_IF_LOWER_BIT
#define BLIS_BITVAL_BUFFER_FOR_A_BLOCK        0x0
#define BLIS_BITVAL_BUFFER_FOR_B_PANEL      ( 0x1 << BLIS_PACK_BUFFER_SHIFT )
#define BLIS_BITVAL_BUFFER_FOR_C_PANEL      ( 0x2 << BLIS_PACK_BUFFER_SHIFT )
#define BLIS_BITVAL_BUFFER_FOR_GEN_USE      ( 0x3 << BLIS_PACK_BUFFER_SHIFT )
#define BLIS_BITVAL_GENERAL                   0x0
#define BLIS_BITVAL_HERMITIAN               ( 0x1 << BLIS_STRUC_SHIFT )
#define BLIS_BITVAL_SYMMETRIC               ( 0x2 << BLIS_STRUC_SHIFT )
#define BLIS_BITVAL_TRIANGULAR              ( 0x3 << BLIS_STRUC_SHIFT )


//
// -- BLIS enumerated type definitions -----------------------------------------
//

// -- Operational parameter types --

typedef enum
{
	BLIS_NO_TRANSPOSE      = 0x0,
	BLIS_TRANSPOSE         = BLIS_BITVAL_TRANS,
	BLIS_CONJ_NO_TRANSPOSE = BLIS_BITVAL_CONJ,
	BLIS_CONJ_TRANSPOSE    = BLIS_BITVAL_CONJ_TRANS
} trans_t;

typedef enum
{
	BLIS_NO_CONJUGATE      = 0x0,
	BLIS_CONJUGATE         = BLIS_BITVAL_CONJ
} conj_t;

typedef enum
{
	BLIS_ZEROS             = BLIS_BITVAL_ZEROS,
	BLIS_LOWER             = BLIS_BITVAL_LOWER,
	BLIS_UPPER             = BLIS_BITVAL_UPPER,
	BLIS_DENSE             = BLIS_BITVAL_DENSE
} uplo_t;

typedef enum
{
	BLIS_LEFT              = 0x0,
	BLIS_RIGHT
} side_t;

typedef enum
{
	BLIS_NONUNIT_DIAG      = 0x0,
	BLIS_UNIT_DIAG         = BLIS_BITVAL_UNIT_DIAG
} diag_t;

typedef enum
{
	BLIS_NO_INVERT_DIAG    = 0x0,
	BLIS_INVERT_DIAG       = BLIS_BITVAL_INVERT_DIAG
} invdiag_t;

typedef enum
{
	BLIS_GENERAL           = BLIS_BITVAL_GENERAL,
	BLIS_HERMITIAN         = BLIS_BITVAL_HERMITIAN,
	BLIS_SYMMETRIC         = BLIS_BITVAL_SYMMETRIC,
	BLIS_TRIANGULAR        = BLIS_BITVAL_TRIANGULAR
} struc_t;


// -- Data type --

typedef enum
{
	BLIS_FLOAT             = BLIS_BITVAL_FLOAT_TYPE,
	BLIS_DOUBLE            = BLIS_BITVAL_DOUBLE_TYPE,
	BLIS_SCOMPLEX          = BLIS_BITVAL_SCOMPLEX_TYPE,
	BLIS_DCOMPLEX          = BLIS_BITVAL_DCOMPLEX_TYPE,
	BLIS_INT               = BLIS_BITVAL_INT_TYPE,
	BLIS_CONSTANT          = BLIS_BITVAL_CONST_TYPE,
	BLIS_DT_LO             = BLIS_FLOAT,
	BLIS_DT_HI             = BLIS_DCOMPLEX
} num_t;

typedef enum
{
	BLIS_REAL              = BLIS_BITVAL_REAL,
	BLIS_COMPLEX           = BLIS_BITVAL_COMPLEX
} dom_t;

typedef enum
{
	BLIS_SINGLE_PREC       = BLIS_BITVAL_SINGLE_PREC,
	BLIS_DOUBLE_PREC       = BLIS_BITVAL_DOUBLE_PREC
} prec_t;

#endif
// end bli_type_defs.h

// End extern "C" construct block.
#ifdef __cplusplus
}
#endif

#endif


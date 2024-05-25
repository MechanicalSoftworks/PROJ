/******************************************************************************
 * Project:  PROJ.4
 * Purpose:  Internal plumbing for the PROJ.4 library.
 *
 * Author:   Lucas Zadrozny
 *
 ******************************************************************************
 * Copyright (c) 2016, 2017, Thomas Knudsen / SDFE
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, including without limitation
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,
 * and/or sell copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO COORD SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *****************************************************************************/

#ifndef PROJ_INTERNAL_SHARED_FUNCTION_IDS_H
#define PROJ_INTERNAL_SHARED_FUNCTION_IDS_H

 /*
  * Generate all function IDs.
  * These are consistent among all files that include this header.
  */
enum { PJ_FUNCTION_BASE_ID = __COUNTER__ + 1 };

// Function IDs start at 1. 0 is reserved for invalid.
#define DEFINE_FUNCTION_ID(name)	name##_id = __COUNTER__ - PJ_FUNCTION_BASE_ID + 1

#define DEFINE_COROUTINE_FUNCTION_ID(name) DEFINE_FUNCTION_ID(name)

/******************************************************************************
 * Coroutine IDs.
 *****************************************************************************/
#define PROJ_COROUTINE(name) DEFINE_FUNCTION_ID(name),
#define PROJ_FWD_2D(name)
#define PROJ_INV_2D(name)
#define PROJ_FWD_3D(name)
#define PROJ_INV_3D(name)
#define PROJ_OPERATOR(name)

typedef enum
{
	PJ_INVALID_COROUTINE,

#	include "pj_function_list_shared.h"
#	include "pj_function_list_host.h"

	PJ_COROUTINE_COUNT
} PJ_COROUTINE_ID;

#undef PROJ_COROUTINE
#undef PROJ_FWD_2D
#undef PROJ_INV_2D
#undef PROJ_FWD_3D
#undef PROJ_INV_3D
#undef PROJ_OPERATOR

/******************************************************************************
 * FWD_2D IDs.
 *****************************************************************************/
#define PROJ_COROUTINE(name)
#define PROJ_FWD_2D(name) DEFINE_FUNCTION_ID(name),
#define PROJ_INV_2D(name)
#define PROJ_FWD_3D(name)
#define PROJ_INV_3D(name)
#define PROJ_OPERATOR(name)

typedef enum
{
	PJ_INVALID_FWD_2D,

#	include "pj_function_list_shared.h"
#	include "pj_function_list_host.h"

	PJ_FWD_2D_COUNT
} PJ_FWD_2D_ID;

#undef PROJ_COROUTINE
#undef PROJ_FWD_2D
#undef PROJ_INV_2D
#undef PROJ_FWD_3D
#undef PROJ_INV_3D
#undef PROJ_OPERATOR

/******************************************************************************
 * INV_2D IDs.
 *****************************************************************************/
#define PROJ_COROUTINE(name)
#define PROJ_FWD_2D(name)
#define PROJ_INV_2D(name) DEFINE_FUNCTION_ID(name),
#define PROJ_FWD_3D(name)
#define PROJ_INV_3D(name)
#define PROJ_OPERATOR(name)

typedef enum
{
	PJ_INVALID_INV_2D,

#	include "pj_function_list_shared.h"
#	include "pj_function_list_host.h"

	PJ_INV_2D_COUNT
} PJ_INV_2D_ID;

#undef PROJ_COROUTINE
#undef PROJ_FWD_2D
#undef PROJ_INV_2D
#undef PROJ_FWD_3D
#undef PROJ_INV_3D
#undef PROJ_OPERATOR

/******************************************************************************
 * FWD_3D IDs.
 *****************************************************************************/
#define PROJ_COROUTINE(name)
#define PROJ_FWD_2D(name)
#define PROJ_INV_2D(name)
#define PROJ_FWD_3D(name) DEFINE_FUNCTION_ID(name),
#define PROJ_INV_3D(name)
#define PROJ_OPERATOR(name)

typedef enum
{
	PJ_INVALID_FWD_3D,

#	include "pj_function_list_shared.h"
#	include "pj_function_list_host.h"

	PJ_FWD_3D_COUNT
} PJ_FWD_3D_ID;

#undef PROJ_COROUTINE
#undef PROJ_FWD_2D
#undef PROJ_INV_2D
#undef PROJ_FWD_3D
#undef PROJ_INV_3D
#undef PROJ_OPERATOR

/******************************************************************************
 * INV_3D IDs.
 *****************************************************************************/
#define PROJ_COROUTINE(name)
#define PROJ_FWD_2D(name)
#define PROJ_INV_2D(name)
#define PROJ_FWD_3D(name)
#define PROJ_INV_3D(name) DEFINE_FUNCTION_ID(name),
#define PROJ_OPERATOR(name)

typedef enum
{
	PJ_INVALID_INV_3D,

#	include "pj_function_list_shared.h"
#	include "pj_function_list_host.h"

	PJ_INV_3D_COUNT
} PJ_INV_3D_ID;

#undef PROJ_COROUTINE
#undef PROJ_FWD_2D
#undef PROJ_INV_2D
#undef PROJ_FWD_3D
#undef PROJ_INV_3D
#undef PROJ_OPERATOR

/******************************************************************************
 * 4D IDs.
 *****************************************************************************/
#define PROJ_COROUTINE(name)
#define PROJ_FWD_2D(name)
#define PROJ_INV_2D(name)
#define PROJ_FWD_3D(name)
#define PROJ_INV_3D(name)
#define PROJ_OPERATOR(name) DEFINE_FUNCTION_ID(name),

typedef enum
{
	PJ_INVALID_OPERATOR,

#	include "pj_function_list_shared.h"
#	include "pj_function_list_host.h"

	PJ_OPERATOR_COUNT
} PJ_OPERATOR_ID;

#undef PROJ_COROUTINE
#undef PROJ_FWD_2D
#undef PROJ_INV_2D
#undef PROJ_FWD_3D
#undef PROJ_INV_3D
#undef PROJ_OPERATOR

#endif // !PROJ_INTERNAL_SHARED_FUNCTION_IDS_H

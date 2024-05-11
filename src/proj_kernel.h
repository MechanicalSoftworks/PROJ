/******************************************************************************
 * Project:  PROJ.4
 * Purpose:  Internal plumbing for the PROJ.4 library.
 *
 * Author:   Thomas Knudsen, <thokn@sdfe.dk>
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

#ifndef PROJ_KERNEL_H
#define PROJ_KERNEL_H

#include "proj_shared.h"

// The constructors need access to the real header.
#ifdef PROJ_OPENCL_DEVICE
#	include "proj_internal_device.h"
#else
#   include "proj_internal.h"
#endif

/*
 * Generate all function IDs.
 * These are consistent among all files that include this header.
 */
enum { PJ_FUNCTION_BASE_ID = __COUNTER__ };
#define DEFINE_FUNCTION_ID(name)	static constexpr int name##_id = __COUNTER__ - PJ_FUNCTION_BASE_ID

#define DEFINE_COROUTINE_FUNCTION_ID(name) DEFINE_FUNCTION_ID(name);
#define DEFINE_FWD_2D_FUNCTION_ID(name) DEFINE_FUNCTION_ID(name)
#define DEFINE_INV_2D_FUNCTION_ID(name) DEFINE_FUNCTION_ID(name)
#define DEFINE_FWD_3D_FUNCTION_ID(name) DEFINE_FUNCTION_ID(name)
#define DEFINE_INV_3D_FUNCTION_ID(name) DEFINE_FUNCTION_ID(name)
#define DEFINE_FWD_4D_FUNCTION_ID(name) DEFINE_FUNCTION_ID(name)
#define DEFINE_INV_4D_FUNCTION_ID(name) DEFINE_FUNCTION_ID(name)

#define PROJ_COROUTINE(name) DEFINE_COROUTINE_FUNCTION_ID(name);
#define PROJ_FWD_2D(name)    DEFINE_FWD_2D_FUNCTION_ID(name);
#define PROJ_INV_2D(name)    DEFINE_INV_2D_FUNCTION_ID(name);
#define PROJ_FWD_3D(name)    DEFINE_FWD_3D_FUNCTION_ID(name);
#define PROJ_INV_3D(name)    DEFINE_INV_3D_FUNCTION_ID(name);
#define PROJ_FWD_4D(name)    DEFINE_FWD_4D_FUNCTION_ID(name);
#define PROJ_INV_4D(name)    DEFINE_INV_4D_FUNCTION_ID(name);
#include "pj_function_list.h"
#undef PROJ_COROUTINE
#undef PROJ_FWD_2D
#undef PROJ_INV_2D
#undef PROJ_FWD_3D
#undef PROJ_INV_3D
#undef PROJ_FWD_4D
#undef PROJ_INV_4D

#endif // !PROJ_KERNEL_H

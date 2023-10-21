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

#ifndef PROJ_INTERNAL_DEVICE_H
#define PROJ_INTERNAL_DEVICE_H

#include "proj_internal_shared.h"

#define proj_log_error(P, fmt, ...)

PJ_COORD proj_coord_error(void);
int proj_errno_set(const PJ* P, int err);

static inline int pj_streq(const char* a, __constant char* b)
{
	while (*a && *b)
	{
		if (*a != *b)
		{
			return 0;
		}
		
		++a;
		++b;
	}
	return 1;
}

void stack_push(cl_local PJstack_t* stack, PJ_COROUTINE_ID fn, PJ* P, PJ_COORD coo);

PJ_COORD proj_trans(cl_local PJstack_t* stack, PJ* P, PJ_DIRECTION direction, PJ_COORD coord);

PJcoroutine_code_t proj_dispatch_coroutine(PJ_COROUTINE_ID fn, cl_local struct PJstack_s* stack, cl_local struct PJstack_entry_s* e);
PJ_XY proj_dispatch_fwd(PJ_FWD_2D_ID fn, PJ_LP lp, PJ* P);
PJ_LP proj_dispatch_inv(PJ_INV_2D_ID fn, PJ_XY xy, PJ* P);
PJ_XYZ proj_dispatch_fwd3d(PJ_FWD_3D_ID fn, PJ_LPZ lpz, PJ* P);
PJ_LPZ proj_dispatch_inv3d(PJ_INV_3D_ID fn, PJ_XYZ xyz, PJ* P);
PJ_COORD proj_dispatch_fwd4d(PJ_FWD_4D_ID fn, PJ_COORD coo, PJ* P);
PJ_COORD proj_dispatch_inv4d(PJ_INV_4D_ID fn, PJ_COORD coo, PJ* P);

#endif // !PROJ_INTERNAL_DEVICE_H

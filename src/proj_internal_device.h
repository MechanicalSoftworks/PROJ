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

PJ_COORD proj_trans(PJ* P, PJ_DIRECTION direction, PJ_COORD coord);
PJ_COORD pj_approx_2D_trans(PJ* P, PJ_DIRECTION direction, PJ_COORD coo);
PJ_COORD pj_approx_3D_trans(PJ* P, PJ_DIRECTION direction, PJ_COORD coo);


#endif // !PROJ_INTERNAL_DEVICE_H

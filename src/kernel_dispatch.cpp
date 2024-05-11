/******************************************************************************
 * Project:  PROJ.4
 * Purpose:  Builds the OpenCL jump table.
 * Author:   Lucas Zadrozny
 *
 ******************************************************************************
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
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *****************************************************************************/
#include "proj_internal_shared.h"
#include "proj_kernel.h"

/*****************************************************************************/
/*****************************************************************************/
#define PROJ_COROUTINE(name) extern PJcoroutine_code_t name(cl_local PJstack_t* stack, cl_local PJstack_entry_t* e);
#define PROJ_FWD_2D(name)    extern PJ_XY name(PJ_LP lp, PJ *P);
#define PROJ_INV_2D(name)    extern PJ_LP name(PJ_XY xy, PJ *P);
#define PROJ_FWD_3D(name)    extern PJ_XYZ name(PJ_LPZ lp, PJ *P);
#define PROJ_INV_3D(name)    extern PJ_LPZ name(PJ_XYZ xy, PJ *P);
#define PROJ_FWD_4D(name)    extern PJ_COORD name(PJ_COORD lp, PJ *P);
#define PROJ_INV_4D(name)    extern PJ_COORD name(PJ_COORD xy, PJ *P);

#include "pj_function_list.h"

#undef PROJ_COROUTINE
#undef PROJ_FWD_2D
#undef PROJ_INV_2D
#undef PROJ_FWD_3D
#undef PROJ_INV_3D
#undef PROJ_FWD_4D
#undef PROJ_INV_4D

/*****************************************************************************/
/*****************************************************************************/
#define PROJ_COROUTINE(name) case name##_id: return name(stack, e);
#define PROJ_FWD_2D(name)
#define PROJ_INV_2D(name)
#define PROJ_FWD_3D(name)
#define PROJ_INV_3D(name)
#define PROJ_FWD_4D(name)
#define PROJ_INV_4D(name)
PJcoroutine_code_t proj_dispatch_coroutine(PJ_COROUTINE_ID fn, cl_local struct PJstack_s* stack, cl_local struct PJstack_entry_s* e)
{
    switch (fn)
    {
#       include "pj_function_list.h"
    }

    return PJ_CO_ERROR;
}
#undef PROJ_COROUTINE
#undef PROJ_FWD_2D
#undef PROJ_INV_2D
#undef PROJ_FWD_3D
#undef PROJ_INV_3D
#undef PROJ_FWD_4D
#undef PROJ_INV_4D

/*****************************************************************************/
/*****************************************************************************/
#define PROJ_COROUTINE(name)
#define PROJ_FWD_2D(name)       case name##_id: return name(in, P);
#define PROJ_INV_2D(name)
#define PROJ_FWD_3D(name)
#define PROJ_INV_3D(name)
#define PROJ_FWD_4D(name)
#define PROJ_INV_4D(name)
PJ_XY proj_dispatch_fwd(PJ_FWD_2D_ID fn, PJ_LP in, PJ *P)
{
    switch (fn)
    {
#       include "pj_function_list.h"
    }

    PJ_XY x{};
    return x;
}
#undef PROJ_COROUTINE
#undef PROJ_FWD_2D
#undef PROJ_INV_2D
#undef PROJ_FWD_3D
#undef PROJ_INV_3D
#undef PROJ_FWD_4D
#undef PROJ_INV_4D

/*****************************************************************************/
/*****************************************************************************/
#define PROJ_COROUTINE(name)
#define PROJ_FWD_2D(name)
#define PROJ_INV_2D(name)       case name##_id: return name(in, P);
#define PROJ_FWD_3D(name)
#define PROJ_INV_3D(name)
#define PROJ_FWD_4D(name)
#define PROJ_INV_4D(name)
PJ_LP proj_dispatch_inv(PJ_INV_2D_ID fn, PJ_XY in, PJ *P)
{
    switch (fn)
    {
#       include "pj_function_list.h"
    }

    PJ_LP x{};
    return x;
}
#undef PROJ_COROUTINE
#undef PROJ_FWD_2D
#undef PROJ_INV_2D
#undef PROJ_FWD_3D
#undef PROJ_INV_3D
#undef PROJ_FWD_4D
#undef PROJ_INV_4D


/*****************************************************************************/
/*****************************************************************************/
#define PROJ_COROUTINE(name)
#define PROJ_FWD_2D(name)
#define PROJ_INV_2D(name)
#define PROJ_FWD_3D(name)       case name##_id: return name(in, P);
#define PROJ_INV_3D(name)
#define PROJ_FWD_4D(name)
#define PROJ_INV_4D(name)
PJ_XYZ proj_dispatch_fwd3d(PJ_FWD_3D_ID fn, PJ_LPZ in, PJ *P)
{
    switch (fn)
    {
#       include "pj_function_list.h"
    }

    PJ_XYZ x{};
    return x;
}
#undef PROJ_COROUTINE
#undef PROJ_FWD_2D
#undef PROJ_INV_2D
#undef PROJ_FWD_3D
#undef PROJ_INV_3D
#undef PROJ_FWD_4D
#undef PROJ_INV_4D

/*****************************************************************************/
/*****************************************************************************/
#define PROJ_COROUTINE(name)
#define PROJ_FWD_2D(name)
#define PROJ_INV_2D(name)
#define PROJ_FWD_3D(name)
#define PROJ_INV_3D(name)       case name##_id: return name(in, P);
#define PROJ_FWD_4D(name)
#define PROJ_INV_4D(name)
PJ_LPZ proj_dispatch_inv3d(PJ_INV_3D_ID fn, PJ_XYZ in, PJ *P)
{
    switch (fn)
    {
#       include "pj_function_list.h"
    }

    PJ_LPZ x{};
    return x;
}
#undef PROJ_COROUTINE
#undef PROJ_FWD_2D
#undef PROJ_INV_2D
#undef PROJ_FWD_3D
#undef PROJ_INV_3D
#undef PROJ_FWD_4D
#undef PROJ_INV_4D


/*****************************************************************************/
/*****************************************************************************/
#define PROJ_COROUTINE(name)
#define PROJ_FWD_2D(name)
#define PROJ_INV_2D(name)
#define PROJ_FWD_3D(name)
#define PROJ_INV_3D(name)
#define PROJ_FWD_4D(name)       case name##_id: return name(in, P);
#define PROJ_INV_4D(name)
PJ_COORD proj_dispatch_fwd4d(PJ_FWD_4D_ID fn, PJ_COORD in, PJ *P)
{
    switch (fn)
    {
#       include "pj_function_list.h"
    }

    PJ_COORD x{};
    return x;
}
#undef PROJ_COROUTINE
#undef PROJ_FWD_2D
#undef PROJ_INV_2D
#undef PROJ_FWD_3D
#undef PROJ_INV_3D
#undef PROJ_FWD_4D
#undef PROJ_INV_4D

/*****************************************************************************/
/*****************************************************************************/
#define PROJ_COROUTINE(name)
#define PROJ_FWD_4D(name)
#define PROJ_INV_4D(name)       case name##_id: return name(in, P);
#define PROJ_FWD_2D(name)
#define PROJ_INV_2D(name)
#define PROJ_FWD_3D(name)
#define PROJ_INV_3D(name)
PJ_COORD proj_dispatch_inv4d(PJ_INV_4D_ID fn, PJ_COORD in, PJ *P)
{
    switch (fn)
    {
#       include "pj_function_list.h"
    }

    PJ_COORD x{};
    return x;
}
#undef PROJ_COROUTINE
#undef PROJ_FWD_2D
#undef PROJ_INV_2D
#undef PROJ_FWD_3D
#undef PROJ_INV_3D
#undef PROJ_FWD_4D
#undef PROJ_INV_4D
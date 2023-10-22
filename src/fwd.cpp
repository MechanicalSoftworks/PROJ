/******************************************************************************
 * Project:  PROJ.4
 * Purpose:  Forward operation invocation
 * Author:   Thomas Knudsen,  thokn@sdfe.dk,  2018-01-02
 *           Based on material from Gerald Evenden (original pj_fwd)
 *           and Piyush Agram (original pj_fwd3d)
 *
 ******************************************************************************
 * Copyright (c) 2000, Frank Warmerdam
 * Copyright (c) 2018, Thomas Knudsen / SDFE
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
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *****************************************************************************/

#include "proj_kernel.h"

#undef INPUT_UNITS
#undef OUTPUT_UNITS
#define INPUT_UNITS  P->left
#define OUTPUT_UNITS P->right


static PJcoroutine_code_t fwd_prepare_co (cl_local PJstack_t* stack, cl_local PJstack_entry_t* e) {
    PJ*         P = e->P;
    PJ_COORD    coo = e->coo;

    switch (e->step) {
        case 0: break;
        case 1: goto p1;
        case 2: goto p2;
        case 3: goto p3;
        case 4: goto p4;
        case 5: goto p5;
        case 6: goto p6;
        default: goto ABORT;
    }

    if (HUGE_VAL==coo.v[0] || HUGE_VAL==coo.v[1] || HUGE_VAL==coo.v[2])
    {
        goto ABORT;
    }

    /* The helmert datum shift will choke unless it gets a sensible 4D coordinate */
    if (HUGE_VAL==coo.v[2] && P->helmert) coo.v[2] = 0.0;
    if (HUGE_VAL==coo.v[3] && P->helmert) coo.v[3] = 0.0;

    /* Check validity of angular input coordinates */
    if (INPUT_UNITS==PJ_IO_UNITS_RADIANS) {
        double t;

        /* check for latitude or longitude over-range */
        t = (coo.lp.phi < 0  ?  -coo.lp.phi  :  coo.lp.phi) - M_HALFPI;
        if (t > PJ_EPS_LAT)
        {
            proj_log_error(P, _("Invalid latitude"));
            proj_errno_set (P, PROJ_ERR_COORD_TRANSFM_INVALID_COORD);
            goto ABORT;
        }
        if (coo.lp.lam > 10  ||  coo.lp.lam < -10)
        {
            proj_log_error(P, _("Invalid longitude"));
            proj_errno_set (P, PROJ_ERR_COORD_TRANSFM_INVALID_COORD);
            goto ABORT;
        }


        /* Clamp latitude to -90..90 degree range */
        if (coo.lp.phi > M_HALFPI)
            coo.lp.phi = M_HALFPI;
        if (coo.lp.phi < -M_HALFPI)
            coo.lp.phi = -M_HALFPI;

        /* If input latitude is geocentrical, convert to geographical */
        if (P->geoc)
            coo = pj_geocentric_latitude (P, PJ_INV, coo);

        /* Ensure longitude is in the -pi:pi range */
        if (0==P->over)
            coo.lp.lam = adjlon(coo.lp.lam);

        if (P->hgridshift) {
            push_proj_trans(stack, P->hgridshift, PJ_INV, coo);
            PJ_YIELD(e, 1);
        }
        else if (P->helmert || (P->cart_wgs84 != nullptr && P->cart != nullptr)) {
            push_proj_trans (stack, P->cart_wgs84, PJ_FWD, coo); /* Go cartesian in WGS84 frame */
            PJ_YIELD(e, 2);

            if (P->helmert) {
                push_proj_trans(stack, P->helmert, PJ_INV, coo); /* Step into local frame */
                PJ_YIELD(e, 3);
            }
            push_proj_trans (stack, P->cart,       PJ_INV, coo); /* Go back to angular using local ellps */
            PJ_YIELD(e, 4);
        }
        if (coo.lp.lam==HUGE_VAL)
            goto DONE;
        if (P->vgridshift) {
            push_proj_trans(stack, P->vgridshift, PJ_FWD, coo); /* Go orthometric from geometric */
            PJ_YIELD(e, 5);
        }

        /* Distance from central meridian, taking system zero meridian into account */
        coo.lp.lam = (coo.lp.lam - P->from_greenwich) - P->lam0;

        /* Ensure longitude is in the -pi:pi range */
        if (0==P->over)
            coo.lp.lam = adjlon(coo.lp.lam);

        goto DONE;
    }


    /* We do not support gridshifts on cartesian input */
    if (INPUT_UNITS == PJ_IO_UNITS_CARTESIAN && P->helmert) {
        push_proj_trans(stack, P->helmert, PJ_INV, coo);
        PJ_YIELD(e, 6);
    }

DONE:
    e->coo = coo;
    return PJ_CO_DONE;

YIELD:
    e->coo = coo;
    return PJ_CO_YIELD;

ABORT:
    e->coo = proj_coord_error();
    return PJ_CO_ERROR;
}


static PJcoroutine_code_t fwd_finalize_co (cl_local PJstack_t* stack, cl_local PJstack_entry_t* e) {
    PJ*         P = e->P;
    PJ_COORD    coo = e->coo;

    switch (e->step) {
        case 0: break;
        case 1: goto p1;
        case 2: goto p2;
        default: goto ABORT;
    }

    switch (OUTPUT_UNITS) {

    /* Handle false eastings/northings and non-metric linear units */
    case PJ_IO_UNITS_CARTESIAN:

        if (P->is_geocent) {
            push_proj_trans (stack, P->cart, PJ_FWD, coo);
            PJ_YIELD(e, 1);
        }
        coo.xyz.x *= P->fr_meter;
        coo.xyz.y *= P->fr_meter;
        coo.xyz.z *= P->fr_meter;

        break;

    /* Classic proj.4 functions return plane coordinates in units of the semimajor axis */
    case PJ_IO_UNITS_CLASSIC:
        coo.xy.x *= P->a;
        coo.xy.y *= P->a;

    /* Falls through */ /* (<-- GCC warning silencer) */
    /* to continue processing in common with PJ_IO_UNITS_PROJECTED */
    case PJ_IO_UNITS_PROJECTED:
        coo.xyz.x = P->fr_meter  * (coo.xyz.x + P->x0);
        coo.xyz.y = P->fr_meter  * (coo.xyz.y + P->y0);
        coo.xyz.z = P->vfr_meter * (coo.xyz.z + P->z0);
        break;

    case PJ_IO_UNITS_WHATEVER:
        break;

    case PJ_IO_UNITS_DEGREES:
        break;

    case PJ_IO_UNITS_RADIANS:
        coo.lpz.z = P->vfr_meter * (coo.lpz.z + P->z0);

        if( P->is_long_wrap_set ) {
            if( coo.lpz.lam != HUGE_VAL ) {
                coo.lpz.lam  = P->long_wrap_center +
                               adjlon(coo.lpz.lam - P->long_wrap_center);
            }
        }

        break;
    }

    if (P->axisswap) {
        push_proj_trans(stack, P->axisswap, PJ_FWD, coo);
        PJ_YIELD(e, 2);
    }

//DONE:
    e->coo = coo;
    return PJ_CO_DONE;

YIELD:
    e->coo = coo;
    return PJ_CO_YIELD;

ABORT:
    e->coo = proj_coord_error();
    return PJ_CO_ERROR;
}

static void push_fwd_prepare(cl_local PJstack_t* stack, PJ* P, PJ_COORD coo)
{
    return stack_push(stack, PJ_FUNCTION_PTR(fwd_prepare_co), P, coo);
}

static void push_fwd_finalize(cl_local PJstack_t* stack, PJ* P, PJ_COORD coo)
{
    return stack_push(stack, PJ_FUNCTION_PTR(fwd_finalize_co), P, coo);
}

static PJ_COORD push_fwd4d(cl_local PJstack_t* stack, PJ* P, PJ_COORD coo)
{
    if (PJ_GET_COROUTINE(P, co_fwd4d)) {
        // This code path is really only taken by the pipeline.
        stack_push(stack, PJ_GET_COROUTINE(P, co_fwd4d), P, coo);
        return coo;
    }

    // Try to avoid dispatching a new coroutine if possible, because:
    //  a) the coroutine stack is limited in size, and
    //  b) refactoring 100 projections, that don't need to be coroutines, into coroutines, isn't very fun.
    return proj_dispatch_fwd4d(PJ_GET_COROUTINE(P, fwd4d), coo, P);
}

static PJ_COORD push_fwd3d(cl_local PJstack_t* stack, PJ* P, PJ_COORD coo)
{
    if (PJ_GET_COROUTINE(P, co_fwd3d)) {
        // This code path is really only taken by the pipeline.
        stack_push(stack, PJ_GET_COROUTINE(P, co_fwd3d), P, coo);
        return coo;
    }

    PJ_COORD r = { 0, 0, 0, 0 };

    // Try to avoid dispatching a new coroutine if possible, because:
    //  a) the coroutine stack is limited in size, and
    //  b) refactoring 100 projections, that don't need to be coroutines, into coroutines, isn't very fun.
    r.xyz = proj_dispatch_fwd3d(PJ_GET_COROUTINE(P, fwd3d), coo.lpz, P);
    return r;
}

static PJ_COORD push_fwd(cl_local PJstack_t* stack, PJ* P, PJ_COORD coo)
{
    if (PJ_GET_COROUTINE(P, co_fwd)) {
        // This code path is really only taken by the pipeline.
        stack_push(stack, PJ_GET_COROUTINE(P, co_fwd), P, coo);
        return coo;
    }

    PJ_COORD r = { 0, 0, 0, 0 };

    // Try to avoid dispatching a new coroutine if possible, because:
    //  a) the coroutine stack is limited in size, and
    //  b) refactoring 100 projections, that don't need to be coroutines, into coroutines, isn't very fun.
    r.xy = proj_dispatch_fwd(PJ_GET_COROUTINE(P, fwd), coo.lp, P);
    return r;
}

PJ_COORD error_or_coord(PJ *P, PJ_COORD coord, int last_errno) {
    if (P->shared_ctx->last_errno)
        return proj_coord_error();

    P->shared_ctx->last_errno = last_errno;

    return coord;
}

PJcoroutine_code_t pj_fwd_co(cl_local PJstack_t* stack, cl_local PJstack_entry_t* e) {
    int         last_errno = e->last_errno;
    PJ*         P = e->P;
    PJ_COORD    coo = e->coo;

    switch (e->step) {
        case 0: break;
        case 1: goto p1;
        case 2: goto p2;
        case 3: goto p3;
        case 4: goto p4;
        case 5: goto p5;
        default: goto ABORT;
    }

    coo.lp = e->coo.lp;
    last_errno = P->shared_ctx->last_errno;
    P->shared_ctx->last_errno = 0;

    if (!P->skip_fwd_prepare) {
        push_fwd_prepare(stack, P, coo);
        PJ_YIELD(e, 1);
    }
    if (HUGE_VAL==coo.v[0] || HUGE_VAL==coo.v[1])
        goto ABORT;

    /* Do the transformation, using the lowest dimensional transformer available */
    if (PJ_GET_COROUTINE(P, fwd) || PJ_GET_COROUTINE(P, co_fwd))
    {
        coo.xy = push_fwd(stack, P, coo).xy;
        PJ_YIELD(e, 2);
    }
    else if (PJ_GET_COROUTINE(P, fwd3d) || PJ_GET_COROUTINE(P, co_fwd3d))
    {
        coo.xyz = push_fwd3d(stack, P, coo).xyz;
        PJ_YIELD(e, 3);
    }
    else if (PJ_GET_COROUTINE(P, fwd4d) || PJ_GET_COROUTINE(P, co_fwd4d))
    {
        coo = push_fwd4d(stack, P, coo);
        PJ_YIELD(e, 4);
    }
    else {
        proj_errno_set (P, PROJ_ERR_OTHER_NO_INVERSE_OP);
        goto ABORT;
    }
    if (HUGE_VAL==coo.v[0])
        goto ABORT;

    if (!P->skip_fwd_finalize) {
        push_fwd_finalize(stack, P, coo);
        PJ_YIELD(e, 5);
    }

    coo.xy = error_or_coord(P, coo, last_errno).xy;

//DONE:
    e->coo = coo;
    e->last_errno = last_errno;
    return PJ_CO_DONE;

YIELD:
    e->coo = coo;
    e->last_errno = last_errno;
    return PJ_CO_YIELD;

ABORT:
    e->coo.xy = proj_coord_error().xy;
    return PJ_CO_ERROR;
}



PJcoroutine_code_t pj_fwd3d_co (cl_local PJstack_t *stack, cl_local PJstack_entry_t* e) {
    int         last_errno = e->last_errno;
    PJ*         P = e->P;
    PJ_COORD    coo = e->coo;

    switch (e->step) {
        case 0: break;
        case 1: goto p1;
        case 2: goto p2;
        case 3: goto p3;
        case 4: goto p4;
        case 5: goto p5;
        default: goto ABORT;
    }

    coo.lpz = e->coo.lpz;
    last_errno = P->shared_ctx->last_errno;
    P->shared_ctx->last_errno = 0;

    if (!P->skip_fwd_prepare) {
        push_fwd_prepare(stack, P, coo);
        PJ_YIELD(e, 1);
    }
    if (HUGE_VAL==coo.v[0])
        goto ABORT;

    /* Do the transformation, using the lowest dimensional transformer feasible */
    if (PJ_GET_COROUTINE(P, fwd3d) || PJ_GET_COROUTINE(P, co_fwd3d))
    {
        coo.xyz = push_fwd3d(stack, P, coo).xyz;
        PJ_YIELD(e, 2);
    }
    else if (PJ_GET_COROUTINE(P, fwd4d) || PJ_GET_COROUTINE(P, co_fwd4d))
    {
        coo = push_fwd4d(stack, P, coo);
        PJ_YIELD(e, 3);
    }
    else if (PJ_GET_COROUTINE(P, fwd) || PJ_GET_COROUTINE(P, co_fwd))
    {
        coo.xy = push_fwd(stack, P, coo).xy;
        PJ_YIELD(e, 4);
    }
    else {
        proj_errno_set (P, PROJ_ERR_OTHER_NO_INVERSE_OP);
        goto ABORT;
    }
    if (HUGE_VAL==coo.v[0])
        goto ABORT;

    if (!P->skip_fwd_finalize) {
        push_fwd_finalize(stack, P, coo);
        PJ_YIELD(e, 5);
    }

    coo.xyz = error_or_coord(P, coo, last_errno).xyz;

//DONE:
    e->coo = coo;
    e->last_errno = last_errno;
    return PJ_CO_DONE;

YIELD:
    e->coo = coo;
    e->last_errno = last_errno;
    return PJ_CO_YIELD;

ABORT:
    e->coo.xyz = proj_coord_error().xyz;
    return PJ_CO_ERROR;
}

PJcoroutine_code_t pj_fwd4d_co (cl_local PJstack_t *stack, cl_local PJstack_entry_t* e) {
    int         last_errno = e->last_errno;
    PJ*         P = e->P;
    PJ_COORD    coo = e->coo;

    switch (e->step) {
        case 0: break;
        case 1: goto p1;
        case 2: goto p2;
        case 3: goto p3;
        case 4: goto p4;
        case 5: goto p5;
        default: goto ABORT;
    }

    last_errno = P->shared_ctx->last_errno;
    P->shared_ctx->last_errno = 0;

    if (!P->skip_fwd_prepare) {
        push_fwd_prepare(stack, P, coo);
        PJ_YIELD(e, 1);
    }
    if (HUGE_VAL == coo.v[0])
        goto ABORT;

    /* Call the highest dimensional converter available */
    if (PJ_GET_COROUTINE(P, fwd4d) || PJ_GET_COROUTINE(P, co_fwd4d))
    {
        coo = push_fwd4d(stack, P, coo);
        PJ_YIELD(e, 2);
    }
    else if (PJ_GET_COROUTINE(P, fwd3d) || PJ_GET_COROUTINE(P, co_fwd3d))
    {
        coo = push_fwd3d(stack, P, coo);
        PJ_YIELD(e, 3);
    }
    else if (PJ_GET_COROUTINE(P, fwd) || PJ_GET_COROUTINE(P, co_fwd))
    {
        coo = push_fwd(stack, P, coo);
        PJ_YIELD(e, 4);
    }
    else {
        proj_errno_set (P, PROJ_ERR_OTHER_NO_INVERSE_OP);
        goto ABORT;
    }
    if (HUGE_VAL==coo.v[0])
        goto ABORT;

    if (!P->skip_fwd_finalize) {
        push_fwd_finalize(stack, P, coo);
        PJ_YIELD(e, 5);
    }

    coo = error_or_coord(P, coo, last_errno);

//DONE:
    e->coo = coo;
    e->last_errno = last_errno;
    return PJ_CO_DONE;

YIELD:
    e->coo = coo;
    e->last_errno = last_errno;
    return PJ_CO_YIELD;

ABORT:
    e->coo = proj_coord_error();
    return PJ_CO_ERROR;
}

#ifndef PROJ_OPENCL_DEVICE

PJ_COORD pj_fwd4d(PJ_COORD coo, PJ* P) {
    PJstack_t   stack;

    stack_new(&stack);
    stack_push(&stack, PJ_FUNCTION_PTR(pj_fwd4d_co), P, coo);

    return stack_exec(&stack);
}

PJ_XYZ pj_fwd3d(PJ_LPZ lpz, PJ* P) {
    PJstack_t   stack;
    PJ_COORD    coo = { 0 };
    coo.lpz = lpz;

    stack_new(&stack);
    stack_push(&stack, PJ_FUNCTION_PTR(pj_fwd3d_co), P, coo);

    return stack_exec(&stack).xyz;
}

PJ_XY PROJ_DLL pj_fwd(PJ_LP lp, PJ* P) {
    PJstack_t   stack;
    PJ_COORD    coo = { 0 };
    coo.lp = lp;

    stack_new(&stack);

    stack_push(&stack, PJ_FUNCTION_PTR(pj_fwd_co), P, coo);

    return stack_exec(&stack).xy;
}

void pj_scan_fwd(PJscan& s) {
    s.add_co(PJ_MAKE_KERNEL(fwd_prepare_co), __FILE__);
    s.add_co(PJ_MAKE_KERNEL(fwd_finalize_co), __FILE__);
    s.add_co(PJ_MAKE_KERNEL(pj_fwd_co), __FILE__);
    s.add_co(PJ_MAKE_KERNEL(pj_fwd3d_co), __FILE__);
    s.add_co(PJ_MAKE_KERNEL(pj_fwd4d_co), __FILE__);
}

#endif

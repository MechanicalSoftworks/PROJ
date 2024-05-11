/******************************************************************************
 * Project:  PROJ.4
 * Purpose:  Inverse operation invocation
 * Author:   Thomas Knudsen,  thokn@sdfe.dk,  2018-01-02
 *           Based on material from Gerald Evenden (original pj_inv)
 *           and Piyush Agram (original pj_inv3d)
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
#define INPUT_UNITS  P->right
#define OUTPUT_UNITS P->left

PJcoroutine_code_t inv_prepare_co (cl_local PJstack_t* stack, cl_local PJstack_entry_t* e) {
    PJ*         P = e->P;
    PJ_COORD    coo = e->coo;

    switch (e->step) {
        case 0: break;
        case 1: goto p1;
        case 2: goto p2;
        default: goto ABORT;
    }

    if (coo.v[0] == HUGE_VAL || coo.v[1] == HUGE_VAL || coo.v[2] == HUGE_VAL) {
        proj_errno_set (P, PROJ_ERR_COORD_TRANSFM_OUTSIDE_PROJECTION_DOMAIN);
        goto ABORT;
    }

    /* The helmert datum shift will choke unless it gets a sensible 4D coordinate */
    if (HUGE_VAL==coo.v[2] && P->helmert) coo.v[2] = 0.0;
    if (HUGE_VAL==coo.v[3] && P->helmert) coo.v[3] = 0.0;

    if (P->axisswap) {
        push_proj_trans (stack, P->axisswap, PJ_INV, coo);
        PJ_YIELD(e, 1);
    }

    /* Handle remaining possible input types */
    switch (INPUT_UNITS) {
    case PJ_IO_UNITS_WHATEVER:
        break;

    case PJ_IO_UNITS_DEGREES:
        break;

    /* de-scale and de-offset */
    case PJ_IO_UNITS_CARTESIAN:
        coo.xyz.x *= P->to_meter;
        coo.xyz.y *= P->to_meter;
        coo.xyz.z *= P->to_meter;
        if (P->is_geocent) {
            push_proj_trans (stack, P->cart, PJ_INV, coo);
            PJ_YIELD(e, 2);
        }
        break;

    case PJ_IO_UNITS_PROJECTED:
    case PJ_IO_UNITS_CLASSIC:
        coo.xyz.x = P->to_meter  * coo.xyz.x - P->x0;
        coo.xyz.y = P->to_meter  * coo.xyz.y - P->y0;
        coo.xyz.z = P->vto_meter * coo.xyz.z - P->z0;
        if (INPUT_UNITS==PJ_IO_UNITS_PROJECTED)
            goto DONE;

        /* Classic proj.4 functions expect plane coordinates in units of the semimajor axis  */
        /* Multiplying by ra, rather than dividing by a because the CalCOFI projection       */
        /* stomps on a and hence (apparently) depends on this to roundtrip correctly         */
        /* (CalCOFI avoids further scaling by stomping - but a better solution is possible)  */
        coo.xyz.x *= P->ra;
        coo.xyz.y *= P->ra;
        break;

    case PJ_IO_UNITS_RADIANS:
        coo.lpz.z = P->vto_meter * coo.lpz.z - P->z0;
        break;
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



PJcoroutine_code_t inv_finalize_co (cl_local PJstack_t* stack, cl_local PJstack_entry_t* e) {
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

    if (coo.xyz.x == HUGE_VAL) {
        proj_errno_set (P, PROJ_ERR_COORD_TRANSFM_OUTSIDE_PROJECTION_DOMAIN);
        goto ABORT;
    }

    if (OUTPUT_UNITS==PJ_IO_UNITS_RADIANS) {

        /* Distance from central meridian, taking system zero meridian into account */
        coo.lp.lam = coo.lp.lam + P->from_greenwich + P->lam0;

        /* adjust longitude to central meridian */
        if (0==P->over)
            coo.lpz.lam = adjlon(coo.lpz.lam);

        if (P->vgridshift) {
            push_proj_trans (stack, P->vgridshift, PJ_INV, coo); /* Go geometric from orthometric */
            PJ_YIELD(e, 1);
        }
        if (coo.lp.lam==HUGE_VAL)
            goto DONE;
        if (P->hgridshift) {
            push_proj_trans(stack, P->hgridshift, PJ_FWD, coo);
            PJ_YIELD(e, 2);
        }
        else if (P->helmert || (P->cart_wgs84 != nullptr && P->cart != nullptr)) {
            push_proj_trans(stack, P->cart,       PJ_FWD, coo); /* Go cartesian in local frame */
            PJ_YIELD(e, 3);
            if( P->helmert ) {
                push_proj_trans(stack, P->helmert,    PJ_FWD, coo); /* Step into WGS84 */
                PJ_YIELD(e, 4);
            }
            push_proj_trans(stack, P->cart_wgs84, PJ_INV, coo); /* Go back to angular using WGS84 ellps */
            PJ_YIELD(e, 5);
        }
        if (coo.lp.lam==HUGE_VAL)
            goto DONE;

        /* If input latitude was geocentrical, convert back to geocentrical */
        if (P->geoc)
            coo = pj_geocentric_latitude (P, PJ_FWD, coo);
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

extern PJ_COORD error_or_coord(PJ* P, PJ_COORD coord, int last_errno);

static void push_inv_prepare(cl_local PJstack_t* stack, PJ* P, PJ_COORD coo)
{
    return stack_push(stack, PJ_FUNCTION_PTR(inv_prepare_co), P, coo);
}

static void push_inv_finalize(cl_local PJstack_t* stack, PJ* P, PJ_COORD coo)
{
    return stack_push(stack, PJ_FUNCTION_PTR(inv_finalize_co), P, coo);
}

static PJ_COORD push_inv4d(cl_local PJstack_t* stack, PJ* P, PJ_COORD coo)
{
    if (PJ_GET_COROUTINE(P, co_inv4d)) {
        // This code path is really only taken by the pipeline.
        stack_push(stack, PJ_GET_COROUTINE(P, co_inv4d), P, coo);
        return coo;
    }

    // Try to avoid dispatching a new coroutine if possible, because:
    //  a) the coroutine stack is limited in size, and
    //  b) refactoring 100 projections, that don't need to be coroutines, into coroutines, isn't very fun.
    return proj_dispatch_inv4d(PJ_GET_COROUTINE(P, inv4d), coo, P);
}

static PJ_COORD push_inv3d(cl_local PJstack_t* stack, PJ* P, PJ_COORD coo)
{
    if (PJ_GET_COROUTINE(P, co_inv3d)) {
        // This code path is really only taken by the pipeline.
        stack_push(stack, PJ_GET_COROUTINE(P, co_inv3d), P, coo);
        return coo;
    }

    PJ_COORD r = { 0, 0, 0, 0 };

    // Try to avoid dispatching a new coroutine if possible, because:
    //  a) the coroutine stack is limited in size, and
    //  b) refactoring 100 projections, that don't need to be coroutines, into coroutines, isn't very fun.
    r.lpz = proj_dispatch_inv3d(PJ_GET_COROUTINE(P, inv3d), coo.xyz, P);
    return r;
}

static PJ_COORD push_inv(cl_local PJstack_t* stack, PJ* P, PJ_COORD coo)
{
    if (PJ_GET_COROUTINE(P, co_inv)) {
        // This code path is really only taken by the pipeline.
        stack_push(stack, PJ_GET_COROUTINE(P, co_inv), P, coo);
        return coo;
    }

    PJ_COORD r = { 0, 0, 0, 0 };

    // Try to avoid dispatching a new coroutine if possible, because:
    //  a) the coroutine stack is limited in size, and
    //  b) refactoring 100 projections, that don't need to be coroutines, into coroutines, isn't very fun.
    r.lp = proj_dispatch_inv(PJ_GET_COROUTINE(P, inv), coo.xy, P);

    return r;
}

PJcoroutine_code_t pj_inv_co(cl_local PJstack_t* stack, cl_local PJstack_entry_t* e) {
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

    coo.xy = e->coo.xy;
    last_errno = P->shared_ctx->last_errno;
    P->shared_ctx->last_errno = 0;

    if (!P->skip_inv_prepare) {
        push_inv_prepare(stack, P, coo);
        PJ_YIELD(e, 1);
    }
    if (HUGE_VAL==coo.v[0])
        goto ABORT;

    /* Do the transformation, using the lowest dimensional transformer available */
    if (PJ_GET_COROUTINE(P, inv) || PJ_GET_COROUTINE(P, co_inv))
    {
        coo.lp = push_inv(stack, P, coo).lp;
        PJ_YIELD(e, 2);
    }
    else if (PJ_GET_COROUTINE(P, inv3d) || PJ_GET_COROUTINE(P, co_inv3d))
    {
        coo.lpz = push_inv3d(stack, P, coo).lpz;
        PJ_YIELD(e, 3);
    }
    else if (PJ_GET_COROUTINE(P, inv4d) || PJ_GET_COROUTINE(P, co_inv4d))
    {
        coo = push_inv4d(stack, P, coo);
        PJ_YIELD(e, 4);
    }
    else {
        proj_errno_set (P, PROJ_ERR_OTHER_NO_INVERSE_OP);
        goto ABORT;
    }
    if (HUGE_VAL==coo.v[0])
        goto ABORT;

    if (!P->skip_inv_finalize) {
        push_inv_finalize (stack, P, coo);
        PJ_YIELD(e, 5);
    }

    coo.lp = error_or_coord(P, coo, last_errno).lp;

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



PJcoroutine_code_t pj_inv3d_co (cl_local PJstack_t* stack, cl_local PJstack_entry_t* e) {
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

    coo.xyz = e->coo.xyz;
    last_errno = P->shared_ctx->last_errno;
    P->shared_ctx->last_errno = 0;

    if (!P->skip_inv_prepare) {
         push_inv_prepare (stack, P, coo);
         PJ_YIELD(e, 1);
    }
    if (HUGE_VAL==coo.v[0])
        goto ABORT;

    /* Do the transformation, using the lowest dimensional transformer feasible */
    if (PJ_GET_COROUTINE(P, inv3d) || PJ_GET_COROUTINE(P, co_inv3d))
    {
        coo.lpz = push_inv3d(stack, P, coo).lpz;
        PJ_YIELD(e, 2);
    }
    else if (PJ_GET_COROUTINE(P, inv4d) || PJ_GET_COROUTINE(P, co_inv4d))
    {
        coo = push_inv4d(stack, P, coo);
        PJ_YIELD(e, 3);
    }
    else if (PJ_GET_COROUTINE(P, inv) || PJ_GET_COROUTINE(P, co_inv))
    {
        coo.lp = push_inv(stack, P, coo).lp;
        PJ_YIELD(e, 4);
    }
    else {
        proj_errno_set (P, PROJ_ERR_OTHER_NO_INVERSE_OP);
        goto ABORT;
    }
    if (HUGE_VAL==coo.v[0])
        goto ABORT;

    if (!P->skip_inv_finalize) {
        push_inv_finalize(stack, P, coo);
        PJ_YIELD(e, 5);
    }

    coo.lpz = error_or_coord(P, coo, last_errno).lpz;

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



PJcoroutine_code_t pj_inv4d_co (cl_local PJstack_t* stack, cl_local PJstack_entry_t* e) {
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

    if (!P->skip_inv_prepare) {
        push_inv_prepare(stack, P, coo);
        PJ_YIELD(e, 1);
    }
    if (HUGE_VAL==coo.v[0])
        goto ABORT;

    /* Call the highest dimensional converter available */
    if (PJ_GET_COROUTINE(P, inv4d) || PJ_GET_COROUTINE(P, co_inv4d))
    {
        coo = push_inv4d(stack, P, coo);
        PJ_YIELD(e, 2);
    }
    else if (PJ_GET_COROUTINE(P, inv3d) || PJ_GET_COROUTINE(P, co_inv3d))
    {
        coo.lpz = push_inv3d(stack, P, coo).lpz;
        PJ_YIELD(e, 3);
    }
    else if (PJ_GET_COROUTINE(P, inv) || PJ_GET_COROUTINE(P, co_inv))
    {
        coo.lp = push_inv(stack, P, coo).lp;
        PJ_YIELD(e, 4);
    }
    else {
        proj_errno_set (P, PROJ_ERR_OTHER_NO_INVERSE_OP);
        goto ABORT;
    }
    if (HUGE_VAL==coo.v[0])
        goto ABORT;

    if (!P->skip_inv_finalize) {
        push_inv_finalize(stack, P, coo);
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

PJ_COORD pj_inv4d(PJ_COORD coo, PJ* P) {
    PJstack_t   stack;

    stack_new(&stack);
    stack_push(&stack, PJ_FUNCTION_PTR(pj_inv4d_co), P, coo);

    return stack_exec(&stack);
}

PJ_LPZ pj_inv3d(PJ_XYZ xyz, PJ* P) {
    PJstack_t   stack;
    PJ_COORD    coo = { 0 };
    coo.xyz = xyz;

    stack_new(&stack);
    stack_push(&stack, PJ_FUNCTION_PTR(pj_inv4d_co), P, coo);

    return stack_exec(&stack).lpz;
}

PJ_LP PROJ_DLL pj_inv(PJ_XY xy, PJ* P) {
    PJstack_t   stack;
    PJ_COORD    coo = { 0 };
    coo.xy = xy;

    stack_new(&stack);

    stack_push(&stack, PJ_FUNCTION_PTR(pj_inv_co), P, coo);

    return stack_exec(&stack).lp;
}

void pj_scan_inv(PJscan& s) {
    s.add_co(PJ_MAKE_KERNEL(inv_prepare_co), __FILE__);
    s.add_co(PJ_MAKE_KERNEL(inv_finalize_co), __FILE__);
    s.add_co(PJ_MAKE_KERNEL(pj_inv_co), __FILE__);
    s.add_co(PJ_MAKE_KERNEL(pj_inv3d_co), __FILE__);
    s.add_co(PJ_MAKE_KERNEL(pj_inv4d_co), __FILE__);
}

#endif

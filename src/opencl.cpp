/******************************************************************************
 * Project:  PROJ.4
 * Purpose:  OpenCL memory management
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
#include "proj.h"
#include "proj_internal.h"

/*****************************************************************************/
void pj_map_svm_ptrs(PJ* P, bool map)
/*****************************************************************************/
{
    P->host->ctx->allocator->svm_map(P, map);
    P->host->ctx->allocator->svm_map(P->shared_ctx, map);
    P->host->ctx->allocator->svm_map(P->opaque, map);

    // Nice little circular reference here. Enforce the constraint that this operation needs to happen top-down.
    //if (P->parent) {
    //    P->parent->host->map_pj(P->parent, map);
    //}

    if (P->geod) {
        P->host->ctx->allocator->svm_map(P->geod, map);
    }

    if (P->opaque) {
        P->host->ctx->allocator->svm_map(P->opaque, map);
    }

    if (P->axisswap) {
        P->axisswap->host->map_pj(P->axisswap, map);
    }

    if (P->cart) {
        P->cart->host->map_pj(P->cart, map);
    }

    if (P->cart_wgs84) {
        P->cart_wgs84->host->map_pj(P->cart_wgs84, map);
    }

    if (P->helmert) {
        P->helmert->host->map_pj(P->helmert, map);
    }

    if (P->hgridshift) {
        P->hgridshift->host->map_pj(P->hgridshift, map);
    }

    if (P->vgridshift) {
        P->vgridshift->host->map_pj(P->vgridshift, map);
    }
}

void proj_host_acquire_svm(PJ* P) {
    P->host->map_pj(P, true);
}

void proj_host_release_svm(PJ* P) {
    P->host->map_pj(P, false);
}

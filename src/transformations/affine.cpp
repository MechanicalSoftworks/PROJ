/************************************************************************
* Copyright (c) 2018, Even Rouault <even.rouault at spatialys.com>
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
*
***********************************************************************/
#define PJ_LIB__

#include <errno.h>
#include <math.h>

#include "proj.h"
#include "proj_internal.h"

PROJ_HEAD(affine, "Affine transformation");
PROJ_HEAD(geogoffset, "Geographic Offset");

namespace { // anonymous namespace
struct pj_affine_coeffs {
    double s11;
    double s12;
    double s13;
    double s21;
    double s22;
    double s23;
    double s31;
    double s32;
    double s33;
    double tscale;
};
} // anonymous namespace

namespace { // anonymous namespace
struct pj_opaque_affine {
    double xoff;
    double yoff;
    double zoff;
    double toff;
    struct pj_affine_coeffs forward;
    struct pj_affine_coeffs reverse;
};
} // anonymous namespace


PJ_COORD affine_forward_4d(PJ_COORD obs, PJ *P) {
    PJ_COORD newObs;
    const struct pj_opaque_affine *Q = (const struct pj_opaque_affine *) P->opaque;
    const struct pj_affine_coeffs *C = &(Q->forward);
    newObs.xyzt.x = Q->xoff + C->s11 * obs.xyzt.x + C->s12 * obs.xyzt.y + C->s13 * obs.xyzt.z;
    newObs.xyzt.y = Q->yoff + C->s21 * obs.xyzt.x + C->s22 * obs.xyzt.y + C->s23 * obs.xyzt.z;
    newObs.xyzt.z = Q->zoff + C->s31 * obs.xyzt.x + C->s32 * obs.xyzt.y + C->s33 * obs.xyzt.z;
    newObs.xyzt.t = Q->toff + C->tscale * obs.xyzt.t;
    return newObs;
}

PJ_XYZ affine_forward_3d(PJ_LPZ lpz, PJ *P) {
    PJ_COORD point = {{0,0,0,0}};
    point.lpz = lpz;
    return affine_forward_4d(point, P).xyz;
}


PJ_XY affine_forward_2d(PJ_LP lp, PJ *P) {
    PJ_COORD point = {{0,0,0,0}};
    point.lp = lp;
    return affine_forward_4d(point, P).xy;
}


PJ_COORD affine_reverse_4d(PJ_COORD obs, PJ *P) {
    PJ_COORD newObs;
    const struct pj_opaque_affine *Q = (const struct pj_opaque_affine *) P->opaque;
    const struct pj_affine_coeffs *C = &(Q->reverse);
    obs.xyzt.x -= Q->xoff;
    obs.xyzt.y -= Q->yoff;
    obs.xyzt.z -= Q->zoff;
    newObs.xyzt.x = C->s11 * obs.xyzt.x + C->s12 * obs.xyzt.y + C->s13 * obs.xyzt.z;
    newObs.xyzt.y = C->s21 * obs.xyzt.x + C->s22 * obs.xyzt.y + C->s23 * obs.xyzt.z;
    newObs.xyzt.z = C->s31 * obs.xyzt.x + C->s32 * obs.xyzt.y + C->s33 * obs.xyzt.z;
    newObs.xyzt.t = C->tscale * (obs.xyzt.t - Q->toff);
    return newObs;
}

PJ_LPZ affine_reverse_3d(PJ_XYZ xyz, PJ *P) {
    PJ_COORD point = {{0,0,0,0}};
    point.xyz = xyz;
    return affine_reverse_4d(point, P).lpz;
}

PJ_LP affine_reverse_2d(PJ_XY xy, PJ *P) {
    PJ_COORD point = {{0,0,0,0}};
    point.xy = xy;
    return affine_reverse_4d(point, P).lp;
}

static struct pj_opaque_affine * initQ(PJ_CONTEXT *ctx) {
    struct pj_opaque_affine *Q = static_cast<struct pj_opaque_affine *>(ctx->allocator->svm_calloc(1, sizeof(struct pj_opaque_affine)));
    if (nullptr==Q)
        return nullptr;

    /* default values */
    Q->forward.s11 = 1.0;
    Q->forward.s22 = 1.0;
    Q->forward.s33 = 1.0;
    Q->forward.tscale = 1.0;

    Q->reverse.s11 = 1.0;
    Q->reverse.s22 = 1.0;
    Q->reverse.s33 = 1.0;
    Q->reverse.tscale = 1.0;

    return Q;
}

static void computeReverseParameters(PJ* P)
{
    struct pj_opaque_affine *Q = (struct pj_opaque_affine *) P->opaque;

    /* cf https://en.wikipedia.org/wiki/Invertible_matrix#Inversion_of_3_%C3%97_3_matrices */
    const double a = Q->forward.s11;
    const double b = Q->forward.s12;
    const double c = Q->forward.s13;
    const double d = Q->forward.s21;
    const double e = Q->forward.s22;
    const double f = Q->forward.s23;
    const double g = Q->forward.s31;
    const double h = Q->forward.s32;
    const double i = Q->forward.s33;
    const double A = e * i - f * h;
    const double B = -(d * i - f * g);
    const double C = (d * h - e * g);
    const double D = -(b * i - c * h);
    const double E = (a * i - c * g);
    const double F = -(a * h - b * g);
    const double G = b * f - c * e;
    const double H = -(a * f - c * d);
    const double I = a * e - b * d;
    const double det = a * A + b * B + c * C;
    if( det == 0.0 || Q->forward.tscale == 0.0 ) {
        if (proj_log_level(P->host->ctx, PJ_LOG_TELL) >= PJ_LOG_DEBUG) {
            proj_log_debug(P, "matrix non invertible");
        }
        P->inv4d = nullptr;
        P->inv3d = nullptr;
        P->inv = nullptr;
    } else {
        Q->reverse.s11 = A / det;
        Q->reverse.s12 = D / det;
        Q->reverse.s13 = G / det;
        Q->reverse.s21 = B / det;
        Q->reverse.s22 = E / det;
        Q->reverse.s23 = H / det;
        Q->reverse.s31 = C / det;
        Q->reverse.s32 = F / det;
        Q->reverse.s33 = I / det;
        Q->reverse.tscale = 1.0 / Q->forward.tscale;
    }
}

PJ *TRANSFORMATION(affine,0 /* no need for ellipsoid */) {
    struct pj_opaque_affine *Q = initQ(P->host->ctx);
    if (nullptr==Q)
        return pj_default_destructor(P, PROJ_ERR_OTHER /*ENOMEM*/);
    P->opaque = (void *) Q;

    P->fwd4d = PJ_MAKE_KERNEL(affine_forward_4d);
    P->inv4d = PJ_MAKE_KERNEL(affine_reverse_4d);
    P->fwd3d  = PJ_MAKE_KERNEL(affine_forward_3d);
    P->inv3d  = PJ_MAKE_KERNEL(affine_reverse_3d);
    P->fwd    = PJ_MAKE_KERNEL(affine_forward_2d);
    P->inv    = PJ_MAKE_KERNEL(affine_reverse_2d);

    P->left   = PJ_IO_UNITS_WHATEVER;
    P->right  = PJ_IO_UNITS_WHATEVER;

    /* read args */
    Q->xoff = pj_param(P->host->ctx, P->host->params, "dxoff").f;
    Q->yoff = pj_param(P->host->ctx, P->host->params, "dyoff").f;
    Q->zoff = pj_param(P->host->ctx, P->host->params, "dzoff").f;
    Q->toff = pj_param(P->host->ctx, P->host->params, "dtoff").f;

    if(pj_param (P->host->ctx, P->host->params, "ts11").i) {
        Q->forward.s11 = pj_param(P->host->ctx, P->host->params, "ds11").f;
    }
    Q->forward.s12 = pj_param(P->host->ctx, P->host->params, "ds12").f;
    Q->forward.s13 = pj_param(P->host->ctx, P->host->params, "ds13").f;
    Q->forward.s21 = pj_param(P->host->ctx, P->host->params, "ds21").f;
    if(pj_param (P->host->ctx, P->host->params, "ts22").i) {
        Q->forward.s22 = pj_param(P->host->ctx, P->host->params, "ds22").f;
    }
    Q->forward.s23 = pj_param(P->host->ctx, P->host->params, "ds23").f;
    Q->forward.s31 = pj_param(P->host->ctx, P->host->params, "ds31").f;
    Q->forward.s32 = pj_param(P->host->ctx, P->host->params, "ds32").f;
    if(pj_param (P->host->ctx, P->host->params, "ts33").i) {
        Q->forward.s33 = pj_param(P->host->ctx, P->host->params, "ds33").f;
    }
    if(pj_param (P->host->ctx, P->host->params, "ttscale").i) {
        Q->forward.tscale = pj_param(P->host->ctx, P->host->params, "dtscale").f;
    }

    computeReverseParameters(P);

    return P;
}


/* Arcsecond to radians */
#define ARCSEC_TO_RAD (DEG_TO_RAD / 3600.0)


PJ *TRANSFORMATION(geogoffset,0 /* no need for ellipsoid */) {
    struct pj_opaque_affine *Q = initQ(P->host->ctx);
    if (nullptr==Q)
        return pj_default_destructor(P, PROJ_ERR_OTHER /*ENOMEM*/);
    P->opaque = (void *) Q;

    P->fwd4d = PJ_MAKE_KERNEL(affine_forward_4d);
    P->inv4d = PJ_MAKE_KERNEL(affine_reverse_4d);
    P->fwd3d  = PJ_MAKE_KERNEL(affine_forward_3d);
    P->inv3d  = PJ_MAKE_KERNEL(affine_reverse_3d);
    P->fwd    = PJ_MAKE_KERNEL(affine_forward_2d);
    P->inv    = PJ_MAKE_KERNEL(affine_reverse_2d);

    P->left   = PJ_IO_UNITS_RADIANS;
    P->right  = PJ_IO_UNITS_RADIANS;

    /* read args */
    Q->xoff = pj_param(P->host->ctx, P->host->params, "ddlon").f * ARCSEC_TO_RAD;
    Q->yoff = pj_param(P->host->ctx, P->host->params, "ddlat").f * ARCSEC_TO_RAD;
    Q->zoff = pj_param(P->host->ctx, P->host->params, "ddh").f;

    return P;
}

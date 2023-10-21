#define PJ_LIB__

#include "../proj_kernel.h"

struct pj_opaque_eqc {
    double rc;
};

PROJ_HEAD(eqc, "Equidistant Cylindrical (Plate Carree)")
    "\n\tCyl, Sph\n\tlat_ts=[, lat_0=0]";


static PJ_XY eqc_s_forward (PJ_LP lp, PJ *P) {           /* Spheroidal, forward */
    PJ_XY xy = {0.0,0.0};
    struct pj_opaque_eqc *Q = (struct pj_opaque_eqc*)(P->opaque);

    xy.x = Q->rc * lp.lam;
    xy.y = lp.phi - P->phi0;

    return xy;
}


static PJ_LP eqc_s_inverse (PJ_XY xy, PJ *P) {           /* Spheroidal, inverse */
    PJ_LP lp = {0.0,0.0};
    struct pj_opaque_eqc *Q = (struct pj_opaque_eqc*)(P->opaque);

    lp.lam = xy.x / Q->rc;
    lp.phi = xy.y + P->phi0;

    return lp;
}

#ifndef PROJ_OPENCL_DEVICE

PJ *PROJECTION(eqc) {
    struct pj_opaque_eqc *Q = static_cast<struct pj_opaque_eqc*>(P->host->ctx->allocator->svm_calloc (1, sizeof (struct pj_opaque_eqc)));
    if (nullptr==Q)
        return pj_default_destructor (P, PROJ_ERR_OTHER /*ENOMEM*/);
    P->opaque = Q;

    if ((Q->rc = cos(pj_param(P->host->ctx, P->host->params, "rlat_ts").f)) <= 0.)
    {
        proj_log_error(P, _("Invalid value for lat_ts: |lat_ts| should be <= 90°"));
        return pj_default_destructor(P, PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE);
    }
    P->host->inv = PJ_MAKE_KERNEL(eqc_s_inverse);
    P->host->fwd = PJ_MAKE_KERNEL(eqc_s_forward);
    P->es = 0.;

    return P;
}

#endif /* !PROJ_OPENCL_DEVICE */

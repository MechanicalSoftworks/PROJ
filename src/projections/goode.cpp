#define PJ_LIB__

#include <errno.h>
#include <math.h>

#include "proj.h"
#include "proj_internal.h"

PROJ_HEAD(goode, "Goode Homolosine") "\n\tPCyl, Sph";

#define Y_COR   0.05280
#define PHI_LIM 0.71093078197902358062

C_NAMESPACE PJ *pj_sinu(PJ *), *pj_moll(PJ *);

namespace { // anonymous namespace
struct pj_opaque {
    PJ *sinu;
    PJ *moll;
};
} // anonymous namespace


static PJ_XY goode_s_forward (PJ_LP lp, PJ *P) {           /* Spheroidal, forward */
    PJ_XY xy;
    struct pj_opaque *Q = static_cast<struct pj_opaque*>(P->opaque);

    if (fabs(lp.phi) <= PHI_LIM)
        xy = Q->sinu->host->fwd(lp, Q->sinu);
    else {
        xy = Q->moll->host->fwd(lp, Q->moll);
        xy.y -= lp.phi >= 0.0 ? Y_COR : -Y_COR;
    }
    return xy;
}


static PJ_LP goode_s_inverse (PJ_XY xy, PJ *P) {           /* Spheroidal, inverse */
    PJ_LP lp;
    struct pj_opaque *Q = static_cast<struct pj_opaque*>(P->opaque);

    if (fabs(xy.y) <= PHI_LIM)
        lp = Q->sinu->host->inv(xy, Q->sinu);
    else {
        xy.y += xy.y >= 0.0 ? Y_COR : -Y_COR;
        lp = Q->moll->host->inv(xy, Q->moll);
    }
    return lp;
}


static PJ *destructor (PJ *P, int errlev) {              /* Destructor */
    if (nullptr==P)
        return nullptr;
    if (nullptr==P->opaque)
        return pj_default_destructor (P, errlev);
    proj_destroy (static_cast<struct pj_opaque*>(P->opaque)->sinu);
    proj_destroy (static_cast<struct pj_opaque*>(P->opaque)->moll);
    return pj_default_destructor (P, errlev);
}



PJ *PROJECTION(goode) {
    struct pj_opaque *Q = static_cast<struct pj_opaque*>(svm_calloc (P->host->ctx, 1, sizeof (struct pj_opaque)));
    if (nullptr==Q)
        return pj_default_destructor (P, PROJ_ERR_OTHER /*ENOMEM*/);
    P->opaque = Q;
    P->host->destructor = destructor;

    P->es = 0.;
    Q->sinu = pj_sinu(nullptr);
    Q->moll = pj_moll(nullptr);
    if (Q->sinu == nullptr || Q->moll == nullptr)
        return destructor (P, PROJ_ERR_OTHER /*ENOMEM*/);
    Q->sinu->es = 0.;
    Q->sinu->host->ctx = P->host->ctx;
    Q->moll->host->ctx = P->host->ctx;
    Q->sinu = pj_sinu(Q->sinu);
    Q->moll = pj_moll(Q->moll);
    if (Q->sinu == nullptr || Q->moll == nullptr)
        return destructor (P, PROJ_ERR_OTHER /*ENOMEM*/);

    P->host->fwd = PJ_MAKE_KERNEL(goode_s_forward);
    P->host->inv = PJ_MAKE_KERNEL(goode_s_inverse);

    return P;
}

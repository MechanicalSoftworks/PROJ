#define PJ_LIB__

#include "proj_internal.h"

PROJ_HEAD(noop,    "No operation");

static PJ_COORD noop(PJ_COORD coord, PJ *P) {
    (void) P;
    return coord;
}

PJ *CONVERSION(noop, 0) {
    P->host->fwd4d = PJ_MAKE_KERNEL(noop);
    P->host->inv4d = PJ_MAKE_KERNEL(noop);
    P->left  = PJ_IO_UNITS_WHATEVER;
    P->right = PJ_IO_UNITS_WHATEVER;
    return P;
}


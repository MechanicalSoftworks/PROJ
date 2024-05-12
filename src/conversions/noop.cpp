#define PJ_LIB__

#include "proj_internal.h"

PROJ_HEAD(noop,    "No operation");

PJ_COORD noop_operator(PJ_COORD coord, PJ *P) {
    (void) P;
    return coord;
}

PJ *CONVERSION(noop, 0) {
    P->fwd4d = PJ_MAKE_KERNEL(noop_operator);
    P->inv4d = PJ_MAKE_KERNEL(noop_operator);
    P->left  = PJ_IO_UNITS_WHATEVER;
    P->right = PJ_IO_UNITS_WHATEVER;
    return P;
}


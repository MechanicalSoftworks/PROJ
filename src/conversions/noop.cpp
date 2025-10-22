#define PJ_LIB__
#include "../proj_kernel.h"

PROJ_HEAD(noop,    "No operation");

PROJ_NOINLINE PJ_COORD noop_operator(PJ_COORD coord, __global PJ *P) {
    (void) P;
    return coord;
}

#ifndef PROJ_OPENCL_DEVICE

PJ *CONVERSION(noop, 0) {
    P->fwd4d = PJ_MAKE_KERNEL(noop_operator);
    P->inv4d = PJ_MAKE_KERNEL(noop_operator);
    P->left  = PJ_IO_UNITS_WHATEVER;
    P->right = PJ_IO_UNITS_WHATEVER;
    return P;
}

#endif /* !PROJ_OPENCL_DEVICE */

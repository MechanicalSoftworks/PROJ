/******************************************************************************
 * Project:  PROJ.4
 * Purpose:  Internal plumbing for the PROJ.4 library.
 *
 * Author:   Thomas Knudsen, <thokn@sdfe.dk>
 *
 ******************************************************************************
 * Copyright (c) 2016, 2017, Thomas Knudsen / SDFE
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
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO COORD SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *****************************************************************************/

#ifndef PROJ_INTERNAL_SHARED_H
#define PROJ_INTERNAL_SHARED_H

#include "proj_shared.h"

// The constructors need access to the real header.
#ifndef PROJ_OPENCL_DEVICE
#   include "proj_internal.h"
#endif

struct PJconsts;

#define STATIC_ASSERT(COND) ((void)sizeof(char[(COND) ? 1 : -1]))

#ifndef PJ_TODEG
#define PJ_TODEG(rad)  ((rad)*180.0/M_PI)
#endif
#ifndef PJ_TORAD
#define PJ_TORAD(deg)  ((deg)*M_PI/180.0)
#endif

/* Maximum latitudinal overshoot accepted */
#define PJ_EPS_LAT 1e-12

#ifndef NULL
#  define NULL 0
#endif

#ifndef FALSE
#  define FALSE 0
#endif

#ifndef TRUE
#  define TRUE  1
#endif

#ifndef MAX
#  define MIN(a,b)      ((a<b) ? a : b)
#  define MAX(a,b)      ((a>b) ? a : b)
#endif

#ifndef ABS
#  define ABS(x)        ((x<0) ? (-1*(x)) : x)
#endif

#if INT_MAX == 2147483647
typedef int pj_int32;
#elif LONG_MAX == 2147483647
typedef long pj_int32;
#else
#warning It seems no 32 - bit integer type is available
#endif

/* If we still haven't got M_PI*, we rely on our own defines.
 * For example, this is necessary when compiling with gcc and
 * the -ansi flag.
 */
#ifndef M_PI
#define M_PI            3.14159265358979323846
#define M_PI_2          1.57079632679489661923
#define M_PI_4          0.78539816339744830962
#define M_2_PI          0.63661977236758134308
#endif

 /* M_SQRT2 might be missing */
#ifndef M_SQRT2
#define M_SQRT2         1.41421356237309504880
#endif

/* some more useful math constants and aliases */
#define M_FORTPI        M_PI_4                   /* pi/4 */
#define M_HALFPI        M_PI_2                   /* pi/2 */
#define M_PI_HALFPI     4.71238898038468985769   /* 1.5*pi */
#define M_TWOPI         6.28318530717958647693   /* 2*pi */
#define M_TWO_D_PI      M_2_PI                   /* 2/pi */
#define M_TWOPI_HALFPI  7.85398163397448309616   /* 2.5*pi */

typedef enum
{
    PJ_CO_YIELD,  // Return control to the dispatcher.
    PJ_CO_DONE,   // No more work, dispatcher should pop the stack.
                  // The 'coo' entry is also transferred to the next top.
    PJ_CO_ERROR
} PJcoroutine_code_t;

#define PJ_YIELD(e, s)  \
    e->step = s;        \
    goto YIELD;         \
p##s:

struct PJstack_entry_s;

// These are generated using the results of PJscan. Transmitted to the kernel at compile time.
typedef int PJ_COROUTINE_ID;
typedef int PJ_FWD_2D_ID;
typedef int PJ_INV_2D_ID;
typedef int PJ_FWD_3D_ID;
typedef int PJ_INV_3D_ID;
typedef int PJ_FWD_4D_ID;
typedef int PJ_INV_4D_ID;

#ifdef PROJ_OPENCL_DEVICE

// Initialised on the host.
#   define PJ_FIELD(type, name, value) type name
#   define PJ_HOST_MUTABLE

#   define PJ_FUNCTION_PTR(name)        name##_id
#   define PJ_GET_COROUTINE(P, type)    P->type

#   define lround round
#   define nullptr 0

#   define cl_constant  __constant
#   define cl_local     __local

#   define PROJ_DLL

#else

#   define PJ_FIELD(type, name, value) type name = value
#   define PJ_HOST_MUTABLE mutable

#   define PJ_FUNCTION_PTR(name)       &name
#   define PJ_GET_COROUTINE(P, type)    P->host->type.fn

#   define cl_constant
#   define cl_local

#   include <math.h>

    typedef PJcoroutine_code_t(*PJ_COROUTINE)(cl_local struct PJstack_s* stack, cl_local struct PJstack_entry_s* e);

#endif

enum pj_io_units {
    PJ_IO_UNITS_WHATEVER = 0,  /* Doesn't matter (or depends on pipeline neighbours) */
    PJ_IO_UNITS_CLASSIC = 1,  /* Scaled meters (right), projected system */
    PJ_IO_UNITS_PROJECTED = 2,  /* Meters, projected system */
    PJ_IO_UNITS_CARTESIAN = 3,  /* Meters, 3D cartesian system */
    PJ_IO_UNITS_RADIANS = 4,  /* Radians */
    PJ_IO_UNITS_DEGREES = 5,  /* Degrees */

};
enum pj_io_units pj_left(struct PJconsts* P);
enum pj_io_units pj_right(struct PJconsts* P);

/* datum_type values */
#define PJD_UNKNOWN   0
#define PJD_3PARAM    1
#define PJD_7PARAM    2
#define PJD_GRIDSHIFT 3
#define PJD_WGS84     4   /* WGS84 (or anything considered equivalent) */

typedef struct PJstack_entry_s
{
#ifdef PROJ_OPENCL_DEVICE
    PJ_COROUTINE_ID fn;
#else
    PJ_COROUTINE    fn;
#endif

    // State.
    PJ_COORD        coo;
    PJ*             P;
    int             step;
    union {
        int         i;          // Pipeline.
        int         last_errno; // pj_fwd+pj_inv.
    };
} PJstack_entry_t;

#define PR_CO_STACK_SIZE 16
typedef struct PJstack_s
{
    PJstack_entry_t s[PR_CO_STACK_SIZE];
    int             n;
} PJstack_t;

/* Context data that's shared between the OpenCL host and evice */
struct pj_ctx_shared {
    PJ_FIELD(int, last_errno, 0);
};

/* base projection data structure */
struct PJconsts {

    /*************************************************************************************

                         G E N E R A L   P A R A M E T E R   S T R U C T

    **************************************************************************************

        TODO: Need some description here - especially about the thread context...
        This is the struct behind the PJ typedef

    **************************************************************************************/

    PJ_FIELD(struct PJhost*, host, nullptr);

    PJ_FIELD(struct pj_ctx_shared*, shared_ctx, nullptr);

    PJ_FIELD(struct PJconsts*, parent, nullptr);                   /* Parent PJ of pipeline steps - nullptr if not a pipeline step */

    PJ_FIELD(struct geod_geodesic*, geod, nullptr);         /* For geodesic computations */
    PJ_FIELD(void*, opaque, nullptr);                       /* Projection specific parameters, Defined in PJ_*.c */
    PJ_FIELD(int, inverted, 0);                             /* Tell high level API functions to swap inv/fwd */
    int pad0;

    /*************************************************************************************

                          F U N C T I O N    P O I N T E R S

    **************************************************************************************

        For projection xxx, these are pointers to functions in the corresponding
        PJ_xxx.c file.

        pj_init() delegates the setup of these to pj_projection_specific_setup_xxx(),
        a name which is currently hidden behind the magic curtain of the PROJECTION
        macro.

    **************************************************************************************/

    PJ_FIELD(PJ_COROUTINE_ID, co_fwd, 0);
    PJ_FIELD(PJ_COROUTINE_ID, co_inv, 0);
    PJ_FIELD(PJ_COROUTINE_ID, co_fwd3d, 0);
    PJ_FIELD(PJ_COROUTINE_ID, co_inv3d, 0);
    PJ_FIELD(PJ_COROUTINE_ID, co_fwd4d, 0);
    PJ_FIELD(PJ_COROUTINE_ID, co_inv4d, 0);

    PJ_FIELD(PJ_FWD_2D_ID, fwd, 0);
    PJ_FIELD(PJ_INV_2D_ID, inv, 0);
    PJ_FIELD(PJ_FWD_3D_ID, fwd3d, 0);
    PJ_FIELD(PJ_INV_3D_ID, inv3d, 0);
    PJ_FIELD(PJ_FWD_4D_ID, fwd4d, 0);
    PJ_FIELD(PJ_INV_4D_ID, inv4d, 0);


    /*************************************************************************************

                          E L L I P S O I D     P A R A M E T E R S

    **************************************************************************************

        Despite YAGNI, we add a large number of ellipsoidal shape parameters, which
        are not yet set up in pj_init. They are, however, inexpensive to compute,
        compared to the overall time taken for setting up the complex PJ object
        (cf. e.g. https://en.wikipedia.org/wiki/Angular_eccentricity).

        But during single point projections it will often be a useful thing to have
        these readily available without having to recompute at every pj_fwd / pj_inv
        call.

        With this wide selection, we should be ready for quite a number of geodetic
        algorithms, without having to incur further ABI breakage.

    **************************************************************************************/

    /* The linear parameters */

    PJ_FIELD(double, a, 0.0);          /* semimajor axis (radius if eccentricity==0) */
    PJ_FIELD(double, b, 0.0);          /* semiminor axis */
    PJ_FIELD(double, ra, 0.0);         /* 1/a */
    PJ_FIELD(double, rb, 0.0);         /* 1/b */

    /* The eccentricities */

    PJ_FIELD(double, alpha, 0.0);               /* angular eccentricity */
    PJ_FIELD(double, e, 0.0);                   /* first  eccentricity */
    PJ_FIELD(double, es, 0.0);                  /* first  eccentricity squared */
    PJ_FIELD(double, e2, 0.0);                  /* second eccentricity */
    PJ_FIELD(double, e2s, 0.0);                 /* second eccentricity squared */
    PJ_FIELD(double, e3, 0.0);                  /* third  eccentricity */
    PJ_FIELD(double, e3s, 0.0);                 /* third  eccentricity squared */
    PJ_FIELD(double, one_es, 0.0);              /* 1 - e^2 */
    PJ_FIELD(double, rone_es, 0.0);             /* 1/one_es */


    /* The flattenings */
    PJ_FIELD(double, f, 0.0);                   /* first  flattening */
    PJ_FIELD(double, f2, 0.0);                  /* second flattening */
    PJ_FIELD(double, n, 0.0);                   /* third  flattening */
    PJ_FIELD(double, rf, 0.0);                  /* 1/f  */
    PJ_FIELD(double, rf2, 0.0);                 /* 1/f2 */
    PJ_FIELD(double, rn, 0.0);                   /* 1/n  */

    /* This one's for GRS80 */
    PJ_FIELD(double, J, 0.0);                   /* "Dynamic form factor" */

    PJ_FIELD(double, es_orig, 0.0);    /* es and a before any +proj related adjustment */
    PJ_FIELD(double, a_orig, 0.0);


    /*************************************************************************************

                          C O O R D I N A T E   H A N D L I N G

    **************************************************************************************/

    PJ_FIELD(int, over, 0);                  /* Over-range flag */
    PJ_FIELD(int, geoc, 0);                  /* Geocentric latitude flag */
    PJ_FIELD(int, is_latlong, 0);            /* proj=latlong ... not really a projection at all */
    PJ_FIELD(int, is_geocent, 0);            /* proj=geocent ... not really a projection at all */
    PJ_FIELD(int, need_ellps, 0);            /* 0 for operations that are purely cartesian */
    PJ_FIELD(int, skip_fwd_prepare, 0);
    PJ_FIELD(int, skip_fwd_finalize, 0);
    PJ_FIELD(int, skip_inv_prepare, 0);
    PJ_FIELD(int, skip_inv_finalize, 0);
    int pad1;

    PJ_FIELD(int, left,  PJ_IO_UNITS_WHATEVER); /* Flags for input/output coordinate types */
    PJ_FIELD(int, right, PJ_IO_UNITS_WHATEVER);
    
    /* These PJs are used for implementing cs2cs style coordinate handling in the 4D API */
    PJ_FIELD(struct PJconsts*, axisswap, nullptr);
    PJ_FIELD(struct PJconsts*, cart, nullptr);
    PJ_FIELD(struct PJconsts*, cart_wgs84, nullptr);
    PJ_FIELD(struct PJconsts*, helmert, nullptr);
    PJ_FIELD(struct PJconsts*, hgridshift, nullptr);
    PJ_FIELD(struct PJconsts*, vgridshift, nullptr);


    /*************************************************************************************

                       C A R T O G R A P H I C       O F F S E T S

    **************************************************************************************/

    PJ_FIELD(double, lam0, 0.0);    /* central meridian */
    PJ_FIELD(double, phi0, 0.0);    /* central parallel */
    PJ_FIELD(double, x0, 0.0);      /* false easting */
    PJ_FIELD(double, y0, 0.0);      /* false northing  */
    PJ_FIELD(double, z0, 0.0);      /* height origin */
    PJ_FIELD(double, t0, 0.0);      /* time origin */


    /*************************************************************************************

                                    S C A L I N G

    **************************************************************************************/

    PJ_FIELD(double, k0, 0.0);                  /* General scaling factor - e.g. the 0.9996 of UTM */
    PJ_FIELD(double, to_meter, 0.0);            /* Plane coordinate scaling. Internal unit [m] */
    PJ_FIELD(double, fr_meter, 0.0);
    PJ_FIELD(double, vto_meter, 0.0);           /* Vertical scaling. Internal unit [m] */
    PJ_FIELD(double, vfr_meter, 0.0);


    /*************************************************************************************

                  D A T U M S   A N D   H E I G H T   S Y S T E M S

    **************************************************************************************

        It may be possible, and meaningful, to move the list parts of this up to the
        PJ_CONTEXT level.

    **************************************************************************************/

    PJ_FIELD(double, datum_params[7], { 0 });   /* Parameters for 3PARAM and 7PARAM */
    PJ_FIELD(int,    datum_type, PJD_UNKNOWN);  /* PJD_UNKNOWN/3PARAM/7PARAM/GRIDSHIFT/WGS84 */

    PJ_FIELD(int  ,  has_geoid_vgrids, 0);      /* used by legacy transform.cpp */
    PJ_FIELD(void*,  hgrids_legacy, nullptr);   /* used by legacy transform.cpp. Is a pointer to a ListOfHGrids* */ 
    PJ_FIELD(void*,  vgrids_legacy, nullptr);   /* used by legacy transform.cpp. Is a pointer to a ListOfVGrids* */ 

    PJ_FIELD(double, from_greenwich, 0.0);      /* prime meridian offset (in radians) */
    PJ_FIELD(double, long_wrap_center, 0.0);    /* 0.0 for -180 to 180, actually in radians*/
    PJ_FIELD(int   , is_long_wrap_set, 0);
    PJ_FIELD(char  , axis[4], { 0 });           /* Axis order, pj_transform/pj_adjust_axis */

    /*************************************************************************************
     ISO-19111 interface
    **************************************************************************************/

    // cache pj_get_type() result to help for repeated calls to proj_factors()
    PJ_HOST_MUTABLE PJ_FIELD(int, type, PJ_TYPE_UNKNOWN);

    /*************************************************************************************
     proj_create_crs_to_crs() alternative coordinate operations
    **************************************************************************************/
    PJ_FIELD(int, iCurCoordOp, -1);

    /*************************************************************************************

                 E N D   O F    G E N E R A L   P A R A M E T E R   S T R U C T

    **************************************************************************************/

#ifndef PROJ_OPENCL_DEVICE
    PJconsts();
    PJconsts(const PJconsts &) = delete;
    PJconsts &operator=(const PJconsts &) = delete;
#endif
};

/* Geographical to geocentric latitude - another of the "simple, but useful" */
PJ_COORD pj_geocentric_latitude(const PJ* P, PJ_DIRECTION direction, PJ_COORD coord);

double PROJ_DLL adjlon(double);
void proj_context_errno_set(struct pj_ctx_shared* ctx, int err);
struct pj_ctx_shared* pj_get_ctx_shared(const PJ*);

// Coroutines.
void stack_new(cl_local PJstack_t* stack);
PJ_COORD stack_exec(cl_local PJstack_t* stack);

void push_proj_trans(cl_local PJstack_t* stack, PJ* P, PJ_DIRECTION direction, PJ_COORD coord);
void push_approx_3D_trans(cl_local PJstack_t* stack, PJ* P, PJ_DIRECTION direction, PJ_COORD coord);
void push_approx_2D_trans(cl_local PJstack_t* stack, PJ* P, PJ_DIRECTION direction, PJ_COORD coord);

PJcoroutine_code_t pj_fwd4d_co(cl_local PJstack_t* stack, cl_local PJstack_entry_t* e);
PJcoroutine_code_t pj_inv4d_co(cl_local PJstack_t* stack, cl_local PJstack_entry_t* e);

PJcoroutine_code_t pj_fwd3d_co(cl_local PJstack_t* stack, cl_local PJstack_entry_t* e);
PJcoroutine_code_t pj_inv3d_co(cl_local PJstack_t* stack, cl_local PJstack_entry_t* e);

PJcoroutine_code_t pj_fwd_co(cl_local PJstack_t* stack, cl_local PJstack_entry_t* e);
PJcoroutine_code_t pj_inv_co(cl_local PJstack_t* stack, cl_local PJstack_entry_t* e);

#endif // !PROJ_INTERNAL_SHARED_H

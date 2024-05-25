/******************************************************************************
 * Project:  PROJ.4
 * Purpose:  Definitions that are common between OpenCL host and device.
 *
 ******************************************************************************
 * Copyright (c) 2016, 2017, Thomas Knudsen / SDFE
 * Copyright (c) 2018, Even Rouault
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

#ifndef PROJ_SHARED_H
#define PROJ_SHARED_H

#ifndef PROJ_OPENCL_DEVICE
#   define __global
#   define __constant
#   define __local
#endif

/* first forward declare everything needed */

/* Data type for generic geodetic 3D data plus epoch information */
union PJ_COORD;
typedef union PJ_COORD PJ_COORD;

struct PJ_AREA;
typedef struct PJ_AREA PJ_AREA;

struct P5_FACTORS {                  /* Common designation */
    double meridional_scale;               /* h */
    double parallel_scale;                 /* k */
    double areal_scale;                    /* s */

    double angular_distortion;             /* omega */
    double meridian_parallel_angle;        /* theta-prime */
    double meridian_convergence;           /* alpha */

    double tissot_semimajor;               /* a */
    double tissot_semiminor;               /* b */

    double dx_dlam, dx_dphi;
    double dy_dlam, dy_dphi;
};
typedef struct P5_FACTORS PJ_FACTORS;

/* The compute type - properly namespaced synonym for pj_compute */
struct pj_allocator;
typedef struct pj_allocator PJ_ALLOCATOR;

/* The context type - properly namespaced synonym for pj_ctx */
struct pj_ctx;
typedef struct pj_ctx PJ_CONTEXT;

/* Data type for projection/transformation information */
struct PJconsts;
typedef struct PJconsts PJ;         /* the PJ object herself */

/* Data type for library level information */
struct PJ_INFO;
typedef struct PJ_INFO PJ_INFO;

struct PJ_PROJ_INFO;
typedef struct PJ_PROJ_INFO PJ_PROJ_INFO;

struct PJ_GRID_INFO;
typedef struct PJ_GRID_INFO PJ_GRID_INFO;

struct PJ_INIT_INFO;
typedef struct PJ_INIT_INFO PJ_INIT_INFO;


/* Geodetic, mostly spatiotemporal coordinate types */
typedef struct { double   x,   y,  z, t; }  PJ_XYZT;
typedef struct { double   u,   v,  w, t; }  PJ_UVWT;
typedef struct { double lam, phi,  z, t; }  PJ_LPZT;
typedef struct { double o, p, k; }          PJ_OPK;  /* Rotations: omega, phi, kappa */
typedef struct { double e, n, u; }          PJ_ENU;  /* East, North, Up */
typedef struct { double s, a1, a2; }        PJ_GEOD; /* Geodesic length, fwd azi, rev azi */

/* Classic proj.4 pair/triplet types - moved into the PJ_ name space */
typedef struct { double   u,   v; }  PJ_UV;
typedef struct { double   x,   y; }  PJ_XY;
typedef struct { double lam, phi; }  PJ_LP;

typedef struct { double   x,   y,  z; }  PJ_XYZ;
typedef struct { double   u,   v,  w; }  PJ_UVW;
typedef struct { double lam, phi,  z; }  PJ_LPZ;


/* Avoid preprocessor renaming and implicit type-punning: Use a union to make it explicit */
union PJ_COORD {
     double v[4];   /* First and foremost, it really is "just 4 numbers in a vector" */
    PJ_XYZT xyzt;
    PJ_UVWT uvwt;
    PJ_LPZT lpzt;
    PJ_GEOD geod;
     PJ_OPK opk;
     PJ_ENU enu;
     PJ_XYZ xyz;
     PJ_UVW uvw;
     PJ_LPZ lpz;
      PJ_XY xy;
      PJ_UV uv;
      PJ_LP lp;
};

/* PROJ error codes */

/** Error codes typically related to coordinate operation initialization
 * Note: some of them can also be emitted during coordinate transformation,
 * like PROJ_ERR_INVALID_OP_FILE_NOT_FOUND_OR_INVALID in case the resource loading
 * is deferred until it is really needed.
 */
#define PROJ_ERR_INVALID_OP                           1024                        /* other/unspecified error related to coordinate operation initialization */
#define PROJ_ERR_INVALID_OP_WRONG_SYNTAX              (PROJ_ERR_INVALID_OP+1)     /* invalid pipeline structure, missing +proj argument, etc */
#define PROJ_ERR_INVALID_OP_MISSING_ARG               (PROJ_ERR_INVALID_OP+2)     /* missing required operation parameter */
#define PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE         (PROJ_ERR_INVALID_OP+3)     /* one of the operation parameter has an illegal value */
#define PROJ_ERR_INVALID_OP_MUTUALLY_EXCLUSIVE_ARGS   (PROJ_ERR_INVALID_OP+4)     /* mutually exclusive arguments */
#define PROJ_ERR_INVALID_OP_FILE_NOT_FOUND_OR_INVALID (PROJ_ERR_INVALID_OP+5)     /* file not found (particular case of PROJ_ERR_INVALID_OP_ILLEGAL_ARG_VALUE) */

/** Error codes related to transformation on a specific coordinate */
#define PROJ_ERR_COORD_TRANSFM                           2048                           /* other error related to coordinate transformation */
#define PROJ_ERR_COORD_TRANSFM_INVALID_COORD             (PROJ_ERR_COORD_TRANSFM+1)     /* for e.g lat > 90deg */
#define PROJ_ERR_COORD_TRANSFM_OUTSIDE_PROJECTION_DOMAIN (PROJ_ERR_COORD_TRANSFM+2)     /* coordinate is outside of the projection domain. e.g approximate mercator with |longitude - lon_0| > 90deg, or iterative convergence method failed */
#define PROJ_ERR_COORD_TRANSFM_NO_OPERATION              (PROJ_ERR_COORD_TRANSFM+3)     /* no operation found, e.g if no match the required accuracy, or if ballpark transformations were asked to not be used and they would be only such candidate */
#define PROJ_ERR_COORD_TRANSFM_OUTSIDE_GRID              (PROJ_ERR_COORD_TRANSFM+4)     /* point to transform falls outside grid or subgrid */
#define PROJ_ERR_COORD_TRANSFM_GRID_AT_NODATA            (PROJ_ERR_COORD_TRANSFM+5)     /* point to transform falls in a grid cell that evaluates to nodata */

/** Other type of errors */
#define PROJ_ERR_OTHER                                   4096
#define PROJ_ERR_OTHER_API_MISUSE                        (PROJ_ERR_OTHER+1)             /* error related to a misuse of PROJ API */
#define PROJ_ERR_OTHER_NO_INVERSE_OP                     (PROJ_ERR_OTHER+2)             /* no inverse method available */
#define PROJ_ERR_OTHER_NETWORK_ERROR                     (PROJ_ERR_OTHER+3)             /* failure when accessing a network resource */

/** \brief Object type. */
typedef enum
{
    PJ_TYPE_UNKNOWN,

    PJ_TYPE_ELLIPSOID,

    PJ_TYPE_PRIME_MERIDIAN,

    PJ_TYPE_GEODETIC_REFERENCE_FRAME,
    PJ_TYPE_DYNAMIC_GEODETIC_REFERENCE_FRAME,
    PJ_TYPE_VERTICAL_REFERENCE_FRAME,
    PJ_TYPE_DYNAMIC_VERTICAL_REFERENCE_FRAME,
    PJ_TYPE_DATUM_ENSEMBLE,

    /** Abstract type, not returned by proj_get_type() */
    PJ_TYPE_CRS,

    PJ_TYPE_GEODETIC_CRS,
    PJ_TYPE_GEOCENTRIC_CRS,

    /** proj_get_type() will never return that type, but
     * PJ_TYPE_GEOGRAPHIC_2D_CRS or PJ_TYPE_GEOGRAPHIC_3D_CRS. */
    PJ_TYPE_GEOGRAPHIC_CRS,

    PJ_TYPE_GEOGRAPHIC_2D_CRS,
    PJ_TYPE_GEOGRAPHIC_3D_CRS,
    PJ_TYPE_VERTICAL_CRS,
    PJ_TYPE_PROJECTED_CRS,
    PJ_TYPE_COMPOUND_CRS,
    PJ_TYPE_TEMPORAL_CRS,
    PJ_TYPE_ENGINEERING_CRS,
    PJ_TYPE_BOUND_CRS,
    PJ_TYPE_OTHER_CRS,

    PJ_TYPE_CONVERSION,
    PJ_TYPE_TRANSFORMATION,
    PJ_TYPE_CONCATENATED_OPERATION,
    PJ_TYPE_OTHER_COORDINATE_OPERATION,

    PJ_TYPE_TEMPORAL_DATUM,
    PJ_TYPE_ENGINEERING_DATUM,
    PJ_TYPE_PARAMETRIC_DATUM,
} PJ_TYPE;

/* Apply transformation to observation - in forward or inverse direction */
enum PJ_DIRECTION {
    PJ_FWD = 1,   /* Forward    */
    PJ_IDENT = 0,   /* Do nothing */
    PJ_INV = -1    /* Inverse    */
};
typedef enum PJ_DIRECTION PJ_DIRECTION;

typedef struct PJstack_entry_s
{
    PJ_COORD        coo;
    __global PJ*    P;

    union {
        int         i;              // Pipeline.
        int         last_errno;     // pj_fwd+pj_inv.
    };

    unsigned short  coroutine_id;   // PJ_COROUTINE_ID
    unsigned short  step;
} PJstack_entry_t;

#define PJ_CO_STACK_SIZE 16
typedef struct PJstack_s
{
    PJstack_entry_t s[PJ_CO_STACK_SIZE];
    int             n;
} PJstack_t;

#ifdef PJ_LIB__

#   define PROJ_HEAD(name, desc) static const char des_##name [] = desc

#endif /* PJ_LIB__ */

#endif /* ndef PROJ_SHARED_H */

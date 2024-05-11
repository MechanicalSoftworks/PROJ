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

#ifndef PROJ_INTERNAL_H
#define PROJ_INTERNAL_H

#ifndef __cplusplus
#error "proj_internal.h can only be included from a C++ file"
#endif

#ifdef _MSC_VER
#  ifndef _CRT_SECURE_NO_DEPRECATE
#    define _CRT_SECURE_NO_DEPRECATE
#  endif
#  ifndef _CRT_NONSTDC_NO_DEPRECATE
#    define _CRT_NONSTDC_NO_DEPRECATE
#  endif
/* enable predefined math constants M_* for MS Visual Studio workaround */
#  ifndef _USE_MATH_DEFINES
#     define _USE_MATH_DEFINES
#  endif
#endif

/* standard inclusions */
#include <limits.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "proj/common.hpp"
#include "proj/coordinateoperation.hpp"

#include <string>
#include <vector>
#include <set>
#include <map>
#include <new>

#include "proj.h"

#ifdef PROJ_RENAME_SYMBOLS
#include "proj_symbol_rename.h"
#endif

#define C_NAMESPACE extern "C"
#define C_NAMESPACE_VAR extern "C"

/* maximum path/filename */
#ifndef MAX_PATH_FILENAME
#define MAX_PATH_FILENAME 1024
#endif


/* maximum tag id length for +init and default files */
#ifndef ID_TAG_MAX
#define ID_TAG_MAX 50
#endif

/* Use WIN32 as a standard windows 32 bit declaration */
#if defined(_WIN32) && !defined(WIN32)
#  define WIN32
#endif

#if defined(_WINDOWS) && !defined(WIN32)
#  define WIN32
#endif

/* directory delimiter for DOS support */
#ifdef WIN32
#define DIR_CHAR '\\'
#else
#define DIR_CHAR '/'
#endif

#include "proj_internal_shared.h"

PJ_COORD PROJ_DLL proj_coord_error (void);

void PROJ_DLL proj_context_set (PJ *P, PJ_CONTEXT *ctx);
void proj_context_inherit (PJ *parent, PJ *child);

struct projCppContext;
/* not sure why we need to export it, but mingw needs it */
void PROJ_DLL proj_context_delete_cpp_context(struct projCppContext* cppContext);

PJ_COORD pj_fwd4d (PJ_COORD coo, PJ *P);
PJ_COORD pj_inv4d (PJ_COORD coo, PJ *P);

PJ_COORD PROJ_DLL pj_approx_2D_trans (PJ *P, PJ_DIRECTION direction, PJ_COORD coo);
PJ_COORD PROJ_DLL pj_approx_3D_trans (PJ *P, PJ_DIRECTION direction, PJ_COORD coo);

/* Provision for gettext translatable strings */
#define _(str) (str)

void PROJ_DLL proj_log_error (const PJ *P, const char *fmt, ...);
void proj_log_debug (PJ *P, const char *fmt, ...);
void proj_log_trace (PJ *P, const char *fmt, ...);

void proj_context_log_debug (PJ_CONTEXT *ctx, const char *fmt, ...);

int pj_ellipsoid (PJ *);
void pj_inherit_ellipsoid_def (const PJ *src, PJ *dst);
int pj_calc_ellipsoid_params (PJ *P, double a, double es);

char  PROJ_DLL *pj_chomp (char *c);
char  PROJ_DLL *pj_shrink (char *c);
size_t pj_trim_argc (char *args);
char **pj_trim_argv (size_t argc, char *args);
char  *pj_make_args (size_t argc, char **argv);

typedef struct { double r, i; } COMPLEX;

/* Forward declarations and typedefs for stuff needed inside the PJ object */
struct PJconsts;

union  PJ_COORD;
struct geod_geodesic;
struct ARG_list;
struct PJ_REGION_S;
typedef struct PJ_REGION_S  PJ_Region;
typedef struct ARG_list paralist;   /* parameter list */

#ifndef PROJ_H
typedef struct PJconsts PJ;         /* the PJ object herself */
typedef union  PJ_COORD PJ_COORD;
#endif

struct PJ_REGION_S {
    double ll_long;        /* lower left corner coordinates (radians) */
    double ll_lat;
    double ur_long;        /* upper right corner coordinates (radians) */
    double ur_lat;
};

struct PJ_AREA {
    int     bbox_set;
    double  west_lon_degree;
    double  south_lat_degree;
    double  east_lon_degree;
    double  north_lat_degree;
};


/*****************************************************************************

    Some function types that are especially useful when working with PJs

******************************************************************************

PJ_CONSTRUCTOR:

    A function taking a pointer-to-PJ as arg, and returning a pointer-to-PJ.
    Historically called twice: First with a 0 argument, to allocate memory,
    second with the first return value as argument, for actual setup.

PJ_DESTRUCTOR:

    A function taking a pointer-to-PJ and an integer as args, then first
    handling the deallocation of the PJ, afterwards handing the integer over
    to the error reporting subsystem, and finally returning a null pointer in
    support of the "return free (P)" (aka "get the hell out of here") idiom.

PJ_OPERATOR:

    A function taking a PJ_COORD and a pointer-to-PJ as args, applying the
    PJ to the PJ_COORD, and returning the resulting PJ_COORD.

*****************************************************************************/
typedef    PJ       *(* PJ_CONSTRUCTOR) (PJ *, PJ_CONTEXT *);
typedef    PJ       *(* PJ_DESTRUCTOR)  (PJ *, int);
typedef    PJ_COORD  (* PJ_OPERATOR)    (PJ_COORD, PJ *);
/****************************************************************************/


struct PJCoordOperation
{
    int idxInOriginalList;
    double minxSrc = 0.0;
    double minySrc = 0.0;
    double maxxSrc = 0.0;
    double maxySrc = 0.0;
    double minxDst = 0.0;
    double minyDst = 0.0;
    double maxxDst = 0.0;
    double maxyDst = 0.0;
    PJ* pj = nullptr;
    std::string name{};
    double accuracy = -1.0;
    bool isOffshore = false;

    PJCoordOperation(int idxInOriginalListIn,
                   double minxSrcIn, double minySrcIn, double maxxSrcIn, double maxySrcIn,
                   double minxDstIn, double minyDstIn, double maxxDstIn, double maxyDstIn,
                   PJ* pjIn, const std::string& nameIn, double accuracyIn, bool isOffshoreIn):
        idxInOriginalList(idxInOriginalListIn),
        minxSrc(minxSrcIn), minySrc(minySrcIn), maxxSrc(maxxSrcIn), maxySrc(maxySrcIn),
        minxDst(minxDstIn), minyDst(minyDstIn), maxxDst(maxxDstIn), maxyDst(maxyDstIn),
        pj(pjIn), name(nameIn),
        accuracy(accuracyIn),
        isOffshore(isOffshoreIn)
    {
    }

    PJCoordOperation(const PJCoordOperation&) = delete;

    PJCoordOperation(PJ_CONTEXT* ctx, const PJCoordOperation& other):
        idxInOriginalList(other.idxInOriginalList),
        minxSrc(other.minxSrc), minySrc(other.minySrc), maxxSrc(other.maxxSrc), maxySrc(other.maxySrc),
        minxDst(other.minxDst), minyDst(other.minyDst), maxxDst(other.maxxDst), maxyDst(other.maxyDst),
        pj(proj_clone(ctx, other.pj)),
        name(std::move(other.name)),
        accuracy(other.accuracy),
        isOffshore(other.isOffshore)
    {
    }

    PJCoordOperation(PJCoordOperation&& other):
        idxInOriginalList(other.idxInOriginalList),
        minxSrc(other.minxSrc), minySrc(other.minySrc), maxxSrc(other.maxxSrc), maxySrc(other.maxySrc),
        minxDst(other.minxDst), minyDst(other.minyDst), maxxDst(other.maxxDst), maxyDst(other.maxyDst),
        name(std::move(other.name)),
        accuracy(other.accuracy),
        isOffshore(other.isOffshore) {
        pj = other.pj;
        other.pj = nullptr;
    }

    PJCoordOperation& operator=(const PJCoordOperation&) = delete;

    bool operator == (const PJCoordOperation& other) const {
        return idxInOriginalList == other.idxInOriginalList &&
               minxSrc == other.minxSrc &&
               minySrc == other.minySrc &&
               maxxSrc == other.maxxSrc &&
               maxySrc == other.maxySrc &&
               minxDst == other.minxDst &&
               minyDst == other.minyDst &&
               maxxDst == other.maxxDst &&
               maxyDst == other.maxyDst &&
               name == other.name &&
               proj_is_equivalent_to(pj, other.pj, PJ_COMP_STRICT) &&
               accuracy == other.accuracy &&
               isOffshore == other.isOffshore;
    }

    bool operator != (const PJCoordOperation& other) const {
        return !(operator==(other));
    }

    ~PJCoordOperation()
    {
        proj_destroy(pj);
    }
};

enum class TMercAlgo
{
    AUTO, // Poder/Engsager if far from central meridian, otherwise Evenden/Snyder
    EVENDEN_SNYDER,
    PODER_ENGSAGER,
};

#define PROJ_STR_HELPER(x) #x
#define PROJ_STR(x) PROJ_STR_HELPER(x)

#define PJ_MAKE_KERNEL(name) PJ_FUNCTION_PTR(name)

struct PJhost
{
    /*************************************************************************************

                         G E N E R A L   P A R A M E T E R   S T R U C T

    **************************************************************************************/

    PJ_CONTEXT *ctx = nullptr;
    
    const char *short_name = nullptr; /* From pj_list.h */
    const char *descr = nullptr;   /* From pj_list.h or individual PJ_*.c file */
    paralist *params = nullptr;    /* Parameter list */
    char *def_full = nullptr;      /* Full textual definition (usually 0 - set by proj_pj_info) */

    /* For debugging / logging purposes */
    char *def_size = nullptr;      /* Shape and size parameters extracted from params */
    char *def_shape = nullptr;
    char *def_spherification = nullptr;
    char *def_ellps = nullptr;

    const char *file = nullptr;

    /*************************************************************************************

                          F U N C T I O N    P O I N T E R S

    **************************************************************************************/

    PJ_DESTRUCTOR destructor = nullptr;
    void   (*reassign_context)(PJ*, PJ_CONTEXT*) = nullptr;

    void   (*map_pj)(PJ* P, bool map) = nullptr;
    void   map_svm(void* ptr, bool map);

    /*************************************************************************************
     ISO-19111 interface
    **************************************************************************************/

    NS_PROJ::common::IdentifiedObjectPtr iso_obj{};

    // cached results
    mutable std::string lastWKT{};
    mutable std::string lastPROJString{};
    mutable std::string lastJSONString{};
    mutable bool gridsNeededAsked = false;
    mutable std::vector<NS_PROJ::operation::GridDescription> gridsNeeded{};

    /*************************************************************************************
     proj_create_crs_to_crs() alternative coordinate operations
    **************************************************************************************/
    std::vector<PJCoordOperation> alternativeCoordinateOperations{};

    /*************************************************************************************

                 E N D   O F    G E N E R A L   P A R A M E T E R   S T R U C T

    **************************************************************************************/

    PJhost();
    PJhost(const PJhost&) = delete;
    PJhost&operator=(const PJhost&) = delete;
};

/* Parameter list (a copy of the +proj=... etc. parameters) */
struct ARG_list {
    paralist *next;
    char used;
#if (defined(__GNUC__) && __GNUC__ >= 8) || (defined(__clang__) && __clang_major__ >= 9)
    char param[]; /* variable-length member */
    /* Safer to use [] for gcc 8. See https://github.com/OSGeo/proj.4/pull/1087 */
    /* and https://gcc.gnu.org/bugzilla/show_bug.cgi?id=86914 */
#else
    char param[1]; /* variable-length member */
#endif
};


typedef union { double  f; int  i; char *s; } PROJVALUE;


struct PJ_DATUMS {
    const char    *id;           /* datum keyword */
    const char    *defn;         /* ie. "to_wgs84=..." */
    const char    *ellipse_id;   /* ie from ellipse table */
    const char    *comments;     /* EPSG code, etc */
};


struct DERIVS {
    double x_l, x_p;       /* derivatives of x for lambda-phi */
    double y_l, y_p;       /* derivatives of y for lambda-phi */
};

struct FACTORS {
    struct DERIVS der;
    double h, k;           /* meridional, parallel scales */
    double omega, thetap;  /* angular distortion, theta prime */
    double conv;           /* convergence */
    double s;              /* areal scale factor */
    double a, b;           /* max-min scale error */
    int    code;           /* always 0 */
};

// Legacy
struct projFileAPI_t;

struct projCppContext;

struct projNetworkCallbacksAndData
{
    bool enabled = false;
    bool enabled_env_variable_checked = false; // whereas we have checked PROJ_NETWORK env variable
    proj_network_open_cbk_type open = nullptr;
    proj_network_close_cbk_type close = nullptr;
    proj_network_get_header_value_cbk_type get_header_value = nullptr;
    proj_network_read_range_type read_range = nullptr;
    void* user_data = nullptr;
};

struct projGridChunkCache
{
    bool enabled = true;
    std::string filename{};
    long long max_size = 300 * 1024 * 1024;
    int ttl = 86400; // 1 day
};

struct projFileApiCallbackAndData
{
    PROJ_FILE_HANDLE* (*open_cbk)(PJ_CONTEXT *ctx, const char *filename, PROJ_OPEN_ACCESS access, void* user_data) = nullptr;
    size_t           (*read_cbk)(PJ_CONTEXT *ctx, PROJ_FILE_HANDLE*, void* buffer, size_t size, void* user_data) = nullptr;
    size_t           (*write_cbk)(PJ_CONTEXT *ctx, PROJ_FILE_HANDLE*, const void* buffer, size_t size, void* user_data) = nullptr;
    int              (*seek_cbk)(PJ_CONTEXT *ctx, PROJ_FILE_HANDLE*, long long offset, int whence, void* user_data) = nullptr;
    unsigned long long (*tell_cbk)(PJ_CONTEXT *ctx, PROJ_FILE_HANDLE*, void* user_data) = nullptr;
    void             (*close_cbk)(PJ_CONTEXT *ctx, PROJ_FILE_HANDLE*, void* user_data) = nullptr;

    int (*exists_cbk)(PJ_CONTEXT *ctx, const char *filename, void* user_data) = nullptr;
    int (*mkdir_cbk)(PJ_CONTEXT *ctx, const char *filename, void* user_data) = nullptr;
    int (*unlink_cbk)(PJ_CONTEXT *ctx, const char *filename, void* user_data) = nullptr;
    int (*rename_cbk)(PJ_CONTEXT *ctx, const char *oldPath, const char *newPath, void* user_data) = nullptr;

    void*            user_data = nullptr;
};

struct pj_allocator
{
    pj_allocator(void* u, PROJ_SVM_MALLOC_FUNCTION m, PROJ_SVM_CALLOC_FUNCTION c, PROJ_SVM_FREE_FUNCTION f, PROJ_SVM_UPDATE_FUNCTION map)
        : user_data(u), m_malloc(m), m_calloc(c), m_free(f), m_map(map)
    {}

    void*   user_data = nullptr;
   
    void* svm_malloc(size_t sz)           { return m_malloc(user_data, sz); }
    void* svm_calloc(size_t n, size_t sz) { return m_calloc(user_data, n, sz); }
    void svm_free(void* ptr)              { m_free(user_data, ptr); }
    void svm_map(void* ptr, bool map)     { m_map(user_data, ptr, map); }

    template<typename T>
    T* svm_new()
    {
        auto* p = svm_malloc(sizeof(T));
        if (p)
        {
            new (p) T;
        }
        return reinterpret_cast<T*>(p);
    }

    template<typename T>
    void svm_delete(T* p)
    {
        if (p)
        {
            p->~T();
        }

        svm_free(p);
    }

private:
    const PROJ_SVM_MALLOC_FUNCTION m_malloc;
    const PROJ_SVM_CALLOC_FUNCTION m_calloc;
    const PROJ_SVM_FREE_FUNCTION   m_free;
    const PROJ_SVM_UPDATE_FUNCTION m_map;
};

/* proj thread context */
struct pj_ctx{
    struct pj_ctx_shared *shared = nullptr;
    struct pj_allocator  *allocator = nullptr;

    std::string lastFullErrorMessage{}; // used by proj_context_errno_string
    int     debug_level = PJ_LOG_ERROR;
    void    (*logger)(void *, int, const char *) = nullptr;
    void    *logger_app_data = nullptr;
    struct projCppContext* cpp_context = nullptr; /* internal context for C++ code */
    int     use_proj4_init_rules = -1; /* -1 = unknown, 0 = no, 1 = yes */
    int     epsg_file_exists = -1; /* -1 = unknown, 0 = no, 1 = yes */
    std::string ca_bundle_path{};

    std::string env_var_proj_lib{}; // content of PROJ_LIB environment variable. Use Filemanager::getProjLibEnvVar() to access
    std::vector<std::string> search_paths{};
    const char **c_compat_paths = nullptr; // same, but for projinfo usage

    const char* (*file_finder_legacy) (const char*) = nullptr; // Only for proj_api compat. To remove once it is removed
    const char* (*file_finder) (PJ_CONTEXT *, const char*, void* user_data) = nullptr;
    void* file_finder_user_data = nullptr;

    bool defer_grid_opening = false; // set transiently by pj_obj_create()

    projFileApiCallbackAndData fileApi{};
    std::string custom_sqlite3_vfs_name{};
    std::string user_writable_directory{};

    // BEGIN ini file settings
    bool iniFileLoaded = false;
    std::string endpoint{};
    projNetworkCallbacksAndData networking{};
    projGridChunkCache gridChunkCache{};
    TMercAlgo defaultTmercAlgo = TMercAlgo::PODER_ENGSAGER; // can be overridden by content of proj.ini
    // END ini file settings

    int projStringParserCreateFromPROJStringRecursionCounter = 0; // to avoid potential infinite recursion in PROJStringParser::createFromPROJString()
    int pipelineInitRecursiongCounter = 0; // to avoid potential infinite recursion in pipeline.cpp


    pj_ctx() = delete;
    pj_ctx(pj_allocator *a);
    pj_ctx(const pj_ctx&);
    pj_ctx(const pj_ctx&, pj_allocator*);
    ~pj_ctx();

    pj_ctx& operator= (const pj_ctx&) = delete;

    projCppContext* get_cpp_context();
    void set_search_paths(const std::vector<std::string>& search_paths_in);
    void set_ca_bundle_path(const std::string& ca_bundle_path_in);

    static pj_ctx createDefault(pj_allocator *allocator);
};

/* Generate pj_list external or make list from include file */
#ifndef PJ_DATUMS__
C_NAMESPACE_VAR struct PJ_DATUMS pj_datums[];
#endif





#ifdef PJ_LIB__

#define OPERATION(name, NEED_ELLPS)                          \
                                                             \
pj_projection_specific_setup_##name (PJ *P);                 \
C_NAMESPACE PJ *pj_##name (PJ *P, PJ_CONTEXT *ctx);          \
                                                             \
C_NAMESPACE_VAR const char * const pj_s_##name = des_##name; \
                                                             \
C_NAMESPACE PJ *pj_##name (PJ *P, PJ_CONTEXT *ctx) {         \
    if (P)                                                   \
        return pj_projection_specific_setup_##name (P);      \
    P = pj_new(ctx);                                         \
    if (nullptr==P)                                          \
        return nullptr;                                      \
    P->host->short_name = #name;                             \
    P->host->descr = des_##name;                             \
    P->need_ellps = NEED_ELLPS;                              \
    P->left  = PJ_IO_UNITS_RADIANS;                          \
    P->right = PJ_IO_UNITS_CLASSIC;                          \
    P->host->file = __FILE__;                                \
    return P;                                                \
}                                                            \
                                                             \
PJ *pj_projection_specific_setup_##name (PJ *P)

/* In ISO19000 lingo, an operation is either a conversion or a transformation */
#define CONVERSION(name, need_ellps)      OPERATION (name, need_ellps)
#define TRANSFORMATION(name, need_ellps)  OPERATION (name, need_ellps)

/* In PROJ.4 a projection is a conversion taking angular input and giving scaled linear output */
#define PROJECTION(name) CONVERSION (name, 1)

#endif /* def PJ_LIB__ */

/* procedure prototypes */
double PROJ_DLL dmstor(const char *, char **);
double dmstor_ctx(PJ_CONTEXT *ctx, const char *, char **);
void   PROJ_DLL set_rtodms(int, int);
char  PROJ_DLL *rtodms(char *, double, int, int);
double aacos(pj_ctx_shared*,double);
double aasin(pj_ctx_shared *,double);
double asqrt(double);
double aatan2(double, double);

PROJVALUE PROJ_DLL pj_param(PJ_CONTEXT *ctx, paralist *, const char *);
paralist PROJ_DLL *pj_param_exists (paralist *list, const char *parameter);
paralist PROJ_DLL *pj_mkparam(const char *);
paralist *pj_mkparam_ws (const char *str, const char **next_str);


int PROJ_DLL pj_ell_set(PJ_CONTEXT *ctx, paralist *, double *, double *);
int pj_datum_set(PJ_CONTEXT *,paralist *, PJ *);
int pj_angular_units_set(paralist *, PJ *);

paralist *pj_clone_paralist( const paralist* );
paralist *pj_search_initcache( const char *filekey );
void      pj_insert_initcache( const char *filekey, const paralist *list);
paralist *pj_expand_init(PJ_CONTEXT *ctx, paralist *init);

void     *free_params (PJ_CONTEXT *ctx, paralist *start, int errlev);


double *pj_enfn(PJ_CONTEXT *, double);
void    pj_free_en(PJ_CONTEXT *, double* en);
double  pj_mlfn(double, double, double, const double *);
double  pj_inv_mlfn(pj_ctx_shared *, double, double, const double *);
double  pj_qsfn(double, double, double);
double  pj_tsfn(double, double, double);
double  pj_msfn(double, double, double);
double  pj_phi2(pj_ctx_shared *, const double, const double);
double  pj_sinhpsi2tanphi(pj_ctx_shared *, const double, const double);
double  pj_qsfn_(double, PJ *);
double *pj_authset(PJ_CONTEXT*, double);
void    pj_free_authset(PJ_CONTEXT*, double* apa);
double  pj_authlat(double, double *);

COMPLEX pj_zpoly1(COMPLEX, const COMPLEX *, int);
COMPLEX pj_zpolyd1(COMPLEX, const COMPLEX *, int, COMPLEX *);

int pj_deriv(PJ_LP, double, const PJ *, struct DERIVS *);
int pj_factors(PJ_LP, const PJ *, double, struct FACTORS *);

void  *proj_mdist_ini(double);
double proj_mdist(double, double, double, const void *);
double proj_inv_mdist(pj_ctx_shared *ctx, double, const void *);
void  *pj_gauss_ini(double, double, double *,double *);
PJ_LP     pj_gauss(pj_ctx_shared*, PJ_LP, const void *);
PJ_LP     pj_inv_gauss(pj_ctx_shared*, PJ_LP, const void *);

struct PJ_DATUMS           PROJ_DLL *pj_get_datums_ref( void );

PJ *pj_new(PJ_CONTEXT* ctx);
PJ *pj_default_destructor (PJ *P, int errlev);

double PROJ_DLL pj_atof( const char* nptr );
double pj_strtod( const char *nptr, char **endptr );
void   pj_freeup_plain (PJ *P);

PJ* pj_init_ctx_with_allow_init_epsg( PJ_CONTEXT *ctx, int argc, char **argv, int allow_init_epsg );

std::string PROJ_DLL pj_add_type_crs_if_needed(const std::string& str);
std::string pj_double_quote_string_param_if_needed(const std::string& str);

PJ *pj_create_internal (PJ_CONTEXT *ctx, const char *definition);
PJ *pj_create_argv_internal (PJ_CONTEXT *ctx, int argc, char **argv);

void pj_map_svm_ptrs(PJ* P, bool map);

// For use by projinfo
void pj_load_ini(PJ_CONTEXT* ctx);

// Exported for testing purposes only
std::string PROJ_DLL pj_context_get_grid_cache_filename(PJ_CONTEXT *ctx);

// For use by projsync
void PROJ_DLL pj_context_set_user_writable_directory(PJ_CONTEXT* ctx, const std::string& path);
std::string PROJ_DLL pj_get_relative_share_proj(PJ_CONTEXT *ctx);

std::vector<PJCoordOperation> pj_create_prepared_operations(PJ_CONTEXT *ctx,
                                                     const PJ *source_crs,
                                                     const PJ *target_crs,
                                                     PJ_OBJ_LIST* op_list);

int pj_get_suggested_operation(PJ_CONTEXT *ctx,
                               const std::vector<PJCoordOperation>& opList,
                               const int iExcluded[2],
                               PJ_DIRECTION direction,
                               PJ_COORD coord);

const PJ_UNITS *pj_list_linear_units();
const PJ_UNITS *pj_list_angular_units();

void pj_clear_hgridshift_knowngrids_cache();
void pj_clear_vgridshift_knowngrids_cache();

void pj_clear_sqlite_cache();

PJ_LP pj_generic_inverse_2d(PJ_XY xy, PJ *P, PJ_LP lpInitial);




/*****************************************************************************/
/*                                                                           */
/*                              proj_api.h                                   */
/*                                                                           */
/*    The rest of this header file includes what used to be "proj_api.h"     */
/*                                                                           */
/*****************************************************************************/

/* pj_init() and similar functions can be used with a non-C locale */
/* Can be detected too at runtime if the symbol pj_atof exists */
#define PJ_LOCALE_SAFE 1

#define RAD_TO_DEG    57.295779513082321
#define DEG_TO_RAD   .017453292519943296




extern char const PROJ_DLL pj_release[]; /* global release id string */

/* procedure prototypes */

PJ_CONTEXT PROJ_DLL *pj_get_default_ctx();
PJ_CONTEXT *pj_get_ctx( PJ *);

PJ_XY PROJ_DLL pj_fwd(PJ_LP, PJ *);
PJ_LP PROJ_DLL pj_inv(PJ_XY, PJ *);

PJ_XYZ pj_fwd3d(PJ_LPZ, PJ *);
PJ_LPZ pj_inv3d(PJ_XYZ, PJ *);


void pj_clear_initcache(void);
void PROJ_DLL pj_pr_list(PJ *); /* used by proj.cpp */
char *pj_get_def(PJ *, int);
int pj_has_inverse(PJ *);


char *pj_strdup(const char *str);
const char PROJ_DLL *pj_get_release(void);
void pj_acquire_lock(void);
void pj_release_lock(void);
void pj_cleanup_lock(void);

bool pj_log_active( PJ_CONTEXT *ctx, int level );
void pj_log( PJ_CONTEXT * ctx, int level, const char *fmt, ... );
void pj_stderr_logger( void *, int, const char * );

int pj_find_file(PJ_CONTEXT * ctx, const char *short_filename,
                 char* out_full_filename, size_t out_full_filename_size);

#endif /* ndef PROJ_INTERNAL_H */

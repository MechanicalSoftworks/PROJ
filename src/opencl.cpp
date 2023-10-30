#include "proj.h"
#include "proj_internal.h"

#include <sstream>

/*****************************************************************************/
void pj_scan_recursive(PJ* P, PJscan& s) {
/*****************************************************************************
    Finds all possible functions used by this projection.
******************************************************************************/
    pj_scan_local(P, s);

    if (P->axisswap) {
        P->axisswap->host->scan(P->axisswap, s);
    }

    if (P->cart) {
        P->cart->host->scan(P->cart, s);
    }

    if (P->cart_wgs84) {
        P->cart_wgs84->host->scan(P->cart_wgs84, s);
    }

    if (P->helmert) {
        P->helmert->host->scan(P->helmert, s);
    }

    if (P->hgridshift) {
        P->hgridshift->host->scan(P->hgridshift, s);
    }

    if (P->vgridshift) {
        P->vgridshift->host->scan(P->vgridshift, s);
    }
}

/*****************************************************************************/
void pj_scan_local(PJ* P, PJscan& s) {
/*****************************************************************************
    Finds all the functions used by this operation, and records them so we can
    populate the OpenCL double dispatch table for this projection.
******************************************************************************/
    P->co_fwd = s.add_co(P->host->co_fwd, P->host->file);
    P->co_inv = s.add_co(P->host->co_inv, P->host->file);
    P->co_fwd3d = s.add_co(P->host->co_fwd3d, P->host->file);
    P->co_inv3d = s.add_co(P->host->co_inv3d, P->host->file);
    P->co_fwd4d = s.add_co(P->host->co_fwd4d, P->host->file);
    P->co_inv4d = s.add_co(P->host->co_inv4d, P->host->file);

    P->fwd = s.add_fwd(P->host->fwd, P->host->file);
    P->inv = s.add_inv(P->host->inv, P->host->file);
    P->fwd3d = s.add_fwd(P->host->fwd3d, P->host->file);
    P->inv3d = s.add_inv(P->host->inv3d, P->host->file);
    P->fwd4d = s.add_fwd4d(P->host->fwd4d, P->host->file);
    P->inv4d = s.add_inv4d(P->host->inv4d, P->host->file);
}

/*****************************************************************************/
void pj_scan_nop(PJ* P, PJscan& s) {
/*****************************************************************************/
    (void) P;
    (void) s;
}

static std::string generate_header_list(PJscan& s) {
    std::ostringstream oss;

    for (const auto& p : s.files)
    {
        std::string path(p);
        const auto i = path.rfind("src");

        // We're scanning the result of the __FILE__ macro.
        // We should always be able to find 'src' in the build path because that's how this repo is set up!
        if (i == std::string::npos)
        {
            // This is here...just in case.
            continue;
        }

        oss << "#include <PROJ/" << path.substr(i) << ">" << std::endl;
    }

    return oss.str();
}

static std::string generate_co_dispatch_table(const PJfunction_to_id& fns) {
    std::ostringstream oss;

    oss << "PJcoroutine_code_t proj_dispatch_coroutine(PJ_COROUTINE_ID fn, cl_local struct PJstack_s* stack, cl_local struct PJstack_entry_s* e) { " << std::endl;
    oss << "    switch (fn) {" << std::endl;
    oss << "        default: break;" << std::endl;
    for (const auto& fn : fns)
    {
        oss << "        case " << fn.second << ": return " << fn.first << "(stack, e);" << std::endl;
    }
    oss << "    }" << std::endl;
    oss << "    return PJ_CO_ERROR;" << std::endl;
    oss << "}" << std::endl;
    oss << std::endl;

    return oss.str();
}

static std::string generate_dispatch_table(const PJfunction_to_id& fns, const char *suffix, const char* fn_type, const char *return_type, const char *parm_type) {
    std::ostringstream oss;

    oss << return_type << " proj_dispatch_" << suffix << "(" << fn_type << " fn, " << parm_type << " in, PJ *P) {" << std::endl;
    oss << "    switch (fn) {" << std::endl;
    oss << "        default: break;" << std::endl;
    for (const auto& fn : fns)
    {
        oss << "        case " << fn.second << ": return " << fn.first << "(in, P);" << std::endl;
    }
    oss << "    }" << std::endl;
    oss << "    " << return_type << " tmp = { 0 };" << std::endl;
    oss << "    return tmp;" << std::endl;
    oss << "}" << std::endl;
    oss << std::endl;

    return oss.str();
}

/*****************************************************************************/
std::string pj_create_opencl_source_from_scan(PJscan& s) {
/*****************************************************************************/
    std::ostringstream oss;

    oss << "#ifndef PROJ_DISPATCH" << std::endl;
    oss << "#define PROJ_DISPATCH" << std::endl;

    oss << "#define PJ_LIB__" << std::endl;

    oss << "#include <PROJ/src/proj_internal_device.h>" << std::endl;
    oss << generate_header_list(s) << std::endl;

    oss << generate_co_dispatch_table(s.co) << std::endl;
    oss << generate_dispatch_table(s.fwd, "fwd", "PJ_FWD_2D_ID", "PJ_XY", "PJ_LP") << std::endl;
    oss << generate_dispatch_table(s.inv, "inv", "PJ_INV_2D_ID", "PJ_LP", "PJ_XY") << std::endl;
    oss << generate_dispatch_table(s.fwd3d, "fwd3d", "PJ_FWD_3D_ID", "PJ_XYZ", "PJ_LPZ") << std::endl;
    oss << generate_dispatch_table(s.inv3d, "inv3d", "PJ_INV_3D_ID", "PJ_LPZ", "PJ_XYZ") << std::endl;
    oss << generate_dispatch_table(s.fwd4d, "fwd4d", "PJ_FWD_4D_ID", "PJ_COORD", "PJ_COORD") << std::endl;
    oss << generate_dispatch_table(s.inv4d, "inv4d", "PJ_INV_4D_ID", "PJ_COORD", "PJ_COORD") << std::endl;
    
    oss << "#include <PROJ/src/shared.cpp>" << std::endl;

    oss << "#endif // !PROJ_DISPATCH" << std::endl;

    return oss.str();
}

std::string PJscan::definitions_for(const PJfunction_to_id& m)
{
    std::ostringstream oss;

    for (const auto& kv : m)
    {
        oss << "-D" << kv.first << "_id=" << kv.second << " ";
    }

    return oss.str();
}

/*****************************************************************************/
std::string pj_create_opencl_definitions_from_scan(PJscan& s) {
/*****************************************************************************/
    std::ostringstream oss;

    oss << PJscan::definitions_for(s.co)
        << PJscan::definitions_for(s.fwd)
        << PJscan::definitions_for(s.inv)
        << PJscan::definitions_for(s.fwd3d)
        << PJscan::definitions_for(s.inv3d)
        << PJscan::definitions_for(s.fwd4d)
        << PJscan::definitions_for(s.inv4d);

    return oss.str();
}

extern void pj_scan_fwd(PJscan& s);
extern void pj_scan_inv(PJscan& s);

/*****************************************************************************/
PJ_SCAN PROJ_DLL* proj_create_scan() {
/*****************************************************************************/
    return new PJscan;
}

/*****************************************************************************/
void PROJ_DLL proj_scan(PJ_SCAN* s, PJ* P) {
/*****************************************************************************/
    P->host->scan(P, *s);
    pj_scan_fwd(*s);
    pj_scan_inv(*s);
}

/*****************************************************************************/
void PROJ_DLL proj_destroy_scan(PJ_SCAN* scan) {
/*****************************************************************************/
    delete scan;
}

/*****************************************************************************/
int proj_create_opencl(PJ_SCAN* s, PJopencl* o) {
/*****************************************************************************/
    const auto src = pj_create_opencl_source_from_scan(*s);
    const auto defs = pj_create_opencl_definitions_from_scan(*s);

    o->kernel_source = (char*)malloc(src.length() + 1);
    o->kernel_definitions = (char*)malloc(defs.length() + 1);

    if (nullptr==o->kernel_source || nullptr==o->kernel_definitions) {
        free(o->kernel_source);
        free(o->kernel_definitions);
        return false;
    }

    memcpy(o->kernel_source, src.data(), src.length());
    o->kernel_source[src.length()] = 0;

    memcpy(o->kernel_definitions, defs.data(), defs.length());
    o->kernel_definitions[defs.length()] = 0;

    return true;
}

/*****************************************************************************/
void proj_free_opencl(PJopencl* o) {
/*****************************************************************************/
    if (o)
    {
        free(o->kernel_source);
        free(o->kernel_definitions);
    }
}

void pj_map_svm_ptrs(PJ* P, bool map)
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

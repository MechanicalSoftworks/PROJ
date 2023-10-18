#include "proj.h"
#include "proj_internal.h"

#include <sstream>

//
// TODO: Probably best to commit here. Next change is very risky.
//          * Repurpose pj_scan_recursive to build a list of functions that need to be executed.
//          * This eliminates the need for the indirection table.
//              * And introduces the need for an indirection queue!
//              * Or maybe this indirection queue is a specifically generated version of proj_trans?
//                  * Each PJ gets its own proj_trans_[name]. Fixes recursion, but increases register pressue?
//          * Probably need dedicated "scan" versions of pj_fwd and pj_inv.
//              * Attach a scan context to P->host->ctx?
//              * Presense of a scan context means
//

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
    bool add_file = false;

    if (P->fwd[0])
    {
        s.fwd.insert(P->fwd);
        add_file = true;
    }

    if (P->inv[0])
    {
        s.inv.insert(P->inv);
        add_file = true;
    }

    if (P->fwd3d[0])
    {
        s.fwd3d.insert(P->fwd3d);
        add_file = true;
    }

    if (P->inv3d[0])
    {
        s.inv3d.insert(P->inv3d);
        add_file = true;
    }

    if (P->fwd4d[0])
    {
        s.fwd4d.insert(P->fwd4d);
        add_file = true;
    }

    if (P->inv4d[0])
    {
        s.inv4d.insert(P->inv4d);
        add_file = true;
    }

    if (add_file)
    {
        s.files.insert(P->host->file);
    }
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

static std::string generate_dispatch_table(std::set<const char*>& fns, const char *suffix, const char *return_type, const char *parm_type) {
    std::ostringstream oss;

    oss << return_type << " pj_double_dispatch_" << suffix << "(" << parm_type << " in, PJ *P, const char *name) {" << std::endl;
    for (const auto& fn : fns)
    {
        oss << "    if (pj_streq(name, \"" << fn << "\")) {" << std::endl;
        oss << "        return " << fn << "(in, P);" << std::endl;
        oss << "    }" << std::endl;
    }
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

    oss << "#define PJ_LIB__" << std::endl;
    oss << "#include <PROJ/src/proj_internal_device.h>" << std::endl;
    oss << generate_header_list(s) << std::endl;

    oss << "#ifndef PROJ_DOUBLE_DISPATCH" << std::endl;
    oss << "#define PROJ_DOUBLE_DISPATCH" << std::endl;
    oss << generate_dispatch_table(s.fwd, "fwd", "PJ_XY", "PJ_LP") << std::endl;
    oss << generate_dispatch_table(s.inv, "inv", "PJ_LP", "PJ_XY") << std::endl;
    oss << generate_dispatch_table(s.fwd3d, "fwd3d", "PJ_XYZ", "PJ_LPZ") << std::endl;
    oss << generate_dispatch_table(s.inv3d, "inv3d", "PJ_LPZ", "PJ_XYZ") << std::endl;
    oss << generate_dispatch_table(s.fwd4d, "fwd4d", "PJ_COORD", "PJ_COORD") << std::endl;
    oss << generate_dispatch_table(s.inv4d, "inv4d", "PJ_COORD", "PJ_COORD") << std::endl;
    oss << "#endif // !PROJ_DOUBLE_DISPATCH" << std::endl;

    oss << "#include <PROJ/src/fwd.cpp>" << std::endl;
    oss << "#include <PROJ/src/inv.cpp>" << std::endl;
    oss << "#include <PROJ/src/shared.cpp>" << std::endl;

    return oss.str();
}

/*****************************************************************************/
char* proj_create_opencl_source(PJ* P) {
/*****************************************************************************/
    PJscan s;

    P->host->scan(P, s);

    const auto str = pj_create_opencl_source_from_scan(s);
    char* src = (char*)malloc(str.length() + 1);
    if (nullptr==src) {
        return nullptr;
    }

    memcpy(src, str.data(), str.length());
    src[str.length()] = 0;
    return src;
}

/*****************************************************************************/
void proj_free_opencl_source(char* src) {
/*****************************************************************************/
    free(src);
}

void pj_map_svm_ptrs(PJ* P, bool map)
{
    P->host->ctx->allocator->svm_map(P, map);
    P->host->ctx->allocator->svm_map(P->shared_ctx, map);

    // Nice little circular reference here. Enforce the constraint that this operation needs to happen top-down.
    //if (P->parent) {
    //    P->parent->host->map_svm(P->parent, map);
    //}

    if (P->geod) {
        P->host->ctx->allocator->svm_map(P->geod, map);
    }

    if (P->opaque) {
        P->host->ctx->allocator->svm_map(P->opaque, map);
    }

    if (P->axisswap) {
        P->axisswap->host->map_svm(P->axisswap, map);
    }

    if (P->cart) {
        P->cart->host->map_svm(P->cart, map);
    }

    if (P->cart_wgs84) {
        P->cart_wgs84->host->map_svm(P->cart_wgs84, map);
    }

    if (P->helmert) {
        P->helmert->host->map_svm(P->helmert, map);
    }

    if (P->hgridshift) {
        P->hgridshift->host->map_svm(P->hgridshift, map);
    }

    if (P->vgridshift) {
        P->vgridshift->host->map_svm(P->vgridshift, map);
    }
}

void proj_host_acquire_svm(PJ* P) {
    P->host->map_svm(P, true);
}

void proj_host_release_svm(PJ* P) {
    P->host->map_svm(P, false);
}

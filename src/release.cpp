/* <<< Release Notice for library >>> */

#include "proj.h"
#include "proj_internal.h"

char const pj_release[] =
    "Rel. "
    PROJ_STR(PROJ_VERSION_MAJOR)"."
    PROJ_STR(PROJ_VERSION_MINOR)"."
    PROJ_STR(PROJ_VERSION_PATCH)", "
    "January 1st, 2022";

const char *pj_get_release() {
    return pj_release;
}

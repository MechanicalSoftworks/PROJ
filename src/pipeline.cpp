/*******************************************************************************

                       Transformation pipeline manager

                    Thomas Knudsen, 2016-05-20/2016-11-20

********************************************************************************

    Geodetic transformations are typically organized in a number of
    steps. For example, a datum shift could be carried out through
    these steps:

    1. Convert (latitude, longitude, ellipsoidal height) to
       3D geocentric cartesian coordinates (X, Y, Z)
    2. Transform the (X, Y, Z) coordinates to the new datum, using a
       7 parameter Helmert transformation.
    3. Convert (X, Y, Z) back to (latitude, longitude, ellipsoidal height)

    If the height system used is orthometric, rather than ellipsoidal,
    another step is needed at each end of the process:

    1. Add the local geoid undulation (N) to the orthometric height
       to obtain the ellipsoidal (i.e. geometric) height.
    2. Convert (latitude, longitude, ellipsoidal height) to
       3D geocentric cartesian coordinates (X, Y, Z)
    3. Transform the (X, Y, Z) coordinates to the new datum, using a
       7 parameter Helmert transformation.
    4. Convert (X, Y, Z) back to (latitude, longitude, ellipsoidal height)
    5. Subtract the local geoid undulation (N) from the ellipsoidal height
       to obtain the orthometric height.

    Additional steps can be added for e.g. change of vertical datum, so the
    list can grow fairly long. None of the steps are, however, particularly
    complex, and data flow is strictly from top to bottom.

    Hence, in principle, the first example above could be implemented using
    Unix pipelines:

    cat my_coordinates | geographic_to_xyz | helmert | xyz_to_geographic > my_transformed_coordinates

    in the grand tradition of Software Tools [1].

    The proj pipeline driver implements a similar concept: Stringing together
    a number of steps, feeding the output of one step to the input of the next.

    It is a very powerful concept, that increases the range of relevance of the
    proj.4 system substantially. It is, however, not a particularly intrusive
    addition to the PROJ.4 code base: The implementation is by and large completed
    by adding an extra projection called "pipeline" (i.e. this file), which
    handles all business, and a small amount of added functionality in the
    pj_init code, implementing support for multilevel, embedded pipelines.

    Syntactically, the pipeline system introduces the "+step" keyword (which
    indicates the start of each transformation step), and reintroduces the +inv
    keyword (indicating that a given transformation step should run in reverse, i.e.
    forward, when the pipeline is executed in inverse direction, and vice versa).

    Hence, the first transformation example above, can be implemented as:

    +proj=pipeline +step proj=cart +step proj=helmert <ARGS> +step proj=cart +inv

    Where <ARGS> indicate the Helmert arguments: 3 translations (+x=..., +y=...,
    +z=...), 3 rotations (+rx=..., +ry=..., +rz=...) and a scale factor (+s=...).
    Following geodetic conventions, the rotations are given in arcseconds,
    and the scale factor is given as parts-per-million.

    [1] B. W. Kernighan & P. J. Plauger: Software tools.
        Reading, Massachusetts, Addison-Wesley, 1976, 338 pp.

********************************************************************************

Thomas Knudsen, thokn@sdfe.dk, 2016-05-20

********************************************************************************
* Copyright (c) 2016, 2017, 2018 Thomas Knudsen / SDFE
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
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
* THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
* FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
********************************************************************************/

#define PJ_LIB__

#ifndef PROJ_OPENCL_DEVICE
#   include "geodesic.h"
#endif

#include "proj_kernel.h"

PROJ_HEAD(pipeline,         "Transformation pipeline manager");
PROJ_HEAD(pop, "Retrieve coordinate value from pipeline stack");
PROJ_HEAD(push, "Save coordinate value on pipeline stack");

/* Projection specific elements for the PJ object */
namespace { // anonymous namespace
struct PipelineStep {
    PJ_FIELD(__global PJ*, pj,       nullptr);
    PJ_FIELD(        bool, omit_fwd, false);
    PJ_FIELD(        bool, omit_inv, false);

#ifdef __cplusplus
    PipelineStep(__global PJ* pjIn, bool omitFwdIn, bool omitInvIn) :
        pj(pjIn), omit_fwd(omitFwdIn), omit_inv(omitInvIn) {}
#endif
};

struct PipelineStack {
    PJ_FIELD(__global double*, ptr,      nullptr);
    PJ_FIELD(          size_t, size,     0);
    PJ_FIELD(          size_t, max_size, 0);
};

struct Pipeline {
    PJ_FIELD(char**, argv,         nullptr);
    PJ_FIELD(char**, current_argv, nullptr);

    PJ_FIELD(__global struct PipelineStep*, steps,      nullptr);
    PJ_FIELD(                       size_t, step_count, 0);

    PipelineStack stack1;
    PipelineStack stack2;
    PipelineStack stack3;
    PipelineStack stack4;
};

struct PushPop {
    bool v1;
    bool v2;
    bool v3;
    bool v4;
};
} // anonymous namespace

#ifndef PROJ_OPENCL_DEVICE

static void PipelineStackNew(struct pj_allocator* allocator, struct PipelineStack* stack, size_t n)
{
    stack->size = 0;
    stack->max_size = n;
    stack->ptr = (double*)allocator->svm_calloc(n, sizeof(stack->ptr[0]));
}

static void PipelineStackDelete(struct pj_allocator* allocator, struct PipelineStack* stack)
{
    stack->size = 0;
    stack->max_size = 0;
    allocator->svm_free(stack->ptr);
    stack->ptr = nullptr;
}

#endif

static double PipelineStackTop(__global struct PipelineStack* stack)
{
    return stack->size
        ? stack->ptr[stack->size - 1]
        : 0.0;
}

static int PipelineStackEmpty(__global struct PipelineStack* stack)
{
    return !stack->size;
}

static void PipelineStackPush(__global struct PipelineStack* stack, double d)
{
    if (stack->size < stack->max_size)
    {
        stack->ptr[stack->size++] = d;
    }
}

static void PipelineStackPop(__global struct PipelineStack* stack)
{
    if (stack->size)
    {
        --stack->size;
    }
}


#ifndef PROJ_OPENCL_DEVICE

static void pipeline_reassign_context( PJ* P, PJ_CONTEXT* ctx )
{
    auto pipeline = static_cast<struct Pipeline*>(P->opaque);
    for (size_t i = 0; i < pipeline->step_count; ++i) {
        auto& step = pipeline->steps[i];
        proj_assign_context(step.pj, ctx);
    }
}

static void pipeline_map_svm_ptrs(PJ* P, bool map)
{
    pj_map_svm_ptrs(P, map);

    struct Pipeline* pipeline = (struct Pipeline*)(P->opaque);
    P->host->map_svm(pipeline->steps, map);
    P->host->map_svm(&pipeline->stack1, map);
    for (size_t i = 0; i < pipeline->step_count; ++i)
    {
        struct PipelineStep* step = pipeline->steps + i;
        if (!step->omit_fwd || !step->omit_inv)
        {
            step->pj->host->map_pj(step->pj, map);
        }
    }
}

#endif

PJcoroutine_code_t pipeline_forward_4d_co (__local PJstack_t* stack, __local PJstack_entry_t* e) {
    PJ*                     P = e->P;
    __global Pipeline*      pipeline = static_cast<__global struct Pipeline*>(P->opaque);
    PJ_COORD                point = e->coo;
    size_t                  i = e->i;
    __global PipelineStep*  step = nullptr;

    switch (e->step) {
        case 0: break;
        case 1: goto p1;
        default: goto ABORT;
    }

    i = 0;

    for(; i < pipeline->step_count; ++i)
    {
        step = pipeline->steps + i;
        if( !step->omit_fwd )
        {
            push_proj_trans (stack, step->pj, PJ_FWD, point);
            PJ_YIELD(e, 1);
            if( point.xyzt.x == HUGE_VAL ) {
                break;
            }
        }
    }

//DONE:
    e->coo = point;
    e->i = (int)i;
    return PJ_CO_DONE;

YIELD:
    e->coo = point;
    e->i = (int)i;
    return PJ_CO_YIELD;

ABORT:
    e->coo = proj_coord_error();
    return PJ_CO_ERROR;
}


PJcoroutine_code_t pipeline_reverse_4d_co(__local PJstack_t* stack, __local PJstack_entry_t* e) {
    PJ*                     P = e->P;
    auto                    pipeline = static_cast<struct Pipeline*>(P->opaque);
    PJ_COORD                point = e->coo;
    size_t                  i = e->i;
    struct PipelineStep*    step = nullptr;

    switch (e->step) {
        case 0: break;
        case 1: goto p1;
        default: goto ABORT;
    }

    i = pipeline->step_count;

    for(; i > 0; --i)
    {
        step = pipeline->steps + i - 1;
        if( !step->omit_inv )
        {
            push_proj_trans (stack, step->pj, PJ_INV, point);
            PJ_YIELD(e, 1);
            if( point.xyzt.x == HUGE_VAL ) {
                break;
            }
        }
    }

//DONE:
    e->coo = point;
    e->i = (int)i;
    return PJ_CO_DONE;

YIELD:
    e->coo = point;
    e->i = (int)i;
    return PJ_CO_YIELD;

ABORT:
    e->coo = proj_coord_error();
    return PJ_CO_ERROR;
}




PJcoroutine_code_t pipeline_forward_3d_co(__local PJstack_t* stack, __local PJstack_entry_t* e) {
    PJ*                     P = e->P;
    auto                    pipeline = static_cast<struct Pipeline*>(P->opaque);
    PJ_COORD                point = e->coo;
    size_t                  i = e->i;
    struct PipelineStep*    step = nullptr;

    switch (e->step) {
        case 0: break;
        case 1: goto p1;
        default: goto ABORT;
    }

    i = 0;
    point.lpzt.t = 0;

    for(; i < pipeline->step_count; ++i)
    {
        step = pipeline->steps + i;
        if( !step->omit_fwd )
        {
            push_approx_3D_trans (stack, step->pj, PJ_FWD, point);
            PJ_YIELD(e, 1);
            if( point.xyzt.x == HUGE_VAL ) {
                break;
            }
        }
    }

//DONE:
    e->coo = point;
    e->i = (int)i;
    return PJ_CO_DONE;

YIELD:
    e->coo = point;
    e->i = (int)i;
    return PJ_CO_YIELD;

ABORT:
    e->coo = proj_coord_error();
    return PJ_CO_ERROR;
}


PJcoroutine_code_t pipeline_reverse_3d_co(__local PJstack_t* stack, __local PJstack_entry_t* e) {
    PJ*                     P = e->P;
    auto                    pipeline = static_cast<struct Pipeline*>(P->opaque);
    PJ_COORD                point = e->coo;
    size_t                  i = e->i;
    struct PipelineStep*    step = nullptr;

    switch (e->step) {
        case 0: break;
        case 1: goto p1;
        default: goto ABORT;
    }

    i = pipeline->step_count;
    point.xyzt.t = 0;

    for(; i > 0; --i)
    {
        step = pipeline->steps + i - 1;
        if( !step->omit_inv )
        {
            push_proj_trans (stack, step->pj, PJ_INV, point);
            PJ_YIELD(e, 1);
            if( point.xyzt.x == HUGE_VAL ) {
                break;
            }
        }
    }

//DONE:
    e->coo = point;
    e->i = (int)i;
    return PJ_CO_DONE;

YIELD:
    e->coo = point;
    e->i = (int)i;
    return PJ_CO_YIELD;

ABORT:
    e->coo = proj_coord_error();
    return PJ_CO_ERROR;
}




PJcoroutine_code_t pipeline_forward_co(__local PJstack_t* stack, __local PJstack_entry_t* e) {
    PJ*                     P = e->P;
    auto                    pipeline = static_cast<struct Pipeline*>(P->opaque);
    PJ_COORD                point = e->coo;
    size_t                  i = e->i;
    struct PipelineStep*    step = nullptr;

    switch (e->step) {
    case 0: break;
    case 1: goto p1;
    default: goto ABORT;
    }

    i = 0;
    point.lpzt.z = point.lpzt.t = 0;

    for (; i < pipeline->step_count; ++i)
    {
        step = pipeline->steps + i;
        if( !step->omit_fwd )
        {
            push_approx_2D_trans (stack, step->pj, PJ_FWD, point);
            PJ_YIELD(e, 1);
            if( point.xyzt.x == HUGE_VAL ) {
                break;
            }
        }
    }

//DONE:
    e->coo = point;
    e->i = (int)i;
    return PJ_CO_DONE;

YIELD:
    e->coo = point;
    e->i = (int)i;
    return PJ_CO_YIELD;

ABORT:
    e->coo = proj_coord_error();
    return PJ_CO_ERROR;
}


PJcoroutine_code_t pipeline_reverse_co(__local PJstack_t* stack, __local PJstack_entry_t* e) {
    PJ*                     P = e->P;
    auto                    pipeline = static_cast<struct Pipeline*>(P->opaque);
    PJ_COORD                point = e->coo;
    size_t                  i = e->i;
    struct PipelineStep*    step = nullptr;

    switch (e->step) {
        case 0: break;
        case 1: goto p1;
        default: goto ABORT;
    }

    i = pipeline->step_count;
    point.xyz.z = point.xyzt.t = 0;

    for(; i > 0; --i)
    {
        step = pipeline->steps + i - 1;
        if( !step->omit_inv )
        {
            push_approx_2D_trans (stack, step->pj, PJ_INV, point);
            PJ_YIELD(e, 1);
            if( point.xyzt.x == HUGE_VAL ) {
                break;
            }
        }
    }

//DONE:
    e->coo = point;
    e->i = (int)i;
    return PJ_CO_DONE;

YIELD:
    e->coo = point;
    e->i = (int)i;
    return PJ_CO_YIELD;

ABORT:
    e->coo = proj_coord_error();
    return PJ_CO_ERROR;
}


#ifndef PROJ_OPENCL_DEVICE

static PJ *destructor (PJ *P, int errlev) {
    if (nullptr==P)
        return nullptr;

    if (nullptr==P->opaque)
        return pj_default_destructor (P, errlev);

    auto pipeline = static_cast<struct Pipeline*>(P->opaque);
    auto stack = &pipeline->stack1;

    free (pipeline->argv);
    free (pipeline->current_argv);

    for (size_t i = 0; i < pipeline->step_count; ++i)
    {
        proj_destroy(pipeline->steps[i].pj);
    }
    P->host->ctx->allocator->svm_free(pipeline->steps);

    for (int i = 0; i < 4; ++i)
    {
        PipelineStackDelete(P->host->ctx->allocator, stack + i);
    }

    for (int i = 0; i < 4; ++i)
    {
        P->host->ctx->allocator->svm_free(stack[i].ptr);
    }

    P->host->ctx->allocator->svm_delete(pipeline);
    P->opaque = nullptr;

    return pj_default_destructor(P, errlev);
}


/* count the number of args in pipeline definition, and mark all args as used */
static size_t argc_params (paralist *params) {
    size_t argc = 0;
    for (; params != nullptr; params = params->next) {
        argc++;
        params->used = 1;
    }
    return ++argc;  /* one extra for the sentinel */
}

/* Sentinel for argument list */
static const char *argv_sentinel = "step";

/* turn paralist into argc/argv style argument list */
static char **argv_params (paralist *params, size_t argc) {
    char **argv;
    size_t i = 0;
    argv = static_cast<char**>(calloc (argc, sizeof (char *)));
    if (nullptr==argv)
        return nullptr;
    for (; params != nullptr; params = params->next)
        argv[i++] = params->param;
    argv[i++] = const_cast<char*>(argv_sentinel);
    return argv;
}




/* Being the special operator that the pipeline is, we have to handle the    */
/* ellipsoid differently than usual. In general, the pipeline operation does */
/* not need an ellipsoid, but in some cases it is beneficial nonetheless.    */
/* Unfortunately we can't use the normal ellipsoid setter in pj_init, since  */
/* it adds a +ellps parameter to the global args if nothing else is specified*/
/* This is problematic since that ellipsoid spec is then passed on to the    */
/* pipeline children. This is rarely what we want, so here we implement our  */
/* own logic instead. If an ellipsoid is set in the global args, it is used  */
/* as the pipeline ellipsoid. Otherwise we use GRS80 parameters as default.  */
/* At last we calculate the rest of the ellipsoid parameters and             */
/* re-initialize P->geod.                                                    */
static void set_ellipsoid(PJ *P) {
    paralist *cur, *attachment;
    int err = proj_errno_reset (P);

    /* Break the linked list after the global args */
    attachment = nullptr;
    for (cur = P->host->params; cur != nullptr; cur = cur->next)
        /* cur->next will always be non 0 given argv_sentinel presence, */
        /* but this is far from being obvious for a static analyzer */
        if (cur->next != nullptr && strcmp(argv_sentinel, cur->next->param) == 0) {
            attachment = cur->next;
            cur->next = nullptr;
            break;
        }

    /* Check if there's any ellipsoid specification in the global params. */
    /* If not, use GRS80 as default                                       */
    if (0 != pj_ellipsoid (P)) {
        P->a  = 6378137.0;
        P->f = 1.0 / 298.257222101;
        P->es = 2*P->f - P->f*P->f;

        /* reset an "unerror": In this special use case, the errno is    */
        /* not an error signal, but just a reply from pj_ellipsoid,      */
        /* telling us that "No - there was no ellipsoid definition in    */
        /* the PJ you provided".                                         */
        proj_errno_reset (P);
    }
    P->a_orig = P->a;
    P->es_orig = P->es;

    pj_calc_ellipsoid_params (P, P->a, P->es);

    geod_init(P->geod, P->a,  P->es / (1 + sqrt(P->one_es)));

    /* Re-attach the dangling list */
    /* Note: cur will always be non 0 given argv_sentinel presence, */
    /* but this is far from being obvious for a static analyzer */
    if( cur != nullptr )
        cur->next = attachment;
    proj_errno_restore (P, err);
}


PJ *OPERATION(pipeline,0) {
    int i, nsteps = 0, argc;
    int i_pipeline = -1, i_first_step = -1, i_current_step;
    char **argv, **current_argv;

    if( P->host->ctx->pipelineInitRecursiongCounter == 5 )
    {
        // Can happen for a string like:
        // proj=pipeline step "x="""," u=" proj=pipeline step ste=""[" u=" proj=pipeline step ste="[" u=" proj=pipeline step ste="[" u=" proj=pipeline step ste="[" u=" proj=pipeline step ste="[" u=" proj=pipeline step ste="[" u=" proj=pipeline step ste="[" u=" proj=pipeline step ste="[" u=" proj=pipeline p step ste="[" u=" proj=pipeline step ste="[" u=" proj=pipeline step ste="[" u=" proj=pipeline step ste="[" u=" proj=pipeline step ""x="""""""""""
        // Probably an issue with the quoting handling code
        // But doesn't hurt to add an extra safety check
        proj_log_error (P, _("Pipeline: too deep recursion"));
        return destructor (P, PROJ_ERR_INVALID_OP_WRONG_SYNTAX); /* ERROR: nested pipelines */
    }

    P->co_fwd4d  =  PJ_MAKE_KERNEL(pipeline_forward_4d_co);
    P->co_inv4d  =  PJ_MAKE_KERNEL(pipeline_reverse_4d_co);
    P->co_fwd3d  =  PJ_MAKE_KERNEL(pipeline_forward_3d_co);
    P->co_inv3d  =  PJ_MAKE_KERNEL(pipeline_reverse_3d_co);
    P->co_fwd    =  PJ_MAKE_KERNEL(pipeline_forward_co);
    P->co_inv    =  PJ_MAKE_KERNEL(pipeline_reverse_co);
    P->host->destructor  =  destructor;
    P->host->reassign_context = pipeline_reassign_context;
    P->host->map_pj = pipeline_map_svm_ptrs;

    /* Currently, the pipeline driver is a raw bit mover, enabling other operations */
    /* to collaborate efficiently. All prep/fin stuff is done at the step levels. */
    P->skip_fwd_prepare  = 1;
    P->skip_fwd_finalize = 1;
    P->skip_inv_prepare  = 1;
    P->skip_inv_finalize = 1;


    P->opaque  = P->host->ctx->allocator->svm_new<Pipeline>();
    if (nullptr==P->opaque)
        return destructor(P, PROJ_ERR_INVALID_OP /* ENOMEM */);

    argc = (int)argc_params (P->host->params);
    auto pipeline = static_cast<struct Pipeline*>(P->opaque);
    auto stack = &pipeline->stack1;
    pipeline->argv = argv = argv_params (P->host->params, argc);
    if (nullptr==argv)
        return destructor (P, PROJ_ERR_INVALID_OP /* ENOMEM */);

    pipeline->current_argv = current_argv = static_cast<char**>(calloc (argc, sizeof (char *)));
    if (nullptr==current_argv)
        return destructor (P, PROJ_ERR_OTHER /*ENOMEM*/);

    /* Do some syntactical sanity checking */
    for (i = 0;  i < argc && argv[i] != nullptr;  i++) {
        if ( 0==strcmp (argv_sentinel, argv[i])) {
            if (-1==i_pipeline) {
                proj_log_error (P, _("Pipeline: +step before +proj=pipeline"));
                return destructor (P, PROJ_ERR_INVALID_OP_WRONG_SYNTAX);
            }
            if (0==nsteps)
                i_first_step = i;
            nsteps++;
            continue;
        }

        if (0==strcmp ("proj=pipeline", argv[i])) {
            if (-1 != i_pipeline) {
                proj_log_error (P, _("Pipeline: Nesting only allowed when child pipelines are wrapped in '+init's"));
                return destructor (P, PROJ_ERR_INVALID_OP_WRONG_SYNTAX); /* ERROR: nested pipelines */
            }
            i_pipeline = i;
        }
    }
    nsteps--; /* Last instance of +step is just a sentinel */

    if (-1==i_pipeline)
        return destructor (P, PROJ_ERR_INVALID_OP_WRONG_SYNTAX); /* ERROR: no pipeline def */

    if (0==nsteps)
        return destructor (P, PROJ_ERR_INVALID_OP_WRONG_SYNTAX); /* ERROR: no pipeline def */

    set_ellipsoid(P);

    /* Now loop over all steps, building a new set of arguments for each init */
    i_current_step = i_first_step;
    std::vector<PipelineStep> steps;
    for (i = 0;  i < nsteps;  i++) {
        int j;
        int  current_argc = 0;
        int  err;
        PJ     *next_step = nullptr;

        /* Build a set of setup args for the current step */
        proj_log_trace (P, "Pipeline: Building arg list for step no. %d", i);

        /* First add the step specific args */
        for (j = i_current_step + 1;  0 != strcmp ("step", argv[j]); j++)
            current_argv[current_argc++] = argv[j];

        i_current_step = j;

        /* Then add the global args */
        for (j = i_pipeline + 1;  0 != strcmp ("step", argv[j]); j++)
            current_argv[current_argc++] = argv[j];

        proj_log_trace (P, "Pipeline: init - %s, %d", current_argv[0], current_argc);
        for (j = 1;  j < current_argc; j++)
            proj_log_trace (P, "    %s", current_argv[j]);

        err = proj_errno_reset (P);

        P->host->ctx->pipelineInitRecursiongCounter ++;
        next_step = pj_create_argv_internal (P->host->ctx, current_argc, current_argv);
        P->host->ctx->pipelineInitRecursiongCounter --;
        proj_log_trace (P, "Pipeline: PipelineStep %d (%s) at %p", i, current_argv[0], next_step);

        if (nullptr==next_step) {
            /* The step init failed, but possibly without setting errno. If so, we say "malformed" */
            int err_to_report = proj_errno(P);
            if (0==err_to_report)
                err_to_report = PROJ_ERR_INVALID_OP_WRONG_SYNTAX;
            proj_log_error (P, _("Pipeline: Bad step definition: %s (%s)"), current_argv[0], proj_context_errno_string (P->host->ctx, err_to_report));
            return destructor (P, err_to_report); /* ERROR: bad pipeline def */
        }
        next_step->parent = P;

        proj_errno_restore (P, err);

        /* Is this step inverted? */
        for (j = 0;  j < current_argc; j++) {
            if (0==strcmp("inv", current_argv[j])) {
                /* if +inv exists in both global and local args the forward operation should be used */
                next_step->inverted = next_step->inverted == 0 ? 1 : 0;
            }
        }

        bool omit_fwd = pj_param(P->host->ctx, next_step->host->params, "bomit_fwd").i != 0;
        bool omit_inv = pj_param(P->host->ctx, next_step->host->params, "bomit_inv").i != 0;
        steps.emplace_back(next_step, omit_fwd, omit_inv);

        proj_log_trace (P, "Pipeline at [%p]:    step at [%p] (%s) done", P, next_step, current_argv[0]);
    }

    /* Require a forward path through the pipeline */
    for( auto& step: steps) {
        PJ *Q = step.pj;
        if ( ( Q->inverted && (Q->inv || Q->inv3d || Q->fwd4d) ) ||
             (!Q->inverted && (Q->fwd || Q->fwd3d || Q->fwd4d) ) ) {
            continue;
        } else {
            proj_log_error (P, _("Pipeline: A forward operation couldn't be constructed"));
            return destructor (P, PROJ_ERR_INVALID_OP_WRONG_SYNTAX);
        }
    }

    /* determine if an inverse operation is possible */
    for( auto& step: steps) {
        PJ *Q = step.pj;
        if ( pj_has_inverse(Q) ) {
            continue;
        } else {
            P->inv   = nullptr;
            P->inv3d = nullptr;
            P->inv4d = nullptr;
            break;
        }
    }


    /* Replace PJ_IO_UNITS_WHATEVER with input/output units of neighbouring steps where */
    /* it make sense. It does in most cases but not always, for instance                */
    /*      proj=pipeline step proj=unitconvert xy_in=deg xy_out=rad step ...           */
    /* where the left-hand side units of the first step shouldn't be changed to RADIANS */
    /* as it will result in deg->rad conversions in cs2cs and other applications.       */

    for (i=nsteps-2; i>=0; --i) {
        auto pj = steps[i].pj;
        if (pj_left(pj) == PJ_IO_UNITS_WHATEVER && pj_right(pj) == PJ_IO_UNITS_WHATEVER) {
            const auto right_pj = steps[i+1].pj;
            const auto right_pj_left = pj_left(right_pj);
            const auto right_pj_right = pj_right(right_pj);
            if (right_pj_left != right_pj_right || right_pj_left != PJ_IO_UNITS_WHATEVER )
            {
                pj->left = right_pj_left;
                pj->right = right_pj_left;
            }
        }
    }

    for (i=1; i<nsteps; i++) {
        auto pj = steps[i].pj;
        if (pj_left(pj) == PJ_IO_UNITS_WHATEVER && pj_right(pj) == PJ_IO_UNITS_WHATEVER) {
            const auto left_pj = steps[i-1].pj;
            const auto left_pj_left = pj_left(left_pj);
            const auto left_pj_right = pj_right(left_pj);
            if (left_pj_left != left_pj_right || left_pj_right != PJ_IO_UNITS_WHATEVER )
            {
                pj->left = left_pj_right;
                pj->right = left_pj_right;
            }
        }
    }

    /* Check that units between each steps match each other, fail if they don't */
    for (i = 0; i + 1 < nsteps; i++) {
        enum pj_io_units curr_step_output = pj_right (steps[i].pj);
        enum pj_io_units next_step_input  = pj_left  (steps[i+1].pj);

        if ( curr_step_output == PJ_IO_UNITS_WHATEVER || next_step_input == PJ_IO_UNITS_WHATEVER )
            continue;

        if ( curr_step_output != next_step_input ) {
            proj_log_error (P, _("Pipeline: Mismatched units between step %d and %d"), i+1, i+2);
            return destructor (P, PROJ_ERR_INVALID_OP_WRONG_SYNTAX);
        }
    }

    proj_log_trace (P, "Pipeline: %d steps built. Determining i/o characteristics", nsteps);

    /* Determine forward input (= reverse output) data type */
    P->left = pj_left (steps.front().pj);

    /* Now, correspondingly determine forward output (= reverse input) data type */
    P->right = pj_right (steps.back().pj);

    pipeline->step_count = steps.size();
    pipeline->steps = (PipelineStep*)P->host->ctx->allocator->svm_calloc(pipeline->step_count, sizeof(PipelineStep));
    memcpy(pipeline->steps, steps.data(), sizeof(steps[0]) * pipeline->step_count);

    for (i = 0; i < 4; ++i)
    {
        PipelineStackNew(P->host->ctx->allocator, stack + i, 128);
    }

    return P;
}

#endif // !PROJ_OPENCL_DEVICE

PJ_COORD pipeline_push(PJ_COORD point, __global PJ *P) {
    if (P->parent == nullptr)
        return point;

    __global Pipeline* pipeline = static_cast<__global Pipeline*>(P->parent->opaque);
    __global PushPop* pushpop = static_cast<__global PushPop*>(P->opaque);
    __global PipelineStack* stack = &pipeline->stack1;

    if (pushpop->v1)
        PipelineStackPush(stack + 0, point.v[0]);
    if (pushpop->v2)
        PipelineStackPush(stack + 1, point.v[1]);
    if (pushpop->v3)
        PipelineStackPush(stack + 2, point.v[2]);
    if (pushpop->v4)
        PipelineStackPush(stack + 3, point.v[3]);

    return point;
}

PJ_COORD pipeline_pop(PJ_COORD point, __global PJ *P) {
    if (P->parent == nullptr)
        return point;

    __global Pipeline* pipeline = static_cast<__global Pipeline*>(P->parent->opaque);
    __global PushPop* pushpop = static_cast<__global PushPop*>(P->opaque);
    __global PipelineStack* stack = &pipeline->stack1;

    if (pushpop->v1 && !PipelineStackEmpty(stack + 0)) {
            point.v[0] = PipelineStackTop(stack + 0);
            PipelineStackPop(stack + 0);
    }

    if (pushpop->v2 && !PipelineStackEmpty(stack + 1)) {
            point.v[1] = PipelineStackTop(stack + 1);
            PipelineStackPop(stack + 1);
    }

    if (pushpop->v3 && !PipelineStackEmpty(stack + 2)) {
            point.v[2] = PipelineStackTop(stack + 2);
            PipelineStackPop(stack + 2);
    }

    if (pushpop->v4 && !PipelineStackEmpty(stack + 3)) {
            point.v[3] = PipelineStackTop(stack + 3);
            PipelineStackPop(stack + 3);
    }

    return point;
}

#ifndef PROJ_OPENCL_DEVICE

static PJ *setup_pushpop(PJ *P) {
    auto pushpop = static_cast<struct PushPop*>(P->host->ctx->allocator->svm_calloc (1, sizeof(struct PushPop)));
    P->opaque = pushpop;
    if (nullptr==P->opaque)
        return destructor(P, PROJ_ERR_OTHER /*ENOMEM*/);

    if (pj_param_exists(P->host->params, "v_1"))
        pushpop->v1 = true;

    if (pj_param_exists(P->host->params, "v_2"))
        pushpop->v2 = true;

    if (pj_param_exists(P->host->params, "v_3"))
        pushpop->v3 = true;

    if (pj_param_exists(P->host->params, "v_4"))
        pushpop->v4 = true;

    P->left  = PJ_IO_UNITS_WHATEVER;
    P->right = PJ_IO_UNITS_WHATEVER;

    return P;
}


PJ *OPERATION(push, 0) {
    P->fwd4d = PJ_MAKE_KERNEL(pipeline_push);
    P->inv4d = PJ_MAKE_KERNEL(pipeline_pop);

    return setup_pushpop(P);
}

PJ *OPERATION(pop, 0) {
    P->inv4d = PJ_MAKE_KERNEL(pipeline_push);
    P->fwd4d = PJ_MAKE_KERNEL(pipeline_pop);

    return setup_pushpop(P);
}

#endif // !PROJ_OPENCL_DEVICE

#include "proj_internal_shared.h"

#ifndef PROJ_OPENCL_DEVICE
#define FROM_PROJ_CPP
#include "proj/internal/io_internal.hpp"
#endif

/* Work around non-constness of MSVC HUGE_VAL by providing functions rather than constants */
PJ_COORD proj_coord_error(void) {
    PJ_COORD c;
    c.v[0] = c.v[1] = c.v[2] = c.v[3] = HUGE_VAL;
    return c;
}

/*************************************************************************************/
PJ_COORD pj_geocentric_latitude (const PJ *P, PJ_DIRECTION direction, PJ_COORD coord) {
/**************************************************************************************
    Convert geographical latitude to geocentric (or the other way round if
    direction = PJ_INV)

    The conversion involves a call to the tangent function, which goes through the
    roof at the poles, so very close (the last centimeter) to the poles no
    conversion takes place and the input latitude is copied directly to the output.

    Fortunately, the geocentric latitude converges to the geographical at the
    poles, so the difference is negligible.

    For the spherical case, the geographical latitude equals the geocentric, and
    consequently, the input is copied directly to the output.
**************************************************************************************/
    const double limit = M_HALFPI - 1e-9;
    PJ_COORD res = coord;
    if ((coord.lp.phi > limit) || (coord.lp.phi < -limit) || (P->es==0))
        return res;
    if (direction==PJ_FWD)
        res.lp.phi = atan (P->one_es * tan (coord.lp.phi) );
    else
        res.lp.phi = atan (P->rone_es * tan (coord.lp.phi) );

    return res;
}

/*****************************************************************************/
int proj_errno_set (const PJ *P, int err) {
/******************************************************************************
    Set context-errno, bubble it up to the thread local errno, return err
******************************************************************************/
    /* Use proj_errno_reset to explicitly clear the error status */
    if (0==err)
        return 0;

    /* For P==0 err goes to the default context */
    proj_context_errno_set (pj_get_ctx_shared(P), err);
#ifndef PROJ_OPENCL_DEVICE
    errno = err;
#endif

    return err;
}

/*****************************************************************************/
void proj_context_errno_set (struct pj_ctx_shared *ctx, int err) {
/******************************************************************************
Raise an error directly on a context, without going through a PJ belonging
to that context.
******************************************************************************/
#ifndef PROJ_OPENCL_DEVICE
    if (nullptr==ctx)
        ctx = pj_get_default_ctx()->shared;
#endif
    ctx->last_errno = err;
#ifndef PROJ_OPENCL_DEVICE
    if( err == 0 )
        return;
    errno = err;
#endif
}

/* reduce argument to range +/- PI */
double adjlon (double lon) {
    /* Let lon slightly overshoot, to avoid spurious sign switching at the date line */
    if (fabs (lon) < M_PI + 1e-12)
        return lon;

    /* adjust to 0..2pi range */
    lon += M_PI;

    /* remove integral # of 'revolutions'*/
    lon -= M_TWOPI * floor(lon / M_TWOPI);

    /* adjust back to -pi..pi range */
    lon -= M_PI;

    return lon;
}

PJ_DIRECTION opposite_direction(PJ_DIRECTION dir) {
    return (PJ_DIRECTION)(-dir);
}

struct pj_ctx_shared* pj_get_ctx_shared(const PJ* P)
{
#ifdef PROJ_OPENCL_DEVICE
    return P->shared_ctx;
#else
    return pj_get_ctx((PJ*)P)->shared;
#endif
}

/**************************************************************************************/
PJ_COORD proj_trans (PJ *P, PJ_DIRECTION direction, PJ_COORD coord) {
/***************************************************************************************
Apply the transformation P to the coordinate coord, preferring the 4D interfaces if
available.

See also pj_approx_2D_trans and pj_approx_3D_trans in pj_internal.c, which work
similarly, but prefers the 2D resp. 3D interfaces if available.
***************************************************************************************/
    PJstack_t stack;

    stack_new(&stack);

    push_proj_trans(&stack, P, direction, coord);

    return stack_exec(&stack);

#if 0
    if( !P->host->alternativeCoordinateOperations.empty() ) {
        constexpr int N_MAX_RETRY = 2;
        int iExcluded[N_MAX_RETRY] = {-1, -1};

        const int nOperations = static_cast<int>(
            P->host->alternativeCoordinateOperations.size());

        // We may need several attempts. For example the point at
        // lon=-111.5 lat=45.26 falls into the bounding box of the Canadian
        // ntv2_0.gsb grid, except that it is not in any of the subgrids, being
        // in the US. We thus need another retry that will select the conus
        // grid.
        for( int iRetry = 0; iRetry <= N_MAX_RETRY; iRetry++ )
        {
            // Do a first pass and select the operations that match the area of use
            // and has the best accuracy.
            int iBest = pj_get_suggested_operation(P->host->ctx,
                                                   P->host->alternativeCoordinateOperations,
                                                   iExcluded,
                                                   direction,
                                                   coord);
            if( iBest < 0 ) {
                break;
            }
            if( iRetry > 0 ) {
                const int oldErrno = proj_errno_reset(P);
                if (proj_log_level(P->host->ctx, PJ_LOG_TELL) >= PJ_LOG_DEBUG) {
                    pj_log(P->host->ctx, PJ_LOG_DEBUG, proj_context_errno_string(P->host->ctx, oldErrno));
                }
                pj_log(P->host->ctx, PJ_LOG_DEBUG,
                    "Did not result in valid result. "
                    "Attempting a retry with another operation.");
            }

            const auto& alt = P->host->alternativeCoordinateOperations[iBest];
            if( P->iCurCoordOp != iBest ) {
                if (proj_log_level(P->host->ctx, PJ_LOG_TELL) >= PJ_LOG_DEBUG) {
                    std::string msg("Using coordinate operation ");
                    msg += alt.name;
                    pj_log(P->host->ctx, PJ_LOG_DEBUG, msg.c_str());
                }
                P->iCurCoordOp = iBest;
            }
            PJ_COORD res = direction == PJ_FWD ?
                        pj_fwd4d( coord, alt.pj ) : pj_inv4d( coord, alt.pj );
            if( proj_errno(alt.pj) == PROJ_ERR_OTHER_NETWORK_ERROR ) {
                return proj_coord_error ();
            }
            if( res.xyzt.x != HUGE_VAL ) {
                return res;
            }
            if( iRetry == N_MAX_RETRY ) {
                break;
            }
            iExcluded[iRetry] = iBest;
        }

        // In case we did not find an operation whose area of use is compatible
        // with the input coordinate, then goes through again the list, and
        // use the first operation that does not require grids.
        NS_PROJ::io::DatabaseContextPtr dbContext;
        try
        {
            if( P->host->ctx->cpp_context ) {
                dbContext = P->host->ctx->cpp_context->getDatabaseContext().as_nullable();
            }
        }
        catch( const std::exception& ) {}
        for( int i = 0; i < nOperations; i++ ) {
            const auto &alt = P->host->alternativeCoordinateOperations[i];
            auto coordOperation = dynamic_cast<
            NS_PROJ::operation::CoordinateOperation*>(alt.pj->host->iso_obj.get());
            if( coordOperation ) {
                if( coordOperation->gridsNeeded(dbContext, true).empty() ) {
                    if( P->iCurCoordOp != i ) {
                        if (proj_log_level(P->host->ctx, PJ_LOG_TELL) >= PJ_LOG_DEBUG) {
                            std::string msg("Using coordinate operation ");
                            msg += alt.name;
                            msg += " as a fallback due to lack of more "
                                   "appropriate operations";
                            pj_log(P->host->ctx, PJ_LOG_DEBUG, msg.c_str());
                        }
                        P->iCurCoordOp = i;
                    }
                    if( direction == PJ_FWD ) {
                        return pj_fwd4d( coord, alt.pj );
                    }
                    else {
                        return pj_inv4d( coord, alt.pj );
                    }
                }
            }
        }

        proj_errno_set (P, PROJ_ERR_COORD_TRANSFM_NO_OPERATION);
        return proj_coord_error ();
    }
#endif
}

/**************************************************************************************/
PJ_COORD pj_approx_2D_trans (PJ *P, PJ_DIRECTION direction, PJ_COORD coo) {
/***************************************************************************************
Behave mostly as proj_trans, but attempt to use 2D interfaces only.
Used in gie.c, to enforce testing 2D code, and by PJ_pipeline.c to implement
chained calls starting out with a call to its 2D interface.
***************************************************************************************/
    PJstack_t stack;
    if (nullptr==P || direction == PJ_IDENT)
        return coo;

    stack_new(&stack);

    push_approx_3D_trans(&stack, P, direction, coo);

    coo = stack_exec(&stack);
    coo.xyzt.z = coo.xyzt.t = 0.0;

    return coo;
}

/**************************************************************************************/
PJ_COORD pj_approx_3D_trans (PJ *P, PJ_DIRECTION direction, PJ_COORD coo) {
/***************************************************************************************
Companion to pj_approx_2D_trans.

Behave mostly as proj_trans, but attempt to use 3D interfaces only.
Used in gie.c, to enforce testing 3D code, and by PJ_pipeline.c to implement
chained calls starting out with a call to its 3D interface.
***************************************************************************************/
    PJstack_t stack;
    if (nullptr== P || direction == PJ_IDENT)
        return coo;

    stack_new(&stack);

    push_approx_3D_trans(&stack, P, direction, coo);

    coo = stack_exec(&stack);
    coo.xyzt.t = 0.0;

    return coo;
}

void stack_new(cl_local PJstack_t* stack)
{
    stack->n = 0;
}

PJ_COORD stack_exec(cl_local PJstack_t* stack)
{
    // Default to an error value.
    PJ_COORD coo = { -1, -1, -1, -1 };

    while (stack->n > 0)
    {
        cl_local PJstack_entry_t* top = stack->s + stack->n - 1;

        switch (top->fn(stack, top))
        {
            case PJ_CO_YIELD:
            {
                // Nothing to do. We have a new top we need to run!
                break;
            }

            case PJ_CO_DONE:
            {
                --stack->n;
                
                // Transfer state back to the caller.
                if (!stack->n)
                {
                    coo = top->coo;
                }
                else
                {
                    (top - 1)->coo = top->coo;
                }

                break;
            }

            case PJ_CO_ERROR:
            {
                stack->n = 0;
                break;
            }
        }
    }

    return coo;
}

void stack_push(cl_local PJstack_t* stack, PJ_COROUTINE fn, PJ* P, PJ_COORD coo)
{
    cl_local PJstack_entry_t* e = stack->s + stack->n;
    
    if (stack->n == PR_CO_STACK_SIZE)
    {
        stack->n = 0;
        return;
    }

    e->fn = fn;
    e->coo = coo;
    e->P = P;

    e->step = 0;
    e->i = 0;
    e->last_errno = 0;

    ++stack->n;
}

void push_proj_trans(cl_local PJstack_t* stack, PJ* P, PJ_DIRECTION direction, PJ_COORD coord)
{
    if (nullptr == P || direction == PJ_IDENT)
        return;
    if (P->inverted)
        direction = opposite_direction(direction);

    stack_push(stack, direction == PJ_FWD ? pj_fwd4d_co : pj_inv4d_co, P, coord);
}

void push_approx_2D_trans(cl_local PJstack_t* stack, PJ* P, PJ_DIRECTION direction, PJ_COORD coord)
{
    if (nullptr == P || direction == PJ_IDENT)
        return;
    if (P->inverted)
        direction = opposite_direction(direction);

    stack_push(stack, direction == PJ_FWD ? pj_fwd_co : pj_inv_co, P, coord);
}

void push_approx_3D_trans(cl_local PJstack_t* stack, PJ* P, PJ_DIRECTION direction, PJ_COORD coord)
{
    if (nullptr == P || direction == PJ_IDENT)
        return;
    if (P->inverted)
        direction = opposite_direction(direction);

    stack_push(stack, direction == PJ_FWD ? pj_fwd3d_co : pj_inv3d_co, P, coord);
}
/* These are functions that are setup to run on host and device. */
PROJ_COROUTINE(fwd_prepare_co)
PROJ_COROUTINE(fwd_finalize_co)
PROJ_COROUTINE(pj_fwd_co)
PROJ_COROUTINE(pj_fwd3d_co)
PROJ_COROUTINE(pj_fwd4d_co)

PROJ_COROUTINE(inv_prepare_co)
PROJ_COROUTINE(inv_finalize_co)
PROJ_COROUTINE(pj_inv_co)
PROJ_COROUTINE(pj_inv3d_co)
PROJ_COROUTINE(pj_inv4d_co)

PROJ_COROUTINE(pipeline_forward_4d_co)
PROJ_COROUTINE(pipeline_reverse_4d_co)
PROJ_COROUTINE(pipeline_forward_3d_co)
PROJ_COROUTINE(pipeline_reverse_3d_co)
PROJ_COROUTINE(pipeline_forward_co)
PROJ_COROUTINE(pipeline_reverse_co)

PROJ_FWD_2D(eqc_s_forward)
PROJ_INV_2D(eqc_s_inverse)

PROJ_OPERATOR(pipeline_push)
PROJ_OPERATOR(pipeline_pop)

PROJ_FWD_2D(qsc_e_forward)
PROJ_INV_2D(qsc_e_inverse)

PROJ_OPERATOR(unitconvert_forward_4d)
PROJ_OPERATOR(unitconvert_reverse_4d)
PROJ_FWD_3D(unitconvert_forward_3d)
PROJ_INV_3D(unitconvert_reverse_3d)
PROJ_FWD_2D(unitconvert_forward_2d)
PROJ_INV_2D(unitconvert_reverse_2d)
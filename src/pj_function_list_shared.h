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

PROJ_OPERATOR(affine_forward_4d)
PROJ_OPERATOR(affine_reverse_4d)
PROJ_FWD_3D(affine_forward_3d)
PROJ_INV_3D(affine_reverse_3d)
PROJ_FWD_2D(affine_forward_2d)
PROJ_INV_2D(affine_reverse_2d)

PROJ_OPERATOR(axisswap_forward_4d)
PROJ_OPERATOR(axisswap_reverse_4d)
PROJ_FWD_3D(axisswap_forward_3d)
PROJ_INV_3D(axisswap_reverse_3d)
PROJ_FWD_2D(axisswap_forward_2d)
PROJ_INV_2D(axisswap_reverse_2d)

PROJ_OPERATOR(noop_operator)

PROJ_OPERATOR(set_fwd_inv)

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

PROJ_OPERATOR(latlong_forward_4d)
PROJ_OPERATOR(latlong_inverse_4d)
PROJ_FWD_3D(latlong_forward_3d)
PROJ_INV_3D(latlong_inverse_3d)
PROJ_FWD_2D(latlong_forward)
PROJ_INV_2D(latlong_inverse)

PROJ_FWD_2D(merc_e_forward)
PROJ_INV_2D(merc_e_inverse)
PROJ_FWD_2D(merc_s_forward)
PROJ_INV_2D(merc_s_inverse)

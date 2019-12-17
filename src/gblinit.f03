! initialise global variables
SUBROUTINE globals_init()
use mdl_param
use mdl_eop
use mdl_num

        ! from module EOP
        EOP_MJD0_glb = 0.d0

        ! from module NUM
        GMsun_glb = 0.d0
        GMmoon_glb = 0.d0

        ! from module param
        POD_MODE_glb = 0
        ORBPRED_ARC_glb = 0.d0
        ORBEXT_glb = 0
        SRP_MOD_glb = 0
        ECOM_PARAM_glb = 0
        ECOM_Bias_glb = 0
        ECOM_CPR_glb = 0
        EMP_param_glb = 0
        EMP_Bias_glb = 0
        Bias_accel_glb = 0.d0
        EMP_CPR_glb = 0
        EMP_nCPR_glb = 0
        CPR_CS_glb = 0.d0
        Frame_EmpiricalForces_glb = 0
        VEQ_integration_glb = 0
        NPARAM_glb = 0
        ESTIM_mode_glb = 0
        ESTIM_iter_glb = 0
        SATblock_glb = 0
        BDSorbtype_glb = "12345"
        return
end

! initialise & finalise global variables
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
        SRP_MOD_arp = 0
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
        ! only require month parameter (3rd arg) when ERM (1st arg) is 2)
        call BOXWINGINIT(1, 1, 1)
        return
end

subroutine globals_fini()
use mdl_param
use mdl_tides
use mdl_planets
use mdl_eop

        integer(kind = prec_int2) DeallocateStatus

        ! from module param
        if (allocated(ECOM_accel_glb)) Deallocate(ECOM_accel_glb, Stat=DeallocateStatus)
        if (allocated(IDAT)) Deallocate(IDAT, Stat=DeallocateStatus)
        if (allocated(DATS)) Deallocate(DATS, Stat=DeallocateStatus)
        if (allocated(GFM_Cnm)) Deallocate(GFM_Cnm, Stat=DeallocateStatus)
        if (allocated(GFM_Snm)) Deallocate(GFM_Snm, Stat=DeallocateStatus)
        if (allocated(pseudobs_ICRF)) Deallocate(pseudobs_ICRF, Stat=DeallocateStatus)
        if (allocated(pseudobs_ITRF)) Deallocate(pseudobs_ITRF, Stat=DeallocateStatus)
        if (allocated(orbext_ICRF)) Deallocate(orbext_ICRF, Stat=DeallocateStatus)
        if (allocated(orbext_ITRF)) Deallocate(orbext_ITRF, Stat=DeallocateStatus)

        ! from module tides
        if (allocated(Doodson_mult_glb)) Deallocate(Doodson_mult_glb, Stat=DeallocateStatus)
        if (allocated(Delaunay_FES)) Deallocate(Delaunay_FES, Stat=DeallocateStatus)
        if (allocated(dCnm_p)) Deallocate(dCnm_p, Stat=DeallocateStatus)
        if (allocated(dSnm_p)) Deallocate(dSnm_p, Stat=DeallocateStatus)
        if (allocated(dCnm_m)) Deallocate(dCnm_m, Stat=DeallocateStatus)
        if (allocated(dSnm_m)) Deallocate(dSnm_m, Stat=DeallocateStatus)

        ! from module planets
        if (allocated(CVAL_2)) Deallocate(CVAL_2, Stat=DeallocateStatus)
        if (allocated(DB_array)) Deallocate(DB_array, Stat=DeallocateStatus)
        
        ! from module eop
        if (allocated(EOP_day_glb)) Deallocate(EOP_day_glb, Stat=DeallocateStatus)

        return
end


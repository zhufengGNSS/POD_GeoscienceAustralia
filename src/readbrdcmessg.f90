SUBROUTINE  readbrdcmessg (UNIT_IN,ISAT,IREC,EPH,CLK)


! ----------------------------------------------------------------------
! MODULE: readbrdcmessg.f90 
! ----------------------------------------------------------------------
! Purpose: Read the orbit and clock parameters from broadcast ephemeris 
!          denotes the year 
! ----------------------------------------------------------------------
! Author :      Dr. Tzupang Tseng
!
! Copyright: GEOSCIENCE AUSTRALIA, AUSTRALIA
!
! Created:      31-05-2019
! ----------------------------------------------------------------------
      USE mdl_precision
      USE mdl_num
      USE mdl_planets
      USE m_shadow
      USE mdl_param
      USE m_reformbrdc
      IMPLICIT NONE


!-------------------------------------------------------------------
      INTEGER (KIND = prec_int4) :: prnnum
      INTEGER (KIND = 4)         :: satsvn
      REAL (KIND = prec_d), INTENT(IN) :: mjd
      REAL (KIND = prec_d), DIMENSION(3), INTENT(IN) :: rsat, vsat
      REAL (KIND = prec_q), INTENT(OUT) :: beta
      REAL (KIND = prec_q) :: lambda, del_u, yaw
!--------------------------------------------------------------------
      DOUBLE PRECISION  JD, Zbody(6)
      INTEGER (KIND = 4) :: i, j
      REAL (KIND = prec_q), DIMENSION(3) :: rbody
      REAL (KIND = prec_q), DIMENSION(3) :: rSun, rMoon
      REAL (KIND = prec_q), DIMENSION(3) :: r_sun1, r_sun2
      REAL (KIND = prec_q) :: u_sun, yaw2
      INTEGER  NTARG_SUN, NTARG_MOON, NCTR
      REAL (KIND = prec_q), DIMENSION(3) :: ed, ey, eb, ex, en, ev, ez, et
      REAL (KIND = prec_q), DIMENSION(3) :: yy
      REAL (KIND = prec_q), DIMENSION(9) :: kepler
      REAL (KIND = prec_q), DIMENSION(4) :: cosang
      REAL (KIND = prec_q) :: R11(3,3),R33(3,3)
      REAL (KIND = prec_q) :: Pi, Ps, AU, sclfa
      REAL (KIND = prec_d) :: GM, Ds
      REAL (KIND = prec_d) :: u_sat, i_sat, omega_sat
      INTEGER (KIND = prec_int2) :: AllocateStatus
      INTEGER              :: PD_Param_ID



K=0
! Read broadcast message 
!-------------------------------

CALL brdcinfo (UNIT_IN, VERSION, ISVN, EPHDAT)

CALL reformbrdc (EPHDAT,EPHV3,CLKV3)

K=K+1

DO I=1,20
EPH(I,ISVN)=EPHV3(I)
CLOCK(I,ISVN)=CLKV3(I)
END DO


END SUBROUTINE  

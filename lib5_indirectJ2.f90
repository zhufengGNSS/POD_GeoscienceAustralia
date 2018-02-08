      program lib6_indirectJ2


! ----------------------------------------------------------------------
! Program: lib6_indirectJ2.f90
! ----------------------------------------------------------------------
! Purpose:
!  Orbital perturbations due to the indirect J2 effect
! ----------------------------------------------------------------------
! Called subroutines:
!
! ----------------------------------------------------------------------
! Remark 1:
!  Lines 149-158 : 
!  The transformation matrix (GCRS to ITRS) will be provided by the subroutines
!  of the TRS function which is currently under development. Therefore,
!  this will be fulfilled in a future version which will include TRS. 
! ----------------------------------------------------------------------
! Dr. Thomas Papanikolaou, Geoscience Australia            November 2015
! ----------------------------------------------------------------------

	  
      USE mdl_precision
      USE mdl_num
      USE mdl_planets
      IMPLICIT NONE

! ----------------------------------------------------------------------
! Variables declaration
! ----------------------------------------------------------------------
! "ET, R(6)" variables enter F77 code and thus, these must be declared in Double Precision.
! ----------------------------------------------------------------------
!      DOUBLE PRECISION  ET
!      DOUBLE PRECISION  R(6), R_Horizon(6), dR(6)
      REAL (KIND = prec_d) :: ET, R(6), R_Horizon(6), dR(6)
      INTEGER  NTARG, NCTR
! ----------------------------------------------------------------------
      CHARACTER (LEN=100) :: fname_header,fname_data,fname_out				
      REAL (KIND = prec_q) :: xsat,ysat,zsat, GMbody, GMearth
      REAL (KIND = prec_q), DIMENSION(3) :: rsat, rSun,rMoon,rMoon_ITRS,rSun_ITRS		
      REAL (KIND = prec_q) :: C20, Re, GM_moon,GM_sun
      REAL (KIND = prec_q), DIMENSION(3) :: a_iJ2
      INTEGER (KIND = prec_int4) :: DeAllocateStatus
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! INPUT
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! C20 spherical harmonic coefficient of geopotential 
      C20 = -4.841694552725D-04
! Earth Radius
      Re = Earth_radius
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Satellite position vector
      xsat = -2.747681940000000D5
      ysat =  5.812939376999999D6
      zsat = -3.213088276000000D6
	  rsat = (/ xsat,ysat,zsat /)
      PRINT *, "rsat", rsat
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Epoch
! ----------------------------------------------------------------------
! Julian Day Number
!430  2015.02.01 2457054.5  1  8  4       -0.02577866640782599000
      ET =  2457054.5D0
!      ET =  2457054.78123467D0
      PRINT *, "ET:", ET	   
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Planetary and Lunar ephemeris data
! ----------------------------------------------------------------------
      fname_out = 'DE.430' 
      fname_header = 'header.430_229'
      fname_data = 'ascp1950.430'
      CALL CATfile (fname_header,fname_data,fname_out)
      CALL asc2eph (fname_out)
      CALL GM_de
! ----------------------------------------------------------------------
! STATE.f, Line 157:    KM or AU units for coordinates is set here 
! ----------------------------------------------------------------------
! FSIZER3.f sets KSIZE (e.g. DE430: KSIZE = 2036) 
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! End of INPUT
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Sun and Moon position vectors in celestial reference frame
! ----------------------------------------------------------------------
! Center celestial body
      NCTR = 3 ! Earth
! ----------------------------------------------------------------------
! Target celestial body : NTARG
! ----------------------------------------------------------------------
! Sun
      NTARG = 11 ! Sun
      CALL PLEPH ( ET, NTARG, NCTR, R )
	  rSun = (/ R(1),R(2),R(3) /) * 1000D0 ! KM to M
      PRINT *,"rSun", rSun  
      GM_sun = GMconst(NTARG)
      PRINT *,"GM_sun", GM_sun
! ----------------------------------------------------------------------
! Numerical comparison with JPL Horizon web calculation service
R_Horizon = (/ -9.782087484269591D+07,  1.011533940461814D+08,  4.385217774571174D+07, &
               -2.275539856197723D+01, -1.824081034278090D+01, -7.907465604900494D+00 /)
dR = R + R_Horizon
!PRINT *,"dR", dR

! ----------------------------------------------------------------------
! Moon
      NTARG = 10
      CALL PLEPH ( ET, NTARG, NCTR, R )
	  rMoon = (/ R(1),R(2),R(3) /) * 1000D0 ! KM to M
      PRINT *,"rMoon", rMoon  
      GM_moon = GMconst(NTARG)
      PRINT *,"GM_moon", GM_moon
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
      DEALLOCATE (CVAL_2,   STAT = DeAllocateStatus)
      DEALLOCATE (DB_array, STAT = DeAllocateStatus)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
!      CALL force_gm3rd (rbody,rsat,GMbody , a_perturb)
!      PRINT *,"a_perturb", a_perturb
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Sun and Moon position vectors in terrestrial reference frame 
! ----------------------------------------------------------------------
! Transformation GCRS to ITRS is not applied in this version : see Remark 1
! ----------------------------------------------------------------------
      rMoon_ITRS = rMoon
      rSun_ITRS = rSun
!      PRINT *,"rMoon_ITRS", rMoon_ITRS
!      PRINT *,"rSun_ITRS", rSun_ITRS
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Indirect J2 effect of Sun and Moon
! ----------------------------------------------------------------------
      CALL indirectJ2(C20,Re,GM_moon,rMoon_ITRS,GM_sun,rSun_ITRS , a_iJ2)
      PRINT *,"----------------------------------------------------------------------"											
      PRINT *,"a_iJ2", a_iJ2
      PRINT *,"----------------------------------------------------------------------"											
! ----------------------------------------------------------------------

	  
      end

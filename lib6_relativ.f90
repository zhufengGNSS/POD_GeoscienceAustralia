      program lib6_relativ


! ----------------------------------------------------------------------
! Program: lib6_relativ.f90
! ----------------------------------------------------------------------
! Purpose:
!  General relativistic effects to satellite orbits
!  Computation of the satellite’s acceleration corrections due to general
!  relativistic effects.
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
      REAL (KIND = prec_q) :: xsat,ysat,zsat, GMbody, GMearth, GM_moon,GM_sun
      REAL (KIND = prec_q), DIMENSION(3) :: rsat, rSun,rMoon,rMoon_ITRS,rSun_ITRS, vSun,vsat	
      INTEGER (KIND = prec_int4) :: DeAllocateStatus
! ----------------------------------------------------------------------
      REAL (KIND = prec_q), DIMENSION(6) :: Zsat_GCRS,Zearth_GCRS
      REAL (KIND = prec_q), DIMENSION(3) :: a_Schwarzschild, a_LenseThirring, a_deSitter
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! INPUT
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Satellite state vector

! -2.747435450000000e+05 5.812884181000001e+06 -3.213190333000000e+06 1.618282138000000e+03 -3.632305352000000e+03 -6.716404161700000e+03
!      xsat = -2.747681940000000D5
!      ysat =  5.812939376999999D6
!      zsat = -3.213088276000000D6
!	  rsat = (/ xsat,ysat,zsat /)
!	  vsat = (/ 1.618282138000000D3, -3.632305352000000D3 , -6.716404161700000D3 /)  

! Satellite state vector
! Position vector
! -0.77666873443201953E+06    0.22052406037103664E+08   -0.14272932415913470E+08   -0.28153818042302833E+04    0.14750243023567887E+04    0.22860949538585546E+04 	  
      xsat = -0.77666873443201953D6
      ysat =  0.22052406037103664D8
      ysat = -0.14272932415913470D8	  
	  rsat = (/ xsat,ysat,zsat /)
! Velocity vector
	  vsat = (/ -0.28153818042302833D4, 0.14750243023567887D4 , 0.22860949538585546D4 /)  
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
      NTARG = 11
      CALL PLEPH ( ET, NTARG, NCTR, R )
	  rSun = (/ R(1),R(2),R(3) /) * 1000D0 ! KM to M
	  vSun = (/ R(4),R(5),R(6) /) * 1000D0 ! KM/sec to m/sec
      GM_sun = GMconst(NTARG)
      PRINT *,"rSun", rSun  
      PRINT *,"vSun", vSun  
      PRINT *,"GM_sun", GM_sun
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Moon
      NTARG = 10
      CALL PLEPH ( ET, NTARG, NCTR, R )
	  rMoon = (/ R(1),R(2),R(3) /) * 1000D0 ! KM to M
      GM_moon = GMconst(NTARG)
      PRINT *,"rMoon", rMoon  
      PRINT *,"GM_moon", GM_moon
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
      DEALLOCATE (CVAL_2,   STAT = DeAllocateStatus)
      DEALLOCATE (DB_array, STAT = DeAllocateStatus)
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
      GMsun_glb = GM_sun
	  GM_global = GM_gfm
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Satellite State Vector in GCRS
      Zsat_GCRS = (/rsat(1), rsat(2), rsat(3), vsat(1), vsat(2), vsat(3) /) 
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Earth state vector w.r.t. Sun
! ----------------------------------------------------------------------
      Zearth_GCRS = -1.0D0 * (/ rSun(1),rSun(2),rSun(3) , vSun(1),vSun(2),vSun(3) /)
      PRINT *,"Zearth_GCRS",Zearth_GCRS
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Relativistic Effects
! ----------------------------------------------------------------------
      CALL rel_schwarzschild (Zsat_GCRS,Zearth_GCRS , a_Schwarzschild)
      CALL rel_LenseThirring (Zsat_GCRS,Zearth_GCRS , a_LenseThirring)
      CALL rel_deSitter (Zsat_GCRS,Zearth_GCRS , a_deSitter)
! ----------------------------------------------------------------------
      PRINT *,"----------------------------------------------------------------------"											
      PRINT *,"a_Schwarzschild", a_Schwarzschild
      PRINT *,"a_LenseThirring", a_LenseThirring
      PRINT *,"a_deSitter", a_deSitter
      PRINT *,"----------------------------------------------------------------------"											
! ----------------------------------------------------------------------
      PRINT *,"a_Schwarzschild", sqrt(a_Schwarzschild(1)**2 + a_Schwarzschild(2)**2 + a_Schwarzschild(3)**2)
      PRINT *,"a_LenseThirring", sqrt(a_LenseThirring(1)**2 + a_LenseThirring(2)**2 + a_LenseThirring(3)**2)
      PRINT *,"a_deSitter", sqrt(a_deSitter(1)**2 + a_deSitter(2)**2 + a_deSitter(3)**2)

	  
      end

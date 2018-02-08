      program lib3_planets


! ------------------------------------------------------------------------------
! Purpose:
!  Computation of the 3rd body orbit perturbations due to Sun, Moon and Planets
! ------------------------------------------------------------------------------
	  
! ------------------------------------------------------------------------------
! Remarks:
!  Modified version of the JPL programs for the DE data processing is used 
! ------------------------------------------------------------------------------
! Major modifications:
! ------------------------------------------------------------------------------
! asc2eph.f 
! - Modified for replacing the construction of the JPLEPH binary data file
! - program procedure has been changed to subroutine procedure
! - Revised for reading directly the DE ephemeris file (ascii) instead of reading through the screen input.
! - Subroutine CATfile.f90 has written for supporting the previous purpose  
! - Revised to Fortran 90 for applying dynamic memory allocation
! - Module mdl_planets.f90 has been written for introducing the allocatable arrays which replace the JPLEPH binary data
! ------------------------------------------------------------------------------
! STATE.f
! - Modified for removing the use of the JPLEPH binary data file
! - Revised to Fortran 90 for reading the module mdl_planets.f90 where the required allocatable arrays are defined
! ------------------------------------------------------------------------------
! Dr. Thomas D. Papanikolaou, Geoscience Australia               27 October 2015
! ------------------------------------------------------------------------------
	  

      USE mdl_planets
      USE mdl_precision
      IMPLICIT NONE

	  
! ----------------------------------------------------------------------
! Variables declaration
! ----------------------------------------------------------------------
  
! ----------------------------------------------------------------------
! ATTENTION: The variables ET and R(6) enter F77 code.
! Thus, these must be declared as DP. Do not declare them as QP.
! ----------------------------------------------------------------------
!      DOUBLE PRECISION  ET
!      DOUBLE PRECISION  R(6), R_Horizon(6), dR(6)
      REAL (KIND = prec_d) :: ET, R(6), R_Horizon(6), dR(6)
      INTEGER  NTARG, NCTR
! ----------------------------------------------------------------------
      CHARACTER (LEN=100) :: fname_header,fname_data,fname_out
      REAL (KIND = prec_q), DIMENSION(3) :: rbody,rsat, a_perturb
      REAL (KIND = prec_q) :: xsat,ysat,zsat, GMbody, GMearth, radius_sat, radius_body
      INTEGER (KIND = prec_int4) :: N_CVAL, I, DeAllocateStatus
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! INPUT
! ----------------------------------------------------------------------
      fname_out = 'DE.430' 
      fname_header = 'header.430_229'
      fname_data = 'ascp1950.430'
      CALL CATfile (fname_header,fname_data,fname_out)
      CALL asc2eph (fname_out)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! STATE.f, Line 157:    KM or AU units for coordinates is set
! ----------------------------------------------------------------------
! FSIZER3.f sets KSIZE (e.g. DE430: KSIZE = 2036) 
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! PLEPH.f
!     NTARG = INTEGER NUMBER OF 'TARGET' POINT.
!
!     NCENT = INTEGER NUMBER OF CENTER POINT.
!
!            THE NUMBERING CONVENTION FOR 'NTARG' AND 'NCENT' IS:
!
!                1 = MERCURY           8 = NEPTUNE
!                2 = VENUS             9 = PLUTO
!                3 = EARTH            10 = MOON
!                4 = MARS             11 = SUN
!                5 = JUPITER          12 = SOLAR-SYSTEM BARYCENTER
!                6 = SATURN           13 = EARTH-MOON BARYCENTER
!                7 = URANUS           14 = NUTATIONS (LONGITUDE AND OBLIQ)
!                            15 = LIBRATIONS, IF ON EPH FILE
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Center celestial body
      NCTR = 3 ! Earth
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Satellite position vector
      xsat = -2.747681940000000D5
      ysat =  5.812939376999999D6
      zsat = -3.213088276000000D6
	  rsat(1) = xsat
	  rsat(2) = ysat
	  rsat(3) = zsat
      radius_sat = sqrt(rsat(1)**2 + rsat(2)**2 + rsat(3)**2)
!      PRINT *,"radius_sat",radius_sat  
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Julian Day Number of the input epoch
!430  2015.02.01 2457054.5  1  8  4       -0.02577866640782599000
      ET =  2457054.5D0
!      ET =  2457054.78123467D0
      PRINT *, "ET:", ET	   
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Target celestial body
      NTARG = 11 ! Sun
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! End of INPUT
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Celestial body's Cartesian coordinates : Sun
! Target body's (NTARG) Cartesian coordinates w.r.t. Center body (NCTR)
! ----------------------------------------------------------------------
      CALL PLEPH ( ET, NTARG, NCTR, R )
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! JPL Horizon web calculation service
! ----------------------------------------------------------------------
R_Horizon = (/ -9.782087484269591D+07,  1.011533940461814D+08,  4.385217774571174D+07, &
               -2.275539856197723D+01, -1.824081034278090D+01, -7.907465604900494D+00 /)
! ----------------------------------------------------------------------
! Numerical difference between source code and Horizon system
dR = R + R_Horizon
PRINT *,"dR", dR
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Celestial body's Cartesian coordinates
! Target body's (NTARG) Cartesian coordinates w.r.t. Center body (NCTR)
! ----------------------------------------------------------------------
!      NTARG = 11 ! Sun
!      NTARG = 10 ! GEOCENTRIC MOON
!      NTARG = 5 ! Jupiter
!      NTARG = 2 ! Venus
!      NTARG = 1 ! Mercury 
!      NTARG = 4 ! MARS 
!      NTARG = 6 ! SATURN 
!      NTARG = 7 ! URANUS 
!      NTARG = 8 ! NEPTUNE
! ----------------------------------------------------------------------
      CALL  PLEPH ( ET, NTARG, NCTR, R )
! ----------------------------------------------------------------------
      PRINT *, "r(XYZ):",R(1),R(2),R(3)	   
      PRINT *, "r(XYZ):",R(4),R(5),R(6)   

	  
! ----------------------------------------------------------------------
! GM gravity constants of the solar system bodies
! ----------------------------------------------------------------------
      CALL GM_de
      GMbody = GMconst(NTARG)
! ----------------------------------------------------------------------
      PRINT *,"GMbody", GMbody


! ----------------------------------------------------------------------
! Cartesian coordinates of the celestial body in meters
! ----------------------------------------------------------------------
! KM to M
	  rbody(1) = R(1) * 1000D0
	  rbody(2) = R(2) * 1000D0
	  rbody(3) = R(3) * 1000D0
radius_body = sqrt(rbody(1)**2 + rbody(2)**2 + rbody(3)**2)
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Point-mass perturbations vector due to the selected celestial body
! ----------------------------------------------------------------------
      CALL force_gm3rd (rbody,rsat,GMbody , a_perturb)
! ----------------------------------------------------------------------
      PRINT *,"a_perturb", a_perturb


	  
! ----------------------------------------------------------------------
      DEALLOCATE (CVAL_2,   STAT = DeAllocateStatus)
      DEALLOCATE (DB_array, STAT = DeAllocateStatus)
! ----------------------------------------------------------------------

	  
      end
  

      program lib3_planets


! ------------------------------------------------------------------------------
! Program: lib3_planets.f90
! ------------------------------------------------------------------------------
! Purpose:
!  Computation of the 3rd body orbit perturbations due to Sun, Moon and Planets
! ------------------------------------------------------------------------------
	  
! ------------------------------------------------------------------------------
! Remarks:
!  A modified version of the JPL programs for the DE data processing has been 
!  implemented and used here 
!
! Major modifications:
! asc2eph.f 
! - Modified for replacing the construction of the JPLEPH binary data file
! - program procedure has been changed to subroutine procedure
! - Revised for reading directly the DE ephemeris file (ascii) instead of reading through the screen input.
! - Subroutine CATfile.f90 has written for supporting the previous purpose  
! - Revised to Fortran 90 for applying dynamic memory allocation
! - Module mdl_planets.f90 has been written for introducing the allocatable arrays which replace the JPLEPH binary data
! STATE.f
! - Modified for removing the use of the JPLEPH binary data file
! - Revised to Fortran 90 for reading the module mdl_planets.f90 where the required allocatable arrays are defined
! ------------------------------------------------------------------------------
! Dr. Thomas D. Papanikolaou, Geoscience Australia               27 October 2015
! ------------------------------------------------------------------------------
	  

      USE mdl_planets
      USE mdl_precision
      USE mdl_num
      IMPLICIT NONE

	  
! ----------------------------------------------------------------------
! Variables declaration
! ----------------------------------------------------------------------
      DOUBLE PRECISION  ET, R(6)
      INTEGER  NTARG, NCTR, NTARG_body
! ----------------------------------------------------------------------
      CHARACTER (LEN=100) :: fname_header,fname_data,fname_out
      REAL (KIND = prec_q), DIMENSION(3) :: rbody,rsat, a_perturb
      REAL (KIND = prec_q) :: xsat,ysat,zsat, GMbody, GMearth
      INTEGER (KIND = prec_int4) :: N_CVAL, I, DeAllocateStatus
! ----------------------------------------------------------------------
      REAL (KIND = prec_q), DIMENSION(3) :: rSun, rMoon, rMoon_ITRS, rSun_ITRS		
      REAL (KIND = prec_q) :: C20, Re, GM_moon,GM_sun
      REAL (KIND = prec_q), DIMENSION(3) :: a_iJ2

	  


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
! Center celestial body
      NCTR = 3 ! Earth
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Satellite position vector
      xsat = 17835.809256D0 * 1.0D3
      ysat =  4999.483252D0 * 1.0D3
      ysat = 19053.904501D0 * 1.0D3
	  
	  rsat(1) = xsat
	  rsat(2) = ysat
	  rsat(3) = zsat
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Julian Day Number of the input epoch
!430  2015.02.01 2457054.5  1  8  4       -0.02577866640782599000
      ET =  2457054.5D0
!      ET =  2457054.78123467D0
      PRINT *, "ET:", ET	   
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! C20 spherical harmonic coefficient of geopotential 
      C20 = -4.841694552725D-04
! Earth Radius
      Re = Earth_radius
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! End of INPUT
! ----------------------------------------------------------------------




      DO NTARG_body = 1 , 11
	    IF (NTARG_body /= 3) THEN 
! ----------------------------------------------------------------------
! Celestial body's Cartesian coordinates
! Target body's (NTARG) Cartesian coordinates w.r.t. Center body (NCTR)
      NTARG = NTARG_body
      CALL  PLEPH ( ET, NTARG, NCTR, R )
! ----------------------------------------------------------------------
      PRINT *, "NTARG:",NTARG	   
      !PRINT *, "r_body:",R(1),R(2),R(3)	   
      !PRINT *, "v_body:",R(4),R(5),R(6)   

	  
! ----------------------------------------------------------------------
! Cartesian coordinates of the celestial body in meters
! ----------------------------------------------------------------------
! KM to M
	  rbody(1) = R(1) * 1000D0
	  rbody(2) = R(2) * 1000D0
	  rbody(3) = R(3) * 1000D0
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! GM gravity constants of the solar system bodies
! ----------------------------------------------------------------------
      CALL GM_de
      GMbody = GMconst(NTARG)
! ----------------------------------------------------------------------
      !PRINT *,"GM_body", GMbody

  
! ----------------------------------------------------------------------
! Point-mass perturbations vector due to the selected celestial body
! ----------------------------------------------------------------------
      CALL force_gm3rd (rbody,rsat,GMbody , a_perturb)
! ----------------------------------------------------------------------
      PRINT *,"a_perturb", a_perturb



! ----------------------------------------------------------------------
! Sun
      if (NTARG == 11) then
	  rSun = rbody
      GM_sun = GMbody     
! ----------------------------------------------------------------------
! Moon
      else if (NTARG == 10) then
	  rMoon = rbody
      GM_moon = GMbody
      end if
! ----------------------------------------------------------------------
	  
	     END IF
      END DO
	  
! ----------------------------------------------------------------------
      DEALLOCATE (CVAL_2,   STAT = DeAllocateStatus)
      DEALLOCATE (DB_array, STAT = DeAllocateStatus)
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Sun and Moon position vectors in terrestrial reference frame 
! ----------------------------------------------------------------------
! Transformation GCRS to ITRS is not applied in this version:
      rMoon_ITRS = rMoon
      rSun_ITRS = rSun
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Indirect J2 effect of Sun and Moon
! ----------------------------------------------------------------------
      CALL indirectJ2(C20,Re,GM_moon,rMoon_ITRS,GM_sun,rSun_ITRS , a_iJ2)
! ----------------------------------------------------------------------

      PRINT *,"rSun", rSun  
      PRINT *,"GM_sun", GM_sun
      PRINT *,"rMoon", rMoon  
      PRINT *,"GM_moon", GM_moon
      PRINT *,"a_iJ2", a_iJ2

	  
      end
  

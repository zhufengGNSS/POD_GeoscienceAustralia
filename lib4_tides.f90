      program lib4_tides


! ----------------------------------------------------------------------
! Program: lib4_tides.f90
! ----------------------------------------------------------------------
! Purpose:
!  Tidal effects to orbits including Solid Earth tides, Ocean Tides and Pole Tide
!  Computation of the gravitational components due to the tidal effects
! ----------------------------------------------------------------------
! Called subroutines:
!
! ----------------------------------------------------------------------
! Remark 1:
!  Lines 149-158 : 
!  The transformation matrix (GCRS to ITRS) will be provided by the subroutines
!  of the TRS function which is currently under development. Therefore,
!  this will be fulfilled in a future version which will include TRS. 
!
! Remark 2:
!  Line 182 : UT1-UTC time difference has been set temporarily to zero
!  Line 214 : Polar motion xp,yp have been set temporarily to zero
!  xp,yp,UT1-UTC are computed by the subroutines of the TRS function which 
!  is currently under development. Therefore, this will be fulfilled in a 
!  future version which will include TRS.  
! ----------------------------------------------------------------------
! Dr. Thomas Papanikolaou, Geoscience Australia            November 2015
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_planets
      USE mdl_tides
      USE mdl_write
      IMPLICIT NONE


! ----------------------------------------------------------------------
! Variables declaration
! ----------------------------------------------------------------------
! "ET, R(6)" variables enter F77 code and thus, these must be declared in Double Precision.
! ----------------------------------------------------------------------
      DOUBLE PRECISION  ET
      DOUBLE PRECISION  R(6), R_Horizon(6), dR(6)
      INTEGER  NTARG, NCTR
! ----------------------------------------------------------------------
      CHARACTER (LEN=100) :: fname_header,fname_data,fname_out, FESxxfname, filename				
      REAL (KIND = prec_q) :: xsat,ysat,zsat, GMbody, GMearth
      REAL (KIND = prec_q), DIMENSION(3) :: rsat, rSun,rMoon,rMoon_ITRS,rSun_ITRS		
      REAL (KIND = prec_q), DIMENSION(5,5) :: dCnm_solid1,dSnm_solid1
      REAL (KIND = prec_q), DIMENSION(3,3) :: dCnm_solid2,dSnm_solid2	
      REAL (KIND = prec_q) :: dC20_perm
      INTEGER (KIND = prec_int8) :: ocean_Nmax, ocean_Mmax, sz_tides, i_C, j_C
      INTEGER (KIND = prec_int2) :: AllocateStatus, DeAllocateStatus, tide_system
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: fx_solid1, fy_solid1, fz_solid1
      REAL (KIND = prec_q) :: fx_solid2, fy_solid2, fz_solid2
      REAL (KIND = prec_q) :: fx_ocean, fy_ocean, fz_ocean
      REAL (KIND = prec_q) :: fx_tides1, fy_tides1, fz_tides1
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: epoch, xp_mean, yp_mean, xp, yp
      !REAL*8 epoch, xp_mean,yp_mean
      INTEGER error,version
      REAL (KIND = prec_q) :: dC21_pse, dS21_pse, dC21_poc, dS21_poc
      REAL (KIND = prec_q) :: mjd, ut1_utc
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: ax,ay,az
      REAL (KIND = prec_q), DIMENSION(3) :: a_tides, a_solid1, a_solid2, a_ocean, a_pse, a_poc 
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! INPUT
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
      PRINT *,"----------------------------------------------------------------------"											
      PRINT *,"see Remark 1: Transformation matrix from GCRF to ITRF has been excluded"											
      PRINT *,"see Remark 2: Temporary set in variables {xp, yp , UT1-UTC, epoch}"											
      PRINT *,"----------------------------------------------------------------------"											
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Satellite position vector
      xsat = 17835.809256D0 * 1.0D3
      ysat =  4999.483252D0 * 1.0D3
      ysat = 19053.904501D0 * 1.0D3

	  rsat = (/ xsat,ysat,zsat /)
      !PRINT *, "rsat", rsat
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Input epoch
! ----------------------------------------------------------------------
! Julian Day Number
!430  2015.02.01 2457054.5  1  8  4       -0.02577866640782599000
      ET =  2457054.5D0
!      ET =  2457054.78123467D0
      PRINT *, "ET:", ET
! ----------------------------------------------------------------------
! Modified Julian Number	  
      mjd = ET - 2400000.5D0
      PRINT *,"MJD", mjd
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
! Ocean Tides model
	  FESxxfname = 'fes2004_Cnm-Snm.dat'
      ocean_Nmax = 20
	  ocean_Mmax = ocean_Nmax
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! 1. Tide-free: tide_system = 0
! 2. Zero-tide: tide_system = 1
      tide_system = 0
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
! Sun coordinates
      NTARG = 11
      CALL PLEPH ( ET, NTARG, NCTR, R )
	  rSun = (/ R(1),R(2),R(3) /) * 1000D0 ! KM to M
      !PRINT *,"rSun", rSun  
! Moon coordinates
      NTARG = 10
      CALL PLEPH ( ET, NTARG, NCTR, R )
	  rMoon = (/ R(1),R(2),R(3) /) * 1000D0 ! KM to M
      !PRINT *,"rMoon", rMoon  
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
! Tides
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
      ut1_utc = 0.0D0														! Temporary set to 0.0D0 : to be upgraded by using TRS subroutines
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Solid Earth Tides
! ----------------------------------------------------------------------
! Frequency-independent : Step 1 (IERS Conventions 2010)
      CALL tides_solid1(rMoon_ITRS,rSun_ITRS , dCnm_solid1,dSnm_solid1)
      !PRINT *,"dCnm_solid1",dCnm_solid1
      !PRINT *,"dSnm_solid1",dSnm_solid1
! ----------------------------------------------------------------------
! Frequency-dependent : Step 2 (IERS Conventions 2010)
      CALL tides_solid2(mjd, ut1_utc , dCnm_solid2,dSnm_solid2)
      !PRINT *,"dCnm_solid2",dCnm_solid2
      !PRINT *,"dSnm_solid2",dSnm_solid2
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Ocean Tides
! ----------------------------------------------------------------------
      CALL tides_fes2004(FESxxfname)
      CALL tides_ocean(ocean_Nmax, ocean_Mmax, mjd, ut1_utc)
      !PRINT *,"dCnm_ocean",dCnm_ocean
      !PRINT *,"dSnm_ocean",dSnm_ocean
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Pole Tide
! ----------------------------------------------------------------------
! Polar motion (EOP data) :
      xp = 0.0D0															! xp,yp : Temporary set to 0.0D0
      yp = 0.0D0
! ----------------------------------------------------------------------
! Solid Earth pole tide
      CALL tide_pole_se (mjd,xp,yp , dC21_pse, dS21_pse)
! ----------------------------------------------------------------------
! Ocean pole tide
      CALL tide_pole_oc (mjd,xp,yp , dC21_poc, dS21_poc)
! ----------------------------------------------------------------------
!      PRINT *,"dC21_pse, dS21_pse", dC21_pse, dS21_pse
!      PRINT *,"dC21_poc, dS21_poc", dC21_poc, dS21_poc


! ----------------------------------------------------------------------
! Zero Tide term
      CALL tide_perm (dC20_perm)
!      PRINT *, "dC20_perm"  , dC20_perm

      if (tide_system == 1) then
          !dCnm_solid1(2+1,0+1) = dCnm_solid1(2+1,0+1) - dC20_perm
          dCnm_solid2(2+1,0+1) = dCnm_solid2(2+1,0+1) - dC20_perm
      end if	  
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Acceleration cartesian components
! ----------------------------------------------------------------------
      PRINT *,"----------------------------------------------------------------------"											

! ----------------------------------------------------------------------
! Solid Earth Tides : Step 1
      sz_tides = SIZE (dCnm_solid1,DIM=1)
      ALLOCATE (dCnm_tides(sz_tides,sz_tides), STAT = AllocateStatus)
      ALLOCATE (dSnm_tides(sz_tides,sz_tides), STAT = AllocateStatus)
	  dCnm_tides = dCnm_solid1
	  dSnm_tides = dSnm_solid1
      CALL force_tides(rsat,sz_tides-1,sz_tides-1 , ax,ay,az)
      a_solid1 (1) = ax
      a_solid1 (2) = ay
      a_solid1 (3) = az
      PRINT *,"a_solid1", a_solid1
      !PRINT *,"a_solid1 norm", sqrt(a_solid1(1)**2 + a_solid1(2)**2 + a_solid1(3)**2)
      DEALLOCATE (dCnm_tides,   STAT = DeAllocateStatus)
      DEALLOCATE (dSnm_tides,   STAT = DeAllocateStatus)
      sz_tides = 0
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Solid Earth Tides : Step 2
      sz_tides = SIZE (dCnm_solid2,DIM=1)
      ALLOCATE (dCnm_tides(sz_tides,sz_tides), STAT = AllocateStatus)
      ALLOCATE (dSnm_tides(sz_tides,sz_tides), STAT = AllocateStatus)
	  dCnm_tides = dCnm_solid2
	  dSnm_tides = dSnm_solid2
      CALL force_tides(rsat,sz_tides-1,sz_tides-1 , ax,ay,az)
      a_solid2 (1) = ax
      a_solid2 (2) = ay
      a_solid2 (3) = az
      PRINT *,"a_solid2", a_solid2
      !PRINT *,"a_solid2 norm", sqrt(a_solid2(1)**2 + a_solid2(2)**2 + a_solid2(3)**2)
      DEALLOCATE (dCnm_tides,   STAT = DeAllocateStatus)
      DEALLOCATE (dSnm_tides,   STAT = DeAllocateStatus)
      sz_tides = 0
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Ocean Tides
      sz_tides = SIZE (dCnm_ocean,DIM=1)
      ALLOCATE (dCnm_tides(sz_tides,sz_tides), STAT = AllocateStatus)
      ALLOCATE (dSnm_tides(sz_tides,sz_tides), STAT = AllocateStatus)
	  dCnm_tides = dCnm_ocean
	  dSnm_tides = dSnm_ocean
      CALL force_tides(rsat,sz_tides-1,sz_tides-1 , ax,ay,az)
      a_ocean (1) = ax
      a_ocean (2) = ay
      a_ocean (3) = az
      PRINT *,"a_ocean", a_ocean
      !PRINT *,"a_ocean norm", sqrt(a_ocean(1)**2 + a_ocean(2)**2 + a_ocean(3)**2)
      DEALLOCATE (dCnm_tides,   STAT = DeAllocateStatus)
      DEALLOCATE (dSnm_tides,   STAT = DeAllocateStatus)
      sz_tides = 0
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Solid Earth pole tide
      sz_tides = 2 + 1
      ALLOCATE (dCnm_tides(sz_tides,sz_tides), STAT = AllocateStatus)
      ALLOCATE (dSnm_tides(sz_tides,sz_tides), STAT = AllocateStatus)
      DO i_C = 1 , sz_tides
         DO j_C = 1 , sz_tides
            IF (i_C == 2+1 .AND. j_C == 1+1) THEN
		        dCnm_tides (i_C , j_C) = dC21_pse
		        dSnm_tides (i_C , j_C) = dS21_pse
            ELSE
		        dCnm_tides (i_C , j_C) = 0.0D0
		        dSnm_tides (i_C , j_C) = 0.0D0
            END IF
		    !PRINT *,"dCnm_tides ij", i_C, j_C, dCnm_tides (i_C , j_C)
         END DO
      END DO
      !PRINT *,"dCnm_tides", dCnm_tides
      !PRINT *,"dSnm_tides", dSnm_tides
      CALL force_tides (rsat,sz_tides-1,sz_tides-1 , ax,ay,az)
      a_pse (1) = ax
      a_pse (2) = ay
      a_pse (3) = az
      PRINT *,"a_poletide_se:", ax,ay,az
      !PRINT *,"a_poletide_se norm", sqrt(a_pse(1)**2 + a_pse(2)**2 + a_pse(3)**2)
      DEALLOCATE (dCnm_tides,   STAT = DeAllocateStatus)
      DEALLOCATE (dSnm_tides,   STAT = DeAllocateStatus)
      sz_tides = 0
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Ocean pole tide
      sz_tides = 2 + 1
      ALLOCATE (dCnm_tides(sz_tides,sz_tides), STAT = AllocateStatus)
      ALLOCATE (dSnm_tides(sz_tides,sz_tides), STAT = AllocateStatus)
      DO i_C = 1 , sz_tides
         DO j_C = 1 , sz_tides
            IF (i_C == 2+1 .AND. j_C == 1+1) THEN
		        dCnm_tides (i_C , j_C) = dC21_poc
		        dSnm_tides (i_C , j_C) = dS21_poc
            ELSE
		        dCnm_tides (i_C , j_C) = 0.0D0
		        dSnm_tides (i_C , j_C) = 0.0D0
            END IF
		    !PRINT *,"dCnm_tides ij", i_C, j_C, dCnm_tides (i_C , j_C)
         END DO
      END DO
      !PRINT *,"dCnm_tides", dCnm_tides
      !PRINT *,"dSnm_tides", dSnm_tides
      CALL force_tides (rsat,sz_tides-1,sz_tides-1 , ax,ay,az)
      a_poc (1) = ax
      a_poc (2) = ay
      a_poc (3) = az
      PRINT *,"a_poletide_oc:", ax,ay,az
      !PRINT *,"a_poletide_oc norm", sqrt(a_poc(1)**2 + a_poc(2)**2 + a_poc(3)**2)
      DEALLOCATE (dCnm_tides,   STAT = DeAllocateStatus)
      DEALLOCATE (dSnm_tides,   STAT = DeAllocateStatus)
      sz_tides = 0
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Tides acceleration cartesian components : Overall effects
! ----------------------------------------------------------------------
      a_tides = a_solid1 + a_solid2 + a_ocean + a_pse + a_poc
      PRINT *,"a_tides", a_tides
! ----------------------------------------------------------------------

      PRINT *,"----------------------------------------------------------------------"											



! ----------------------------------------------------------------------
! Write dCnm,dSnm arrays to files
! ----------------------------------------------------------------------
      ALLOCATE (wrtArray(5,5), STAT = AllocateStatus)
	  wrtArray = dCnm_solid1
      filename = 'dCnm_solid1.wrt'
      CALL write_array (filename)	  
	  wrtArray = dSnm_solid1
      filename = 'dSnm_solid1.wrt'
      CALL write_array (filename)	
      DEALLOCATE (wrtArray,   STAT = DeAllocateStatus)
! ----------------------------------------------------------------------
      ALLOCATE (wrtArray(3,3), STAT = AllocateStatus)
	  wrtArray = dCnm_solid2
      filename = 'dCnm_solid2.wrt'
      CALL write_array (filename)	  
	  wrtArray = dSnm_solid2
      filename = 'dSnm_solid2.wrt'
      CALL write_array (filename)	  
      DEALLOCATE (wrtArray,   STAT = DeAllocateStatus)
! ----------------------------------------------------------------------
      ALLOCATE (wrtArray(ocean_Nmax+1,ocean_Nmax+1), STAT = AllocateStatus)
	  wrtArray = dCnm_ocean
      filename = 'dCnm_ocean.wrt'
      CALL write_array (filename)	  
	  wrtArray = dSnm_ocean
      filename = 'dSnm_ocean.wrt'
      CALL write_array (filename)	  
      DEALLOCATE (wrtArray,   STAT = DeAllocateStatus)
! ----------------------------------------------------------------------
	

	
      end

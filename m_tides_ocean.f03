MODULE m_tides_ocean


! ----------------------------------------------------------------------
! MODULE: m_tides_ocean.f03
! ----------------------------------------------------------------------
! Purpose:
!  Module for calling the tides_ocean subroutine 
! 
! Subroutines contained within the module:
! - tides_ocean
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
!			Cooperative Research Centre for Spatial Information, Australia
! Created:	16 November 2017
! ----------------------------------------------------------------------


      IMPLICIT NONE
      !SAVE 			
  
	  
Contains


SUBROUTINE tides_ocean(n_max, m_max , mjd, ut1_utc, dCnm_ocean, dSnm_ocean)


! ----------------------------------------------------------------------
! SUBROUTINE: tides_ocean
! ----------------------------------------------------------------------
! Purpose:
!  Ocean Tides perturbations according to the formula of the IERS Conv. 2010
!  Computation of the variations in the Stokes coefficients are based on
!  an ocean tide model.
!
! ----------------------------------------------------------------------
! Input arguments:
! - n_max :			Degree limit of the ocean tides coefficients
! - m_max :			Order limit of the ocean tides coefficients
! - mjd   :			Modified Julian Day in Terrestrial Time (TT) scale including the fraction 
!             		of the day 
! - ut1_utc:   		Time difference between UT1 and UTC (seconds)
!
! Allocatable arrays called via module : mdl_tides.f90
! - Delaunay_FES : Delaunay variables multipliers for FES2004 tidal waves
! - dCnm_p		 : Spherical harmonic coefficients of ocean tide model to
!   dSnm_p         be used in Eq. 6.15 of IERS Conventions 2010
!   dCnm_m         Coefficients are stored in four 3D arrays
!   dSnm_m
!
! Output arguments:
! - dCnm_ocean :	Geopotential Cnm corrections matrix due to Ocean Tides
! - dSnm_ocean :	Geopotential Snm corrections matrix due to Ocean Tides
! ----------------------------------------------------------------------
! Remark 1:
!  Computed dCnm and dSnm are formed into lower triangular matrices.
!  Coefficient dCnm corresponds to matrix element dCnm(n+1,m+1).
!
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou, Cooperative Research Centre for Spatial Information, Australia
! Created:	12 November 2015
! ----------------------------------------------------------------------
! Last modified:
! - Dr. Thomas Papanikolaou, 27 October 2017:
!   The subroutine has been modified in order to call the subroutine gmst_iers.f03 
!   for the computation of the Greenwich mean sidereal time
! - Dr. Thomas Papanikolaou, 16 November 2017:
!	Upgraded from Fortran 90 to Fortran 2003 for considering the
!	Fortran 2003 advantages in dynamic memory allocation
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      USE mdl_tides
      IMPLICIT NONE

	  
! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
! IN
      INTEGER (KIND = prec_int8), INTENT(IN) :: n_max, m_max
      REAL (KIND = prec_q), INTENT(IN) :: mjd, ut1_utc
	  
!	  "Delaunay_FES, dCnm_p,dSnm_p, dCnm_m,dSnm_m" variables are declared in module "mdl_tides.f90"

! OUT
! 	  "dCnm_ocean, dSnm_ocean" variables are declared as allocatable arrays in module "mdl_tides.f90"
      REAL (KIND = prec_q), INTENT(OUT), DIMENSION(:,:), ALLOCATABLE :: dCnm_ocean, dSnm_ocean		
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: thetag, thetaf
      REAL (KIND = prec_q) :: F1,F2,F3,F4,F5
      REAL (KIND = prec_d), DIMENSION(5) :: delaunay_arr
      INTEGER (KIND = prec_int4) :: Nfrq, i, j, ifrq, AllocateStatus
      REAL (KIND = prec_q) :: dCnm_fi,dSnm_fi, dCnm_f,dSnm_f, matrix_mult
      INTEGER (KIND = prec_int8) :: n, m, m_limit
      REAL (KIND = prec_q) :: pi
! ----------------------------------------------------------------------


pi = PI_global																	


! ----------------------------------------------------------------------
! Delaunay variables (in radians)
CALL delaunay (mjd , F1,F2,F3,F4,F5)
delaunay_arr = (/ F1 , F2 , F3 , F4 , F5 /)							
! ----------------------------------------------------------------------
!PRINT *,"Delaunay",delaunay_arr


! ----------------------------------------------------------------------  
! Greenwich Mean Sidereal Time (GMST) in radians
Call gmst_iers(mjd, ut1_utc , thetag)
! ----------------------------------------------------------------------  
!PRINT *,"thetag",thetag


! ----------------------------------------------------------------------
Nfrq = SIZE (Delaunay_FES,DIM=1)
!sz2 = SIZE (Delaunay_FES,DIM=2)
!PRINT *,"tides_ocean.f90 Nfrq", Nfrq


! ---------------------------------------------------------------------------
! Allocate Arrays	  
! ---------------------------------------------------------------------------
      ALLOCATE (dCnm_ocean(n_max+1,n_max+1), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT *, "Error: Not enough memory"
         PRINT *, "Error: SUBROUTINE tides_ocean.f90"
         PRINT *, "Error: Allocatable Array: dCnm_ocean | Nmax:", n_max
!         STOP "*** Not enough memory ***"
      END IF  

      ALLOCATE (dSnm_ocean(n_max+1,n_max+1), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT *, "Error: Not enough memory"
         PRINT *, "Error: SUBROUTINE tides_ocean.f90"
         PRINT *, "Error: Allocatable Array: dSnm_ocean | Nmax:", n_max
!         STOP "*** Not enough memory ***"
      END IF  
! ---------------------------------------------------------------------------


! ---------------------------------------------------------------------------
DO n = 2 , n_max
    IF (n > m_max) THEN
        m_limit = m_max
    ELSE
        m_limit = n
    END IF
    DO m = 0 , m_limit
        dCnm_f = 0.0D0
        dSnm_f = 0.0D0
        DO ifrq = 1 , Nfrq
			! ifrq.NE.5 .and. 
			IF (ifrq.NE.9 .and. ifrq.NE.10 .and. ifrq.NE.11 .and. ifrq.NE.12 .and. ifrq.NE.13 .and. ifrq.NE.Nfrq) then
            ! thetaf (in radians)
            matrix_mult = 0.0D0
            DO j = 1 , 5
               matrix_mult = matrix_mult + ( Delaunay_FES(ifrq,j+1) * delaunay_arr(j) )
            END DO
			thetaf = m * (thetag + pi) - 1.0D0 * matrix_mult														 		

            dCnm_fi = (dCnm_p(n+1,m+1,ifrq) + dCnm_m(n+1,m+1,ifrq)) * cos(thetaf) &
			        + (dSnm_p(n+1,m+1,ifrq) + dSnm_m(n+1,m+1,ifrq)) * sin(thetaf)
            dSnm_fi = (dSnm_p(n+1,m+1,ifrq) - dSnm_m(n+1,m+1,ifrq)) * cos(thetaf) &
			        - (dCnm_p(n+1,m+1,ifrq) - dCnm_m(n+1,m+1,ifrq)) * sin(thetaf)
            dCnm_f = dCnm_f + dCnm_fi
            dSnm_f = dSnm_f + dSnm_fi	
			END IF
        END DO
		
        ! IERS Conventions update 10/08/2012 Section 6.3.2
        IF (m == 0) THEN
          dSnm_f = 0
        END IF
		
        dCnm_ocean(n+1,m+1) = dCnm_f
        dSnm_ocean(n+1,m+1) = dSnm_f
    END DO
END DO
! ---------------------------------------------------------------------------


END SUBROUTINE



END Module

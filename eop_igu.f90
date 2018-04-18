SUBROUTINE eop_igu (mjd, ERP_fname, EOP_fname, EOP_int)


! ----------------------------------------------------------------------
! Subroutine:  eop_igu.f90
! ----------------------------------------------------------------------
! Purpose:
!  EOP data reading and processing by using:
!  - ERP (Earth Rotation Parameters) data from the ultra-rapid products 
!    provided by the IGS (International GNSS Service).
!  - Corrections to the precession-nutation model are obtained from the 
!    daily solutions (finals2000A.daily) provided by the International 
!    Earth Rotation Service and Reference Systems (IERS) 
!    Rapid Service/Prediction Center (RS/PC) 
! ----------------------------------------------------------------------
! Input arguments:
! - mjd:			Modified Julian Day number at the required epoch
!					(including fraction of the day)
! - ERP_fname:		IGS ultra-rapid ERP data file name e.g. igu18861_00.erp
! - EOP_fname:		IERS RS/PC EOP data file name e.g. finals2000A.daily
!
! Output arguments:
! - eop_int:		EOP data array at the input epoch
!   				eop_int = [MJD xp yp UT1_UTC LOD dX dY] 
!   				MJD:     MJD at the input epoch (including fraction of the day)
!   				x,y:     Polar motion coordinates (arcsec) 
!   				UT1_UTC: Difference between UT1 and UTC (sec)
!					dX,dY:   Corrections to Precession-Nutation model (arcsec)
! ----------------------------------------------------------------------
! Dr. Thomas Papanikolaou, Geoscience Australia               March 2016
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      IMPLICIT NONE

! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
! IN
      REAL (KIND = prec_d), INTENT(IN) :: mjd
      CHARACTER (LEN=50), INTENT(IN) :: ERP_fname, EOP_fname
! OUT
      REAL (KIND = prec_d), INTENT(OUT) :: EOP_int(7)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      REAL (KIND = prec_d) :: ERP_igu_data(2,5), EOP_data(7), ERP_int(5)
      LOGICAL :: igu_flag
      INTEGER (KIND = prec_int8) :: mjd_UTC_day
      REAL (KIND = prec_d) :: mjd_ar(2), Xpole_ar(2), Ypole_ar(2), UT1UTC_ar(2), LOD_ar(2)
      REAL (KIND = prec_d) :: mjd_int, Xpole_int, Ypole_int, UT1UTC_int, LOD_int
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! ERP data reading 
      CALL erp_igu (ERP_fname, mjd, ERP_igu_data, igu_flag)
      !if (igu_flag == .FALSE.) then
      if (igu_flag .EQV. .FALSE.) then
         PRINT *,"--------------------------------------------------------"
         PRINT *, "Warning error: Subroutine erp_igu.f90"
         PRINT *, "Input epoch is out of the range covered by the IGS ultra-rapid ERP file" 
         PRINT *, "Check the input iguwwwwd_hh.erp file"
         PRINT *,"--------------------------------------------------------"
         STOP  ! END PROGRAM
      end if		 
! ----------------------------------------------------------------------

 
! ----------------------------------------------------------------------
! ERP interpolation
      mjd_int = mjd
      mjd_ar    = ERP_igu_data(1:2,1)
      Xpole_ar  = ERP_igu_data(1:2,2)
      Ypole_ar  = ERP_igu_data(1:2,3)
      UT1UTC_ar = ERP_igu_data(1:2,4)
      LOD_ar    = ERP_igu_data(1:2,5)
      CALL interp_lin (mjd_ar, Xpole_ar , mjd_int, Xpole_int)
      CALL interp_lin (mjd_ar, Ypole_ar , mjd_int, Ypole_int)
      CALL interp_lin (mjd_ar, UT1UTC_ar, mjd_int, UT1UTC_int)
      CALL interp_lin (mjd_ar, LOD_ar   , mjd_int, LOD_int)
      ERP_int (1) = mjd_int
      ERP_int (2) = Xpole_int
      ERP_int (3) = Ypole_int
      ERP_int (4) = UT1UTC_int
      ERP_int (5) = LOD_int
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! dX,dY : Corrections w.r.t Precession-Nutation model
      mjd_UTC_day = INT (mjd)
      CALL eop_finals2000A (EOP_fname, mjd_UTC_day , EOP_data)
      !dX = EOP_data(6) 
      !dY = EOP_data(7)  
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
      EOP_int (1) = mjd_int
      EOP_int (2) = Xpole_int
      EOP_int (3) = Ypole_int
      EOP_int (4) = UT1UTC_int
      EOP_int (5) = LOD_int
	  ! dX,dY : Precession-Nutation model corrections
      EOP_int (6) = EOP_data (6) 
      EOP_int (7) = EOP_data (7) 
! ----------------------------------------------------------------------


!      PRINT *,"--------------------------------------------------------"
!      PRINT *, "ERP_igu_data"
!      PRINT *, ERP_igu_data
!      PRINT *, "ERP_int"
!      PRINT *, ERP_int
!      PRINT *,"--------------------------------------------------------"


	  
END

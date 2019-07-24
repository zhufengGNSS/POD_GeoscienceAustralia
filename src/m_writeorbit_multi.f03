MODULE m_writeorbit_multi


! ----------------------------------------------------------------------
! MODULE: m_writeorbit_multi.f03
! ----------------------------------------------------------------------
! Purpose:
!  Module for write orbit array data to output (ascii) files 
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
!			Geoscience Australia, Frontier-SI
! Created:	21 March 2019
! ----------------------------------------------------------------------


      IMPLICIT NONE
      !SAVE 			

	  
Contains


SUBROUTINE writeorbit_multi (orbitsmatrix_crf,orbitsmatrix_trf, orbits_ics_icrf, PRN_array, filename)

! ----------------------------------------------------------------------
! SUBROUTINE: writeorbit_multi 
! ----------------------------------------------------------------------
! Purpose:
!  Write orbit and partial derivatives matrices to an output ascii file
!
!   writeorbit.f03 subroutine has been modified in order to write the 
!	estimated orbits and partial derivatives (solution of the Variational Equations) 
!   to output file (ascii) based on an internal adopted orbit format:
!   {MJD Sec_00h r(XYZ) v(XYZ) State_Transition_Matrix Sensitivity_matrix} 
! ----------------------------------------------------------------------
! Input arguments:
! - wrtArray:       Input allocatable array
! - filename:       file name to be used for writing array data
!
! Output arguments:
!
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
!			Geoscience Australia, Frontier-SI
! Created:	21 March 2019
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      USE mdl_config
      USE mdl_param
      IMPLICIT NONE
	  
! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
! IN
      REAL (KIND = prec_q), INTENT(IN), DIMENSION(:,:,:), ALLOCATABLE :: orbitsmatrix_crf
      REAL (KIND = prec_q), INTENT(IN), DIMENSION(:,:,:), ALLOCATABLE :: orbitsmatrix_trf
      REAL (KIND = prec_q), INTENT(IN), DIMENSION(:,:), ALLOCATABLE :: orbits_ics_icrf
      CHARACTER (LEN=3), ALLOCATABLE :: PRN_array(:)
      CHARACTER (LEN=100), INTENT(IN) :: filename
! OUT
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: i, i_write
      INTEGER (KIND = prec_int2) :: UNIT_IN, ios, ios_ith
      INTEGER (KIND = prec_int8) :: sz1, sz2, sz3
      INTEGER (KIND = prec_int2) :: wrt_opt
      INTEGER (KIND = prec_int2) :: FMT_opt
! ----------------------------------------------------------------------
      CHARACTER (LEN=1) :: RealT
      INTEGER (KIND = prec_int2) :: RealW, RealD
      CHARACTER (LEN=70) :: fmt_wrt, fmt_wrt0, fmt_sz2
      REAL (KIND = prec_q) :: wrtArrayLN 
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: Nepochs, Nparam, Nsat 
      INTEGER (KIND = prec_int8) :: i_epoch, i_sat
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
UNIT_IN = 7  												
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Orbit arrays dimensions
sz1 = SIZE (orbitsmatrix_crf,DIM=1)
sz2 = SIZE (orbitsmatrix_crf,DIM=2)
sz3 = SIZE (orbitsmatrix_crf,DIM=3)

Nepochs = sz1
Nparam  = sz2
Nsat    = sz3
   
! PRN
Nsat = SIZE (PRN_array,DIM=1)
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Format definition
! ----------------------------------------------------------------------
! Orbit format: {MJD Sec_00h r(XYZ) v(XYZ)} 
!fmt_wrt = '(F25.15,F25.9,3F25.3,3F25.7)'
! Orbit-VEQ format: {PRN MJD Sec_00h r(XYZ) v(XYZ) VEQ-Z VEQ-P} 
!fmt_wrt = '(5A3,F25.12,F25.9,3F25.4,3F25.9, F25)'
fmt_wrt = '(A3,A1,F25.12,F25.9,3F25.4,3F25.9, A)'
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Open file
      OPEN (UNIT=UNIT_IN,FILE=filename,ACTION="WRITE",POSITION="REWIND", IOSTAT=ios)
      IF (ios /= 0) THEN
         PRINT *, "Error in opening file:", filename
         PRINT *, "OPEN IOSTAT=", ios
      END IF
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Thomas - Please complete adding the following header information to the output tabular file
! Write file header information
! ----------------------------------------------------------------------
WRITE (UNIT=UNIT_IN,FMT='(a)'              ,IOSTAT=ios_ith) '#INFO    POD tabular ephemeris file: '
WRITE (UNIT=UNIT_IN,FMT='(a,a)'            ,IOSTAT=ios_ith) '#INFO    Generated from: ',trim(pseudobs_orbit_filename_cfg)
WRITE (UNIT=UNIT_IN,FMT='(a,f19.12,f13.6)' ,IOSTAT=ios_ith) '#INFO    Epoch initial conditions:   ',orbitsmatrix_crf(1,1:2,1)
WRITE (UNIT=UNIT_IN,FMT='(a,f19.12,f13.6)' ,IOSTAT=ios_ith) '#INFO    Epoch Start:                ',orbitsmatrix_crf(1,1:2,1)
WRITE (UNIT=UNIT_IN,FMT='(a,f19.12,f13.6)' ,IOSTAT=ios_ith) '#INFO    Epoch End:                  ',orbitsmatrix_crf(Nepochs,1:2,1)
WRITE (UNIT=UNIT_IN,FMT='(a,i5)'           ,IOSTAT=ios_ith) '#INFO    Tabular interval (sec):     ',900
WRITE (UNIT=UNIT_IN,FMT='(a,i5)'           ,IOSTAT=ios_ith) '#INFO    Number of Epochs:           ',Nepochs
WRITE (UNIT=UNIT_IN,FMT='(a)'              ,IOSTAT=ios_ith) '#INFO    Model information:           [TIME_SYS] [GRAV_MODEL]&
                                           &[PN_MODEL] [EPH_MODEL] [ALBEDO_MODEL]'
WRITE (UNIT=UNIT_IN,FMT='(a,a)'            ,IOSTAT=ios_ith) '#INFO    EOP file: ',trim(EOP_fname_cfg)
WRITE (UNIT=UNIT_IN,FMT='(a,i3)'           ,IOSTAT=ios_ith) '#INFO    Number of Satellites:       ',Nsat
WRITE (UNIT=UNIT_IN,FMT='(a,i4)'           ,IOSTAT=ios_ith) '#INFO    Number of Partials:         ',Nparam-8
WRITE (UNIT=UNIT_IN,FMT='(a)'              ,IOSTAT=ios_ith) '#INFO    Satellite ICS:              '
DO i_sat = 1 , Nsat
   WRITE (UNIT=UNIT_IN,FMT='(a,a3,a,i3,a)' ,IOSTAT=ios_ith) '#IC_INFO ',PRN_array(i_sat),' [SVN] [BLK_TYP] [ANT_TH] &
                                           &[ECOM1] ',(Nparam-8)/6, '[IC PARAM LIST - X Y Z XD YD ZD DRAD .....]'
   WRITE (UNIT=UNIT_IN,FMT='(a,a3,a,1x,f14.4,f14.6,1x,15(d17.10,1x))',IOSTAT=ios_ith) '#IC_XYZ  ',PRN_array(i_sat),&
                                           &' [SVN] [BLK_TYP] ICRF ',orbits_ics_icrf(:,i_sat)   
! orbitsmatrix_crf(1,3:8,i_sat), ' DR YR BR DC DS YC YS BC BS''(a3,1x,f14.4,f14.6,1x,15(d17.10,1x))'
END DO 
WRITE (UNIT=UNIT_IN,FMT='(a)'              ,IOSTAT=ios_ith) '#INFO PRN MJD SOD, ICRF [X Y Z ZD YD ZD], ITRF [X Y Z XD YD ZD], &
                                           &Partials [dx/dX dx/dY dx/dZ dx/dXD dx/dYD dx/dZD dY/dX dY/dY dY/dZ &
                                           &dY/dXD dY/dYD dY/dZD ... dx/dRAD1,dx/dRAD1,dx/dRAD3,dx/dRAD3 dx/dRAD4 &
                                           &... dx/dRADN ... dy/dRAD1,dx/dRAD2 ... dy/dRADN ... ]'   
                                                       
! ----------------------------------------------------------------------
! Write data to file | Line by line	  
! ----------------------------------------------------------------------
! Write orbit-partials matrix per epoch per satellite
! ----------------------------------------------------------------------
DO i_epoch = 1 , Nepochs
	DO i_sat = 1 , Nsat

!print *,"PRN_array(i_sat)", PRN_array(i_sat)	
!print *,"orbitsmatrix_crf(i_epoch,:,i_sat)", orbitsmatrix_crf(i_epoch,:,i_sat)
	
! Based on the format definition by fmt_wrt 
!WRITE (UNIT=UNIT_IN,FMT=fmt_wrt,IOSTAT=ios_ith) PRN_array(i_sat),' ', orbitsmatrix_crf(i_epoch,:,i_sat)
WRITE (UNIT=UNIT_IN,FMT= * ,IOSTAT=ios_ith) PRN_array(i_sat),' ', orbitsmatrix_crf(i_epoch,1:8,i_sat), &
                                            orbitsmatrix_trf(i_epoch,3:8,i_sat),orbitsmatrix_crf(i_epoch,9:NParam,i_sat)
		 
IF (ios_ith /= 0) THEN
   PRINT *, "Error in writing to file: ", TRIM (filename)
   PRINT *, "WRITE IOSTAT=", ios_ith
END IF
	
	END DO
END DO


ENDFILE (UNIT = UNIT_IN) 
CLOSE (UNIT = UNIT_IN)
! ----------------------------------------------------------------------



END SUBROUTINE



End Module


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


SUBROUTINE writeorbit_multi (orbitsmatrix, PRN_array, filename)

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
      IMPLICIT NONE
	  
! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
! IN
      REAL (KIND = prec_q), INTENT(IN), DIMENSION(:,:,:), ALLOCATABLE :: orbitsmatrix 
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
sz1 = SIZE (orbitsmatrix,DIM=1)
sz2 = SIZE (orbitsmatrix,DIM=2)
sz3 = SIZE (orbitsmatrix,DIM=3)

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
! Write data to file | Line by line	  
! ----------------------------------------------------------------------
! Write orbit-partials matrix per epoch per satellite
! ----------------------------------------------------------------------
DO i_epoch = 1 , Nepochs
	DO i_sat = 1 , Nsat

!print *,"PRN_array(i_sat)", PRN_array(i_sat)	
!print *,"orbitsmatrix(i_epoch,:,i_sat)", orbitsmatrix(i_epoch,:,i_sat)
	
! Based on the format definition by fmt_wrt 
!WRITE (UNIT=UNIT_IN,FMT=fmt_wrt,IOSTAT=ios_ith) PRN_array(i_sat),' ', orbitsmatrix(i_epoch,:,i_sat)
WRITE (UNIT=UNIT_IN,FMT= * ,IOSTAT=ios_ith) PRN_array(i_sat),' ', orbitsmatrix(i_epoch,:,i_sat)
		 
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


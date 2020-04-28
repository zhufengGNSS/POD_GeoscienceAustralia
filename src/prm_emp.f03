SUBROUTINE prm_emp (PRMfname)


! ----------------------------------------------------------------------
! SUBROUTINE: prm_emp.f03
! ----------------------------------------------------------------------
! Purpose:
!  Read the configaration of the orbit emprical parameters 
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
!			Geoscience Australia, CRC-SI
! Created:	3 September 2018
! ----------------------------------------------------------------------
	  
	  
      USE mdl_precision
      USE mdl_num
      USE mdl_param
      USE mdl_config
      IMPLICIT NONE


! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
! IN
      CHARACTER (LEN=100), INTENT(IN) :: PRMfname				
! OUT

! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Variables declaration
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int2) :: ForceMod  
! ----------------------------------------------------------------------
      CHARACTER (LEN=100) :: filename				
      INTEGER (KIND = prec_int2) :: AllocateStatus,DeAllocateStatus
      INTEGER (KIND = prec_int8) :: i, read_i
      INTEGER (KIND = prec_int2) :: UNIT_IN, ios
      INTEGER (KIND = prec_int2) :: ios_line, ios_key, ios_data
      INTEGER (KIND = prec_int2) :: space_i
      CHARACTER (LEN=7) :: Format1, Format2, Format3
      CHARACTER (LEN=500) :: line_ith	  
      CHARACTER (LEN=150) :: word1_ln, word_i, t0	  
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int2) :: bias_r, bias_t, bias_n
      INTEGER (KIND = prec_int2) :: cpr_r, cpr_t, cpr_n, cpr_freq
      REAL (KIND = prec_q) :: Bias_radial, Bias_along, Bias_cross
      REAL (KIND = prec_q) :: Cterm, Sterm



! ----------------------------------------------------------------------
! Orbit parameterization INPUT file read:
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
      UNIT_IN = 9  												
      Format1 = '(A)'
      Format2 = '(F)'
      Format3 = '(I100)'
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Open .in file
      OPEN (UNIT = UNIT_IN, FILE = TRIM (PRMfname), IOSTAT = ios)
      IF (ios /= 0) THEN
         PRINT *, "Error in opening file:", PRMfname
         PRINT *, "OPEN IOSTAT=", ios
      END IF
! ----------------------------------------------------------------------
! Read EMP parameterizations
! ----------------------------------------------------------------------
! Read input file
i = 0
DO

READ (UNIT=UNIT_IN,FMT=Format1,IOSTAT=ios_line) line_ith
i = i + 1
!PRINT *, "READ Line (i,ios):", i, ios_line

! ----------------------------------------------------------------------
! End of file
         IF (ios_line < 0) THEN
!            PRINT *, "End of file, i=", i
            EXIT		
         END IF
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! 1st Word of Line ith
READ (line_ith, * , IOSTAT=ios_data) word1_ln  ! 1st word
!READ (line_ith, * , IOSTAT=ios_data) word1_ln, charN 
! ----------------------------------------------------------------------
!PRINT *, "word1_ln: ", word1_ln


! ----------------------------------------------------------------------
! Parameters Keywords read 
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Empirical parameters related to orbit force modelling
! ----------------------------------------------------------------------
!IF (word1_ln == "EMP_param") THEN
!   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, EMP_param_glb 
!END IF
! ----------------------------------------------------------------------
! Empirical parameters: Bias accelerations 
! 
! EMP_Bias_glb(1) : Bias acceleration in Radial direction
! EMP_Bias_glb(2) : Bias acceleration in Tangential direction
! EMP_Bias_glb(3) : Bias acceleration in Normal direction
! 
! 1. EMP_Bias_glb(i)=1 : Effect is considered
! 2. EMP_Bias_glb(i)=0 : Effect is neglected 
!
! Radial direction
IF (word1_ln == "bias_r") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, bias_r 
	EMP_Bias_glb(1) = bias_r
END IF
! Tangential direction
IF (word1_ln == "bias_t") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, bias_t 
	EMP_Bias_glb(2) = bias_t
END IF
! Normal direction
IF (word1_ln == "bias_n") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, bias_n 
	EMP_Bias_glb(3) = bias_n
END IF
! ----------------------------------------------------------------------
! Empirical parameters: Cycle per revolution accelerations
! Number of cycles per revolution
IF (word1_ln == "cpr_no") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, cpr_freq 
	EMP_nCPR_glb = cpr_freq
END IF
! 
! EMP_CPR_glb(1) : Bias acceleration in Radial direction
! EMP_CPR_glb(2) : Bias acceleration in Tangential direction
! EMP_CPR_glb(3) : Bias acceleration in Normal direction
! 
! 1. EMP_CPR_glb(i)=1 : Effect is considered
! 2. EMP_CPR_glb(i)=0 : Effect is neglected 
!
! Radial direction
IF (word1_ln == "cpr_r") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, cpr_r 
	EMP_CPR_glb(1) = cpr_r
END IF
! Tangential direction
IF (word1_ln == "cpr_t") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, cpr_t 
	EMP_CPR_glb(2) = cpr_t
END IF
! Normal direction
IF (word1_ln == "cpr_n") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, cpr_n 
	EMP_CPR_glb(3) = cpr_n
END IF
! ----------------------------------------------------------------------

END DO
CLOSE (UNIT=UNIT_IN)
! Close of input parameterization file
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Empirical parameters :: Bias acceleration :: Read input file
! ----------------------------------------------------------------------
! Bias acceleration vector in orbital frame
! Read bias input file
filename = 'emp_bias.in'
! ----------------------------------------------------------------------
! Open .in file
OPEN (UNIT = UNIT_IN, FILE = TRIM (filename), IOSTAT = ios)
IF (ios /= 0) THEN
   PRINT *, "Error in opening file:", filename
   PRINT *, "OPEN IOSTAT=", ios
END IF

! Read input file
i = 0
DO
READ (UNIT=UNIT_IN,FMT=Format1,IOSTAT=ios_line) line_ith
i = i + 1
!PRINT *, "READ Line (i,ios):", i, ios_line

! End of file
IF (ios_line < 0) THEN
	! PRINT *, "End of file, i=", i
   EXIT		
END IF

! 1st Word of Line ith
READ (line_ith, * , IOSTAT=ios_data) word1_ln  ! 1st word

! Empirical parameters :: Bias accelerations
IF (word1_ln == "bias_rtn") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, Bias_radial, Bias_along, Bias_cross  
   Bias_accel_glb(1) = Bias_radial
   Bias_accel_glb(2) = Bias_along
   Bias_accel_glb(3) = Bias_cross
   !print *,"Bias_accel_glb", Bias_accel_glb
END IF

END DO
CLOSE (UNIT=UNIT_IN)
! Close of input file
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Empirical parameters :: Cycle-per-revolution :: Read input file
! ----------------------------------------------------------------------
! Cycle-per-revolution acceleration vector in orbital frame
! Read input file
filename = 'emp_cpr.in'
! ----------------------------------------------------------------------
! Open .in file
OPEN (UNIT = UNIT_IN, FILE = TRIM (filename), IOSTAT = ios)
IF (ios /= 0) THEN
   PRINT *, "Error in opening file:", filename
   PRINT *, "OPEN IOSTAT=", ios
END IF

! Read input file
i = 0
DO
READ (UNIT=UNIT_IN,FMT=Format1,IOSTAT=ios_line) line_ith
i = i + 1
!PRINT *, "READ Line (i,ios):", i, ios_line

! End of file
IF (ios_line < 0) THEN
	! PRINT *, "End of file, i=", i
   EXIT		
END IF

! 1st Word of Line ith
READ (line_ith, * , IOSTAT=ios_data) word1_ln  ! 1st word

! Empirical parameters :: CPR coefficients
IF (word1_ln == "cpr_radial") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, Cterm, Sterm
   CPR_CS_glb(1,1) = Cterm
   CPR_CS_glb(1,2) = Sterm
END IF
IF (word1_ln == "cpr_along") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, Cterm, Sterm  
   CPR_CS_glb(2,1) = Cterm
   CPR_CS_glb(2,2) = Sterm
END IF
IF (word1_ln == "cpr_cross") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, Cterm, Sterm  
   CPR_CS_glb(3,1) = Cterm
   CPR_CS_glb(3,2) = Sterm
END IF


END DO
CLOSE (UNIT=UNIT_IN)
! Close of input file
! ----------------------------------------------------------------------


END

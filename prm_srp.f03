SUBROUTINE prm_srp (PRMfname)


! ----------------------------------------------------------------------
! SUBROUTINE: prm_srp.f03
! ----------------------------------------------------------------------
! Purpose:
!  Read the configaration of the SRP parameters 
! ----------------------------------------------------------------------
! Author :	Tzupang Tseng,  Thomas Papanikolaou
!
! Copyright:	GEOSCIENCE AUSTRALIA 
!
! Created:	13 DEC 2018
!
! Changes:      13-12-2018 Dr. Tzupang Tseng : added ECOM-based options in the configuration file
!                                              for ACS POD
!
! ----------------------------------------------------------------------
	  
	  
      USE mdl_precision
      USE mdl_num
      USE mdl_param
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
      INTEGER (KIND = prec_int2) :: PD_Param_ID
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
      INTEGER (KIND = prec_int2) :: bias_D, bias_Y, bias_B
      INTEGER (KIND = prec_int2) :: cpr_D, cpr_Y, cpr_B, cpr_freq
      INTEGER (KIND = prec_int2) :: cpr_D2, cpr_D4
      REAL (KIND = prec_q) :: Bias_radial, Bias_along, Bias_cross
      REAL (KIND = prec_q) :: Cterm, Sterm
      REAL (KIND = prec_q) :: D0, Y0, B0, DC, DS, YC, YS, BC, BS
      REAL (KIND = prec_q) :: D2C, D2S, D4C, D4S

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
! ECOM parameters related to orbit force modelling
! ----------------------------------------------------------------------
IF (word1_ln == "ECOM_param") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, ECOM_param_glb 
END IF
! ----------------------------------------------------------------------
! Empirical parameters: Bias accelerations 
! 
! ECOM_Bias_glb(1) : D0 acceleration 
! ECOM_Bias_glb(2) : Y0 acceleration 
! ECOM_Bias_glb(3) : B0 acceleration 
! 
! 1. ECOM_Bias_glb(i)=1 : Effect is considered
! 2. ECOM_Bias_glb(i)=0 : Effect is neglected 
!
! ECOM1 is applied
! **********************************************************************
IF (ECOM_param_glb == 1) THEN

! D direction
IF (word1_ln == "bias_D") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, bias_D 
	ECOM_Bias_glb(1) = bias_D
END IF
! Y direction
IF (word1_ln == "bias_Y") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, bias_Y 
	ECOM_Bias_glb(2) = bias_Y
END IF
! B direction
IF (word1_ln == "bias_B") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, bias_B 
	ECOM_Bias_glb(3) = bias_B
END IF
! ----------------------------------------------------------------------
! ECOM_CPR_glb(1) : CPR acceleration in D direction
! ECOM_CPR_glb(2) : CPR acceleration in Y direction
! ECOM_CPR_glb(3) : CPR acceleration in B direction
! 
! 1. ECOM_CPR_glb(i)=1 : Effect is considered
! 2. ECOM_CPR_glb(i)=0 : Effect is neglected 
!
! D direction
IF (word1_ln == "cpr_D") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, cpr_D 
	ECOM_CPR_glb(1) = cpr_D
END IF
! Y direction
IF (word1_ln == "cpr_Y") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, cpr_Y 
	ECOM_CPR_glb(2) = cpr_Y
END IF
! B direction
IF (word1_ln == "cpr_B") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, cpr_B 
	ECOM_CPR_glb(3) = cpr_B
END IF
! ----------------------------------------------------------------------
! ECOM2 is applied
! **********************************************************************
ELSE IF (ECOM_param_glb == 2) THEN

IF (word1_ln == "bias_D") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, bias_D
        ECOM_Bias_glb(1) = bias_D
END IF
! Y direction
IF (word1_ln == "bias_Y") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, bias_Y
        ECOM_Bias_glb(2) = bias_Y
END IF
! B direction
IF (word1_ln == "bias_B") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, bias_B
        ECOM_Bias_glb(3) = bias_B
END IF
! CPR terms
IF (word1_ln == "cpr_D") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, cpr_D
        ECOM_CPR_glb(1) = cpr_D
END IF
! Y direction
IF (word1_ln == "cpr_Y") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, cpr_Y
        cpr_D4 = cpr_Y
        ECOM_CPR_glb(2) = cpr_D4
END IF
! B direction
IF (word1_ln == "cpr_B") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, cpr_B
        ECOM_CPR_glb(3) = cpr_B
END IF

END IF

END DO
CLOSE (UNIT=UNIT_IN)
! Close of input parameterization file
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Read the estimation values of the SRP model 
! ----------------------------------------------------------------------
IF (ECOM_param_glb == 1) THEN

! Read the input file
filename = 'ECOM1_srp.in'

ELSE IF (ECOM_param_glb == 2) THEN

filename = 'ECOM2_srp.in'

END IF
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

PD_Param_ID = 0
If (ECOM_Bias_glb(1) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
End IF
If (ECOM_Bias_glb(2) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
End IF
If (ECOM_Bias_glb(3) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
End IF
If (ECOM_CPR_glb(1) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
! S term
        PD_Param_ID = PD_Param_ID + 1
End IF
If (ECOM_CPR_glb(2) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
! S term
        PD_Param_ID = PD_Param_ID + 1
End IF
If (ECOM_CPR_glb(3) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
! S term
        PD_Param_ID = PD_Param_ID + 1
End If

ALLOCATE (ECOM_accel_glb(PD_Param_ID), STAT = AllocateStatus)

NPARAM_glb = PD_Param_ID

! 1st Word of Line ith
READ (line_ith, * , IOSTAT=ios_data) word1_ln  ! 1st word

! ECOM parameters :: Bias accelerations
IF (word1_ln == "ECOM1") THEN

PD_Param_ID = 0

   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, &
   D0, Y0, B0, DC, DS, YC, YS, BC, BS  
If (ECOM_Bias_glb(1) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = D0
End IF
If (ECOM_Bias_glb(2) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = Y0
End IF
If (ECOM_Bias_glb(3) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = B0
End IF
If (ECOM_CPR_glb(1) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = DC
! S term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = DS
End IF
If (ECOM_CPR_glb(2) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = YC
! S term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = YS
End IF
If (ECOM_CPR_glb(3) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = BC
! S term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = BS
End If

ELSE IF (word1_ln == "ECOM2") THEN

PD_Param_ID = 0

   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, &
   D0, Y0, B0, D2C, D2S, D4C, D4S, BC, BS
If (ECOM_Bias_glb(1) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = D0
End IF
If (ECOM_Bias_glb(2) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = Y0
End IF
If (ECOM_Bias_glb(3) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = B0
End IF
If (ECOM_CPR_glb(1) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = D2C
! S term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = D2S
End IF
If (ECOM_CPR_glb(2) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = D4C
! S term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = D4S
End IF
If (ECOM_CPR_glb(3) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = BC
! S term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = BS
End If


END IF

END DO

CLOSE (UNIT=UNIT_IN)
! Close of input file
! ----------------------------------------------------------------------



END

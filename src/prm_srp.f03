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
! Changes:      13-12-2018 Tzupang Tseng : added ECOM-based options in the configuration file
!                                              for ACS POD
!               21-02-2019 Tzupang Tseng : added a function for switching on and
!                                          off some coefficients in ECOM models   
!               03-12-2019 Tzupang Tseng : added a function of estimating
!                                          parameters in simple box wing model
!
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
      REAL(KIND = prec_q), DIMENSION(:), ALLOCATABLE :: SRP_PARA

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
! Read SRP parameterizations
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
!
IF(ECOM_param_glb/=0 .and.ECOM_param_glb <= 2 .or. ECOM_param_glb == 12) THEN

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

   IF (word1_ln == "cpr_D2") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, cpr_D2
        ECOM_CPR_glb(4) = cpr_D2
   END IF

   IF (word1_ln == "cpr_D4") THEN
   READ ( line_ith, FMT = * , IOSTAT=ios_key ) word_i, cpr_D4
        ECOM_CPR_glb(5) = cpr_D4
   END IF

! ----------------------------------------------------------------------
! SBOXW is applied
! **********************************************************************
ELSE IF (ECOM_param_glb == 3) THEN

   ECOM_Bias_glb = 0.d0
   ECOM_CPR_glb  = 0.d0

ELSE IF(ECOM_param_glb == 0) THEN

!PRINT*,'ECOM NOT ACTIVATED'

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

ELSE IF (ECOM_param_glb == 3) THEN

filename = 'SBOXW_srp.in'

ELSE IF (ECOM_param_glb == 12) THEN

filename = 'ECOM12_srp.in'

END IF

IF (ECOM_param_glb/=0 .and.ECOM_param_glb <= 2 .or. ECOM_param_glb == 12) THEN
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
    End IF
    If (ECOM_CPR_glb(4) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
! S term
        PD_Param_ID = PD_Param_ID + 1
    End IF
    If (ECOM_CPR_glb(5) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
! S term
        PD_Param_ID = PD_Param_ID + 1
    End IF

END IF

IF (ECOM_param_glb == 3) PD_Param_ID = 7

!PRINT*,'Reading config,  PD_Param_ID =', PD_Param_ID

ALLOCATE (ECOM_accel_glb(PD_Param_ID), STAT = AllocateStatus)

IF(ECOM_param_glb==0)THEN

RETURN

ELSE

ALLOCATE (SRP_PARA(PD_Param_ID), STAT = AllocateStatus)
SRP_PARA = 0.d0

END IF

! ----------------------------------------------------------------------
! Open .in file if ECOM activated
IF (ECOM_param_glb > 0 ) THEN

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

! ECOM parameters :: Bias accelerations
    IF (word1_ln=="ECOM1".OR.word1_ln=='ECOM2'.OR.word1_ln=='ECOM12') THEN

      PD_Param_ID = 0

      READ (line_ith, FMT = * , IOSTAT=ios_key) word_i, SRP_PARA(:)  
!print*,"SRP_PARA=",SRP_PARA(:)
      If (ECOM_Bias_glb(1) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = SRP_PARA(PD_Param_ID)
!print*,'D0=',PD_Param_ID
      End IF

      If (ECOM_Bias_glb(2) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = SRP_PARA(PD_Param_ID)
!print*,'Y0=',PD_Param_ID
      End IF

      If (ECOM_Bias_glb(3) == 1) Then
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = SRP_PARA(PD_Param_ID)
!print*,'B0=',PD_Param_ID
      End IF

      If (ECOM_CPR_glb(1) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = SRP_PARA(PD_Param_ID)
!print*,'DC=',PD_Param_ID
! S term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = SRP_PARA(PD_Param_ID)
!print*,'DS=',PD_Param_ID
      End IF

      If (ECOM_CPR_glb(2) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = SRP_PARA(PD_Param_ID)
!print*,'YC=',PD_Param_ID
! S term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = SRP_PARA(PD_Param_ID)
!print*,'YS=',PD_Param_ID
      End IF

      If (ECOM_CPR_glb(3) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = SRP_PARA(PD_Param_ID)
!print*,'BC=',PD_Param_ID
! S term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = SRP_PARA(PD_Param_ID)
!print*,'BS=',PD_Param_ID
      End If

      If (ECOM_CPR_glb(4) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = SRP_PARA(PD_Param_ID)
!print*,'D2C=',PD_Param_ID
! S term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = SRP_PARA(PD_Param_ID)
!print*,'D2S=',PD_Param_ID
      End If

      If (ECOM_CPR_glb(5) == 1) THEN
! C term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = SRP_PARA(PD_Param_ID)
!print*,'D4C=',PD_Param_ID
! S term
        PD_Param_ID = PD_Param_ID + 1
        ECOM_accel_glb(PD_Param_ID) = SRP_PARA(PD_Param_ID)
!print*,'D4S=',PD_Param_ID
      End If


    ELSE IF (word1_ln == "SBOXW") THEN
      READ (line_ith, FMT = * , IOSTAT=ios_key) word_i, SRP_PARA(:)
!print*,"SRP_PARA=",SRP_PARA(:)
      DO PD_Param_ID =1, 7 
        ECOM_accel_glb(PD_Param_ID)= SRP_PARA(PD_Param_ID)
      END DO
    END IF

  END DO

CLOSE (UNIT=UNIT_IN)
! Close of input file
! ----------------------------------------------------------------------
ENDIF

END

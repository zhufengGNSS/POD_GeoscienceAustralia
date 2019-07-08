PROGRAM  main_brdcorbit 

! ----------------------------------------------------------------------
! PROGRAM: brdcorbit.f90
! ----------------------------------------------------------------------
! Purpose: Create a sp3 file from the broadcast ephemeris *.yyn file, 
!          where yy denotes the year 
! ----------------------------------------------------------------------
! Author :      Dr. Tzupang Tseng
!
! Copyright: GEOSCIENCE AUSTRALIA, AUSTRALIA
!
! Created:      31-05-2019
!
! Changes:
!
! Remarks:      08-07-2019 Tzupang Tseng: CUrrently, the program only works for
!                                         GPS-only.
! ----------------------------------------------------------------------
      USE mdl_precision
      USE mdl_num
      USE mdl_param
      USE m_write_brd2sp3
      IMPLICIT NONE

!--------------------------------------------------------------------
      INTEGER (KIND = 4) :: I, J, K, II, JJ, INEW
      INTEGER (KIND = prec_int2) :: AllocateStatus
      
      CHARACTER (LEN=100) :: INFILENAME, OUTFILENAME
      CHARACTER (LEN=100) :: NEW
      INTEGER (KIND = prec_int4) :: UNIT_IN1,UNIT_IN2,ios
      INTEGER (KIND = prec_int4) :: VERSION
      INTEGER (KIND = prec_int4) :: LEAP, NWKUTC
      INTEGER (KIND = prec_int4) :: ISAT
      INTEGER (KIND = prec_int4) :: IYEAR4,MONTH,IDAY
      REAL (KIND = prec_d),DIMENSION(4)::IONALPHA, IONBETA
      REAL (KIND = prec_d),DIMENSION(:,:,:),ALLOCATABLE ::EPH, CLK
      REAL (KIND = prec_d),DIMENSION(:,:,:),ALLOCATABLE ::EPHNEW
      REAL (KIND = prec_d),DIMENSION(:,:,:),ALLOCATABLE ::EPHNEW2
      REAL (KIND = prec_d),DIMENSION(:,:,:),ALLOCATABLE ::ECEFPOS
      REAL (KIND = prec_d),DIMENSION(:,:),ALLOCATABLE ::NEWEPOCH
      REAL (KIND = prec_d),DIMENSION(:,:),ALLOCATABLE ::NEWEPOCH2
      REAL (KIND = prec_d) :: A0UTC,A1UTC,ITUTC
      REAL (KIND = prec_d) :: CPU_t0, CPU_t1
      REAL (KIND = prec_d) :: DT, DA, DE, DI
      
      INTEGER (KIND = 4) :: MAXNSAT, MAXNPAR, MAXEPO
      INTEGER (KIND = 4) :: sz1, sz2, sz3, IREC
      INTEGER (KIND = 4) :: IG, IR, IE, IC
!------------------------------------------------------------------
      REAL (KIND = prec_d) :: DELT, DELT2, PI, DWEEK
      REAL (KIND = prec_d) :: XM021, DM0
      REAL (KIND = prec_d) :: M01, M02
      REAL (KIND = prec_q) :: MJD0, MJD
      REAL (KIND = prec_q) :: SEC00
      REAL (KIND = prec_d)       :: GPS_day
      INTEGER (KIND = prec_int8) :: GPS_week, GPSweek_mod1024

      INTEGER (KIND = 4) :: SAMPLE
      INTEGER (KIND = prec_d),DIMENSION(:),ALLOCATABLE :: IPRN
      INTEGER (KIND = prec_d),DIMENSION(:,:),ALLOCATABLE :: IBAD
! ----------------------------------------------------------------
PI = 4*ATAN(1.D0)        

! START CPU COUNTER
! -----------------
CALL cpu_time (CPU_t0)

! LOCAL VARIABLES
! ---------------
MAXNSAT = 200 ! MAXIMUM NUMBER OF SATELLITE RECORD IN BRDC ORBIT
MAXNPAR = 32  ! MAXIMUM NUMBER OF BRDC ORBI PARAMETERS
IREC = 2      ! SAMPLING RATE OF GPS BROADCAST ORBIT
MAXEPO  =MAXNSAT* 24/IREC ! MAXIMUM NUMBER OF EPOCH
SAMPLE = 900    ! SAMPLING RATE FOR SP3 OUTPUT 

ALLOCATE(CLK(MAXNPAR,MAXEPO,MAXNSAT), STAT = AllocateStatus)
ALLOCATE(EPH(MAXNPAR,MAXEPO,MAXNSAT), STAT = AllocateStatus)
ALLOCATE(EPHNEW(MAXNPAR,15,MAXNSAT), STAT = AllocateStatus)
ALLOCATE(EPHNEW2(MAXNPAR,15,MAXNSAT), STAT = AllocateStatus)
ALLOCATE(IPRN(MAXNSAT), STAT = AllocateStatus)
ALLOCATE(IBAD(MAXEPO,MAXNSAT), STAT = AllocateStatus)
ALLOCATE(ECEFPOS(86400/SAMPLE,3,MAXNSAT), STAT = AllocateStatus)
ALLOCATE(NEWEPOCH(86400/SAMPLE,MAXNSAT), STAT = AllocateStatus)
ALLOCATE(NEWEPOCH2(86400/SAMPLE,MAXNSAT), STAT = AllocateStatus)

EPHNEW = 0.d0
ECEFPOS = 0.d0
NEWEPOCH= 0.d0
CLK = 0.d0
EPH = 0.d0
IPRN= 0

INFILENAME = 'brdc3210.16n'
UNIT_IN1 = 66
! ----------------------------------------------------------------------
! Open *.yyn broadcast ephemeris file
!      OPEN (UNIT = UNIT_IN, FILE = TRIM (PRMfname), IOSTAT = ios)
OPEN (UNIT = UNIT_IN1, FILE = INFILENAME, IOSTAT = ios)
IF (ios /= 0) THEN
! PRINT *, "Error in opening file:", PRMfname
PRINT *, "Error in opening file:", INFILENAME
PRINT *, "OPEN IOSTAT=", ios
END IF

! Create an output file
OUTFILENAME = 'brdc3210.sp3'
!UNIT_IN2 = 77
!OPEN (UNIT = UNIT_IN2, FILE = OUTFILENAME)

! Read broadcast ephemeris file (GPS-only)
!-----------------------------------------
CALL readbrdcheader (UNIT_IN1, VERSION, IONALPHA, IONBETA, A0UTC, A1UTC, &
                     ITUTC, NWKUTC, LEAP)
PRINT*,'VERSION',VERSION

CALL readbrdcmessg (UNIT_IN1,VERSION,MAXNSAT,MAXNPAR,MAXEPO,IYEAR4,MONTH,IDAY,ISAT,EPH,CLK)


CLOSE (UNIT=UNIT_IN1)
PRINT*,'ALL SATELLITES HAVE BEEN RECORD !'


! Assign a new array for each satellite with the epoch recounting (K value)
! Convert the orbital elements to the ecef frame
JJ = 0
K = 0
DO  I=1, MAXNSAT
    DO  J =1,MAXEPO
! Check the epoch is effective
        IF (EPH(3,J,I) .GT. 0.d0) THEN 
            K= K+1
!print*,'PRN SAT', I
            EPHNEW(1:20,K,I) = EPH(1:20,J,I)
          IF (K .EQ. 1) THEN 
              DT = 0.d0
            DO II=1, 7
               JJ = JJ + 1
               DELT = SAMPLE*(II-1) + EPHNEW(2,K,I)
               NEWEPOCH(JJ,I) = DELT 
            END DO
    
          ELSE
                DT =  EPHNEW(2,K,I)-EPHNEW(2,K-1,I)
                CALL chkbrdc (EPHNEW(1:20,K-1,I), EPHNEW(1:20,K,I))
               IF (DT .LT. 0.d0) EXIT
               IF (DT .LT. 1800.d0) THEN
                   PRINT*,'PREVIOUS EPOCH =',EPHNEW(2,K-1,I)
                   PRINT*,'CURRENT  EPOCH =',EPHNEW(2,K,I)
                   PRINT*,'BAD RECORD IN TIME TOE',DT, 'PRN SAT =', I , &
                          'EPOCH NUMBER =', K
               IBAD(K,I) = 1
               ELSE 
               IBAD(K,I) = 0
               END IF

               
                  DO II=1, NINT(DT/SAMPLE+1)-1
                  JJ = JJ + 1
                  DELT = SAMPLE*(II-1) + EPHNEW(2,K,I)
                  NEWEPOCH(JJ,I) = DELT 
                  END DO


          END IF
              
        END IF
     
    END DO

     K=0 
     JJ = 0            
END DO




! The leap second needs to be checked!! (by using dat.for)
! -----------------------------------------
!         LEAPTST=DGPSUT(TFIRST)
!         IF (LEAP.EQ.0 .OR. LEAPTST.NE.LEAP) LEAP=LEAPTST
!         NEPO=IDINT((TLAST-TFIRST)/DTTAB*86400.D0)+1



IG = 0 ! GPS
IR = 0 ! GLONASS
IE = 0 ! GALILEO
IC = 0 ! BDS
! Count how many satellites are recorded
! --------------------------------------
DO I = 1,MAXNSAT 
   IF (I < 100 .AND. EPHNEW(3,1,I) > 0.d0) THEN
        IG = IG + 1 
        IPRN(IG) = I
!print*,'GPS PRN =', IG, IPRN(IG)
   ELSE IF (I < 200 .AND. I > 100 .AND. EPHNEW(3,1,I) > 0.d0) THEN
       IR = IR + 1
       IPRN(IR+100) = I 
!print*,'GLONASS PRN =', IR, IPRN(IR+100)
   ELSE IF (I < 300 .AND. I > 200 .AND. EPHNEW(3,1,I) > 0.d0) THEN
       IE = IE + 1
       IPRN(IE+200) = I 
!print*,'GALILEO PRN =', IE, IPRN(IE+200)
   ELSE IF (I < 400 .AND. I > 300 .AND. EPHNEW(3,1,I) > 0.d0) THEN
       IC = IC + 1
       IPRN(IC+300) = I 
!print*,'BDS PRN =', IC, IPRN(IC+300)
   END IF
END DO





JJ=0
INEW = 0
DO I = 1, IG
!print*,'SAT PRN =', I
   DO K = 1, 15
      IF (IBAD(K,I) /= 1)THEN
        INEW = INEW + 1
!print*,'NEW EPOCH =', INEW
        EPHNEW2(1:20,INEW,I) = EPHNEW(1:20,K,I)
        IF (INEW == 1) THEN
           DT = 0
        DO II=1,7 
           JJ = JJ + 1
!print*,'NEW EPOCH =', JJ, 'NEW TOE RECORD =', INEW, 'PRN SAT =', I
           DELT2 = SAMPLE*(II-1) + EPHNEW2(2,INEW,I)
           NEWEPOCH2(JJ,I) = DELT2
           CALL brdc2ecef(DELT2,EPHNEW2(1:20,INEW,I),ECEFPOS(JJ,1:3,I))
!IF (I == 5) PRINT*,'ECEFPOS =', ECEFPOS(JJ,1:3,I)
        END DO
        ELSE
           DT =  EPHNEW2(2,INEW,I)-EPHNEW2(2,INEW-1,I)
           IF (DT .LT. 0.d0) EXIT
           IF (DT .LT. 7200.d0) EPHNEW2(2,INEW,I) = 2*3600.d0 + EPHNEW2(2,INEW-1,I) 
           IF (DT .GT. 7200.d0) EPHNEW2(2,INEW,I) = 2*3600.d0 + EPHNEW2(2,INEW-1,I)
!CALL chkbrdc (EPHNEW2(1:20,INEW-1,I), EPHNEW2(1:20,INEW,I))
!print*,'NEW TOE RECORD =', INEW, 'PRN SAT =', I
       ! END IF
        DO II = 1,NINT(DT/SAMPLE) 
           JJ = JJ + 1
!print*,'NEW EPOCH =', JJ, 'NEW TOE RECORD=', INEW, 'PRN SAT =', I
!CALL chkbrdc (EPHNEW2(1:20,INEW-1,I), EPHNEW2(1:20,INEW,I))
           DELT2 = SAMPLE*(II-1) + EPHNEW2(2,INEW,I)
           NEWEPOCH2(JJ,I) = DELT2
!print*,'NEWEPOCH2(JJ,I)=',DELT2,'JJ=',JJ,'EPHNEW2(2,INEW,I)=',EPHNEW2(2,INEW,I)
           CALL brdc2ecef(DELT2,EPHNEW2(1:20,INEW,I),ECEFPOS(JJ,1:3,I))
!IF (I == 5) PRINT*,'ECEFPOS =', ECEFPOS(JJ,1:3,I)
        END DO
        END IF
      END IF
   END DO
INEW = 0
JJ = 0
END DO


! WRITE THE ECEF POSITIONS IN A SP3 FORMAT
! ----------------------------------------
CALL write_brd2sp3 (IG,SAMPLE,IYEAR4,MONTH,IDAY,NEWEPOCH2, ECEFPOS, OUTFILENAME, 0)


! SATELLITE CLOCK CORRRECTION
! ----------------------------------
!DT=TEPO-GPSMJD(CLKNEW(L0+11,JSAT),IDNINT(CLKNEW(L0+1,JSAT)))
!DT=DT*86400.D0
!DTSATC(ISAT)=CLKNEW(L0+14,JSAT)+CLKNEW(L0+13,JSAT)*DT+CLKNEW(L0+12,JSAT)*DT**2


CALL cpu_time (CPU_t1)

PRINT *,"CPU Time (sec)", CPU_t1-CPU_t0

END PROGRAM 

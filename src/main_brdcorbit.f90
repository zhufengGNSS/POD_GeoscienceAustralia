PROGRAM  main_brdcorbit 

! ----------------------------------------------------------------------
! PROGRAM: brdcorbit.f90
! ----------------------------------------------------------------------
! Purpose: Create a sp3 file from the multi-GNSS broadcast ephemeris *.yyn file, 
!          where yy denotes the year or *.rnx file
! ----------------------------------------------------------------------
! Author :      Dr. Tzupang Tseng
!
! Copyright: GEOSCIENCE AUSTRALIA, AUSTRALIA
!
! Created:   31-05-2019
!
! Changes:   23-07-2019 Tzupang Tseng: The multi-GNSS option has been implemented 
!                                      in the program, except for GLONASS   
!
! Remarks:   08-07-2019 Tzupang Tseng: CUrrently, the program only works for
!                                      GPS-only.
! ----------------------------------------------------------------------
      USE mdl_precision
      USE mdl_num
      USE mdl_param
      USE mdl_brdconfig
      USE m_write_brd2sp3
      IMPLICIT NONE

!--------------------------------------------------------------------
      INTEGER (KIND = 4) :: I, J, K, II, JJ, INEW
      INTEGER (KIND = prec_int2) :: AllocateStatus, DeAllocateStatus
       
      !CHARACTER (LEN=100) :: INFILENAME, OUTFILENAME
      CHARACTER (LEN=100) :: NEW,INPUT,OUTPUT
      CHARACTER (LEN=3)   :: SATTYPE
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
      REAL (KIND = prec_d) :: DTGAL,XT, KKK
      
      INTEGER (KIND = 4) :: MAXNSAT, MAXNPAR, MAXEPO
      INTEGER (KIND = 4) :: sz1, sz2, sz3, IREC, NDT
      INTEGER (KIND = 4) :: IG, IR, IE, IC, IJ, TOTG
      INTEGER (KIND = 4) :: IGAL, IBDS, IQZSS
      INTEGER (KIND = 4) :: MAXNGPS, MAXNGLN, MAXNGAL, MAXNBDS, MAXNQZSS
!------------------------------------------------------------------
      REAL (KIND = prec_d) :: DELT, DELT2, PI, DWEEK
      REAL (KIND = prec_d) :: XM021, DM0
      REAL (KIND = prec_d) :: M01, M02
      REAL (KIND = prec_q) :: MJD0, MJD
      REAL (KIND = prec_q) :: SEC00
      REAL (KIND = prec_q) :: MJD_ti
      REAL (KIND = prec_d)       :: GPS_day, GPS_wsec
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
MAXNGPS = 50
MAXNGLN = 30
MAXNGAL = 50
MAXNBDS = 50
MAXNQZSS = 10

MAXNSAT = 500 ! MAXIMUM NUMBER OF SATELLITE RECORD IN BRDC ORBIT
MAXNPAR = 32  ! MAXIMUM NUMBER OF BRDC ORBI PARAMETERS
IREC = 2      ! SAMPLING RATE OF GPS BROADCAST ORBIT
MAXEPO = MAXNSAT* 24/IREC ! MAXIMUM NUMBER OF EPOCH
SAMPLE = 900    ! SAMPLING RATE FOR SP3 OUTPUT 

ALLOCATE(CLK(MAXNPAR,MAXEPO,MAXNSAT), STAT = AllocateStatus)
ALLOCATE(EPH(MAXNPAR,MAXEPO,MAXNSAT), STAT = AllocateStatus)
ALLOCATE(EPHNEW(MAXNPAR,MAXEPO,MAXNSAT), STAT = AllocateStatus)
ALLOCATE(EPHNEW2(MAXNPAR,MAXEPO,MAXNSAT), STAT = AllocateStatus)
ALLOCATE(IPRN(MAXNSAT), STAT = AllocateStatus)
ALLOCATE(IBAD(MAXEPO,MAXNSAT), STAT = AllocateStatus)
ALLOCATE(ECEFPOS(86400/SAMPLE,3,MAXNSAT), STAT = AllocateStatus)
ALLOCATE(NEWEPOCH(86400/SAMPLE,MAXNSAT), STAT = AllocateStatus)
ALLOCATE(NEWEPOCH2(86400/SAMPLE,MAXNSAT), STAT = AllocateStatus)

EPHNEW = 0.d0
ECEFPOS = 0.d0
NEWEPOCH= 0.d0
!CLK = 0.d0
!EPH = 0.d0
IPRN= 0

CALL brdc_cmdline
INPUT = TRIM(INFILENAME)
UNIT_IN1 = 66
! ----------------------------------------------------------------------
! Open *.yyn broadcast ephemeris file
OPEN (UNIT = UNIT_IN1, FILE = INPUT, IOSTAT = ios)
IF (ios /= 0) THEN
! PRINT *, "Error in opening file:", PRMfname
PRINT *, "Error in opening file:", INFILENAME
PRINT *, "OPEN IOSTAT=", ios
END IF

! Create an output file
OUTPUT = TRIM(OUTFILENAME)

! Select the GNSS constellation to be processed
SATTYPE = TRIM(GNSSTYPE)

! Read broadcast ephemeris file 
!-----------------------------------------
CALL readbrdcheader (UNIT_IN1, VERSION, IONALPHA, IONBETA, A0UTC, A1UTC, &
                     ITUTC, NWKUTC, LEAP)
!PRINT*,'RINEX VERSION =',VERSION

CALL readbrdcmessg (UNIT_IN1,VERSION,MAXNSAT,MAXNPAR,MAXEPO,IYEAR4,MONTH,IDAY,ISAT,EPH,CLK)

CLOSE (UNIT=UNIT_IN1)
PRINT*,'ALL SATELLITES HAVE BEEN RECORD !'
PRINT*,


! Create the reference epoch of the broadcast file
! ------------------------------------------------
CALL iau_CAL2JD ( IYEAR4, MONTH, IDAY, MJD0, MJD_ti, J )

! The seconds of the GPS week
! ----------------------------
CALL time_GPSweek (MJD_ti , GPS_week, GPS_wsec, GPSweek_mod1024)



! Assign a new array for each satellite with the epoch recounting (K value)
! and Detect the bad records in TOE by using an index of IBAD
! -------------------------------------------------------------------------

! Check all records (Here the daily brdc file is used, so the first epoch should
! be the UTC 0h.)
! ------------------------------------------------------------------------------
JJ = 0
K = 0
DO I=1, MAXNSAT

IF (I < 100 .OR. I < 300 .AND. I > 200 .OR. &
    I < 400 .AND. I > 300 .OR. I > 400) THEN

    DO  J =1,MAXEPO
        IF (EPH(3,J,I) .GT. 0.d0) THEN 
            K= K+1
            EPHNEW(1:20,K,I) = EPH(1:20,J,I)
            IF (K .EQ. 1) THEN 
               DT = 0.d0
               DO II=1, 7
               JJ = JJ + 1
               DELT = SAMPLE*(II-1) + EPHNEW(2,K,I)
               NEWEPOCH(JJ,I) = DELT 
!print*,'PRN SAT', I, 'DELT =', DELT,  'JJ =', JJ
               END DO
            ELSE
               DT =  EPHNEW(2,K,I)-EPHNEW(2,K-1,I)
               IF (DT .LT. 0.d0) EXIT
               ! For GPS and BDS case
               IF (I < 100 .OR.  I < 400 .AND. I > 300) THEN
                  IF (DT .LT. 1800.d0) THEN
!                  PRINT*,'PREVIOUS EPOCH =',EPHNEW(2,K-1,I)
!                  PRINT*,'CURRENT  EPOCH =',EPHNEW(2,K,I)
!                  PRINT*,'BAD RECORD IN TIME TOE',DT, 'PRN SAT =', I , &
!                         'EPOCH NUMBER =', K
                  IBAD(K,I) = 1
                  ELSE 
                  IBAD(K,I) = 0
                  END IF
               ! For Galileo case
               ELSE IF (I < 300 .AND. I > 200 ) THEN
                  IF (DT .LT. 600.d0) THEN
!                  PRINT*,'PREVIOUS EPOCH =',EPHNEW(2,K-1,I)
!                  PRINT*,'CURRENT  EPOCH =',EPHNEW(2,K,I)
!                  PRINT*,'BAD RECORD IN TIME TOE',DT, 'PRN SAT =', I , &
!                         'EPOCH NUMBER =', K
                  IBAD(K,I) = 1
                  ELSE
                  IBAD(K,I) = 0
                  END IF
              ! For QZSS case
               ELSE IF (I > 400 ) THEN
                  IF (DT .LT. 800.d0) THEN
!                  PRINT*,'PREVIOUS EPOCH =',EPHNEW(2,K-1,I)
!                  PRINT*,'CURRENT  EPOCH =',EPHNEW(2,K,I)
!                  PRINT*,'BAD RECORD IN TIME TOE',DT, 'PRN SAT =', I , &
!                         'EPOCH NUMBER =', K
                  IBAD(K,I) = 1
                  ELSE
                  IBAD(K,I) = 0
                  END IF
               END IF
               
               DO II=1, NINT(DT/SAMPLE+1)-1
               JJ = JJ + 1
               DELT = SAMPLE*(II-1) + EPHNEW(2,K,I)
               NEWEPOCH(JJ,I) = DELT 
!print*,'PRN SAT', I, 'DELT =', DELT,  'JJ =', JJ
               END DO
            END IF
              
        END IF
     
    END DO

END IF
K=0 
JJ = 0            
END DO



! Check leap seconds ? (by using dat.for) (To be processed)
! -----------------------------------------
!         LEAPTST=DGPSUT(TFIRST)
!         IF (LEAP.EQ.0 .OR. LEAPTST.NE.LEAP) LEAP=LEAPTST
!         NEPO=IDINT((TLAST-TFIRST)/DTTAB*86400.D0)+1



IG = 0 ! GPS
IR = 0 ! GLONASS
IE = 0 ! GALILEO
IC = 0 ! BDS
IJ = 0 ! QZSS
! Count how many satellites are recorded
! --------------------------------------
DO I = 1,MAXNSAT 
   IF (I <= 100 .AND. EPHNEW(3,1,I) > 0.d0) THEN
        IG = IG + 1 
        IPRN(IG) = I
!print*,'GPS PRN =', IG, IPRN(IG)
   ELSE IF (I <= 200 .AND. I > 100 .AND. EPHNEW(6,1,I) > 0.d0) THEN
       IR = IR + 1
       IPRN(IR+100) = I 
!print*,'GLONASS PRN =', IR, IPRN(IR+100)
   ELSE IF (I <= 300 .AND. I > 200 .AND. EPHNEW(3,1,I) > 0.d0) THEN
       IE = IE + 1
       IPRN(IE+200) = I 
!print*,'GALILEO PRN =', IE, IPRN(IE+200)
   ELSE IF (I <= 400 .AND. I > 300 .AND. EPHNEW(3,1,I) > 0.d0) THEN
       IC = IC + 1
       IPRN(IC+300) = I 
!print*,'BDS PRN =', IC, IPRN(IC+300)
   ELSE IF (I <= 500 .AND. I > 400 .AND. EPHNEW(3,1,I) > 0.d0) THEN
       IJ = IJ + 1
       IPRN(IJ+400) = I
!print*,'OZSS PRN =', IJ, IPRN(IJ+400)
   END IF
END DO


! Process the desire GNSS constellation type
! ------------------------------------------
IF (SATTYPE == 'G')THEN 
TOTG = IG
DO I =1, MAXNSAT
   IF (I > 100) THEN
       EPHNEW(:,:,I) = 0.d0
       IPRN(I) = 0
   END IF
END DO

ELSE IF (SATTYPE == 'R') THEN 
TOTG = IR
DO I =1, MAXNSAT
   IF (I <= 100 .OR. I > 200) THEN
       EPHNEW(:,:,I) = 0.d0
       IPRN(I) = 0
   END IF
END DO

ELSE IF (SATTYPE == 'E') THEN 
TOTG = IE
DO I =1, MAXNSAT
   IF (I <= 200 .OR. I > 300) THEN
       EPHNEW(:,:,I) = 0.d0
       IPRN(I) = 0
   END IF
END DO

ELSE IF (SATTYPE == 'C') THEN 
TOTG = IC
DO I =1, MAXNSAT
   IF (I <= 300 .OR. I > 400) THEN
       EPHNEW(:,:,I) = 0.d0
       IPRN(I) = 0
   END IF
END DO

ELSE IF (SATTYPE == 'J') THEN 
TOTG = IJ
DO I =1, MAXNSAT
   IF (I <= 400 ) THEN 
       EPHNEW(:,:,I) = 0.d0
       IPRN(I) = 0
   END IF
END DO

ELSE IF (SATTYPE == 'A') THEN
TOTG = IG + IE + IC + IJ
END IF

! Remove the bad records and convert the broadcast elements to ECEF coordinate
! system 
! -----------------------------------------------------------------------------
JJ=0

DO I = 1, MAXNSAT
INEW = 0
IGAL = 0
IBDS = 0
IQZSS = 0
KKK = 0
IF (I < 100 .AND. EPHNEW(3,1,I) > 0.d0 ) THEN ! GNSS/GPS sampling rate: 2 hour
   DO K = 1, 15
         IF (IBAD(K,I) /= 1 )THEN
            INEW = INEW + 1
            EPHNEW2(1:20,INEW,I) = EPHNEW(1:20,K,I)
            IF (INEW == 1 .AND. EPHNEW2(2,1,I) >= GPS_wsec) THEN
            DO II=1,8 ! 15-minute resolution 
               JJ = JJ + 1
               DELT2 = SAMPLE*(II-1) + GPS_wsec
               NEWEPOCH2(JJ,I) = DELT2
               CALL brdc2ecef(DELT2,EPHNEW2(1:20,INEW,I),ECEFPOS(JJ,1:3,I))
!print*,'NEW EPOCH =', JJ, 'PRN SAT =', I,'TOE=',EPHNEW2(2,INEW,I), &
!       'DELT2 =', DELT2,'POSITIONS =',ECEFPOS(JJ,1:3,I)

            END DO
           
            ELSE IF (EPHNEW2(2,INEW,I) == GPS_wsec+7200*(INEW-1)) THEN
                 DT = EPHNEW2(2,INEW,I) - GPS_wsec 
                 IF (DT .LT. 0.d0) EXIT
                 IF (NINT(DT) - NINT(KKK) > 7200.d0) THEN
!                 PRINT*,'ADD POINTS FOR INTERPOLATION'
                 NDT =  (NINT(DT) - NINT(KKK))/SAMPLE
                 DO II=1,NDT-8
                 JJ = JJ + 1
                 DELT2 = SAMPLE*(II-1) + NINT(KKK) + GPS_wsec + 7200 
                 NEWEPOCH2(JJ,I) = DELT2
                 CALL brdc2ecef(DELT2,EPHNEW2(1:20,INEW,I),ECEFPOS(JJ,1:3,I))
!print*,'NEW EPOCH =', JJ, 'PRN SAT =', I,'TOE=',EPHNEW2(2,INEW,I), &
!       'DELT2 =', DELT2,'POSITIONS =',ECEFPOS(JJ,1:3,I)

                 END DO
                 END IF
                 KKK = DT
 
            DO II = 1, 8
               JJ = JJ + 1
               DELT2 = SAMPLE*(II-1) + NINT(DT)+ GPS_wsec
               IF (DELT2 > EPHNEW2(2,1,I)+86400-900 ) EXIT
               NEWEPOCH2(JJ,I) = DELT2
               CALL brdc2ecef(DELT2,EPHNEW2(1:20,INEW,I),ECEFPOS(JJ,1:3,I))
!print*,'NEW EPOCH =', JJ, 'PRN SAT =', I,'TOE=',EPHNEW2(2,INEW,I), &
!       'DELT2 =', DELT2,'POSITIONS =',ECEFPOS(JJ,1:3,I)

            END DO
         
            END IF
         END IF
   END DO

ELSE IF (I < 300 .AND. I >= 200 .AND. EPHNEW(3,1,I) > 0.d0) THEN ! GALILEO sampling rate: 10 minutes

   DO K = 1, 150
         IF (IBAD(K,I) /= 1 )THEN
            IGAL = IGAL +1
            EPHNEW2(1:20,IGAL,I) = EPHNEW(1:20,K,I)
        
            IF (IGAL == 1 .AND. EPHNEW2(2,1,I) >= GPS_wsec) THEN
               DO II=1,8
                  JJ = JJ + 1
                  DELT2 = SAMPLE*(II-1) + GPS_wsec
                  NEWEPOCH2(JJ,I) = DELT2
                  CALL brdc2ecef(DELT2,EPHNEW2(1:20,IGAL,I),ECEFPOS(JJ,1:3,I))
!print*,'NEW EPOCH =', JJ, 'PRN SAT =', I,'TOE=',EPHNEW2(2,IGAL,I), &
!       'DELT2 =', DELT2,'POSITIONS =',ECEFPOS(JJ,1:3,I)

               END DO
            ELSE IF (IGAL > 1 .AND. EPHNEW2(2,IGAL,I)-EPHNEW2(2,IGAL-1,I) == 600.d0 &
                     .AND.  MOD(EPHNEW2(2,IGAL,I),7200.d0) == 0) THEN
               DT = EPHNEW2(2,IGAL,I) - GPS_wsec 
!PRINT*,'DT =', NINT(DT), NINT(KKK), IGAL
               IF (NINT(DT)-NINT(KKK) > 7200.d0) THEN
!               print*,'ADD POINTS FOR INTERPOLATION' 
               NDT =  (NINT(DT) - NINT(KKK))/SAMPLE
                       DO II = 1, NDT-8
                          JJ = JJ + 1
                       DELT2 = SAMPLE*(II-1) + NINT(KKK) + GPS_wsec + 7200  
                       NEWEPOCH2(JJ,I) = DELT2
                       CALL brdc2ecef(DELT2,EPHNEW2(1:20,IGAL,I),ECEFPOS(JJ,1:3,I))
!print*,'NEW EPOCH =', JJ, 'PRN SAT =', I,'TOE=',EPHNEW2(2,IGAL,I), &
!       'DELT2 =', DELT2,'POSITIONS =',ECEFPOS(JJ,1:3,I)

                       END DO
               END IF
               KKK = DT
               DO II = 1, 8 
                  JJ = JJ + 1
                  DELT2 = SAMPLE*(II-1) +  NINT(DT) + GPS_wsec 
                  NEWEPOCH2(JJ,I) = DELT2
                  CALL brdc2ecef(DELT2,EPHNEW2(1:20,IGAL,I),ECEFPOS(JJ,1:3,I))
!print*,'NEW EPOCH =', JJ, 'PRN SAT =', I,'TOE=',EPHNEW2(2,IGAL,I), &
!       'DELT2 =', DELT2,'POSITIONS =',ECEFPOS(JJ,1:3,I)

               END DO
            ELSE IF (EPHNEW2(2,IGAL,I) >= 79200.d0 + GPS_wsec ) THEN
               IF (JJ < 96) THEN       
               DO II = 1, 8
                  JJ = JJ + 1
                  DELT2 = SAMPLE*(II-1) +  79200.d0 + GPS_wsec
                  NEWEPOCH2(JJ,I) = DELT2
                  CALL brdc2ecef(DELT2,EPHNEW2(1:20,IGAL,I),ECEFPOS(JJ,1:3,I))
!print*,'NEW EPOCH =', JJ, 'PRN SAT =', I,'TOE=',EPHNEW2(2,IGAL,I), &
!       'DELT2 =', DELT2,'POSITIONS =',ECEFPOS(JJ,1:3,I)
                
               END DO
               ELSE
               EXIT 
               END IF

            END IF 
         END IF
   END DO

ELSE IF (I < 400 .AND. I >= 300 .AND. EPHNEW(3,1,I) > 0.d0) THEN ! BDS sampling rate: 1 hour

   DO K = 1, 24
         IF (IBAD(K,I) /= 1 )THEN
            IBDS = IBDS +1
            EPHNEW2(1:20,IBDS,I) = EPHNEW(1:20,K,I)
            IF (IBDS == 1 .AND. EPHNEW2(2,1,I) >= GPS_wsec) THEN
               DO II=1,8
                  JJ = JJ + 1
                  DELT2 = SAMPLE*(II-1) + GPS_wsec
                  NEWEPOCH2(JJ,I) = DELT2
                  CALL brdc2ecef(DELT2,EPHNEW2(1:20,IBDS,I),ECEFPOS(JJ,1:3,I))
!print*,'NEW EPOCH =', JJ, 'PRN SAT =', I,'TOE=',EPHNEW2(2,IBDS,I), &
!       'DELT2 =', DELT2, 'POSITIONS =',ECEFPOS(JJ,1:3,I)

               END DO
            ELSE IF (IBDS > 1 .AND. EPHNEW2(2,IBDS,I)-EPHNEW2(2,IBDS-1,I) == 3600.d0 &
                     .AND.  MOD(EPHNEW2(2,IBDS,I),7200.d0) == 0) THEN
               DT = EPHNEW2(2,IBDS,I) - GPS_wsec
!PRINT*,'DT =', NINT(DT), NINT(KKK)
               IF(NINT(DT)-NINT(KKK) > 7200.d0) THEN
!               print*,'ADD POINTS FOR INTERPOLATION'
               NDT =  (NINT(DT) - NINT(KKK))/SAMPLE
                       DO II = 1, NDT-8
                          JJ = JJ + 1
                       DELT2 = SAMPLE*(II-1) + NINT(KKK) + GPS_wsec + 7200
                       NEWEPOCH2(JJ,I) = DELT2
                       CALL brdc2ecef(DELT2,EPHNEW2(1:20,IBDS,I),ECEFPOS(JJ,1:3,I))
!print*,'NEW EPOCH =', JJ, 'PRN SAT =', I,'TOE=',EPHNEW2(2,IBDS,I), &
!       'DELT2 =', DELT2,'POSITIONS =',ECEFPOS(JJ,1:3,I)
                       END DO

               END IF
               KKK = DT
               DO II = 1, 8
                  JJ = JJ + 1
                  DELT2 = SAMPLE*(II-1) + NINT(DT) + GPS_wsec 
                  NEWEPOCH2(JJ,I) = DELT2
                  CALL brdc2ecef(DELT2,EPHNEW2(1:20,IBDS,I),ECEFPOS(JJ,1:3,I))
!print*,'NEW EPOCH =', JJ, 'PRN SAT =', I,'TOE=',EPHNEW2(2,IBDS,I), &
!       'DELT2 =', DELT2, 'POSITIONS =',ECEFPOS(JJ,1:3,I)

               END DO
            END IF
         END IF
   END DO

ELSE IF (I >= 400 .AND. EPHNEW(3,1,I) > 0.d0) THEN ! QZSS sampling rate: 15 minutes

   DO K = 1, 100
         IF (IBAD(K,I) /= 1 )THEN
            IQZSS = IQZSS +1
            EPHNEW2(1:20,IQZSS,I) = EPHNEW(1:20,K,I)

            IF (IQZSS == 1 .AND. EPHNEW2(2,1,I) >= GPS_wsec)  THEN
               DO II=1,8
                  JJ = JJ + 1
                  DELT2 = SAMPLE*(II-1) + GPS_wsec 
                  NEWEPOCH2(JJ,I) = DELT2
                  CALL brdc2ecef(DELT2,EPHNEW2(1:20,IQZSS,I),ECEFPOS(JJ,1:3,I))
!print*,'NEW EPOCH =', JJ, 'PRN SAT =', I,'TOE=',EPHNEW2(2,IQZSS,I), &
!       'DELT2 =', DELT2,'POSITIONS =',ECEFPOS(JJ,1:3,I)

               END DO
             
            ELSE IF (IQZSS > 1 .AND. EPHNEW2(2,IQZSS,I)-EPHNEW2(2,IQZSS-1,I) >=500.d0 &
                     .AND.  MOD(EPHNEW2(2,IQZSS,I),7200.d0) == 0) THEN
               DT = EPHNEW2(2,IQZSS,I) - GPS_wsec
!PRINT*,'DT =', NINT(DT), NINT(KKK)
               IF(NINT(DT)-NINT(KKK) > 7200.d0) THEN
!               print*,'ADD POINTS FOR INTERPOLATION'
               NDT = (NINT(DT) - NINT(KKK))/SAMPLE
                       DO II = 1, NDT-8
                          JJ = JJ + 1
                       DELT2 = SAMPLE*(II-1) +  NINT(KKK) + GPS_wsec + 7200
                       NEWEPOCH2(JJ,I) = DELT2
                       CALL brdc2ecef(DELT2,EPHNEW2(1:20,IQZSS,I),ECEFPOS(JJ,1:3,I))
!print*,'NEW EPOCH =', JJ, 'PRN SAT =', I,'TOE=',EPHNEW2(2,IQZSS,I), &
!       'DELT2 =', DELT2,'POSITIONS =',ECEFPOS(JJ,1:3,I)

                       END DO
               END IF
               KKK = DT
               
               DO II = 1, 8 
                  JJ = JJ + 1
                  DELT2 = SAMPLE*(II-1) +  NINT(DT) + GPS_wsec 
                  NEWEPOCH2(JJ,I) = DELT2
                  CALL brdc2ecef(DELT2,EPHNEW2(1:20,IQZSS,I),ECEFPOS(JJ,1:3,I))
!print*,'NEW EPOCH =', JJ, 'PRN SAT =', I,'TOE=',EPHNEW2(2,IQZSS,I), &
!       'DELT2 =', DELT2, 'POSITIONS =',ECEFPOS(JJ,1:3,I)

               END DO
            END IF
         END IF
   END DO

END IF ! GNSS/GPS
  
INEW = 0
JJ = 0
END DO


! WRITE THE ECEF POSITIONS IN A SP3 FORMAT
! ----------------------------------------
CALL write_brd2sp3 (IPRN,TOTG,IG,IR,IE,IC,IJ,SATTYPE, SAMPLE,IYEAR4,MONTH,IDAY,NEWEPOCH2, ECEFPOS, OUTPUT, 0)


! SATELLITE CLOCK CORRRECTION (TO BE PROCESSED)
! ----------------------------------

DEALLOCATE(EPHNEW,  STAT = DeAllocateStatus)
DEALLOCATE(EPHNEW2, STAT = DeAllocateStatus)
DEALLOCATE(CLK, STAT = DeAllocateStatus)
DEALLOCATE(IPRN, STAT = DeAllocateStatus)
DEALLOCATE(IBAD, STAT = DeAllocateStatus)
DEALLOCATE(ECEFPOS,  STAT = DeAllocateStatus)
DEALLOCATE(NEWEPOCH, STAT = DeAllocateStatus)
DEALLOCATE(NEWEPOCH2, STAT = DeAllocateStatus)
DEALLOCATE(EPH, STAT = DeAllocateStatus)

CALL cpu_time (CPU_t1)
PRINT *,"CPU Time (sec)", CPU_t1-CPU_t0


END PROGRAM 

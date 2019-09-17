C*
      SUBROUTINE SRPFBOXW(ERM,ANT,GRD,REFF,YSAT,SUN,KAPPA,MONTH,
     1                    BLKNUM,SVN,MJD,ACCEL)
CC
CC NAME       :  SRPFBOXW
CC
CC PURPOSE    :  COMPUTATION OF EARTH RADIATION PRESSURE ACTING ON A
CC               BOW-WING SATELLITE
CC
CC PARAMETERS :
CC         IN :  ERM       : EARTH RADIATION MODEL
CC                           0 = NONE
CC                           1 = EARTH RADIATION PRESSURE (ANALYTICAL)
CC                           2 = EARTH RADIATION PRESSURE (CERES DATA)
CC               CERES DAT : CERES DATA (EXTERNAL ASCII FILES) -> SET PATH IN DATPATH
CC               ANT       : 0 = NO ANTENNA THRUST
CC                         : 1 = WITH ANTENNA THRUST
CC               GRD       : 1 = 2.5 degree grid
CC                         : 2 = 5.0 degree grid
CC                         : 3 = 10  degree grid
CC               REFF      : ACCELERATION IN REFERENCE FRAME:
CC                           0 = INERTIAL
CC                           1 = BODY FIXED (Z,Y,X)
CC                           2 = SUN FIXED (D,Y,B)
CC                           3 = ORBITAL (RADIAL, ALONG- AND CROSS-TRACK)
CC               YSAT      : SATELLITE POSITION [m] AND VELOCITY [m/s] (INERTIAL), 
CC                           (POSITION,VELOCITY) = (RX,RY,RZ,VX,VY,VZ)
CC               SUN       : SUN POSITION VECTOR [m] (INERTIAL)
CC               KAPPA     : ROTATION MATRIX FROM INERTIAL TO EARTH-FIXED
CC                           REFERENCE FRAME (FOR ERM=2)
CC               MONTH     : MONTH NUMBER (1 ... 12), FOR ERM=2
CC               BLKNUM    : BLOCK NUMBER
CC                             1 = GPS-I
CC                             2 = GPS-II
CC                             3 = GPS-IIA
CC                             4 = GPS-IIR
CC                             5 = GPS-IIR-A
CC                             6 = GPS-IIR-B
CC                             7 = GPS-IIR-M
CC                             8 = GPS-IIF
CC                             9 = GPS-IIIA  (Updated from acc_albedo_propboxw.f)
CC                           101 = GLONASS
CC                           102 = GLONASS-M (Added TAH 190702)
CC                           103 = GLONASS-K (Added TAH 190702)
CC                           201 = Galileo (IOV) (Added from acc_albedo_propboxw.f)
CC                           202 = Galileo (FOC) (Added from acc_albedo_propboxw.f)
CC                           NB: when adding a new block number, check last block in this function!
CC               SVN       : SPACE VEHICLE NUMBER          
CC               MJD       : MODIFIED JULIAN DAY
CC
CC        OUT : ACCEL      : ACCELERATION VECTOR [m/s^2]
CC
CC AUTHOR     : C.J. RODRIGUEZ-SOLANO
CC              rodriguez@bv.tum.de
CC
CC VERSION    : 1.0 (OCT 2010)
CC
CC CREATED    : 2010/10/18
CC
C*
C
      IMPLICIT NONE

C MOD TAH 190722: Introduced max_blk for number of block types
C     that can be stored (was 15 for GPS+GLONASS).
      integer*4 max_blk   ! Max number of block types indices
                          ! 0-10 GPS, 11-20 Glonass, 21-30 Galileo
      parameter ( max_blk = 30 ) 
C
      INTEGER*4 BLKNUM,ERM,ANT,GRD,REFF,SBLK,INDB,MONTH,SVN
      INTEGER*4 IFIRST,II,JJ,K,LFNLOC,IOSTAT,LENPATH
      INTEGER*4 LATK,LONK,LATKMX,LONKMX,GRDCER
      integer*4 i,j
C
      REAL*8 PI,C,AU,S0,TOA,ALB,MASS
      REAL*8 YSAT(6),SUN(3),FORCE(3),ACCEL(3)
C
C     EARTH AND SATELLITE PROPERTIES
      REAL*8 CERES_R(72,144),CERES_E(72,144)
      REAL*8 CERGRE(72,144),CERGEM(72,144)
      REAL*8 D_AREA_ALL(72,144),V_NS_ALL(72,144,3)
      REAL*8 AREA(4,2,max_blk),REFL(4,2,max_blk),DIFU(4,2,max_blk),
     .       ABSP(4,2,max_blk)
      REAL*8 AREA2(4,2,max_blk),REFL2(4,2,max_blk),DIFU2(4,2,max_blk),
     .       ABSP2(4,2,max_blk)
      REAL*8 REFLIR(4,2,max_blk),DIFUIR(4,2,max_blk),ABSPIR(4,2,max_blk)
      REAL*8 AREAS(4,2),REFLS(4,2),DIFUS(4,2),ABSPS(4,2)
      REAL*8 AREA2S(4,2),REFL2S(4,2),DIFU2S(4,2),ABSP2S(4,2)
      REAL*8 REFLIRS(4,2),DIFUIRS(4,2),ABSPIRS(4,2)
C
C     ATTITUDE
      REAL*8 RADVEC(3),ALGVEC(3),CRSVEC(3),ESUN(3)
      REAL*8 RSUN,ABSPOS,ABSCRS,ABSY0V,Y0SAT(3)
      REAL*8 Z_SAT(3),D_SUN(3),Y_SAT(3),B_SUN(3),X_SAT(3)
      REAL*8 ATTSURF(3,4)
C
      REAL*8 ABSNCFVI,ABSNCFIR,ALBFAC,PHASEVI,PHASEIR
      REAL*8 NCFVEC(3)
      REAL*8 D_ANG,GRDANG,GRDNUM,LATIND,LONIND
      REAL*8 D_AREA,LAT_IN,LON_IN,DIST2,ABSDIST
      REAL*8 COSLAT,PSI
      REAL*8 V_NS(3),V_INS(3),V_DIST(3),V_SAT(3)
      REAL*8 COS_IN,COS_RE,PSIDOT,ABSSUN
      REAL*8 REFL_CF,EMIT_CF,E_REFL,E_EMIT
      REAL*8 FORCE_LL(3),FREF(3)
      REAL*8 ANTFORCE,ANTPOW,MJD
      REAL*8 KAPPA(3,3)
C
      CHARACTER*2 F_MONTH(12)
      CHARACTER*100 DATPATH,FILREFL,FILEMIT
      CHARACTER*5000 LINE
      CHARACTER*25 ITEM

C MOD TAH 190722: introduce "known_blks" arrway with know block
C     types for GPS, Glonass and Galileo.  Needs to be extended
C     when other constellations added.
      integer*4 known_blks(max_blk)  ! Block numbers by the type
                  ! we known and have coded.  0 value unknown

      data known_blks /  1,   2,   3,   4,   5,  6,  7,  8,  9,  0,
     .                 101, 102, 103,   0,   0,  0,  0,  0,  0,  0,
     .                 201, 202,   0,   0,   0,  0,  0,  0,  0,  0 /
C
      DATA IFIRST/1/
      DATA (F_MONTH(K),K=1,12)/
     . '01', '02', '03', '04', '05', '06', 
     . '07', '08', '09', '10', '11', '12'/

C DATA FILES PATH AND NAME (CHANGE TO LOCAL PATH)
      DATPATH = './CERES/'
            
C COMPLETE FILE NAMES
      II = LEN(DATPATH)
      DO WHILE (DATPATH(II:II).EQ.' ')
         II = II-1
      ENDDO
      LENPATH = II

      FILREFL = DATPATH(1:LENPATH) // 'REFLMO' // F_MONTH(MONTH)
      FILEMIT = DATPATH(1:LENPATH) // 'EMITMO' // F_MONTH(MONTH)

C Constants needed
      PI = 4D0*DATAN(1D0)
      C  = 299792458D0
      AU = 149597870691D0
C Solar Constant(W/m2)
      S0 = 1367D0
C Top of Atmosphere for CERES
      TOA = 6371000D0 + 30000D0
C Albedo of the Earth
      ALB = 0.3D0
C Initialization of force vector
      FORCE(1) = 0D0
      FORCE(2) = 0D0
      FORCE(3) = 0D0

      ACCEL(1) = 0D0
      ACCEL(2) = 0D0
      ACCEL(3) = 0D0

      LFNLOC = 999

C ----------------------------------------
C LOAD SATELLITE PROPERTIES AND CERES DATA
C ----------------------------------------
      IF ((IFIRST==1).AND.(ERM.GT.0)) THEN

C     PROPERTIES FOR ALL SATELLITES BLOCKS
C MOD TAH 190722: Extended loop to handle galileo (201,202) 
C        DO SBLK = 1,105
         DO INDB = 1,max_blk
            SBLK = known_blks(INDB) 
C           IF((SBLK.LE.10).OR.(SBLK.GT.100))THEN 
            IF( SBLK.ne.0 )THEN 

            CALL PROPBOXW(SBLK,AREAS,REFLS,DIFUS,ABSPS,AREA2S,REFL2S,
     1                    DIFU2S,ABSP2S,REFLIRS,DIFUIRS,ABSPIRS)

C           IF(SBLK.LE.10)THEN
C              INDB = SBLK
C           ELSEIF(SBLK.GT.100)THEN
C              INDB = SBLK - 90
C           ENDIF

            DO II = 1,4
               DO JJ = 1,2
                  AREA(II,JJ,INDB) = AREAS(II,JJ)
                  REFL(II,JJ,INDB) = REFLS(II,JJ)
                  DIFU(II,JJ,INDB) = DIFUS(II,JJ)
                  ABSP(II,JJ,INDB) = ABSPS(II,JJ)

                  AREA2(II,JJ,INDB) = AREA2S(II,JJ)
                  REFL2(II,JJ,INDB) = REFL2S(II,JJ)
                  DIFU2(II,JJ,INDB) = DIFU2S(II,JJ)
                  ABSP2(II,JJ,INDB) = ABSP2S(II,JJ)

                  REFLIR(II,JJ,INDB) = REFLIRS(II,JJ)
                  DIFUIR(II,JJ,INDB) = DIFUIRS(II,JJ)
                  ABSPIR(II,JJ,INDB) = ABSPIRS(II,JJ)
               ENDDO
            ENDDO
            ENDIF
         ENDDO


         IF(ERM.EQ.2)THEN

C     REFLECTIVITY 
            CERES_R = 0
            OPEN(UNIT=LFNLOC,FILE=FILREFL,STATUS='UNKNOWN',
     1           FORM='FORMATTED',IOSTAT=IOSTAT)
            DO II=1,72
               READ(LFNLOC,"(A)")LINE
               DO JJ=1,144
                  ITEM = LINE((JJ-1)*25+1:JJ*25)
                  IF (INDEX(ITEM,'NaN') /= 0) THEN
                     CERES_R(II,JJ) = 0D0
                  ELSE
                     READ(ITEM,*) CERES_R(II,JJ)
                  ENDIF
               ENDDO
            ENDDO
            CLOSE(LFNLOC)
C     EMISSIVITY 
            CERES_E = 0
            OPEN(UNIT=LFNLOC,FILE=FILEMIT,STATUS='UNKNOWN',
     1           FORM='FORMATTED',IOSTAT=IOSTAT)
            DO II=1,72
               READ(LFNLOC,"(A)")LINE
               DO JJ=1,144
                  ITEM = LINE((JJ-1)*25+1:JJ*25)
                  IF (INDEX(ITEM,'NaN') /= 0) THEN
                     CERES_E(II,JJ) = 0D0
                  ELSE
                     READ(ITEM,*) CERES_E(II,JJ)
                  ENDIF
               ENDDO
            ENDDO
            CLOSE(LFNLOC)

C     PRE-INTEGRATION, INDEPENDENT OF SATELLITE POSITION
            GRDANG = 2.5D0
            LATIND = 36.5D0
            LONIND = 72.5D0
            IF(GRD.EQ.1)THEN
               GRDNUM = 1D0
               LATKMX = 72
               LONKMX = 144
            ELSEIF(GRD.EQ.2)THEN
               GRDNUM = 2D0
               GRDCER = 2
               LATKMX = 36
               LONKMX = 72
            ELSEIF(GRD.EQ.3)THEN
               GRDNUM = 4D0
               GRDCER = 4
               LATKMX = 18
               LONKMX = 36
            ENDIF
 
            GRDANG = GRDANG*GRDNUM
            LATIND = (LATIND-0.5D0)/GRDNUM + 0.5D0
            LONIND = (LONIND-0.5D0)/GRDNUM + 0.5D0
              
            D_ANG = (PI*GRDANG/180D0)**2

            DO LATK = 1,LATKMX
               DO LONK = 1,LONKMX

                  LAT_IN = (LATK-LATIND)*GRDANG*(PI/180D0)
                  LON_IN = (LONK-LONIND)*GRDANG*(PI/180D0)

C                 Sphere normal vector and differential of area
                  COSLAT = DCOS(LAT_IN)
                  D_AREA_ALL(LATK,LONK) = (TOA**2)*COSLAT*D_ANG
                  V_NS_ALL(LATK,LONK,1) = COSLAT*DCOS(LON_IN)
                  V_NS_ALL(LATK,LONK,2) = COSLAT*DSIN(LON_IN)
                  V_NS_ALL(LATK,LONK,3) = DSIN(LAT_IN)

C                 New matrix of Reflectivity and Emissivity
                  CERGRE(LATK,LONK) = 0D0
                  CERGEM(LATK,LONK) = 0D0
                  IF(GRD.EQ.1)THEN
                     CERGRE(LATK,LONK) = CERES_R(LATK,LONK)
                     CERGEM(LATK,LONK) = CERES_E(LATK,LONK)
                  ELSEIF((GRD.EQ.2).OR.(GRD.EQ.3))THEN
                     DO II = 0,(GRDCER-1)
                        DO JJ = 0,(GRDCER-1)
                           CERGRE(LATK,LONK) = CERGRE(LATK,LONK) 
     1                     + CERES_R(GRDCER*LATK-II,GRDCER*LONK-JJ)
                           CERGEM(LATK,LONK) = CERGEM(LATK,LONK)
     1                     + CERES_E(GRDCER*LATK-II,GRDCER*LONK-JJ)
                        ENDDO
                     ENDDO
                     CERGRE(LATK,LONK) = CERGRE(LATK,LONK)/(GRDNUM**2)
                     CERGEM(LATK,LONK) = CERGEM(LATK,LONK)/(GRDNUM**2)
                  ENDIF

               ENDDO
            ENDDO
C
          ENDIF
         IFIRST = 0
      ENDIF


C --------------------------
C NOMINAL SATELLITE ATTITUDE
C --------------------------

      ABSPOS = DSQRT(YSAT(1)**2+YSAT(2)**2+YSAT(3)**2)
      DO K=1,3
         RADVEC(K) = YSAT(K)/ABSPOS
      ENDDO
      
      CRSVEC(1) = YSAT(2)*YSAT(6)-YSAT(3)*YSAT(5)
      CRSVEC(2) = YSAT(3)*YSAT(4)-YSAT(1)*YSAT(6)
      CRSVEC(3) = YSAT(1)*YSAT(5)-YSAT(2)*YSAT(4)
      ABSCRS = DSQRT(CRSVEC(1)**2+CRSVEC(2)**2+CRSVEC(3)**2)
      DO K=1,3
         CRSVEC(K) = CRSVEC(K)/ABSCRS
      ENDDO

      ALGVEC(1) = CRSVEC(2)*RADVEC(3)-CRSVEC(3)*RADVEC(2)
      ALGVEC(2) = CRSVEC(3)*RADVEC(1)-CRSVEC(1)*RADVEC(3)
      ALGVEC(3) = CRSVEC(1)*RADVEC(2)-CRSVEC(2)*RADVEC(1) 

      IF(ERM.GT.0)THEN

C        DISTANCE FROM SATELLITE TO SUN
         RSUN = DSQRT((YSAT(1)-SUN(1))**2+(YSAT(2)-SUN(2))**2+
     1           (YSAT(3)-SUN(3))**2)

C        D VECTOR AND Z VECTOR
         DO K=1,3
            D_SUN(K) = (SUN(K)-YSAT(K))/RSUN
            Z_SAT(K) = -YSAT(K)/ABSPOS
         ENDDO

C        Y VECTOR
         Y0SAT(1) = Z_SAT(2)*D_SUN(3)-Z_SAT(3)*D_SUN(2)
         Y0SAT(2) = Z_SAT(3)*D_SUN(1)-Z_SAT(1)*D_SUN(3)
         Y0SAT(3) = Z_SAT(1)*D_SUN(2)-Z_SAT(2)*D_SUN(1)
         ABSY0V = DSQRT(Y0SAT(1)**2 + Y0SAT(2)**2 + Y0SAT(3)**2)
         DO K=1,3
            Y_SAT(K) = Y0SAT(K)/ABSY0V
         ENDDO

C        B VECTOR
         B_SUN(1) = Y_SAT(2)*D_SUN(3) - Y_SAT(3)*D_SUN(2)
         B_SUN(2) = Y_SAT(3)*D_SUN(1) - Y_SAT(1)*D_SUN(3)
         B_SUN(3) = Y_SAT(1)*D_SUN(2) - Y_SAT(2)*D_SUN(1)

C        X VECTOR
         X_SAT(1) = Y_SAT(2)*Z_SAT(3) - Y_SAT(3)*Z_SAT(2)
         X_SAT(2) = Y_SAT(3)*Z_SAT(1) - Y_SAT(1)*Z_SAT(3)
         X_SAT(3) = Y_SAT(1)*Z_SAT(2) - Y_SAT(2)*Z_SAT(1)

         DO K=1,3
            ATTSURF(K,1) = Z_SAT(K)
            ATTSURF(K,2) = Y_SAT(K)
            ATTSURF(K,3) = X_SAT(K)
            ATTSURF(K,4) = D_SUN(K)
         ENDDO
      ENDIF

C ---------------------------- 
C OPTICAL PROPERTIES PER BLOCK
C ----------------------------

      IF(ERM.GT.0)THEN
         IF(BLKNUM.LE.10)THEN
            INDB = BLKNUM
         ELSEIF(BLKNUM.GT.100 .and. BLKNUM.lt.200 )THEN
            INDB = BLKNUM - 90
C MOD TAH 190722: Added Galileo
         ELSEIF(BLKNUM.gt.200) then
            indb = blknum - 180
         ENDIF
         DO II = 1,4
            DO JJ = 1,2
               AREAS(II,JJ) = AREA(II,JJ,INDB)
               REFLS(II,JJ) = REFL(II,JJ,INDB)
               DIFUS(II,JJ) = DIFU(II,JJ,INDB)
               ABSPS(II,JJ) = ABSP(II,JJ,INDB)

               AREA2S(II,JJ) = AREA2(II,JJ,INDB)
               REFL2S(II,JJ) = REFL2(II,JJ,INDB)
               DIFU2S(II,JJ) = DIFU2(II,JJ,INDB)
               ABSP2S(II,JJ) = ABSP2(II,JJ,INDB)

               REFLIRS(II,JJ) = REFLIR(II,JJ,INDB)
               DIFUIRS(II,JJ) = DIFUIR(II,JJ,INDB)
               ABSPIRS(II,JJ) = ABSPIR(II,JJ,INDB)
            ENDDO
         ENDDO
      ENDIF

C ----------------------
C EARTH RADIATION MODELS (This part is modified for SRP modeling)
C ----------------------

      IF(ERM.GT.0)THEN
         ABSSUN = DSQRT(SUN(1)**2 + SUN(2)**2 + SUN(3)**2)
         DO K=1,3
            ESUN(K) = SUN(K)/ABSSUN
         ENDDO

         PSIDOT = ESUN(1)*RADVEC(1)+ESUN(2)*RADVEC(2)+ESUN(3)*RADVEC(3)
         IF(DABS(PSIDOT).GT.(1D0-1D-6))THEN
            PSI = 0D0
         ELSE
            PSI = DACOS(PSIDOT)
         ENDIF
         S0 = S0*(AU/ABSSUN)**2
      ENDIF

C     ANALYTICAL MODEL
      IF(ERM.EQ.1)THEN

c         NCFVEC(1) = RADVEC(1)
c         NCFVEC(2) = RADVEC(2)
c         NCFVEC(3) = RADVEC(3)
c         ALBFAC = (PI*TOA**2)*(S0/C)/(ABSPOS**2)
c         PHASEVI = (2*ALB/(3*PI**2))*((PI-PSI)*DCOS(PSI)+DSIN(PSI))
c         PHASEIR = (1-ALB)/(4*PI)
c         ABSNCFVI = ALBFAC*PHASEVI
c         ABSNCFIR = ALBFAC*PHASEIR

!     CHANGE THE RADIATION DIRECTION TO SUN-SAT
         NCFVEC(1) = D_SUN(1)
         NCFVEC(2) = D_SUN(2)
         NCFVEC(3) = D_SUN(3)
         ALBFAC = (S0/C)
         ABSNCFVI = ALBFAC*1.d0
         ABSNCFIR = ALBFAC*1.d0


         CALL SURFBOXW(AREAS,REFLS,DIFUS,ABSPS,
     1                    AREA2S,REFL2S,DIFU2S,ABSP2S,
     2                    REFLIRS,DIFUIRS,ABSPIRS,
     3                    ABSNCFVI,ABSNCFIR,NCFVEC,ATTSURF,FORCE)
              
cd         print *,'In ERPFBOXW blk ',indb
cd         write(*,*) 'areas ', ((AREAS(i,j),i=1,4),j=1,2)
cd         write(*,*) 'refls ', ((REFLS(i,j),i=1,4),j=1,2)
cd         write(*,*) 'difus ', ((DIFUS(i,j),i=1,4),j=1,2)
cd         write(*,*) 'absps ', ((ABSPS(i,j),i=1,4),j=1,2)                     
cd        write(*,*) 'areas2s ',((AREA2S(i,j),i=1,4),j=1,2)
cd         write(*,*) 'refl2s ', ((REFL2S(i,j),i=1,4),j=1,2)
cd         write(*,*) 'difu2s ',  ((DIFU2S(i,j),i=1,4),j=1,2)
cd         write(*,*) 'absp2s',   ((ABSP2S(i,j),i=1,4),j=1,2)
cd         write(*,*) 'reflirs ', ((REFLIRS(i,j),i=1,4),j=1,2)
cd         write(*,*) 'difuirs ', ((DIFUIRS(i,j),i=1,4),j=1,2)
cd         write(*,*) 'abspirs ', (( ABSPIRS(i,j),i=1,4),j=1,2) 
cd         write(*,*) 'absncfvi ',ABSNCFVI
cd         write(*,*) 'absncfir ' ,ABSNCFIR
cd         write(*,*) 'ncfvec ',  ( NCFVEC(i),i=1,3)        
cd         write(*,*) 'attsurf ', (( ATTSURF(i,j),i=1,4),j=1,4)
cd         write(*,*) 'force ', (FORCE(i),i=1,3)

C     NUMERICAL MODEL (CERES DATA)
      ELSEIF(ERM.EQ.2)THEN

         ABSSUN = DSQRT(SUN(1)**2 + SUN(2)**2 + SUN(3)**2)
         DO K=1,3
            ESUN(K) = SUN(K)/ABSSUN
         ENDDO

         DO LATK = 1,LATKMX
            DO LONK = 1,LONKMX

               D_AREA = D_AREA_ALL(LATK,LONK)
               V_NS(1) = V_NS_ALL(LATK,LONK,1)
               V_NS(2) = V_NS_ALL(LATK,LONK,2)
               V_NS(3) = V_NS_ALL(LATK,LONK,3)

               DO II=1,3
                  V_INS(II)=0D0
                  DO JJ=1,3
                     V_INS(II) = V_INS(II) + KAPPA(JJ,II)*V_NS(JJ)
                  ENDDO
               ENDDO

C              Distance and direction from point in the Earth to satellite
               V_DIST(1) = YSAT(1)-TOA*V_INS(1) 
               V_DIST(2) = YSAT(2)-TOA*V_INS(2)
               V_DIST(3) = YSAT(3)-TOA*V_INS(3)
               DIST2 = V_DIST(1)**2 +V_DIST(2)**2 +V_DIST(3)**2
               ABSDIST = DSQRT(DIST2)
               V_SAT(1) = V_DIST(1)/ABSDIST
               V_SAT(2) = V_DIST(2)/ABSDIST
               V_SAT(3) = V_DIST(3)/ABSDIST

C              Cosine of angles of incident and reflected radiation
               COS_IN = ESUN(1)*V_INS(1) + ESUN(2)*V_INS(2) 
     1                + ESUN(3)*V_INS(3)
               COS_RE =V_SAT(1)*V_INS(1) +V_SAT(2)*V_INS(2)
     1                +V_SAT(3)*V_INS(3)

               IF(COS_RE.GE.0)THEN

C              Reflectivity and emissivity coefficients
                  REFL_CF = CERGRE(LATK,LONK)
                  EMIT_CF = CERGEM(LATK,LONK)

C                 Reflected Irradiance
                  IF(COS_IN.GE.0)THEN
                     E_REFL=(REFL_CF/(PI*DIST2))*COS_RE*COS_IN*S0*D_AREA
                  ELSE
                     E_REFL=0D0
                  ENDIF

C                 Emitted Irradiance
                  E_EMIT = (EMIT_CF/(4*PI*DIST2))*COS_RE*S0*D_AREA

C                 Non-conservative force
                  ABSNCFVI = E_REFL/C
                  ABSNCFIR = E_EMIT/C
                  NCFVEC(1) = V_SAT(1)
                  NCFVEC(2) = V_SAT(2)
                  NCFVEC(3) = V_SAT(3)

                  CALL SURFBOXW(AREAS,REFLS,DIFUS,ABSPS,
     1                    AREA2S,REFL2S,DIFU2S,ABSP2S,
     2                    REFLIRS,DIFUIRS,ABSPIRS,
     3                    ABSNCFVI,ABSNCFIR,NCFVEC,ATTSURF,FORCE_LL)

                  FORCE(1) = FORCE(1) + FORCE_LL(1)
                  FORCE(2) = FORCE(2) + FORCE_LL(2)
                  FORCE(3) = FORCE(3) + FORCE_LL(3)
               ENDIF

            ENDDO
         ENDDO

      ENDIF


C     ANTENNA POWER OF GPS SATELLITES (IN WATTS)
C     IGS MODEL (JIM RAY, 2011)

C     GPS BLOCK IIA (ASSUMED THE SAME FOR BLOCK I AND II) 
      IF(BLKNUM.LE.3)THEN
         ANTPOW = 76D0

C     GPS BLOCK IIR
      ELSEIF((BLKNUM.GE.4).AND.(BLKNUM.LE.6))THEN
         ANTPOW = 85D0

C     GPS BLOCK IIR-M
      ELSEIF(BLKNUM.EQ.7)THEN
         ANTPOW = 198D0

C     GPS BLOCK IIF
      ELSEIF(BLKNUM.EQ.8)THEN
         ANTPOW = 249D0
         IF((SVN.EQ.62).AND.(MJD.GE.55656D0))THEN
C           NO M-CODE FOR SVN62/PRN25 STARTING 05APR2011; NANU 2011026
            ANTPOW = 154D0
         ENDIF
      ENDIF

C     NO ANTENNA POWER INFORMATION FOR GLONASS SATELLITES


C     NAVIGATION ANTENNA THRUST (SIMPLE MODEL)
      IF((ANT.EQ.1).AND.(BLKNUM.LT.100))THEN
         ANTFORCE = ANTPOW/C
         FORCE(1) = FORCE(1) + ANTFORCE*RADVEC(1)
         FORCE(2) = FORCE(2) + ANTFORCE*RADVEC(2)
         FORCE(3) = FORCE(3) + ANTFORCE*RADVEC(3)
      ENDIF

      IF((REFF.GT.0).AND.((ERM.GT.0).OR.(ANT.EQ.1)))THEN
         DO K=1,3
            FREF(K) = 0D0
         ENDDO

C        FORCE IN BODY-FIXED REFERENCE FRAME
         IF(REFF.EQ.1)THEN
            DO K=1,3
               FREF(1) = FREF(1) + FORCE(K)*Z_SAT(K)
               FREF(2) = FREF(2) + FORCE(K)*Y_SAT(K)
               FREF(3) = FREF(3) + FORCE(K)*X_SAT(K)
            ENDDO

C        FORCE IN SUN-FIXED REFERENCE FRAME
         ELSEIF(REFF.EQ.2)THEN
            DO K=1,3
               FREF(1) = FREF(1) + FORCE(K)*D_SUN(K)
               FREF(2) = FREF(2) + FORCE(K)*Y_SAT(K)
               FREF(3) = FREF(3) + FORCE(K)*B_SUN(K)
            ENDDO

C        FORCE IN ORBITAL REFERENCE FRAME
         ELSEIF(REFF.EQ.3)THEN
            DO K=1,3
               FREF(1) = FREF(1) + FORCE(K)*RADVEC(K)
               FREF(2) = FREF(2) + FORCE(K)*ALGVEC(K)
               FREF(3) = FREF(3) + FORCE(K)*CRSVEC(K)
            ENDDO
         ENDIF

         DO K=1,3
            FORCE(K) = FREF(K)
         ENDDO
      ENDIF


C     MASS OF SATELLITES
* MOD TAH 190722: Updated the masses for Block IIF, IIIA, and
*     Galileo.  Made values consistent with igs_metadata.snx
      IF(BLKNUM.EQ.1)THEN
         MASS = 455D0
      ELSEIF(BLKNUM.EQ.2)THEN 
         MASS = 843D0
      ELSEIF(BLKNUM.EQ.3)THEN
         MASS = 930D0
      ELSEIF((BLKNUM.GE.4).AND.(BLKNUM.LE.7))THEN
         MASS = 1080D0
      ELSEIF(BLKNUM.EQ.8)THEN
         MASS = 1633D0
      ELSEIF(BLKNUM.EQ.9)THEN
         MASS = 2161D0
* GLONASS
      ELSEIF((BLKNUM.EQ.101).OR.(BLKNUM.EQ.102))THEN
         MASS = 1415D0
      ELSEIF(BLKNUM.EQ.103)THEN
         MASS = 995D0
* GALILEO
      ELSEIF(BLKNUM.eq.201) then
         mass =  695.D0    ! Average value
      ELSEIF(BLKNUM.eq.202) then
*        Need to treat hi-eccentricity SV with different mass
         if ( svn.ge.201 .and. svn.le.202 ) then
            mass = 660.d0   ! Typical value
         else
            mass = 710.d0   ! Average 19/07/22 == 708.597kg
         endif
      else
         print *,'No MASS for BLKNUM ',BLKNUM
         stop 'NO MASS: ERPFBOXW'
      ENDIF

C     CONVERSION TO ACCELERATION
      DO K=1,3
         ACCEL(K) = FORCE(K)/MASS
      ENDDO
                    
cd      write(*,*) 'force mass accel '
cd     .   ,(force(i),i=1,3),mass,(accel(i),i=1,3)
      END SUBROUTINE

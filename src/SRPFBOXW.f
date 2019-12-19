C*
      SUBROUTINE SRPFBOXW(REFF,YSAT,SUN,BLKNUM,SVN,ACCEL)
CC
CC NAME       :  SRPFBOXW
CC
CC PURPOSE    :  COMPUTATION OF EARTH RADIATION PRESSURE ACTING ON A
CC               BOW-WING SATELLITE
CC
CC PARAMETERS :
CC         IN :  
CC               REFF      : ACCELERATION IN REFERENCE FRAME:
CC                           0 = INERTIAL
CC                           1 = BODY FIXED (Z,Y,X)
CC                           2 = SUN FIXED (D,Y,B)
CC                           3 = ORBITAL (RADIAL, ALONG- AND CROSS-TRACK)
CC               YSAT      : SATELLITE POSITION [m] AND VELOCITY [m/s] (INERTIAL), 
CC                           (POSITION,VELOCITY) = (RX,RY,RZ,VX,VY,VZ)
CC               SUN       : SUN POSITION VECTOR [m] (INERTIAL)
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
      INTEGER*4 BLKNUM,REFF,SBLK,INDB,SVN
      INTEGER*4 JFIRST,II,JJ,K
      INTEGER*4 LATK,LONK,LATKMX,LONKMX
      integer*4 i,j
C
      REAL*8 PI,C,AU,S0,TOA,ALB,MASS
      REAL*8 YSAT(6),SUN(3),FORCE(3),ACCEL(3)
C
C     EARTH AND SATELLITE PROPERTIES
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
      REAL*8 D_ANG,LATIND,LONIND
      REAL*8 D_AREA,LAT_IN,LON_IN,DIST2,ABSDIST
      REAL*8 COSLAT,PSI
      REAL*8 V_NS(3),V_INS(3),V_DIST(3),V_SAT(3)
      REAL*8 COS_IN,COS_RE,PSIDOT,ABSSUN
      REAL*8 REFL_CF,EMIT_CF,E_REFL,E_EMIT
      REAL*8 FORCE_LL(3),FREF(3)
      REAL*8 ANTFORCE,ANTPOW,MJD
C
      
C MOD TAH 190722: introduce "known_blks" arrway with know block
C     types for GPS, Glonass and Galileo.  Needs to be extended
C     when other constellations added.
      integer*4 known_blks(max_blk)  ! Block numbers by the type
                  ! we known and have coded.  0 value unknown

      data known_blks /  1,   2,   3,   4,   5,  6,  7,  8,  9,  0,
     .                 101, 102, 103,   0,   0,  0,  0,  0,  0,  0,
     .                 201, 202,   0,   0,   0,  0,  0,  0,  0,  0 /
C
      DATA JFIRST/1/

C DATA FILES PATH AND NAME (CHANGE TO LOCAL PATH)
            


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

C ----------------------------------------
C LOAD SATELLITE PROPERTIES AND CERES DATA
C ----------------------------------------
      IF ((JFIRST==1)) THEN

C     PROPERTIES FOR ALL SATELLITES BLOCKS
C MOD TAH 190722: Extended loop to handle galileo (201,202) 
C        DO SBLK = 1,105
         DO INDB = 1,max_blk
            SBLK = known_blks(INDB) 
C           IF((SBLK.LE.10).OR.(SBLK.GT.100))THEN 
            IF( SBLK.ne.0 )THEN 

            CALL PROPBOXW(SBLK,AREAS,REFLS,DIFUS,ABSPS,AREA2S,REFL2S,
     1                    DIFU2S,ABSP2S,REFLIRS,DIFUIRS,ABSPIRS)


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

!         JFIRST = 0
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


C ---------------------------- 
C OPTICAL PROPERTIES PER BLOCK
C ----------------------------


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


C ----------------------
C EARTH RADIATION MODELS (This part is modified for SRP modeling)
C ----------------------

c         NCFVEC(1) = RADVEC(1)
c         NCFVEC(2) = RADVEC(2)
c         NCFVEC(3) = RADVEC(3)
c         ALBFAC = (PI*TOA**2)*(S0/C)/(ABSPOS**2)
c         PHASEVI = (2*ALB/(3*PI**2))*((PI-PSI)*DCOS(PSI)+DSIN(PSI))
c         PHASEIR = (1-ALB)/(4*PI)
c         ABSNCFVI = ALBFAC*PHASEVI
c         ABSNCFIR = ALBFAC*PHASEIR

!     CHANGE THE RADIATION DIRECTION TO SUN->SAT
         NCFVEC(1) = (-1.d0)*D_SUN(1)
         NCFVEC(2) = (-1.d0)*D_SUN(2)
         NCFVEC(3) = (-1.d0)*D_SUN(3)
         ALBFAC = (S0/C)
         ABSNCFVI = ALBFAC*1.d0
         ABSNCFIR = ALBFAC*1.d0


         CALL SURFBOXW(AREAS,REFLS,DIFUS,ABSPS,
     1                    AREA2S,REFL2S,DIFU2S,ABSP2S,
     2                    REFLIRS,DIFUIRS,ABSPIRS,
     3                    ABSNCFVI,ABSNCFIR,NCFVEC,ATTSURF,FORCE)
              
cd         print *,'In ERPFBOXW blk ',indb
cd         write(*,*) 'areas ',((AREAS(i,j),i=1,4),j=1,2)
cd         write(*,*) 'refls ', ((REFLS(i,j),i=1,4),j=1,2)
cd         write(*,*) 'difus ', ((DIFUS(i,j),i=1,4),j=1,2)
cd         write(*,*) 'absps ', ((ABSPS(i,j),i=1,4),j=1,2)                     
cd         write(*,*) 'areas2s ',((AREA2S(i,j),i=1,4),j=1,2)
cd         write(*,*) 'albfac, C, S0', ALBFAC, C, S0
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


      IF(REFF.GT.0)THEN
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
!      IF(BLKNUM.EQ.1)THEN
!         MASS = 455D0
!      ELSEIF(BLKNUM.EQ.2)THEN 
!         MASS = 843D0
!      ELSEIF(BLKNUM.EQ.3)THEN
!         MASS = 930D0
!      ELSEIF((BLKNUM.GE.4).AND.(BLKNUM.LE.7))THEN
!         MASS = 1080D0
!      ELSEIF(BLKNUM.EQ.8)THEN
!         MASS = 1633D0
!      ELSEIF(BLKNUM.EQ.9)THEN
!         MASS = 2161D0
* GLONASS
!      ELSEIF((BLKNUM.EQ.101).OR.(BLKNUM.EQ.102))THEN
!         MASS = 1415D0
!      ELSEIF(BLKNUM.EQ.103)THEN
!         MASS = 995D0
* GALILEO
!      ELSEIF(BLKNUM.eq.201) then
!         mass =  695.D0    ! Average value
!      ELSEIF(BLKNUM.eq.202) then
*        Need to treat hi-eccentricity SV with different mass
!         if ( svn.ge.201 .and. svn.le.202 ) then
!            mass = 660.d0   ! Typical value
!         else
!            mass = 710.d0   ! Average 19/07/22 == 708.597kg
!         endif
* BDS
!      ELSEIF(BLKNUM.eq.301) then
!         MASS = 1550D0
!      ELSEIF(BLKNUM.eq.302 .or.BLKNUM.eq.303) then
!         MASS = 1900D0

!      else
!         print *,'No MASS for BLKNUM ',BLKNUM
!         stop 'NO MASS: SRPFBOXW'
!      ENDIF

C     CONVERSION TO ACCELERATION
      DO K=1,3
!         ACCEL(K) = FORCE(K)/MASS
         ACCEL(K) = FORCE(K)
      ENDDO
                    
cd      write(*,*) 'force mass accel '
cd     .   ,(force(i),i=1,3),mass,(accel(i),i=1,3)
      END SUBROUTINE

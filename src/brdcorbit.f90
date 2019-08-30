PROGRAM  brdcorbit 

! ----------------------------------------------------------------------
! PROGRAM: brdcorbit.f90
! ----------------------------------------------------------------------
! Purpose: Create a sp3 format file from the broadcast ephemeris *.yyn file, 
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
! ----------------------------------------------------------------------
      USE mdl_precision
      USE mdl_num
      USE mdl_planets
      USE m_shadow
      USE mdl_param
      IMPLICIT NONE


!-------------------------------------------------------------------
      INTEGER (KIND = prec_int4) :: prnnum
      INTEGER (KIND = 4)         :: satsvn
      REAL (KIND = prec_d), INTENT(IN) :: mjd
      REAL (KIND = prec_d), DIMENSION(3), INTENT(IN) :: rsat, vsat
      REAL (KIND = prec_q), INTENT(OUT) :: beta
      REAL (KIND = prec_q) :: lambda, del_u, yaw
!--------------------------------------------------------------------
      DOUBLE PRECISION  JD, Zbody(6)
      INTEGER (KIND = 4) :: i, j
      REAL (KIND = prec_q), DIMENSION(3) :: rbody
      REAL (KIND = prec_q), DIMENSION(3) :: rSun, rMoon
      REAL (KIND = prec_q), DIMENSION(3) :: r_sun1, r_sun2
      REAL (KIND = prec_q) :: u_sun, yaw2
      INTEGER  NTARG_SUN, NTARG_MOON, NCTR
      REAL (KIND = prec_q), DIMENSION(3) :: ed, ey, eb, ex, en, ev, ez, et
      REAL (KIND = prec_q), DIMENSION(3) :: yy
      REAL (KIND = prec_q), DIMENSION(9) :: kepler
      REAL (KIND = prec_q), DIMENSION(4) :: cosang
      REAL (KIND = prec_q) :: R11(3,3),R33(3,3)
      REAL (KIND = prec_q) :: Pi, Ps, AU, sclfa
      REAL (KIND = prec_d) :: GM, Ds
      REAL (KIND = prec_d) :: u_sat, i_sat, omega_sat
      INTEGER (KIND = prec_int2) :: AllocateStatus

      INTEGER (KIND = prec_int2) :: UNIT_IN
      REAL (KIND = prec_d) :: CPU_t0, CPU_t1
! START CPU COUNTER
! -----------------
CALL cpu_time (CPU_t0)

      UNIT_IN = 66
! ----------------------------------------------------------------------
! Open *.yyn broadcast ephemeris file
      OPEN (UNIT = UNIT_IN, FILE = TRIM (PRMfname), IOSTAT = ios)
      IF (ios /= 0) THEN
         PRINT *, "Error in opening file:", PRMfname
         PRINT *, "OPEN IOSTAT=", ios
      END IF

! Read broadcast ephemeris file (GPS-only)
!-----------------------------------------
CALL readbrdcheader (UNIT_IN, VERSION, IONALPHA, IONBETA, A0UTC, A1UTC, &
                     ITUTC, NWKUTC, LEAP)

CALL readbrdcmessg (UNIT_IN,ISAT,IREC,EPH,CLK)

CLOSE (UNIT=UNIT_IN)


! DEFINE GPS-UTC (DO NOT TRUST RINEX FILES)
! The leap second needs to be checked!! (by using dat.for)
! -----------------------------------------
         LEAPTST=DGPSUT(TFIRST)
         IF (LEAP.EQ.0 .OR. LEAPTST.NE.LEAP) LEAP=LEAPTST
         NEPO=IDINT((TLAST-TFIRST)/DTTAB*86400.D0)+1
! DEFINE ASCENDING ORDER FOR ALL SATELLITE
! ----------------------------------------
         CALL IORDUP(NRSNEW,NSANEW,ISATN)
         CALL IORDUP(NRGASNEW,NGASANEW,IGASATN)
         CALL IORDUP(NRGSAN,NGSANW,ISATG)

! WRITE HEADER OF PRECISE FILE
! ----------------------------
         NSTOT=NSANEW+NGASANEW+NGSANW
         DO ISAT=1,NSANEW
           NRSTOT(ISAT)=NRSNEW(ISATN(ISAT))
         ENDDO
         DO ISAT=1,NGASANEW
           NRSTOT(ISAT+NSANEW)=NRGASNEW(IGASATN(ISAT))
         ENDDO
         DO ISAT=1,NGSANW
           NRSTOT(ISAT+NGASANEW+NSANEW)=NRGSAN(ISATG(ISAT))
         ENDDO

         DATDES = 'U    '
         ORBTYP = 'BRD'
         TIMSYS = 'GPS'

         BASPOS = 1.25D0
         BASCLK = 1.025D0
         VEL    = 0D0
         ACCPOS = 0
         ACCVEL = 0
         EVTFLG = ' '
         SDEVP  = 0D0
         SDEVV  = 0D0
         CORRP  = 0D0
         CORRV  = 0D0

         CALL WTPREH(FILTAB(4,IFIL),LFNPRE,IFRMAT,NSTOT,NRSTOT,SATWGT,&
                TFIRST,NEPO,DTTAB,TITLE,DATDES,COOSYS,ORBTYP, AGENCY,&
                TIMSYS,BASPOS,BASCLK)
CALL brd2ecef(SECEPO,EPHNEW(L0+1,JSAT),POS(1,ISAT))
! SATELLITE CLOCK CORRRECTION
! ----------------------------------
DT=TEPO-GPSMJD(CLKNEW(L0+11,JSAT),IDNINT(CLKNEW(L0+1,JSAT)))
DT=DT*86400.D0
DTSATC(ISAT)=CLKNEW(L0+14,JSAT)+CLKNEW(L0+13,JSAT)*DT+CLKNEW(L0+12,JSAT)*DT**2

! WRITE VALUES OF EPOCH "TEPO" INTO PRECISE FILE
! ----------------------------------------------
         CALL WTPREI(LFNPRE,IFRMAT,(/0,0/),NSTOT,NRSTOT,TMJD, POS,VEL,DTSATC,DDTSAT,ACCPOS,ACCVEL,EVTFLG,&
                    SDEVP,SDEVV,CORRP,CORRV,IRCODE)

CALL cpu_time (CPU_t1)
END PROGRAM 

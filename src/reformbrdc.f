      SUBROUTINE reformbrdc(EPHDAT,EPH,CLOCK)
CC
CC NAME       :  reformbrdc
CC
CC PURPOSE    :  rearrange the broadcast information into orbit-related
CC               and clock related groups (REFER TO BERNESE RXV3BR.f)
CC
CC PARAMETERS :
CC         IN :  EPHDAT(1:28)          : ARRAY CONTAINING THE        R*8
CC                                       RINEX NAVIGATION MESSAGE
CC                 EPHDAT(1)           : TOC
CC                 EPHDAT(2)           : AF0
CC                 EPHDAT(3)           : AF1
CC                 EPHDAT(4)           : AF2
CC                 EPHDAT(5)           : AODE
CC                 EPHDAT(6)           : CRS
CC                 EPHDAT(7)           : DELTA N
CC                 EPHDAT(8)           : M0
CC                 EPHDAT(9)           : CUC
CC                 EPHDAT(10)          : E
CC                 EPHDAT(11)          : CUS
CC                 EPHDAT(12)          : SQRT(A)
CC                 EPHDAT(13)          : TOE
CC                 EPHDAT(14)          : CIC
CC                 EPHDAT(15)          : OMEGA0 (R.A.ASCEND NODE)
CC                 EPHDAT(16)          : CIS
CC                 EPHDAT(17)          : I0
CC                 EPHDAT(18)          : CRC
CC                 EPHDAT(19)          : SMALL OMEGA0 (PERIGEE)
CC                 EPHDAT(20)          : OMEGA DOT
CC                 EPHDAT(21)          : I DOT
CC                 EPHDAT(22)          : L2CODE
CC                 EPHDAT(23)          : GPS WEEK (BELONGS TO TOE)
CC                 EPHDAT(24)          : L2 P DATA FLAG
CC                 EPHDAT(25)          : SV ACCURACY
CC                 EPHDAT(26)          : SV HEALTH (MSB ONLY)
CC                 EPHDAT(27)          : TGD  (SEC)
CC                 EPHDAT(28)          : AODC (SEC)
CC        OUT :  EPH(I),I=1,2,..,20    : ARRAY CONTAINING THE        R*8
CC                                       EPHEMERIS-INFORMATION
CC                 EPH(1)              : GPS-WEEK
CC                 EPH(2)              : T0E
CC                 EPH(3)              : A
CC                 EPH(4)              : E
CC                 EPH(5)              : I
CC                 EPH(6)              : R.A. OF ASCENDING NODE
CC                 EPH(7)              : PERIGEE
CC                 EPH(8)              : MEAN ANOMALY (T0E)
CC                 EPH(9)              : DN (CORRECTION TO MEAN MOTION)
CC                 EPH(10)             : RATE OF NODE
CC                 EPH(11)             : CUS
CC                 EPH(12)             : CUC
CC                 EPH(13)             : CRS
CC                 EPH(14)             : CRC
CC                 EPH(15)             : CIS
CC                 EPH(16)             : CIC
CC                 EPH(17)             : AODE
CC                 EPH(18)             : IDOT
CC                 EPH(I),I=19,20      : NOT USED
CC               CLOCK(I),I=1,2,..,20  : ARRAY CONTAINING THE        R*8
CC                                       SATELLITE CLOCK INFORMATION
CC                 CLOCK(1)            : GPS-WEEK
CC                 CLOCK(2)            : L2 CODE INDICATOR
CC                 CLOCK(3)            : USER RANGE ACCURACY (M)
CC                 CLOCK(4)            : SV HEALTH MSB (NAVIGATION DATA)
CC                 CLOCK(5)            : SV HEALTH LSB'S
CC                 CLOCK(6)            : L2 P DATA FLAG
CC                 CLOCK(7)            : NOT USED
CC                 CLOCK(8)            : NOT USED
CC                 CLOCK(9)            : TGD
CC                 CLOCK(10)           : AODC
CC                 CLOCK(11)           : TOC
CC                 CLOCK(12)           : A2
CC                 CLOCK(13)           : A1
CC                 CLOCK(14)           : A0
CC                 CLOCK(I),I=15,...,20: NOT USED
CC
CC
CC AUTHOR     :  W. GURTNER
CC
CC VERSION    :  3.4  (JAN 93)
CC
CC CREATED    :  89/04/05

C*
      IMPLICIT NONE
C
C DECLARATIONS INSTEAD OF IMPLICIT
C --------------------------------
      INTEGER*4 I
C
CCC       IMPLICIT REAL*8 (A-H,O-Z)
C
      REAL*8      CLOCK(*),EPH(*),EPHDAT(*)
C
      DO 10 I = 1,20
       CLOCK(I) = 0.D0
10    CONTINUE
      EPH(19) = 0.D0
      EPH(20) = 0.D0
C
C COPY EPHEMERIS AND CLOCK CORRECTIONS TO EPH AND CLOCK ARRAYS
C ----------------------------------------------------------------
      EPH(1)  = EPHDAT(23)
C
C Rounding necessary???
C ---------------------
c      EPH(2)  = EPHDAT(13)
      EPH(2)  = DNINT(EPHDAT(13)*1D6)/1D6
      EPH(3)  = EPHDAT(12)**2
      EPH(4)  = EPHDAT(10)
      EPH(5)  = EPHDAT(17)
      EPH(6)  = EPHDAT(15)
      EPH(7)  = EPHDAT(19)
      EPH(8)  = EPHDAT(8)
      EPH(9)  = EPHDAT(7)
      EPH(10) = EPHDAT(20)
      EPH(11) = EPHDAT(11)
      EPH(12) = EPHDAT(9)
      EPH(13) = EPHDAT(6)
      EPH(14) = EPHDAT(18)
      EPH(15) = EPHDAT(16)
      EPH(16) = EPHDAT(14)
      EPH(17) = EPHDAT(5)
      EPH(18) = EPHDAT(21)
C
      CLOCK(1)  = EPHDAT(23)
      CLOCK(2)  = EPHDAT(22)
      CLOCK(3)  = EPHDAT(25)
      CLOCK(4)  = EPHDAT(26)
      CLOCK(6)  = EPHDAT(24)
      CLOCK(9)  = EPHDAT(27)
      CLOCK(10) = EPHDAT(28)
      CLOCK(11) = EPHDAT(1)
      IF(CLOCK(11)-EPH(2).GT.+302400.D0) CLOCK(1)=CLOCK(1)-1
      IF(CLOCK(11)-EPH(2).LT.-302400.D0) CLOCK(1)=CLOCK(1)+1
      CLOCK(12) = EPHDAT(4)
      CLOCK(13) = EPHDAT(3)
      CLOCK(14) = EPHDAT(2)
C
      RETURN
      END SUBROUTINE


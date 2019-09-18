      SUBROUTINE reformbrdc(MAXNPAR,eph1,eph2)
CC
CC NAME       :  reformbrdc
CC
CC PURPOSE    :  rearrange the broadcast information into orbit-related
CC               and clock related groups. 
CC
CC PARAMETERS :
CC         IN :  eph1          : Broadcast dynamic elements 
CC                                       
CC        OUT :  eph2          : the rearranged array
CC
CC
CC Author:  Tzupang Tseng
CC
CC

C*
      USE mdl_precision
      IMPLICIT NONE
C
      INTEGER (KIND = prec_int4) :: MAXNPAR
      REAL    (KIND = prec_q)    :: eph2(MAXNPAR)
      REAL    (KIND = prec_q)    :: eph1(MAXNPAR)

C
C
C Rearrange EPH and CLOCK 
C ------------------------
      eph2(1)  = eph1(23)                ! GPS week              
      eph2(2)  = DNINT(eph1(13)*1D6)/1D6 ! toe
      eph2(3)  = eph1(12)**2             ! semi-major axis
      eph2(4)  = eph1(10)                ! eccentricity
      eph2(5)  = eph1(17)                ! inclination
      eph2(6)  = eph1(15)                ! right ascension of ascending node 
      eph2(7)  = eph1(19)                ! perigee 
      eph2(8)  = eph1(8)                 ! mean anomaly 
      eph2(9)  = eph1(7)                 ! correction to mean motion 
      eph2(10) = eph1(20)                ! rate of node 
      eph2(11) = eph1(11)                ! cus
      eph2(12) = eph1(9)                 ! cuc
      eph2(13) = eph1(6)                 ! crs
      eph2(14) = eph1(18)                ! crc
      eph2(15) = eph1(16)                ! cis
      eph2(16) = eph1(14)                ! cic
      eph2(17) = eph1(5)                 ! aode
      eph2(18) = eph1(21)                ! idot 
      eph2(19) = 0.d0
      eph2(20) = 0.d0

      RETURN
      END SUBROUTINE


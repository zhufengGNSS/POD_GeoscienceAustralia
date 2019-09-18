SUBROUTINE  glnorbint(EPH,TEPO,EPOEPH,POS)

! NAME       :  glnorbint
!
! PURPOSE    :  Computation of Glonass orbit and clock from 
!               the interpolation method specified in GLONASS-ICD document
!
!
! PARAMETERS :
!         IN :  EPH    : Broadcast message for GLONASS 
!               TEPO   : Desire epoch for new positions from interpolation  
!               EPOEPH : reference epoch for interpolation 
!        OUT :  POS      : new positions at the desired epoch 
!                        POS(1)=X
!                        POS(2)=Y
!                        POS(3)=Z
!
!
! Author     :  Tzupang Tseng
!
!
! Created    :  18-09-2019
!
! Copyright  :  Geoscinece Australia


USE mdl_precision
USE mdl_num
USE mdl_param
     
IMPLICIT NONE

REAL (KIND = prec_q) :: EPH(*)
REAL (KIND = prec_q) :: POS(3)
REAL (KIND = prec_q) :: X(3), Y(3), Z(3)
REAL (KIND = prec_q) :: EPOEPH, TEPO, intv, ACC
REAL (KIND = prec_q) :: hstep(500), h
INTEGER (KIND=4)     :: i, NSTEP

REAL (KIND = prec_q) :: XACC2, YACC2, ZACC2
REAL (KIND = prec_q) :: XACC3, YACC3, ZACC3
REAL (KIND = prec_q) :: XACC4, YACC4, ZACC4
REAL (KIND = prec_q) :: XACC5, YACC5, ZACC5
REAL (KIND = prec_q) :: XVEL2, YVEL2, ZVEL2
REAL (KIND = prec_q) :: XVEL3, YVEL3, ZVEL3
REAL (KIND = prec_q) :: XVEL4, YVEL4, ZVEL4
REAL (KIND = prec_q) :: XPOS2, YPOS2, ZPOS2
REAL (KIND = prec_q) :: XPOS3, YPOS3, ZPOS3
REAL (KIND = prec_q) :: XPOS4, YPOS4, ZPOS4
REAL (KIND = prec_q) :: XACCTOTAL, YACCTOTAL, ZACCTOTAL
REAL (KIND = prec_q) :: XPOSFINAL, YPOSFINAL, ZPOSFINAL

! ----------------------------------------------------------------------

! Integration interval (sec)
intv = 10.d0

! X coordinates in meter
X(1) = EPH(5)*1000 ! position
X(2) = EPH(6)*1000 ! velocity
X(3) = EPH(7)*1000 ! acceleration

! Y coordinates in meter
Y(1) = EPH(9)*1000 ! position
Y(2) = EPH(10)*1000 ! velocity
Y(3) = EPH(11)*1000 ! acceleration

! Z coordinates in meter
Z(1) = EPH(13)*1000 ! position
Z(2) = EPH(14)*1000 ! velocity
Z(3) = EPH(15)*1000 ! acceleration

! Numerical integration
! ---------------------
NSTEP = ABS(INT(TEPO-EPOEPH)/intv)+1
do i = 1, NSTEP
hstep(i) = intv
if(TEPO-EPOEPH < 0.d0) hstep(i) = hstep(i)*(-1.d0)
h = hstep(i)

CALL glnacc(1,X(1),Y(1),Z(1),X(2),Y(2),ACC)
XACC2=ACC+X(3)
XVEL2=X(2)+XACC2*h/2
XPOS2=X(1)+XVEL2*h/2

CALL glnacc(2,X(1),Y(1),Z(1),X(2),Y(2),ACC)
YACC2=ACC+Y(3)
YVEL2=Y(2)+YACC2*h/2
YPOS2=Y(1)+YVEL2*h/2

CALL glnacc(3,X(1),Y(1),Z(1),X(2),Y(2),ACC)
ZACC2=ACC+Z(3)
ZVEL2=Z(2)+ZACC2*h/2
ZPOS2=Z(1)+ZVEL2*h/2


CALL glnacc(1,XPOS2,YPOS2,ZPOS2,XVEL2,YVEL2,ACC)
XACC3=ACC+X(3)
XVEL3=X(2)+XACC3*h/2
XPOS3=X(1)+XVEL3*h/2

CALL glnacc(2,XPOS2,YPOS2,ZPOS2,XVEL2,YVEL2,ACC)
YACC3=ACC+Y(3)
YVEL3=Y(2)+YACC3*h/2
YPOS3=Y(1)+YVEL3*h/2

CALL glnacc(3,XPOS2,YPOS2,ZPOS2,XVEL2,YVEL2,ACC)
ZACC3=ACC+Z(3)
ZVEL3=Z(2)+ZACC3*h/2
ZPOS3=Z(1)+ZVEL3*h/2


CALL glnacc(1,XPOS3,YPOS3,ZPOS3,XVEL3,YVEL3,ACC)
XACC4=ACC+X(3)
XVEL4=X(2)+XACC4*h
XPOS4=X(1)+XVEL4*h

CALL glnacc(2,XPOS3,YPOS3,ZPOS3,XVEL3,YVEL3,ACC)
YACC4=ACC+Y(3)
YVEL4=Y(2)+YACC4*h
YPOS4=Y(1)+YVEL4*h

CALL glnacc(3,XPOS3,YPOS3,ZPOS3,XVEL3,YVEL3,ACC)
ZACC4=ACC+Z(3)
ZVEL4=Z(2)+ZACC4*h
ZPOS4=Z(1)+ZVEL4*h



CALL glnacc(1,XPOS4,YPOS4,ZPOS4,XVEL4,YVEL4,ACC)
XACC5=ACC+X(3)
CALL glnacc(2,XPOS4,YPOS4,ZPOS4,XVEL4,YVEL4,ACC)
YACC5=ACC+Y(3)
CALL glnacc(3,XPOS4,YPOS4,ZPOS4,XVEL4,YVEL4,ACC)
ZACC5=ACC+Z(3)

XACCTOTAL = (XACC2+2*XACC3+2*XACC4+XACC5)/6
YACCTOTAL = (YACC2+2*YACC3+2*YACC4+YACC5)/6
ZACCTOTAL = (ZACC2+2*ZACC3+2*ZACC4+ZACC5)/6

XPOSFINAL = X(1) + XACCTOTAL*h
YPOSFINAL = Y(1) + YACCTOTAL*h
ZPOSFINAL = Z(1) + ZACCTOTAL*h


end do

POS(1) = XPOSFINAL
POS(2) = YPOSFINAL
POS(3) = ZPOSFINAL

END

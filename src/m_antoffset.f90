MODULE m_antoffset

      IMPLICIT NONE
      !SAVE
Contains


SUBROUTINE antoffset(R1,R2)


! ----------------------------------------------------------------------
! SUBROUTINE: antoffset.f90
! ----------------------------------------------------------------------
! Purpose:
!  Correct the position of transmitting antenna to the center of mass of satellite
!
! NOTE: If the orbital information from SP3 file is based on the center of mass
! of satellite, then this routine cannot be used for the position correction.
! Conversely, if the satellite position is given in the phase center of transmitting antenna, 
! then this routine is used to correct the poistion from antenna to the center of mass.
! ----------------------------------------------------------------------
! Input:
! - R1 : The position of the transmitting antenna 
!
! Output: 
! - R2 : The position of center of mass of satellite 
!
! ----------------------------------------------------------------------
! Author : Dr. Tzupang Tseng, Geoscience Australia         
!
! Create: 21-08-2019
!
! ----------------------------------------------------------------------

USE mdl_precision
IMPLICIT NONE

! ---------------------------------------------------------------------------
REAL (KIND = prec_q), DIMENSION(3) :: R1, R2
REAL (KIND = prec_q), DIMENSION(3) :: XTF
REAL (KIND = prec_q), DIMENSION(3) :: ANTOFF
REAL (KIND = prec_q), DIMENSION(3) :: UX, UY, UZ
REAL (KIND = prec_q) :: D
INTEGER (KIND = 4)   :: K
! ---------------------------------------------------------------------------
D  = SQRT(R1(1)**2+R1(2)**2+R1(3)**2)
UX = R1(1)/D
UY = R1(2)/D
UZ = R1(3)/D 

! Get antenna offset (the following information is used for GLN testing)
! ----------------------------------------------------------------------
ANTOFF(1) = -0.545d0
ANTOFF(2) =  0.d0
ANTOFF(3) =  2.5d0

DO K=1,3
R2(K)=R1(K)-UX(K)*ANTOFF(1) &
           -UY(K)*ANTOFF(2) & 
           -UZ(K)*ANTOFF(3)
END DO
        
END
END

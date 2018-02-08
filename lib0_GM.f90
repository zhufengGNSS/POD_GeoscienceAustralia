PROGRAM lib0_GM
	  
	  
      USE mdl_precision
      USE mdl_num
      IMPLICIT NONE
  
! ----------------------------------------------------------------------
! Variables declaration
! ----------------------------------------------------------------------
      REAL (KIND = prec_q) :: x,y,z, phi,lamda,radius
      REAL (KIND = prec_q), DIMENSION(3) :: r
      REAL (KIND = prec_q) :: fx,fy,fz
! ----------------------------------------------------------------------  
      REAL (KIND = prec_q) :: a1,a2,a3, atan_argument
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
      x = 36745.6783123456789   !False: State D0 or Q0
      y = 19034.1278123423414Q0
      z = 19034.1278123423414D0
!      z = 54678.9934563411223D0
	  ! r = (x,y,z)
      !DATA r(:,1)/x,y,z/
	  r(1) = x
	  r(2) = y
	  r(3) = z
      print *,"XYZ=", x, y, z
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! arctan subroutine test
! ----------------------------------------------------------------------
      CALL arctan (y,x, lamda)
      print *,"lamda=", lamda
      print *,'Kind: ', kind(lamda)
! ----------------------------------------------------------------------	  
!       atan_argument = abs( y/x )  
!       a1 = atan( abs( y/x ) )
!       a2 = atan2(y,x)
!       a3 = qatan( abs( y/x ) )	  
!       print *,"a=", a1
!       print *,"a=", a2
!       print *,"a=", a3
!       print *,"atan_argument=", atan_argument
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! coord_r2sph subroutine test
! ----------------------------------------------------------------------
      CALL coord_r2sph (r , phi,lamda,radius)	  
      print *,"phi=", phi
      print *,'Kind: ', kind(phi)
      print *,"lamda=", lamda
      print *,'Kind: ', kind(lamda)
      print *,"radius=", radius
      print *,'Kind: ', kind(radius)
! ----------------------------------------------------------------------
      
! ----------------------------------------------------------------------
! force_gm subroutine test
! ----------------------------------------------------------------------
      CALL force_gm ( r, fx,fy,fz )	  
      print *,"fx=", fx
      print *,'Kind: ', kind(fx)
      print *,"fy=", fy
      print *,'Kind: ', kind(fy)
      print *,"fz=", fz
      print *,'Kind: ', kind(fz)
! ----------------------------------------------------------------------	  
	  
END
	  
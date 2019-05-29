MODULE m_integrVEQ


! ----------------------------------------------------------------------
! MODULE: m_integrVEQ.f03
! ----------------------------------------------------------------------
! Purpose:
!  Module for numerical integration of the Variational Equations 
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
!			Cooperative Research Centre for Spatial Information, Australia
! Created:	22 January 2018
! ----------------------------------------------------------------------


      IMPLICIT NONE
      !SAVE 			
	  


Contains


SUBROUTINE integr_VEQ (MJDo, tsec_start, ro, vo, arc, integID, step, Nparam, orbc, Smatrix, Pmatrix)


! ----------------------------------------------------------------------
! SUBROUTINE: integr_VEQ
! ----------------------------------------------------------------------
! Purpose:
!  Variational Equations solution based on Runge-Kutta numerical integration methods
! ----------------------------------------------------------------------
! Input arguments:
! - to:			Initial Epoch: Modified Julian Day (MJD) in Terrestrial Time (including the fraction of the day)
! - ro: 		Satellite position vector (m) in ICRF
! - vo: 		Satellite velocity vector (m/sec) in ICRF
! - arc:		Orbit arc lenth (seconds)
! - integID: 	Numerical integration method ID number
! 				1. RKN7(6)8:	Runge-Kutta-Nystrom 7th order   
! 				2. RK4:			Runge-Kutta 4th order
! 				3. RK8(7)13:	Runge-Kutta 8th order
! - step: 		Numerical Integrator Stepsize (seconds)
! - Nparam:		Number of parameters to be estimated through the sensitivity matrix
!
! Output arguments:
! - orbc: 		Satellite orbit array in ICRF including the following per epoch:
!               - Modified Julian Day number in TT (including the fraction of the day) 
!				- Seconds since 00h in TT
!				- Position vector (m)
!				- Velocity vector (m/sec)
!- Smatrix:		State transition matrix (6*Epochs x 6)
!- Pmatrix:		Sensitivity matrix      (6*Epochs x Np)
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
!			Cooperative Research Centre for Spatial Information, Australia
! Created:	22 January 2018
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      !USE mdl_param
      USE m_veq_rkn768
      IMPLICIT NONE

	  
! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
! IN
      REAL (KIND = prec_d), INTENT(IN) :: MJDo, tsec_start
      REAL (KIND = prec_d), INTENT(IN), DIMENSION(3) :: ro, vo
      REAL (KIND = prec_d), INTENT(IN) :: arc
      INTEGER (KIND = prec_int2), INTENT(IN) :: integID
      REAL (KIND = prec_d), INTENT(IN) :: step
	  INTEGER (KIND = prec_int8), INTENT(IN) :: Nparam
! ----------------------------------------------------------------------
! OUT
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: orbc, Smatrix, Pmatrix  
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: Nepochs 
      INTEGER (KIND = prec_int8) :: i, j, iparam, i1, j1 
      INTEGER (KIND = prec_int2) :: AllocateStatus, DeAllocateStatus
      REAL (KIND = prec_d) :: lamda_h
      REAL (KIND = prec_d) :: to_sec, tmax, TT, TTo, t
      REAL (KIND = prec_d), DIMENSION(6) :: Zq
      REAL (KIND = prec_d), DIMENSION(3) :: er
      REAL (KIND = prec_d) :: mjd_to, mjd_th, mjdn
      REAL (KIND = prec_d), DIMENSION(3) :: rto, vto, rth, vth  
      REAL (KIND = prec_d), DIMENSION(8) :: Zo 
      REAL (KIND = prec_d), DIMENSION(6) :: yo, yn, ey 
! ----------------------------------------------------------------------	  
      REAL (KIND = prec_d), DIMENSION(6,6) :: veqZo, veqZ  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: veqPo, veqP 




! ----------------------------------------------------------------------
! Initial Conditions
! ----------------------------------------------------------------------
! Time conversion to seconds (MJD from days to seconds)

! Initial Epoch's Fraction of the day (in seconds)
to_sec = tsec_start !to_sec = (MJDo - INT(MJDo)) * (24.D0 * 3600.D0)

! Final Epoch
tmax = to_sec + arc

! Number of epochs
Nepochs = INT(arc / step) + 1


! ----------------------------------------------------------------------
! Dynamic allocatable arrays 
! ----------------------------------------------------------------------
ALLOCATE (orbc(Nepochs,8), STAT = AllocateStatus)

!ALLOCATE (Smatrix(Nepochs*6,6), STAT = AllocateStatus)
ALLOCATE (Smatrix(Nepochs,6*6+2), STAT = AllocateStatus)

!ALLOCATE (Pmatrix(Nepochs,Nparam*6), STAT = AllocateStatus)
ALLOCATE (Pmatrix(Nepochs,Nparam*6+2), STAT = AllocateStatus)

If (Nparam /= 0) Then
ALLOCATE (veqPo(6,Nparam), STAT = AllocateStatus)
Else
ALLOCATE (veqPo(6,1), STAT = AllocateStatus)
End IF
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Initial state vector into the orbit array
! ----------------------------------------------------------------------
orbc(1,1) = MJDo
orbc(1,2) = to_sec
orbc(1,3:5) = ro
orbc(1,6:8) = vo
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! VEQ matrices initial conditions
! ----------------------------------------------------------------------
! State Transition Matrix Initial values: veqZo = I (6x6)
veqZo(1,1:6) = (/ 1.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0 /)
veqZo(2,1:6) = (/ 0.D0, 1.D0, 0.D0, 0.D0, 0.D0, 0.D0 /)
veqZo(3,1:6) = (/ 0.D0, 0.D0, 1.D0, 0.D0, 0.D0, 0.D0 /)
veqZo(4,1:6) = (/ 0.D0, 0.D0, 0.D0, 1.D0, 0.D0, 0.D0 /)
veqZo(5,1:6) = (/ 0.D0, 0.D0, 0.D0, 0.D0, 1.D0, 0.D0 /)
veqZo(6,1:6) = (/ 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 1.D0 /)

! Sensitivity Matrix Initial values: veqPo = 0(6xNp)
Do i1 = 1 , 6
   Do j1 = 1 , Nparam
      veqPo(i1,j1) = 0.0D0
   End Do	  
End Do
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! VEQ matrices initial epoch values
! ----------------------------------------------------------------------
! State Transition Matrix
Smatrix(1,1) = MJDo
Smatrix(1,2) = to_sec
Smatrix(1,  3:8) = veqZo(1,1:6)
Smatrix(1, 9:14) = veqZo(2,1:6)
Smatrix(1,15:20) = veqZo(3,1:6)
Smatrix(1,21:26) = veqZo(4,1:6)
Smatrix(1,27:32) = veqZo(5,1:6)
Smatrix(1,33:38) = veqZo(6,1:6)

! Sensitivity matrix
Pmatrix(1,1) = MJDo
Pmatrix(1,2) = to_sec
!i = 1
!Pmatrix(1, 2+(i-1)*6+1 : 2+(i-1)*6+6 ) = (/ veqPo(1,i), veqPo(2,i), veqPo(3,i), veqPo(4,i), veqPo(5,i), veqPo(6,i) /)
!Do iparam = 1 , Nparam
!   Pmatrix(1, 2+(iparam-1)*6+1 : 2+(iparam-1)*6+6 ) = & 
!   (/ veqPo(1,iparam), veqPo(2,iparam), veqPo(3,iparam), veqPo(4,iparam), veqPo(5,iparam), veqPo(6,iparam) /)
!End Do
i1 = 0
iparam = 0	
Do i1 = 1 , 6 
Do iparam = 1 , Nparam
   Pmatrix(1, 2+(i1-1)*Nparam+iparam) = veqPo(i1,iparam) 
   !(/ veqP(1,iparam), veqP(2,iparam), veqP(3,iparam), veqP(4,iparam), veqP(5,iparam), veqP(6,iparam) /)
End Do
End Do
i1 = 0
!print *,"Pmatrix(1,:)", Pmatrix(1,:)
! ----------------------------------------------------------------------



! ----------------------------------------------------------------------
! Numerical integration methods
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Runge-Kutta-Nystrom RKN7(6)-8 method
! ----------------------------------------------------------------------
If (integID == 1) Then

! RKN7(6)-8 lamda parameter for stepsize control
lamda_h = 0.001D0

!i = 1
!Do t = to_sec , step , tmax-1

TT = to_sec
Do j = 1 , Nepochs-1
	
	! Epoch to	
    Zo(1) = orbc(j,1) 
    Zo(2) = orbc(j,2) 
    Zo(3:8) = orbc(j,3:8) 
	
    ! Numerical integration for the next epoch
    !Call integr_rkn768(Zo, step, lamda_h, Zq, er)
	Call veq_rkn768(Zo, veqZo, veqPo, step, lamda_h, Nparam, Zq, er, veqZ, veqP)	
    
	! Next epoch TT (to+h)
    TT = TT + step    

	! Seconds since 00h
	TTo = TT
	If (TT >= 86400.D0) Then
		TTo = TT - INT(TT / 86400.D0) * 86400.D0
	End IF 

	! MJD and Seconds
	orbc(j+1,1) = INT(MJDo) + TT / (24.0D0 * 3600.0D0) 
	orbc(j+1,2) = TTo
		
	! State vector at the next epoch TT (to+h) in the GCRS
	orbc(j+1,3:8) = Zq	
!	print *,"orbc t", orbc(j+1,1:2)
!	print *,"orbc r", orbc(j+1,3:5)

	! VEQ matrices at the next epoch TT (to+h)
	! State Transition Matrix
	Smatrix(j+1,1) = orbc(j+1,1)
	Smatrix(j+1,2) = orbc(j+1,2)
	Smatrix(j+1,  3:8) = veqZ(1,1:6)
	Smatrix(j+1, 9:14) = veqZ(2,1:6)
	Smatrix(j+1,15:20) = veqZ(3,1:6)
	Smatrix(j+1,21:26) = veqZ(4,1:6)
	Smatrix(j+1,27:32) = veqZ(5,1:6)
	Smatrix(j+1,33:38) = veqZ(6,1:6)

	! Sensitivity matrix
	Pmatrix(j+1,1) = orbc(j+1,1)
	Pmatrix(j+1,2) = orbc(j+1,2)
	!Do iparam = 1 , Nparam
	!   Pmatrix(j+1, 2+(iparam-1)*6+1 : 2+(iparam-1)*6+6 ) = & 
	!   (/ veqP(1,iparam), veqP(2,iparam), veqP(3,iparam), veqP(4,iparam), veqP(5,iparam), veqP(6,iparam) /)
	!End Do

	i1 = 0
	iparam = 0	
	Do i1 = 1 , 6 
	Do iparam = 1 , Nparam
	   Pmatrix(j+1, 2+(i1-1)*Nparam+iparam) = veqP(i1,iparam) 
	End Do
	End Do
!print *,"m_integrVEQf03: j,Nepochs",j,Nepochs
!print *,"m_integrVEQf03: Nparam ", Nparam
!print *,"m_integrVEQf03: veqP   ", veqP(1,:)
!print *,"m_integrVEQf03: Pmatrix", Pmatrix(j+1,3:Nparam+2)
!print *,"m_integrVEQf03: veqPo  ", veqPo
!print *,"m_integrVEQf03: veqP   ", veqP
!print *,"m_integrVEQf03: Pmatrix", Pmatrix(j+1,:)

    ! Initial VEQ arrays at next epoch numerical integration 
    veqZo = veqZ
    veqPo = veqP
End DO
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Runge-Kutta RK8(7)-13 method
! ----------------------------------------------------------------------
Else If (integID == 3) Then


TT = to_sec
Do j = 1 , Nepochs-1
	! Epoch to	
    mjd_to = orbc(j,1) 
    yo = orbc(j,3:8) 
	
    ! Numerical integration for computation of the state vector at next epoch
	!Call integr_rk87 (mjd_to, yo, step, mjdn, yn, ey)
   
	! Next epoch TT (to+h)
    TT = TT + step    

	! Seconds since 00h
	TTo = TT
	If (TT >= 86400.D0) Then
		TTo = TT - INT(TT / 86400.D0) * 86400.D0
	End IF 

	! MJD and Seconds
	orbc(j+1,1) = INT(MJDo) + TT / (24.0D0 * 3600.0D0) 
	orbc(j+1,2) = TTo
		
	! State vector at the next epoch TT (to+h) in the GCRS
	orbc(j+1,3:8) = yn
	
	! VEQ matrices at the next epoch TT (to+h)
	! State Transition Matrix
	Smatrix(j+1,1) = orbc(j+1,1)
	Smatrix(j+1,2) = orbc(j+1,2)
	Smatrix(j+1,  3:8) = veqZ(1,1:6)
	Smatrix(j+1, 9:14) = veqZ(2,1:6)
	Smatrix(j+1,15:20) = veqZ(3,1:6)
	Smatrix(j+1,21:26) = veqZ(4,1:6)
	Smatrix(j+1,27:32) = veqZ(5,1:6)
	Smatrix(j+1,33:38) = veqZ(6,1:6)

	! Sensitivity matrix
	Pmatrix(j+1,1) = orbc(j+1,1)
	Pmatrix(j+1,2) = orbc(j+1,2)
	!Do iparam = 1 , Nparam
	!   Pmatrix(j+1, 2+(iparam-1)*6+1 : 2+(iparam-1)*6+6 ) = & 
	!   (/ veqP(1,iparam), veqP(2,iparam), veqP(3,iparam), veqP(4,iparam), veqP(5,iparam), veqP(6,iparam) /)
	!End Do
	i1 = 0
	iparam = 0	
	Do i1 = 1 , 6 
	Do iparam = 1 , Nparam
	   Pmatrix(j+1, 2+(i1-1)*Nparam+iparam) = veqP(i1,iparam) 
	End Do
	End Do
End DO
! ----------------------------------------------------------------------

End If


End Subroutine



End Module



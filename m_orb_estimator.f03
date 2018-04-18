MODULE m_orb_estimator


! ----------------------------------------------------------------------
! MODULE: m_orb_estimator
! ----------------------------------------------------------------------
! Purpose:
!  Module for calling the following subroutines 
! 
! Subroutines contained within the module:
! - orb_estimator
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
!			Cooperative Research Centre for Spatial Information, Australia
! Created:	19 March 2018
! ----------------------------------------------------------------------


      IMPLICIT NONE
      !SAVE 			
  
	  
Contains


SUBROUTINE orb_estimator(orbref, veqZarray, veqParray, orbobs, Xmatrix, Wmatrix, Amatrix)

! ----------------------------------------------------------------------
! SUBROUTINE: orb_estimator
! ----------------------------------------------------------------------
! Purpose:
!  Orbit parameter estimation based on least-squares method and pseudo-observations
! ----------------------------------------------------------------------
! Input arguments:
! - obs_orb:	Pseudo-observations matrix with dimensions N1x5; obs_orb(ti)=[MJD_ti Sec_ti X Y Z]
! - orb_dyn:	Orbit computed from numerical integration solution (dimensions N2x8)
!				orb_dyn(ti)=[MJD_ti Sec_ti X Y Z Vx Vy Vz]
!
! Output arguments:
! - Xmatrix: 	Optimum estimated parameters corrections (initial state vector + additional parameters)
! ----------------------------------------------------------------------
! Dr. Thomas Papanikolaou, Geoscience Australia            19 March 2018
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      USE m_matrixinv
      IMPLICIT NONE

	  
! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
! IN
      REAL (KIND = prec_d), INTENT(IN), DIMENSION(:,:), ALLOCATABLE :: orbref, veqZarray, veqParray, orbobs
! ----------------------------------------------------------------------
! OUT
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: Xmatrix  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: Wmatrix  
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE, INTENT(OUT) :: Amatrix  
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: sz1, sz2, sz3, sz4, sz5, sz6, sz7, sz8 
      INTEGER (KIND = prec_int8) :: Norb, Nobs, Nparam, n_NEQn 
      INTEGER (KIND = prec_int8) :: Nepochs, Ncommon, iobs0, iobs, iref   
      INTEGER (KIND = prec_int2) :: Nobsset   
      REAL (KIND = prec_d) :: tiref, tiobs, dt_limit
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: AmatrixZ
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: AmatrixP
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: Amatrix_T
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: NEQn, NEQu, NEQn_inv
      REAL (KIND = prec_d), DIMENSION(:,:), ALLOCATABLE :: ObsEpochs
      INTEGER (KIND = prec_int2) :: AllocateStatus, DeAllocateStatus  
! ----------------------------------------------------------------------	  
      INTEGER (KIND = prec_int8) :: An
      !REAL (KIND = prec_q) :: Nmatrix(An,An)
      !REAL (KIND = prec_q) :: Nmatrix_inv(An,An)
      REAL (KIND = prec_q) :: Nmatrix(3,3)
      REAL (KIND = prec_q) :: Nmatrix_inv(3,3)


	  
! ----------------------------------------------------------------------
sz1 = size(orbobs, DIM = 1)
sz2 = size(orbobs, DIM = 2)
Nobs = sz1

sz3 = size(orbref, DIM = 1)
sz4 = size(orbref, DIM = 2)
Norb = sz3

sz5 = size(veqZarray, DIM = 1)
sz6 = size(veqZarray, DIM = 2)

sz7 = size(veqParray, DIM = 1)
sz8 = size(veqParray, DIM = 2)
Nparam = (sz8 - 2) / 6
! ----------------------------------------------------------------------

! Temp
Nparam = 0																	! 777 Temp

! ----------------------------------------------------------------------
! Number of common epochs
! Set the numerical differences digits limit for the common epochs
dt_limit = 1.D-08

Nepochs = 0
iref = 0
iobs = 0
iobs0 = 1
Do iref = 1 , Norb
    ! Reference Orbit
    tiref = orbref(iref,1)
    ! Pseduo-Observations (sp3 Orbit)
    Do iobs = iobs0 , Nobs
        tiobs = orbobs(iobs,1)
        ! Common Epochs
        If ( ABS(tiref - tiobs) < dt_limit) Then
            Nepochs = Nepochs + 1
            iobs0 = iobs + 1			
        end IF
    end Do   
end Do
Ncommon = Nepochs
! ----------------------------------------------------------------------
!print *,"Ncommon", Ncommon


! ----------------------------------------------------------------------
! Dynamic Allocatable arrays
ALLOCATE (Xmatrix(6+Nparam,1), STAT = AllocateStatus)
ALLOCATE (Wmatrix(3*Ncommon,1), STAT = AllocateStatus)

ALLOCATE (Amatrix(3*Ncommon,6+Nparam), STAT = AllocateStatus)
ALLOCATE (AmatrixZ(3*Ncommon,6), STAT = AllocateStatus)
ALLOCATE (AmatrixP(3*Ncommon,Nparam), STAT = AllocateStatus)
ALLOCATE (Amatrix_T(6+Nparam,3*Ncommon), STAT = AllocateStatus)

ALLOCATE (NEQn(6+Nparam,6+Nparam), STAT = AllocateStatus)
ALLOCATE (NEQu(6+Nparam,1), STAT = AllocateStatus)
ALLOCATE (NEQn_inv(6+Nparam,6+Nparam), STAT = AllocateStatus)

ALLOCATE (ObsEpochs(Ncommon,2), STAT = AllocateStatus)
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
iobs0 = 1
Nepochs = 0
iref = 0
iobs = 0
Do iref = 1 , Norb
    ! Reference Orbit
    tiref = orbref(iref,1)
    ! Pseduo-Observations (sp3 Orbit)
    Do iobs = iobs0 , Nobs
        tiobs = orbobs(iobs,1)
        ! Common Epochs
        If ( ABS(tiref - tiobs) < dt_limit) Then
            Nepochs = Nepochs + 1
            ObsEpochs(Nepochs,1) = orbref(iref,1)
            ObsEpochs(Nepochs,2) = orbref(iref,2)
! ----------------------------------------------------------------------
			Nobsset = 3
			! Wmatrix
            Wmatrix( (Nepochs-1)*3+1 : Nepochs*3 , 1) = orbobs(iobs,3:5) - orbref(iref,3:5)
            ! Design Matrix values obtained by VEQ arrays values
            ! veqZarray		
			! dX/dZo, dY/dZo, dZ/dZo
			AmatrixZ((Nepochs-1)*Nobsset + 1,:) = veqZarray(iref,3:8)
			AmatrixZ((Nepochs-1)*Nobsset + 2,:) = veqZarray(iref,9:14)
			AmatrixZ((Nepochs-1)*Nobsset + 3,:) = veqZarray(iref,15:20)	
! ----------------------------------------------------------------------
			! veqParray
!            AmatrixP((Nepochs-1)*Nobsset+1, : ) = ...
! ----------------------------------------------------------------------
            iobs0 = iobs + 1
        end IF
    end Do   
end Do

! ----------------------------------------------------------------------
! Design matrix
! Amatrix summarized
!Amatrix = [AmatrixZ AmatrixP]
Amatrix = AmatrixZ
! ----------------------------------------------------------------------
! Normal Equations
Amatrix_T = TRANSPOSE(Amatrix)
NEQn = MATMUL(Amatrix_T , Amatrix) 
NEQu = MATMUL(Amatrix_T , Wmatrix) 

! N matrix inversion
!NEQn_inv = inv(NEQn)
n_NEQn = size(Amatrix, DIM = 2)
Call matrixinv (NEQn, NEQn_inv, n_NEQn)
!Print *,"NEQn", NEQn
!Print *,"NEQn_inv", NEQn_inv

! NEQ matrices test
! NEQ : I test

! ----------------------------------------------------------------------
! Parameters Estimation
Xmatrix = MATMUL(NEQn_inv , NEQu) 
! ----------------------------------------------------------------------


END SUBROUTINE


END Module


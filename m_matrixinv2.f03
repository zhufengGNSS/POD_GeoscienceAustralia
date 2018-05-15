MODULE m_matrixinv


! ----------------------------------------------------------------------
! MODULE: m_matrixinv
! ----------------------------------------------------------------------
! Purpose:
!  Module for calling the following subroutines 
! 
! Subroutines contained within the module:
! - matrixinv
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
!			Cooperative Research Centre for Spatial Information, Australia
! Created:	22 March 2018
! ----------------------------------------------------------------------


      IMPLICIT NONE
      !SAVE 			
  
	  
Contains


SUBROUTINE matrixinv (Amatrix, AmatrixInv)


! ----------------------------------------------------------------------
! SUBROUTINE: matrixinv
! ----------------------------------------------------------------------
! Purpose:
!  Matrix inversion based on LAPACK library 
! ----------------------------------------------------------------------
! Input arguments:
! - R1 : Matrix with dimensions ixj
!
! Output arguments:
! - R3 : Inverse of input matrix R1
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou
! 			Geoscience Australia
!			Cooperative Research Centre for Spatial Information, Australia
! Created:	22 March 2018
! ----------------------------------------------------------------------


      USE mdl_precision
      IMPLICIT NONE

! ---------------------------------------------------------------------------
! Dummy arguments declaration
! ---------------------------------------------------------------------------
! IN
      REAL (KIND = prec_q), INTENT(IN), DIMENSION(:,:), ALLOCATABLE :: Amatrix
! OUT
      REAL (KIND = prec_q), INTENT(OUT), DIMENSION(:,:), ALLOCATABLE :: AmatrixInv
! ---------------------------------------------------------------------------

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: n, m
      INTEGER (KIND = prec_int8) :: info
	  Integer, Dimension(:), ALLOCATABLE :: ipiv   
      REAL (KIND = prec_q), DIMENSION(:), ALLOCATABLE :: work_matrix
      INTEGER (KIND = prec_int2) :: AllocateStatus, DeAllocateStatus	  
! ----------------------------------------------------------------------
	  REAL (KIND = prec_q), dimension(:,:) :: A
  real(dp), dimension(size(A,1),size(A,2)) :: Ainv

  real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
  integer, dimension(size(A,1)) :: ipiv   ! pivot indices
  integer :: n, info
! ----------------------------------------------------------------------
  

! External procedures in LAPACK library
external DGETRF
external DGETRI


n = size(Amatrix, DIM = 1)
m = size(Amatrix, DIM = 2)

ALLOCATE (AmatrixInv(n,n), STAT = AllocateStatus)


ALLOCATE (work_matrix(n), STAT = AllocateStatus)
ALLOCATE (ipiv(n), STAT = AllocateStatus)


! LU factorization using partial pivoting with row interchanges
call DGETRF(n, n, AmatrixInv, n, ipiv, info)

if (info /= 0) then
   stop 'Matrix is numerically singular!'
end if

! DGETRI computes the inverse of an LU-factored general matrix computed by DGETRF.
call DGETRI(n, AmatrixInv, n, ipiv, work_matrix, n, info)

if (info /= 0) then
   stop 'Matrix inversion failed!'
end if



END SUBROUTINE


END Module


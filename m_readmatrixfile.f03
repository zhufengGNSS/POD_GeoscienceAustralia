MODULE m_readmatrixfile


! ----------------------------------------------------------------------
! MODULE: m_readmatrixfile.f03
! ----------------------------------------------------------------------
! Purpose:
!  Module for the subroutine readmatrixfile that read a file 
! 
! Subroutines contained within this module:
! - readmatrixfile.f03 
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou, Geoscience Australia
! Created:	6 February 2019
! ----------------------------------------------------------------------


      IMPLICIT NONE
      !SAVE 			
  
	  
Contains


   

SUBROUTINE readmatrixfile (filename, matrix)

! ----------------------------------------------------------------------
! SUBROUTINE: readmatrixfile
! ----------------------------------------------------------------------
! Purpose:
!  Read gravity field model's .gfc data file in the ICGEM (International
!  Centre for Global Earth Models) format
! ----------------------------------------------------------------------
! Input arguments:
! - filename:		Data file name that includes only columns with numerical values 
! 
! Output arguments:
! - matrix:        	Dynamic allocatable array with the numerical values of the file being read
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou, Geoscience Australia
! Created:	6 February 2019
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      IMPLICIT NONE
	  
! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
! IN
      CHARACTER (LEN=100), INTENT(IN) :: filename
! OUT
      REAL (KIND = prec_q), INTENT(OUT), DIMENSION(:,:), ALLOCATABLE :: matrix
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: i, read_i
      INTEGER (KIND = prec_int2) :: UNIT_IN, ios
      INTEGER (KIND = prec_int2) :: ios_line, ios_key, ios_data
      INTEGER (KIND = prec_int2) :: space_i
      INTEGER (KIND = prec_int2) :: AllocateStatus
      CHARACTER (LEN=7) :: Format1, Format2, Format3, Format4
      CHARACTER (LEN=500) :: line_ith	  
      CHARACTER (LEN=150) :: word1_ln, word_i	  
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: n_matrix, m_matrix
      REAL (KIND = prec_q) :: Omikron
      REAL (KIND = prec_q) :: num_line(5), num_i
      INTEGER (KIND = prec_int8) :: j, j_line, kapa


 !m_matrix = 5
 
! ----------------------------------------------------------------------
      UNIT_IN = 7  												
      Format1 = '(A)'
      Format2 = '(F)'
      Format3 = '(I100)'
! ----------------------------------------------------------------------
      Format4 = '(F25.15)' !'(F)'


! ----------------------------------------------------------------------
! Open file
      OPEN (UNIT = UNIT_IN, FILE = TRIM (filename), IOSTAT = ios)
      IF (ios /= 0) THEN
         PRINT *, "Error in opening file:", filename
         PRINT *, "OPEN IOSTAT=", ios
      END IF
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Read dimensions of the file
! ----------------------------------------------------------------------
      i = 0
      DO
	    i = i + 1
! ----------------------------------------------------------------------
		If (i > 1) Then
		
		READ (UNIT=UNIT_IN,FMT=Format1,IOSTAT=ios_line) line_ith
  	     !PRINT *, "READ Line (i,ios):", i, ios_line, line_ith 

! End of file
         IF (ios_line < 0) THEN
!            PRINT *, "End of file, i=", i
            EXIT		
         END IF
		 
! ----------------------------------------------------------------------
		Else If (i == 1) Then ! if (1<) 0 Then
		
		j = 0
		Do 					
			READ (UNIT=UNIT_IN, FMT=Format4, advance='no', IOSTAT=ios_data) num_i			
			j = j + 1
			print *,"j, num_i", j, num_i, ios_data
			CALL SLEEP(1)
			! End of line
			IF (ios_data .NE. 0) THEN
				PRINT *, "End of line, i,j:", i,j
				EXIT		
			END IF
		End Do
		READ(UNIT=UNIT_IN, FMT=Format4)
		print *,"j, num_i", j, num_i
		m_matrix = j
		
		End IF
! ----------------------------------------------------------------------
		
      END DO
! Close File
      CLOSE (UNIT=UNIT_IN)
! ----------------------------------------------------------------------
n_matrix = i
PRINT *, "n_matrix = ", n_matrix, i		 
PRINT *, "m_matrix = ", m_matrix, j		 



! ----------------------------------------------------------------------
! Allocatable arrays
! ----------------------------------------------------------------------
      ALLOCATE (matrix(n_matrix,m_matrix), STAT = AllocateStatus)
      IF (AllocateStatus /= 0) THEN
         PRINT *, "Error: Not enough memory"
         PRINT *, "Error: m_readmatrixfile.f03"
         PRINT *, "Error: Allocatable Array: matrix"
!         STOP "*** Not enough memory ***"
      END IF   
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Read file for storing the numerical values to the dynamic allocatable array
! ----------------------------------------------------------------------
 ! ----------------------------------------------------------------------
! Open file
      OPEN (UNIT = UNIT_IN, FILE = TRIM (filename), IOSTAT = ios)
      IF (ios /= 0) THEN
         PRINT *, "Error in opening file:", filename
         PRINT *, "OPEN IOSTAT=", ios
      END IF
! ----------------------------------------------------------------------
     i = 0
      DO
	     READ (UNIT=UNIT_IN,FMT=Format1,IOSTAT=ios_line) line_ith
	     i = i + 1
!  	     PRINT *, "READ Line (i,ios):", i, ios_line
! ----------------------------------------------------------------------
! End of file
         IF (ios_line < 0) THEN
!            PRINT *, "End of file, i=", i
            EXIT		
         END IF
! ----------------------------------------------------------------------

		!Do j_line = 1 , m_matrix
  	     !READ (line_ith, FMT= Format4, advance='no' , IOSTAT=ios_data) Omikron
		 !matrix(i,j_line) = Omikron
		!End Do

		Do j_line = 1 , m_matrix					
			READ (UNIT=UNIT_IN, FMT=Format4, advance='no', IOSTAT=ios_data) Omikron			
			matrix(i,j_line) = Omikron
			!print *,"j, num_i", j, num_i, ios_data
			!CALL SLEEP(1)
			! End of line
			IF (ios_data .NE. 0) THEN
				PRINT *, "End of line, i,j:", i,j
				EXIT		
			END IF
		End Do
		READ(UNIT=UNIT_IN, FMT=Format4)
		

	END DO
! ----------------------------------------------------------------------

      CLOSE (UNIT=UNIT_IN)
	  
END SUBROUTINE



End Module


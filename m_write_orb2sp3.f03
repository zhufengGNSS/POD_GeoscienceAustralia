MODULE m_write_orb2sp3


! ----------------------------------------------------------------------
! MODULE: m_write_orb2sp3.f03
! ----------------------------------------------------------------------
! Purpose:
!  Module for write the POD *.orb format to sp3 orbit format 
! ----------------------------------------------------------------------
! Author :	Dr. Thomas Papanikolaou, Geoscience Australia
! Created:	21 February 2019
! ----------------------------------------------------------------------


      IMPLICIT NONE
      !SAVE 			

	  
Contains


SUBROUTINE write_orb2sp3 (ORBmatrix, PRNmatrix, filename, sat_vel)

! ----------------------------------------------------------------------
! SUBROUTINE: writesp3_hd 
! ----------------------------------------------------------------------
! Purpose:
!  Write orbit sp3 format header 
! ----------------------------------------------------------------------
! Input arguments:
! - ORBmatrix:       Input allocatable array
! - filename:       file name to be used for writing array data
!
! Output arguments:
!
! ----------------------------------------------------------------------
! Dr. Thomas Papanikolaou, Geoscience Australia         21 February 2019
! ----------------------------------------------------------------------


      USE mdl_precision
      USE mdl_num
      IMPLICIT NONE
	  
! ----------------------------------------------------------------------
! Dummy arguments declaration
! ----------------------------------------------------------------------
! IN
      REAL (KIND = prec_q), INTENT(IN), DIMENSION(:,:,:), ALLOCATABLE :: ORBmatrix 
	  CHARACTER (LEN=3), ALLOCATABLE, INTENT(IN) :: PRNmatrix(:)
      CHARACTER (LEN=100), INTENT(IN) :: filename
      !CHARACTER (LEN=3), INTENT(IN) :: sat_prn
      INTEGER (KIND = prec_int2), INTENT(IN) :: sat_vel	  
! OUT
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------
! Local variables declaration
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: i, i_write
      INTEGER (KIND = prec_int2) :: UNIT_IN, ios, ios_ith
      INTEGER (KIND = prec_int8) :: sz1, sz2, sz3
      INTEGER (KIND = prec_int2) :: wrt_opt
      INTEGER (KIND = prec_int2) :: FMT_opt
! ----------------------------------------------------------------------
      CHARACTER (LEN=1) :: RealT
      INTEGER (KIND = prec_int2) :: RealW, RealD
      CHARACTER (LEN=70) :: fmt_wrt, fmt_wrt0, fmt_sz2
      REAL (KIND = prec_q) :: wrtArrayLN 
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int2) :: vel
      CHARACTER (LEN=100) :: fmt_epoch, fmt_pos, fmt_vel

      REAL (KIND = prec_q) :: MJD_ti, Sec_00
      CHARACTER (LEN=3) :: PRN_ti
      REAL (KIND = prec_q) :: r_ti(3), v_ti(3), cl_ti, Dcl_ti

      DOUBLE PRECISION DJ1, DJ2
      INTEGER IY, IM, ID
      DOUBLE PRECISION FD
      INTEGER J

      !INTEGER (KIND = prec_int4) :: year, month, day, hour_ti, min_ti 
      INTEGER :: year, month, day, hour_ti, min_ti 
      REAL (KIND = prec_q) :: sec_ti 	  
! ----------------------------------------------------------------------
      INTEGER (KIND = prec_int8) :: Nepochs, Nelem, Nsat
      INTEGER (KIND = prec_int8) :: i_sat
! ----------------------------------------------------------------------


UNIT_IN = 7  												


! ----------------------------------------------------------------------
! Format definition
! ----------------------------------------------------------------------
!*  2010 12 25  0  0  0.0000 
!PG01  -6582.270015  18452.592641 -17946.851511    -67.214659  6  9  8  87       
fmt_epoch = '(A3,I5.4,4I3.2,F11.8)'
!fmt_pos = '(A1,A1,I2.2,4F14.6,A)'
!fmt_vel = '(A1,A1,I2.2,4F14.6,A)'
fmt_pos = '(A1,A3,A1,4F14.6)'
fmt_vel = '(A1,A3,A1,4F14.6)'
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Orbit array dimensions
sz1 = SIZE (ORBmatrix,DIM=1)
sz2 = SIZE (ORBmatrix,DIM=2)
sz3 = SIZE (ORBmatrix,DIM=3)
Nepochs = sz1
Nelem   = sz2
Nsat    = sz3
! ----------------------------------------------------------------------

 

! ----------------------------------------------------------------------
! Writing satellite velocity vector (optional)
vel = sat_vel
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Open file
      OPEN (UNIT=UNIT_IN,FILE=filename,ACTION="WRITE",POSITION="REWIND", IOSTAT=ios)
      IF (ios /= 0) THEN
         PRINT *, "Error in opening file:", filename
         PRINT *, "OPEN IOSTAT=", ios
      END IF
! ----------------------------------------------------------------------


! ----------------------------------------------------------------------
! Write data to file | Line by line	  
! ----------------------------------------------------------------------
! Epochs loop
DO i_write = 1 , Nepochs       


! Satellites loop
DO i_sat = 1 , Nsat       


! MJD of epoch (including fraction of the day)
MJD_ti = ORBmatrix(i_write,1, i_sat)
! Sec of Day (since 0h)
Sec_00 = ORBmatrix(i_write,2, i_sat)

! Satellite PRN number 
!PRN_ti = filename(1:3)
!PRN_ti = sat_prn
PRN_ti = PRNmatrix(i_sat)

! Position vector
r_ti(1) = ORBmatrix(i_write,3, i_sat) !/ 1.D03
r_ti(2) = ORBmatrix(i_write,4, i_sat) !/ 1.D03 
r_ti(3) = ORBmatrix(i_write,5, i_sat) !/ 1.D03
! Unit conversion: m to Km                 
r_ti = r_ti * 1.0D-3

cl_ti = 999999.999999D0

if (vel>0) then
! Velocity vector in dm/sec
v_ti(1) = ORBmatrix(i_write,6, i_sat) 
v_ti(2) = ORBmatrix(i_write,7, i_sat) 
v_ti(3) = ORBmatrix(i_write,8, i_sat) 
! Unit conversion: m/sec to dm/sec                 
v_ti = v_ti * 1.0D1
end if


Dcl_ti = 999999.999999D0


! ----------------------------------------------------------------------
IF (i_sat == 1) THEN
! Write the Epoch line

! MJD of Epoch (including fraction of the day)
DJ1 = 2400000.5D0
DJ2 = MJD_ti
Call iau_JD2CAL ( DJ1, DJ2, IY, IM, ID, FD, J )

year  = IY
month = IM
day   = ID

if (1<0) then
hour_ti = INT(FD * 24.0D0) 
min_ti  = INT(FD * 24.0D0 * 60.0D0)
sec_ti  = (FD * 24.0D0 * 60.0D0 * 60.0D0)
else
hour_ti = INT(Sec_00 / (60.0D0 * 60.0D0))
min_ti  = INT(Sec_00/60.0D0 - hour_ti*60.0D0)  
sec_ti  = (Sec_00 - hour_ti*3600.0D0 - min_ti*60.D0)
end if

!print *,"sec_ti print", sec_ti 
!WRITE (*,FMT='(A6,F17.6)'),"sec_ti", sec_ti
!WRITE (*,FMT='(A6,F17.6)'),"sec_ti", sec_ti

! Epoch line
!*  2010 12 25  0  0  0.0000 
!READ (line_ith, * , IOSTAT=ios_data) char3, year, month, day, hr, minute, sec  
WRITE (UNIT=UNIT_IN,FMT=fmt_epoch,IOSTAT=ios_ith) '*  ', year, month, day, hour_ti, min_ti, sec_ti
!WRITE (UNIT=UNIT_IN,FMT=*,IOSTAT=ios_ith) '*  ', year,' ', month,' ',day,' ',hour_ti,' ',min_ti,' ',sec_ti,' '

END IF
! ----------------------------------------------------------------------


if (1>0) then
! ----------------------------------------------------------------------
! Write Satellite Position and Velocity Vector per PRN 	  
! ----------------------------------------------------------------------
! Position per epoch	  
WRITE (UNIT=UNIT_IN,FMT=fmt_pos,IOSTAT=ios_ith) 'P',PRN_ti,' ', r_ti, cl_ti

! Velocity per epoch	  
if (vel>0) then
WRITE (UNIT=UNIT_IN,FMT=fmt_vel,IOSTAT=ios_ith) 'V',PRN_ti,' ', v_ti, Dcl_ti
end if
! ----------------------------------------------------------------------
end if


IF (ios_ith /= 0) THEN
    PRINT *, "Error in writing to file: ", TRIM (filename)
    PRINT *, "WRITE IOSTAT=", ios_ith
END IF
! ----------------------------------------------------------------------

 
 END DO
 END DO
! ----------------------------------------------------------------------


ENDFILE (UNIT = UNIT_IN) 
CLOSE (UNIT = UNIT_IN)
! ----------------------------------------------------------------------


END SUBROUTINE


End Module

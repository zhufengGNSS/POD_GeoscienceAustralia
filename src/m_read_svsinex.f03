MODULE m_read_svsinex

! ----------------------------------------------------------------------
! MODULE: m_read_svsinex.f03
! ----------------------------------------------------------------------
! Purpose:
!
!  Module for reading IGS metadat SINEX file and stroing info as global
!  variables 
! 
! ----------------------------------------------------------------------
! Author :	Tzupang Tseng, Geoscience Australia, Australia
! Created:	27-09-2019
! ----------------------------------------------------------------------

      IMPLICIT NONE
!      SAVE		
  	  
Contains
  
SUBROUTINE read_svsinex (UNIT_IN,idir,iyr,iday,ihr,imin,gnss,isat, &
                         SVNID,BLKTYP,BLKID,MASS,POWER)

! ----------------------------------------------------------------------
! SUBROUTINE: read_svsinex.f03
! ----------------------------------------------------------------------
! Purpose:
!  Read and store satellite metadata information
! ----------------------------------------------------------------------
! Input arguments:
!
!  UNIT     : logical unit number for read              
!  idir    : idir = 1 SVN to PRN, idir = -1 PRN to SVN  
!  iyr     : 4-digit year                               
!  iday    : 3-digit day of year                        
!  ihr     : 2-digit hr                                 
!  imin    : 2-digit minute                             
!  gnss    : 1-character GNSS code (G R E C J I)        
!  isat    : input SVN or PRN number                    
! 
! Ouptut arguments:
!
!  SVNID   : SVN NUMBER                    
!  BLKTYP  : BLOCK TYPE   
!  MASS    : S/C mass in kg 
!  POWER   : transmitted power in watts                 
!
! ----------------------------------------------------------------------
! Remarks:
!  
      USE mdl_precision
      USE mdl_num
      IMPLICIT NONE

      INTEGER (KIND = prec_int4) :: UNIT_IN
      INTEGER (KIND = prec_int2) :: SVNID, BLKID
      CHARACTER(LEN=1) :: gnss, gnss_tmp
      CHARACTER(LEN=3) :: prn, prn_in 
      CHARACTER(LEN=4) :: svn, svn_in 
      CHARACTER(LEN=10):: cospar_id   
      CHARACTER(LEN=6) :: SatCat 
      CHARACTER(LEN=20):: antbody_in  ! Body type read from metadata snx

      INTEGER(KIND = prec_int2):: idir     ! Direction of PRN->SVN or visa-vers
      INTEGER(KIND = prec_int2):: yr1,  yr2   ! Year read from sinex file
      INTEGER(KIND = prec_int2):: doy1, doy2  ! DOY read from sinex file

      REAL(KIND = prec_d) ::  sod1, sod2  ! Seconds of day read from sinex file

      CHARACTER(LEN=20):: BLKTYP
      CHARACTER(LEN=80):: prog_name
      CHARACTER(LEN=128):: record
      CHARACTER(LEN=256):: message

      INTEGER(KIND = prec_int4):: isat,satid,frqchn,POWER,iyr,iday,ihr,imin,rcpar,ioerr,i

      REAL(KIND = prec_d) :: time_id, time1, time2
      REAL(KIND = prec_d) :: MASS

      logical found  ! Used to indicate that blocks and SVNs have been found.
      logical debug / .false. /   ! Turn on ouput.

! Put the requested time into a scale of year 
time_id = iyr + iday/365.d0 + (ihr*3600.d0+imin*60.d0)/86400.d0/365.d0

! Initialize frquency channel (GLONASS only)
frqchn = 0

! Which conversion is going to be implemented? 
IF(idir==-1) THEN  ! PRN to SVN
   WRITE(prn_in,'(a1,I2.2)') gnss, isat
ELSE                   ! SVN to PRN
   WRITE(svn_in,'(a1,I3.3)') gnss, isat
END IF

! Start to read the SINEX file
found = .false.
REWIND(UNIT_IN)
DO WHILE (.not.found )
   READ(UNIT_IN,'(a)',IOSTAT=ioerr) record
   IF(record(1:3)=='-SAT') THEN 
     print*,'NO SVN CAn BE FOUND ',ioerr
   ELSE IF (record(1:1)==' ' ) THEN
     READ(record,'(1x,a4,2(1x,i4,1x,i3,1x,f5.0),1x,a3)',IOSTAT=ioerr) svn,yr1,doy1,sod1,yr2,doy2,sod2,prn
     IF(yr2==0000) THEN
     yr2   = 2100
     doy2  = 365
     sod2  = 86400
     time2 = yr2 + doy2/365.d0 + sod2/86400.d0/365.d0
     END IF
     time1 = yr1 + doy1/365.d0 + sod1/86400.d0/365.d0
     time2 = yr2 + doy2/365.d0 + sod2/86400.d0/365.d0
     IF(idir==-1) THEN
        IF(prn==prn_in) THEN
           IF(time_id>=time1 .and. time_id<=time2) THEN
              found = .true.
              READ(svn(2:4),'(i3)') satid
              svn_in = svn   
           END IF
        END IF
     ELSE 
        IF( idir==1 ) THEN
        IF(svn==svn_in) THEN
           IF(time_id>=time1 .and. time_id<=time2) THEN
              found = .true.
              READ(prn(2:3),'(i2)') satid
              prn_in = prn    
           END IF
        END IF
        END IF
      END IF
   END IF
END DO


! This only needed if gnss == R for Glonass.
IF(gnss(1:1)=='R') THEN
REWIND(UNIT_IN)
found = .false.
   DO WHILE(.not.found)
   READ(UNIT_IN,'(a)',IOSTAT=ioerr) record
      IF(ioerr/=0 ) THEN
         print*,'Failed to find FREQUENCY_CHANNEL SINEX block',ioerr
      ELSEIF(record(12:15)=='FREQ') THEN
         found = .true.
      END IF
   END DO 
   found = .false.
   DO WHILE (.not.found)
   READ(UNIT_IN,'(a)',IOSTAT=ioerr) record
      IF(record(1:3)=='-SAT'.or.ioerr.ne.0) THEN
         If( ioerr/=0 ) print*,'Failed to find SV in FREQUENCY_CHANNEL SINEX block',ioerr
      ELSE 
         IF(record(1:1)==' ') THEN
         READ(record,'(1x,a4,2(1x,i4,1x,i3,1x,f5.0),1x,i4)',IOSTAT=ioerr) svn,yr1,doy1,sod1,yr2,doy2,sod2,frqchn
         IF(ioerr/=0) print*,'Reading GLONASS frequency',ioerr
         IF(yr2==0000) THEN
         yr2  = 2100
         doy2 = 365
         sod2 = 86400
         time2 = yr2 + doy2/365.d0 + sod2/86400.d0/365.d0
         END IF
         time1 = yr1 + doy1/365.d0 + sod1/86400.d0/365.d0
         time2 = yr2 + doy2/365.d0 + sod2/86400.d0/365.d0
         IF(svn==svn_in) THEN
            IF(time_id>=time1 .and. time_id<=time2) THEN
            found = .true.
            END IF
         END IF
         END IF
      END IF
   END DO 
END IF         ! gnss = R

! Get the SATELLITE/IDENTIFIER block

REWIND(UNIT_IN)
found = .false.
DO WHILE (.not.found)
READ(UNIT_IN,'(a)',iostat=ioerr) record
   IF(ioerr/=0 ) THEN
      print*,'Failed to find IDENTIFIER SINEX block',ioerr
   ELSEIF ( record(12:21)=='IDENTIFIER' ) THEN
      found = .true.
   END IF 
END DO 
found = .false.
DO WHILE (.not.found)
READ(UNIT_IN,'(a)',iostat=ioerr) record
   IF (record(1:3)=='-SAT'.or.ioerr/=0) THEN
      If(ioerr/=0) print*,'Failed to find SV in IDENTIFIER SINEX block',ioerr
   ELSE 
      IF(record(1:1)==' ') THEN
      READ(record,'(1x,a4,1x,a9,1x,A6,1x,a15)',iostat=ioerr)svn,cospar_id,SatCat,antbody_in
      IF(ioerr/=0) print*,'Fail to reading satellite id ',ioerr
         IF(svn==svn_in) THEN
            found = .true.
          BLKTYP = antbody_in
         END IF 
      END IF 
   END IF 
END DO 

! Prepare SVNID and BLKID for the global variables
! ------------------------------------------------
SVNID = satid
IF(BLKTYP=='GPS-I')      BLKID = 1
IF(BLKTYP=='GPS-II')     BLKID = 2
IF(BLKTYP=='GPS-IIA')    BLKID = 3
IF(BLKTYP=='GPS-IIR')    BLKID = 4
IF(BLKTYP=='GPS-IIR-A')  BLKID = 5
IF(BLKTYP=='GPS-IIR-B')  BLKID = 6
IF(BLKTYP=='GPS-IIR-M')  BLKID = 7
IF(BLKTYP=='GPS-IIF')    BLKID = 8
IF(BLKTYP=='GPS-IIIA')   BLKID = 9
IF(BLKTYP=='GLO')        BLKID = 101
IF(BLKTYP=='GLO-M'  .or.BLKTYP == 'GLO-M+')  BLKID = 102
IF(BLKTYP=='GLO-K1A'.or.BLKTYP == 'GLO-K1B') BLKID = 103
IF(BLKTYP=='GAL-1')     BLKID = 201 ! Galileo (IOV)
IF(BLKTYP=='GAL-2')     BLKID = 202 ! Galileo (FOC)
IF(BLKTYP=='BDS-2G'.or.BLKTYP == 'BDS-3G')            BLKID = 301 ! BDS GEO
IF(BLKTYP=='BDS-2I'.or.BLKTYP == 'BDS-3I'.or.&
   BLKTYP=='BDS-3SI-SECM'.or.BLKTYP =='BDS-3SI-CAST') BLKID = 302 ! BDS IGSO
IF(BLKTYP=='BDS-2M'.or.BLKTYP == 'BDS-3M'.or.&
   BLKTYP=='BDS-3M-SECM'.or.BLKTYP =='BDS-3M-CAST')   BLKID = 303 ! BDS MEO


! Get the satellite mass ( SATELLITE/MASS block )
REWIND(UNIT_IN)
found = .false.
DO WHILE (.not.found)
   READ(UNIT_IN,'(a)',IOSTAT=ioerr) record
   IF(ioerr/=0 ) THEN
   print*,'Failed to find MASS SINEX block',ioerr
   ELSEIF( record(12:15)=='MASS' ) THEN
   found = .true.
   END IF 
END DO 
found = .false.
DO WHILE (.not.found)
   READ(UNIT_IN,'(a)',IOSTAT=ioerr) record
   IF(record(1:3)=='-SAT'.or.ioerr/=0 ) THEN
     IF( ioerr/=0 ) print*,'Failed to find SV in MASS SINEX block',ioerr
   ELSE
! Only try to decode record if it is not a comment
   IF(record(1:1)==' ' ) THEN
      READ(record,'(1x,a4,2(1x,i4,1x,i3,1x,f5.0),1x,f9.3)',IOSTAT=ioerr) svn,yr1,doy1,sod1,yr2,doy2,sod2,MASS
      IF(ioerr/=0) print*,'Reading satellite mass ',ioerr
      IF(yr2==0000) THEN
         yr2 = 2100
         doy2 = 365
         sod2 = 86400
         time2 = yr2+doy2/365.d0+sod2/86400.d0/365.d0
      END IF
         time1 = yr1+doy1/365.d0+sod1/86400.d0/365.d0
         time2 = yr2+doy2/365.d0+sod2/86400.d0/365.d0
      IF(svn==svn_in) THEN
        IF(time_id>=time1 .and. time_id<=time2) THEN
        MASS = MASS
        found = .true.
        END IF 
      END IF
    END IF
    END IF 
END DO 

! Get the transmitter power ( SATELLITE/TX_POWER block)

REWIND(UNIT_IN)
found = .false.
DO WHILE (.not.found)
   READ(UNIT_IN,'(a)',IOSTAT=ioerr) record
   IF(ioerr/=0 ) THEN
   print*,'svnav.dat','Failed to find TX_POWER SINEX block',ioerr
   ELSEIF( record(12:15)=='TX_P' ) THEN
   found = .true.
   END IF
END DO
found = .false.
DO WHILE (.not.found)
   READ(UNIT_IN,'(a)',IOSTAT=ioerr) record
   IF(record(1:3)=='-SAT'.or.ioerr/=0) THEN
   IF( ioerr/=0 ) then
   print*,'Failed to find SV in TX_POWER block and the POWER is set to 185 W'
   POWER=185 ! assumed to be consistent with BDS IGSO
   RETURN
   END IF
   ELSE IF(record(1:1)==' ') THEN
      READ(record,'(1x,a4,2(1x,i4,1x,i3,1x,F5.0),1x,i4)',IOSTAT=ioerr) svn,yr1,doy1,sod1,yr2,doy2,sod2,POWER
      IF(yr2==0000) THEN
         yr2 = 2100
        doy2 = 365
        sod2 = 86400
       time2 = yr2 + doy2/365.d0 + sod2/86400.d0/365.d0
      END IF 
      IF(ioerr/=0) print*,'Check metadata file'
       time1 = yr1 + doy1/365.d0 + sod1/86400.d0/365.d0
       time2 = yr2 + doy2/365.d0 + sod2/86400.d0/365.d0
      IF(svn==svn_in) THEN
         IF(time_id>=time1 .and. time_id<=time2) THEN
         found = .true.
         END IF
      END IF
    END IF
END DO

! Prepare SVNID and BLKID for the global variables
! ------------------------------------------------
!SVNID = satid
!IF(BLKTYP=='GPS-I')      BLKID = 1
!IF(BLKTYP=='GPS-II')     BLKID = 2
!IF(BLKTYP=='GPS-IIA')    BLKID = 3
!IF(BLKTYP=='GPS-IIR')    BLKID = 4
!IF(BLKTYP=='GPS-IIR-A')  BLKID = 5
!IF(BLKTYP=='GPS-IIR-B')  BLKID = 6
!IF(BLKTYP=='GPS-IIR-M')  BLKID = 7
!IF(BLKTYP=='GPS-IIF')    BLKID = 8
!IF(BLKTYP=='GPS-IIIA')   BLKID = 9
!IF(BLKTYP=='GLO')        BLKID = 101
!IF(BLKTYP=='GLO-M'  .or.BLKTYP == 'GLO-M+')  BLKID = 102
!IF(BLKTYP=='GLO-K1A'.or.BLKTYP == 'GLO-K1B') BLKID = 103
!IF(BLKTYP=='GLA-1')     BLKID = 201 ! Galileo (IOV)
!IF(BLKTYP=='GLA-2')     BLKID = 202 ! Galileo (FOC)
!IF(BLKTYP=='BDS-2G'.or.BLKTYP == 'BDS-3G')            BLKID = 301 ! BDS GEO
!IF(BLKTYP=='BDS-2I'.or.BLKTYP == 'BDS-3I'.or.&
!   BLKTYP=='BDS-3SI-SECM'.or.BLKTYP =='BDS-3SI-CAST') BLKID = 302 ! BDS IGSO
!IF(BLKTYP=='BDS-2M'.or.BLKTYP == 'BDS-3M'.or.&
!   BLKTYP=='BDS-3M-SECM'.or.BLKTYP =='BDS-3M-CAST')   BLKID = 303 ! BDS MEO

END SUBROUTINE

END

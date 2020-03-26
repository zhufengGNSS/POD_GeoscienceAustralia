MODULE m_read_svsinex

! ----------------------------------------------------------------------
! MODULE: m_read_svsinex.f03
! ----------------------------------------------------------------------
! Purpose:
!
!  Module for reading IGS metadat SINEX file and storing info as global
!  variables 
! 
! ----------------------------------------------------------------------
! Author :	Tzupang Tseng, Geoscience Australia, Australia
! Created:	27-09-2019
! ----------------------------------------------------------------------

      IMPLICIT NONE
!      SAVE		

  
Contains

subroutine QSORT_SINEX_SVN(array)
        use mdl_precision
        use mdl_param
        implicit none
        type (Sinex) :: hold, array(:)
        integer, parameter :: QSORT_THRESHOLD=32
        include "qsort_inline.inc"
        contains
logical function QSORT_COMPARE(ind1, ind2)
        use mdl_precision
        use mdl_param
        implicit none

        integer , intent(in) :: ind1, ind2
        type(Sinex) :: s1, s2

        integer*2   :: snum_s1, snum_s2, g1, g2
        integer*4   :: s1yr, s1doy, s1sod, s2yr, s2doy, s2sod
        character*1   :: gnss_s1, gnss_s2
        integer*2  :: ioerr

        s1 = satellites(ind1)
        s2 = satellites(ind2)
         
        gnss_s1 = s1%SVN(1:1)
        gnss_s2 = s2%SVN(1:1)

        g1 = ichar(gnss_s1)
        g2 = ichar(gnss_s2)

        s1yr = s1%STARTYR
        s1doy = s1%STARTDOY
        s1sod = s1%STARTSOD
        s2yr = s2%STARTYR
        s2doy = s2%STARTDOY
        s2sod = s2%STARTSOD
        QSORT_COMPARE = .false.

        READ (s1%SVN, '(1x,i3)', iostat=ioerr) snum_s1 
        if (ioerr /= 0) then
            return
        end if
        READ (s2%SVN, '(1x,i3)', iostat=ioerr) snum_s2 
        if (ioerr /= 0) then
            return
        end if

        if (g1 == g2) then
            if (snum_s1 == snum_s2) then
                if (s1yr == s2yr ) then
                    if (s1doy == s2doy) then
                        if (s1sod <= s2sod) then
                            QSORT_COMPARE = .true.
                        endif
                    elseif (s1doy < s2doy) then
                        QSORT_COMPARE = .true.
                    endif
                elseif (s1yr < s2yr) then
                    QSORT_COMPARE = .true.
                end if
           elseif (snum_s1 < snum_s2) then
                QSORT_COMPARE = .true.
           end if
        elseif (g1 < g2) then
            QSORT_COMPARE = .true.
        end if
        return
end function QSORT_COMPARE
end subroutine QSORT_SINEX_SVN

subroutine QSORT_SINEX_PRN(array)
        use mdl_precision
        use mdl_param
        implicit none
        type (Sinex) :: hold, array(:)
        integer, parameter :: QSORT_THRESHOLD=32
        include "qsort_inline.inc"
        contains
logical function QSORT_COMPARE(ind1, ind2)
        use mdl_precision
        use mdl_param
        implicit none

        integer , intent(in) :: ind1, ind2
        type(Sinex) :: s1, s2

        integer*2   :: snum_s1, snum_s2, g1, g2
        integer*4   :: s1yr, s1doy, s1sod, s2yr, s2doy, s2sod
        character*1   :: gnss_s1, gnss_s2
        integer*2  :: ioerr

        s1 = satellites(ind1)
        s2 = satellites(ind2)
         
        gnss_s1 = s1%PRN(1:1)
        gnss_s2 = s2%PRN(1:1)

        g1 = ichar(gnss_s1)
        g2 = ichar(gnss_s2)

        s1yr = s1%STARTYR
        s1doy = s1%STARTDOY
        s1sod = s1%STARTSOD
        s2yr = s2%STARTYR
        s2doy = s2%STARTDOY
        s2sod = s2%STARTSOD
        QSORT_COMPARE = .false.

        READ (s1%PRN, '(1x,i2)', iostat=ioerr) snum_s1 
        if (ioerr /= 0) then
            return
        end if
        READ (s2%PRN, '(1x,i2)', iostat=ioerr) snum_s2 
        if (ioerr /= 0) then
            return
        end if

        if (g1 == g2) then
            if (snum_s1 == snum_s2) then
                if (s1yr == s2yr ) then
                    if (s1doy == s2doy) then
                        if (s1sod <= s2sod) then
                            QSORT_COMPARE = .true.
                        endif
                    elseif (s1doy < s2doy) then
                        QSORT_COMPARE = .true.
                    endif
                elseif (s1yr < s2yr) then
                    QSORT_COMPARE = .true.
                end if
           elseif (snum_s1 < snum_s2) then
                QSORT_COMPARE = .true.
           end if
        elseif (g1 < g2) then
            QSORT_COMPARE = .true.
        end if
        return
end function QSORT_COMPARE
end subroutine QSORT_SINEX_PRN

SUBROUTINE lookup_sinex (idir,found, iyr,iday,ihr,imin,gnss,isat, &
                         ls_SVNID,ls_BLKTYP,ls_BLKID,ls_MASS,ls_POWER)

! ----------------------------------------------------------------------
! SUBROUTINE: lookup_sinex
! ----------------------------------------------------------------------
! Purpose:
!  locate satellite metadata information
! ----------------------------------------------------------------------
! Input arguments:
!
!  idir    : idir = 1 SVN to PRN, idir = -1 PRN to SVN  
!  iyr     : 4-digit year                               
!  iday    : 3-digit day of year                        
!  ihr     : 2-digit hr                                 
!  imin    : 2-digit minute                             
!  gnss    : 1-character GNSS code (G R E C J I)        
!  isat    : input SVN or PRN number                    
! 
! Output arguments:
!
!  ls_SVNID   : SVN NUMBER                    
!  ls_BLKTYP  : BLOCK TYPE   
!  ls_BLKID   : BLKID of type
!  ls_MASS    : S/C mass in kg 
!  ls_POWER   : transmitted power in watts                 
!  found   : we found the values in the array
!
! ----------------------------------------------------------------------
! Remarks:
!  
      use mdl_precision
      use mdl_param
      use mdl_config
      integer(kind = prec_int4) idir, isat, iyr, iday, ihr, imin
      logical found
      integer i
      character(LEN=1) :: gnss
      integer(Kind=2) ls_SVNID
      CHARACTER(LEN=20) :: ls_BLKTYP
      integer(kind=2) :: ls_BLKID
      REAL (KIND=prec_q) :: ls_MASS
      integer(kind=4)  ls_POWER
      CHARACTER(LEN=3) :: prn_in 
      CHARACTER(LEN=4) :: svn_in 
      real (kind=prec_d) :: time_id
      CHARACTER(LEN=256) :: message

! Which conversion is going to be implemented? 
IF(idir==-1) THEN  ! PRN to SVN
   WRITE(prn_in,'(a1,I2.2)') gnss, isat
   call QSORT_SINEX_PRN(satellites(1:SAT_COUNT))
ELSE                   ! SVN to PRN
   WRITE(svn_in,'(a1,I3.3)') gnss, isat
   call QSORT_SINEX_SVN(satellites(1:SAT_COUNT))
END IF

! Put the requested time into a scale of year 
time_id = iyr + iday/365.d0 + (ihr*3600.d0+imin*60.d0)/86400.d0/365.d0

found = .false.
ls_BLKTYP = ' '

! scan satellites array for matching value
if (idir==-1) then
  Do i = 1, SAT_COUNT
    if (satellites(i)%PRN == prn_in .and. time_id .ge. satellites(i)%TSTART .and. &
       satellites(i)%TSTOP .ge. time_id) then
      found = .true. !successful search
      read(satellites(i)%SVN(2:4), '(i3)') ls_SVNID
      ls_BLKTYP = TRIM(satellites(i)%BLKTYP)
      ls_BLKID = satellites(i)%BLKID
      ls_MASS = satellites(i)%MASS
      ls_POWER = satellites(i)%POWER
      write(message, '(a, i4, a, i4, a, i4, 3(a1, i4))') &
              'Located satellite '//prn_in//' at index ', i, ' for day ', iyr, ';', iday, ':', ihr, ':', imin
      call report('STATUS', pgrm_name, 'lookup_sinex', &
              message, ' ', 0)
      exit
    endif
  enddo
else
  do i = 1, SAT_COUNT
    if (satellites(i)%SVN == svn_in .and. time_id .ge. satellites(i)%TSTART .and. &
        satellites(i)%TSTOP .ge. time_id) then
      found = .true. !successful search
      read(satellites(i)%PRN(2:3), '(i2)') ls_SVNID
      ls_BLKTYP = TRIM(satellites(i)%BLKTYP)
      ls_BLKID = satellites(i)%BLKID
      ls_MASS = satellites(i)%MASS
      ls_POWER = satellites(i)%POWER
      write(message, '(a, i4, a, i4, a, i4, 3(a1, i4))') &
              'Located satellite '//svn_in//' at index ', i, ' for day ', iyr, ';', iday, ':', ihr, ':', imin
      call report('STATUS', pgrm_name, 'lookup_sinex', &
              message, ' ', 0)
      exit
    endif
  enddo
endif

!otherwise was not found - up to caller to decide what to do next

END SUBROUTINE

SUBROUTINE read_sinex_file (UNIT_IN)

! ----------------------------------------------------------------------
! SUBROUTINE: read_sinex_file
! ----------------------------------------------------------------------
! Purpose:
!  Read and store satellite metadata information
! ----------------------------------------------------------------------
! Input arguments:
!
!  UNIT     : logical unit number for read              
!
! ----------------------------------------------------------------------
! Remarks:
!  
      USE mdl_precision
      USE mdl_num
      use mdl_param
      use mdl_config
      IMPLICIT NONE

      INTEGER (KIND = prec_int4) :: UNIT_IN
      INTEGER (KIND = prec_int2) :: s_SVNID, s_BLKID
      CHARACTER(LEN=1) :: gnss, gnss_tmp, etype
      CHARACTER(LEN=3) :: satprn
      CHARACTER(LEN=4) :: satsvn
      CHARACTER(LEN=10):: cospar_id   
      CHARACTER(LEN=6) :: SatCat 
      CHARACTER(LEN=20):: antbody_in  ! Body type read from metadata snx

      INTEGER(KIND = prec_int2):: yr1,  yr2   ! Year read from sinex file
      INTEGER(KIND = prec_int2):: doy1, doy2  ! DOY read from sinex file

      REAL(KIND = prec_d) ::  sod1, sod2  ! Seconds of day read from sinex file

      CHARACTER(LEN=20):: s_BLKTYP
      CHARACTER(LEN=128):: record
      CHARACTER(LEN=256):: message

      INTEGER(KIND = prec_int4):: isat,satid,frqchn,s_POWER,iyr,iday,ihr,imin,ioerr,i,j,k

      REAL(KIND = prec_d) :: time_id, time1, time2, e_x, e_y, e_z
      REAL(KIND = prec_d) :: s_MASS
      
      LOGICAL found

!init ints to 0
yr1 = 0
yr2 = 0
doy1 = 0
doy2 = 0
sod1 = 0
sod2 = 0

! Initialize frquency channel (GLONASS only)
frqchn = 0

! Start with 0 satellites
SAT_COUNT = 0

! Start to read the SINEX file. We first look for '+SATELLITE/IDENTIFIER'
REWIND(UNIT_IN)
found = .false.
DO WHILE (.not.found)
READ(UNIT_IN,'(a)',iostat=ioerr) record
   IF(ioerr/=0 ) THEN
      call report('FATAL', pgrm_name, 'read_sinex_file', &
              'Failed to find IDENTIFIER SINEX block', ' ', ioerr)
   ELSEIF ( record(12:21)=='IDENTIFIER' ) THEN
      found = .true.
   END IF 
END DO 
found = .false.
DO WHILE (.not.found)
READ(UNIT_IN,'(a)',iostat=ioerr) record
   IF (record(1:4)=='-SAT'.or.ioerr/=0) THEN
      If(ioerr/=0) call report('WARNING', pgrm_name, 'read_sinex_file', &
              'No termination of IDENTIFIER SINEX block', ' ',ioerr)
      found = .true.
   ELSE 
      IF(record(1:1)==' ') THEN
      READ(record,'(1x,a4,1x,a9,1x,A6,1x,a15)',iostat=ioerr)satsvn,cospar_id,SatCat,antbody_in
      IF(ioerr/=0) print*,'Fail to reading satellite id ',ioerr
      SAT_COUNT = SAT_COUNT+1
      if (SAT_COUNT > MAX_SAT) THEN
         call report ('FATAL', pgrm_name, 'm_read_sinex_file', &
                      'too many satellites in sinex file', ' ', 0)
      end if
      s_BLKTYP = antbody_in
      ioerr = 0
      satellites(SAT_COUNT)%SVN = satsvn
      satellites(SAT_COUNT)%BLKTYP=TRIM(antbody_in)
      satellites(SAT_COUNT)%COSPAR=TRIM(cospar_id)
      !don't care about SatCat

      ! Prepare SVNID and BLKID for the global variables
      ! ------------------------------------------------
      s_BLKID = 0
      IF(TRIM(s_BLKTYP)=='GPS-I')      s_BLKID = 1
      IF(TRIM(s_BLKTYP)=='GPS-II')     s_BLKID = 2
      IF(TRIM(s_BLKTYP)=='GPS-IIA')    s_BLKID = 3
      IF(TRIM(s_BLKTYP)=='GPS-IIR')    s_BLKID = 4
      IF(TRIM(s_BLKTYP)=='GPS-IIR-A')  s_BLKID = 5
      IF(TRIM(s_BLKTYP)=='GPS-IIR-B')  s_BLKID = 6
      IF(TRIM(s_BLKTYP)=='GPS-IIR-M')  s_BLKID = 7
      IF(TRIM(s_BLKTYP)=='GPS-IIF')    s_BLKID = 8
      IF(TRIM(s_BLKTYP)=='GPS-IIIA')   s_BLKID = 9
      IF(TRIM(s_BLKTYP)=='GLO')        s_BLKID = 101
      IF(TRIM(s_BLKTYP)=='GLO-M'  .or.TRIM(s_BLKTYP) == 'GLO-M+')  s_BLKID = 102
      IF(TRIM(s_BLKTYP)=='GLO-K1A'.or.TRIM(s_BLKTYP) == 'GLO-K1B') s_BLKID = 103
      IF(TRIM(s_BLKTYP)=='GAL-1')     s_BLKID = 201 ! Galileo (IOV)
      IF(TRIM(s_BLKTYP)=='GAL-2')     s_BLKID = 202 ! Galileo (FOC)
      IF(TRIM(s_BLKTYP)=='BDS-2G'.or.TRIM(s_BLKTYP) == 'BDS-3G')            s_BLKID = 301 ! BDS GEO
      IF(TRIM(s_BLKTYP)=='BDS-2I'.or.TRIM(s_BLKTYP) == 'BDS-3I'.or. &
         TRIM(s_BLKTYP)=='BDS-3SI-SECM'.or.TRIM(s_BLKTYP) =='BDS-3SI-CAST') s_BLKID = 302 ! BDS IGSO
      IF(TRIM(s_BLKTYP)=='BDS-2M'.or.TRIM(s_BLKTYP) == 'BDS-3M'.or. &
         TRIM(s_BLKTYP)=='BDS-3M-SECM'.or.TRIM(s_BLKTYP) =='BDS-3M-CAST')   s_BLKID = 303 ! BDS MEO
      IF(TRIM(s_BLKTYP)=='QZS-1')     s_BLKID = 401 
      IF(TRIM(s_BLKTYP)=='QZS-2I')    s_BLKID = 402 ! QZSS-IGSO
      IF(TRIM(s_BLKTYP)=='QZS-2G')    s_BLKID = 403 ! QZSS-GEO
      satellites(SAT_COUNT)%BLKID = s_BLKID
      ! default all other variables for now
      satellites(SAT_COUNT)%PRN = ' '
      satellites(SAT_COUNT)%TSTART = 0.d0
      satellites(SAT_COUNT)%TSTOP = 0.d0
      satellites(SAT_COUNT)%MASS = 0.d0
      satellites(SAT_COUNT)%POWER = 0
      satellites(SAT_COUNT)%FRQCHN = 0
      satellites(SAT_COUNT)%E_PX = 0.d0
      satellites(SAT_COUNT)%E_PY = 0.d0
      satellites(SAT_COUNT)%E_PZ = 0.d0
      satellites(SAT_COUNT)%E_LX = 0.d0
      satellites(SAT_COUNT)%E_LY = 0.d0
      satellites(SAT_COUNT)%E_LZ = 0.d0
      satellites(SAT_COUNT)%STARTYR = 0
      satellites(SAT_COUNT)%STARTDOY = 0
      satellites(SAT_COUNT)%STARTSOD = 0
      satellites(SAT_COUNT)%STOPYR = 0
      satellites(SAT_COUNT)%STOPDOY = 0
      satellites(SAT_COUNT)%STOPSOD = 0
      END IF 
   END IF 
END DO 

call report('STATUS', pgrm_name, 'read_sinex_file', &
        'Finished reading Sinex Identifiers from file', ' ', 0)

! We have read all the IDENTIFIERS now. SAT_COUNT is now fixed. But because the remaining data
! does not necessarily align with the PRN data (to be read next) we may add more rows -
! one for each time window

! Now for eccentricity block. No times in eccentricity block so no row additions to be made, (SATELLITE/ECCENTRICITY block)
! just update whatever rows are there (P & L for each svn)
REWIND(UNIT_IN)
found = .false.
DO WHILE (.not.found)
   READ(UNIT_IN,'(a)',IOSTAT=ioerr) record
   IF(ioerr/=0 ) THEN
      call report('WARNING', pgrm_name, 'read_sinex_file', &
              'Failed to find ECCENTRICITY SINEX block', ' ', ioerr)
      found = .true.
   ELSEIF( record(12:15)=='ECCE' ) THEN
   found = .true.
   END IF
END DO
! if we didn't find the section skip to next optional segment
if (ioerr /= 0) found = .false.
DO WHILE (.not.found)
   READ(UNIT_IN,'(a)',IOSTAT=ioerr) record
   IF(record(1:4)=='-SAT'.or.ioerr/=0) THEN
      If(ioerr/=0) call report('WARNING', pgrm_name, 'read_sinex_file', &
              'No termination of ECCENTRICITY SINEX block', ' ',ioerr)
   found = .true.
   ELSE IF(record(1:1)==' ') THEN
      READ(record,'(1x,a4,23x,a1,3(1x,f9.4))',IOSTAT=ioerr) satsvn, etype, e_x, e_y, e_z
      IF(ioerr/=0)  call report ('WARNING', pgrm_name, 'read_sinex_file', &
              'Error reading satellite ECCENTRICITY', ' ', ioerr)

      DO i=1,SAT_COUNT
        if (satellites(i)%SVN == satsvn) then
           if (etype == 'P') then
              satellites(i)%E_PX = e_x
              satellites(i)%E_PY = e_y
              satellites(i)%E_PZ = e_z
           else
              satellites(i)%E_LX = e_x
              satellites(i)%E_LY = e_y
              satellites(i)%E_LZ = e_z
           endif
        endif
      ENDDO
  endif
ENDDO
ioerr = 0

! sort array based on svn (first) then start time
call Qsort_Sinex_SVN(satellites(1:SAT_COUNT))
call report('STATUS', pgrm_name, 'read_sinex_file', &
        'Finished reading Sinex Eccentricities from file', ' ', 0)

! now look for '+SATELLITE/PRN'
REWIND(UNIT_IN)
found = .false.
ioerr = 0
DO WHILE (.not.found)
   READ(UNIT_IN,'(a)',IOSTAT=ioerr) record
   ! error if we don't find the line
   IF (ioerr/=0) then
     call report ('FATAL', pgrm_name, 'read_sinex_file', &
                  'Could not find PRN identifier line', ' ', 0)
   else IF(record(1:14)=='+SATELLITE/PRN') THEN 
     !found it. Move to next part
     exit
   end if
end do

!Now read lines until we finish this section, updating the entry for each line
DO WHILE (.not.found)
   READ(UNIT_IN,'(a)',IOSTAT=ioerr) record
   IF (record(1:4)=='-SAT') THEN
       found = .true.
       cycle
   ! better detect if we are reading a valid line of the file!!!
   ELSE IF (record(1:1)==' ' .and. record(11:11)==':') THEN
   !ELSE IF (record(1:1)==' ') THEN
     yr2 = 0
     READ(record,'(1x,a4,2(1x,i4,1x,i3,1x,f5.0),1x,a3)',IOSTAT=ioerr) satsvn,yr1,doy1,sod1,yr2,doy2,sod2,satprn
     ! If no data for end time default to the end of 2100 ..
     IF(yr2==0000) THEN
     yr2   = 2100
     doy2  = 365
     sod2  = 86400
     END IF
     time1 = yr1 + doy1/365.d0 + sod1/86400.d0/365.d0
     time2 = yr2 + doy2/365.d0 + sod2/86400.d0/365.d0
     do i = 1, SAT_COUNT
       if (satellites(i)%SVN == satsvn) then
         if (satellites(i)%PRN == ' ') then
           ! no existing entry. Fill in the blanks
           satellites(i)%PRN = satprn
           satellites(i)%TSTART=time1
           satellites(i)%TSTOP=time2
           satellites(i)%STARTYR=yr1
           satellites(i)%STARTDOY=doy1
           satellites(i)%STARTSOD=sod1
           satellites(i)%STOPYR=yr2
           satellites(i)%STOPDOY=doy2
           satellites(I)%STOPSOD=sod2
         else
           ! we have an entry for this SVN already. Create new entry at end of array
           ! NB We do not check there is no overlap of epochs. Results will be undefined 
           ! if they do
           SAT_COUNT=SAT_COUNT + 1
           if (SAT_COUNT > MAX_SAT) call report('FATAL', pgrm_name, 'read_sinex_file', &
                   'Too many satellite rows in sinex file (increase MAX_SAT in mdl_param.f03)', &
                   ' ', 0)
           satellites(SAT_COUNT)%SVN = satsvn
           satellites(SAT_COUNT)%BLKTYP = satellites(i)%BLKTYP
           satellites(SAT_COUNT)%BLKID = satellites(i)%BLKID
           satellites(SAT_COUNT)%COSPAR = satellites(i)%COSPAR
           satellites(SAT_COUNT)%PRN = satprn
           satellites(SAT_COUNT)%TSTART = time1
           satellites(SAT_COUNT)%TSTOP = time2
           satellites(SAT_COUNT)%MASS = 0.d0
           satellites(SAT_COUNT)%POWER = 0
           satellites(SAT_COUNT)%STARTYR=yr1
           satellites(SAT_COUNT)%STARTDOY=doy1
           satellites(SAT_COUNT)%STARTSOD=sod1
           satellites(SAT_COUNT)%STOPYR=yr2
           satellites(SAT_COUNT)%STOPDOY=doy2
           satellites(SAT_COUNT)%STOPSOD=sod2
           satellites(SAT_COUNT)%FRQCHN = 0
           satellites(SAT_COUNT)%E_PX = satellites(i)%E_PX
           satellites(SAT_COUNT)%E_PY = satellites(i)%E_PY
           satellites(SAT_COUNT)%E_PZ = satellites(i)%E_PZ
           satellites(SAT_COUNT)%E_LX = satellites(i)%E_LX
           satellites(SAT_COUNT)%E_LY = satellites(i)%E_LY
           satellites(SAT_COUNT)%E_LZ = satellites(i)%E_LZ
         end if
         exit ! exit loop - we have inserted the data
       end if
     end do
   END IF
END DO

! sort array based on svn (first) then start time
call Qsort_Sinex_SVN(satellites(1:SAT_COUNT))
call report('STATUS', pgrm_name, 'read_sinex_file', &
        'Finished reading Sinex PRNs from file', ' ', 0)

ioerr = 0

! Get the satellite mass ( SATELLITE/MASS block )
REWIND(UNIT_IN)
found = .false.
ioerr = 0
DO WHILE (.not.found)
   READ(UNIT_IN,'(a)',IOSTAT=ioerr) record
   IF(ioerr/=0 ) THEN
      call report('FATAL', pgrm_name, 'read_sinex_file', &
              'Failed to find MASS SINEX block', ' ', ioerr)
   ELSEIF( record(12:15)=='MASS' ) THEN
   found = .true.
   END IF 
END DO 
found = .false.
ioerr = 0
DO WHILE (.not.found)
   READ(UNIT_IN,'(a)',IOSTAT=ioerr) record
   IF(record(1:4)=='-SAT'.or.ioerr/=0 ) THEN
      If(ioerr/=0) call report('WARNING', pgrm_name, 'read_sinex_file', &
              'No termination of MASS SINEX block', ' ',ioerr)
      found = .true.
! Only try to decode record if it is not a comment
   ELSE IF(record(1:1)==' ' .and. record(11:11)== ':') THEN
      READ(record,'(1x,a4,2(1x,i4,1x,i3,1x,f5.0),1x,f9.3)',IOSTAT=ioerr) satsvn,yr1,doy1,sod1,yr2,doy2,sod2,s_MASS
      IF(yr2==0000) THEN
         yr2 = 2100
         doy2 = 365
         sod2 = 86400
      END IF
      time1 = yr1+doy1/365.d0+sod1/86400.d0/365.d0
      time2 = yr2+doy2/365.d0+sod2/86400.d0/365.d0
      IF(ioerr/=0)  call report ('WARNING', pgrm_name, 'read_sinex_file', &
              'Error reading satellite mass ', ' ', ioerr)
      do i = 1, SAT_COUNT
        if (satellites(i)%SVN /= satsvn .or. satellites(i)%TSTOP < time1 .or. satellites(i)%TSTART > time2) cycle
        j = i
        exit ! found the first one that fits the time window given. FIXME: if outside the window should we not
             ! create a new row?
      enddo
      i = j
        if (satellites(i)%SVN == satsvn) then
          if (satellites(i)%TSTART == time1 .and. satellites(i)%TSTOP == time2) then
            !yippee. It aligns. no splitting of record required. Update MASS, read next record
            satellites(i)%MASS = s_MASS
            cycle
          endif
          if (satellites(i)%TSTART > time1) then
            ! need to append new row for this SVN, from time1 to TSTART with no PRN
            satellites(i)%MASS = s_MASS
            SAT_COUNT = SAT_COUNT + 1
            if (SAT_COUNT > MAX_SAT) call report('FATAL', pgrm_name, 'read_sinex_file', &
                   'Too many satellite rows in sinex file (increase MAX_SAT in mdl_param.f03)', &
                   ' ', 0)
            satellites(SAT_COUNT)%SVN = satsvn
            satellites(SAT_COUNT)%BLKTYP = satellites(i)%BLKTYP
            satellites(SAT_COUNT)%BLKID = satellites(i)%BLKID
            satellites(SAT_COUNT)%COSPAR = satellites(i)%COSPAR
            satellites(SAT_COUNT)%PRN = ' '
            satellites(SAT_COUNT)%TSTART = time1
            satellites(SAT_COUNT)%TSTOP = satellites(i)%TSTART
            satellites(SAT_COUNT)%MASS = s_MASS
            satellites(SAT_COUNT)%POWER = satellites(i)%POWER
            satellites(SAT_COUNT)%FRQCHN = satellites(i)%FRQCHN
            satellites(SAT_COUNT)%E_PX = satellites(i)%E_PX
            satellites(SAT_COUNT)%E_PY = satellites(i)%E_PY
            satellites(SAT_COUNT)%E_PZ = satellites(i)%E_PZ
            satellites(SAT_COUNT)%E_LX = satellites(i)%E_LX
            satellites(SAT_COUNT)%E_LY = satellites(i)%E_LY
            satellites(SAT_COUNT)%E_LZ = satellites(i)%E_LZ
            satellites(SAT_COUNT)%STARTYR=yr1
            satellites(SAT_COUNT)%STARTDOY=doy1
            satellites(SAT_COUNT)%STARTSOD=sod1
            satellites(SAT_COUNT)%STOPYR=satellites(i)%STARTYR
            satellites(SAT_COUNT)%STOPDOY=satellites(i)%STARTDOY
            satellites(SAT_COUNT)%STOPSOD=satellites(i)%STARTSOD
          elseif (satellites(i)%TSTART < time1) then
            ! need to append new row for this SVN, from TSTART to time1 with PRN but no MASS
            ! amend current row to start at time1
            SAT_COUNT = SAT_COUNT + 1
            if (SAT_COUNT > MAX_SAT) call report('FATAL', pgrm_name, 'read_sinex_file', &
                   'Too many satellite rows in sinex file (increase MAX_SAT in mdl_param.f03)', &
                   ' ', 0)
            satellites(SAT_COUNT)%SVN = satsvn
            satellites(SAT_COUNT)%BLKTYP = satellites(i)%BLKTYP
            satellites(SAT_COUNT)%BLKID = satellites(i)%BLKID
            satellites(SAT_COUNT)%COSPAR = satellites(i)%COSPAR
            satellites(SAT_COUNT)%PRN = satellites(i)%PRN
            satellites(SAT_COUNT)%TSTART = satellites(i)%TSTART
            satellites(SAT_COUNT)%STARTYR = satellites(i)%STARTYR
            satellites(SAT_COUNT)%STARTDOY = satellites(i)%STARTDOY
            satellites(SAT_COUNT)%STARTSOD = satellites(i)%STARTSOD
            satellites(SAT_COUNT)%TSTOP = time1
            satellites(SAT_COUNT)%STOPYR = yr1
            satellites(SAT_COUNT)%STOPDOY = doy1
            satellites(SAT_COUNT)%STOPSOD = sod1
            satellites(SAT_COUNT)%MASS = 0.d0
            satellites(SAT_COUNT)%POWER = satellites(i)%POWER
            satellites(SAT_COUNT)%FRQCHN = satellites(i)%FRQCHN
            satellites(SAT_COUNT)%E_PX = satellites(i)%E_PX
            satellites(SAT_COUNT)%E_PY = satellites(i)%E_PY
            satellites(SAT_COUNT)%E_PZ = satellites(i)%E_PZ
            satellites(SAT_COUNT)%E_LX = satellites(i)%E_LX
            satellites(SAT_COUNT)%E_LY = satellites(i)%E_LY
            satellites(SAT_COUNT)%E_LZ = satellites(i)%E_LZ
            satellites(i)%TSTART = time1
            satellites(i)%STARTYR=yr1
            satellites(i)%STARTDOY=doy1
            satellites(i)%STARTSOD=sod1
          endif 
          if (satellites(i)%TSTOP < time2) then
            ! check next row(s) and update mass again, we may also need to create a new row for any trailing period
            ! note that sort ordering means we only need to check forwards
            Do while (satellites(i)%SVN == satsvn .and. satellites(i)%TSTOP <= time2)
                satellites(i)%MASS=s_MASS
                i=i+1
            enddo
          endif
          if (satellites(i)%SVN == satsvn .and. satellites(i)%TSTOP < time2) then
            ! create new row with MASS but no PRN
            satellites(i)%MASS = s_MASS
            SAT_COUNT = SAT_COUNT + 1
            if (SAT_COUNT > MAX_SAT) call report('FATAL', pgrm_name, 'read_sinex_file', &
               'Too many satellite rows in sinex file (increase MAX_SAT in mdl_param.f03)', &
               ' ', 0)
            satellites(SAT_COUNT)%SVN = satsvn
            satellites(SAT_COUNT)%BLKTYP = satellites(i)%BLKTYP
            satellites(SAT_COUNT)%BLKID = satellites(i)%BLKID
            satellites(SAT_COUNT)%COSPAR = satellites(i)%COSPAR
            satellites(SAT_COUNT)%PRN = ' '
            satellites(SAT_COUNT)%TSTART = satellites(i)%TSTOP
            satellites(SAT_COUNT)%STARTYR = satellites(i)%STOPYR
            satellites(SAT_COUNT)%STARTDOY = satellites(i)%STOPDOY
            satellites(SAT_COUNT)%STARTSOD = satellites(i)%STOPSOD
            satellites(SAT_COUNT)%TSTOP = time2
            satellites(SAT_COUNT)%STOPYR=yr2
            satellites(SAT_COUNT)%STOPDOY = doy2
            satellites(SAT_COUNT)%STOPSOD=sod2
            satellites(SAT_COUNT)%MASS = s_MASS
            satellites(SAT_COUNT)%POWER = satellites(i)%POWER
            satellites(SAT_COUNT)%FRQCHN = satellites(i)%FRQCHN
            satellites(SAT_COUNT)%E_PX = satellites(i)%E_PX
            satellites(SAT_COUNT)%E_PY = satellites(i)%E_PY
            satellites(SAT_COUNT)%E_PZ = satellites(i)%E_PZ
            satellites(SAT_COUNT)%E_LX = satellites(i)%E_LX
            satellites(SAT_COUNT)%E_LY = satellites(i)%E_LY
            satellites(SAT_COUNT)%E_LZ = satellites(i)%E_LZ
          elseif (satellites(i)%SVN == satsvn .and. satellites(i)%TSTOP > time2) then
            ! need to append new row for this SVN, from time2 to TSTOP with no MASS
            satellites(i)%MASS = s_MASS
            SAT_COUNT = SAT_COUNT + 1
            if (SAT_COUNT > MAX_SAT) call report('FATAL', pgrm_name, 'read_sinex_file', &
                   'Too many satellite rows in sinex file (increase MAX_SAT in mdl_param.f03)', &
                   ' ', 0)
            satellites(SAT_COUNT)%SVN = satsvn
            satellites(SAT_COUNT)%BLKTYP = satellites(i)%BLKTYP
            satellites(SAT_COUNT)%BLKID = satellites(i)%BLKID
            satellites(SAT_COUNT)%COSPAR = satellites(i)%COSPAR
            satellites(SAT_COUNT)%PRN = satellites(i)%PRN
            satellites(SAT_COUNT)%TSTART = time2
            satellites(SAT_COUNT)%STARTYR = yr2
            satellites(SAT_COUNT)%STARTDOY = doy2
            satellites(SAT_COUNT)%STARTSOD = sod2
            satellites(SAT_COUNT)%TSTOP = satellites(i)%TSTOP
            satellites(SAT_COUNT)%STOPYR = satellites(i)%STOPYR
            satellites(SAT_COUNT)%STOPDOY = satellites(i)%STOPDOY
            satellites(SAT_COUNT)%STOPSOD = satellites(i)%STOPSOD
            satellites(SAT_COUNT)%MASS = 0.d0
            satellites(SAT_COUNT)%POWER = satellites(i)%POWER
            satellites(SAT_COUNT)%FRQCHN = satellites(i)%FRQCHN
            satellites(SAT_COUNT)%E_PX = satellites(i)%E_PX
            satellites(SAT_COUNT)%E_PY = satellites(i)%E_PY
            satellites(SAT_COUNT)%E_PZ = satellites(i)%E_PZ
            satellites(SAT_COUNT)%E_LX = satellites(i)%E_LX
            satellites(SAT_COUNT)%E_LY = satellites(i)%E_LY
            satellites(SAT_COUNT)%E_LZ = satellites(i)%E_LZ
            satellites(i)%TSTOP = time2
            satellites(i)%STOPYR=yr2
            satellites(i)%STOPDOY=doy2
            satellites(i)%STOPSOD=sod2
          endif
        endif
    END IF 
END DO 

! sort array based on svn (first) then start time
call Qsort_Sinex_Svn(satellites(1:SAT_COUNT))
ioerr = 0
call report('STATUS', pgrm_name, 'read_sinex_file', &
        'Finished reading Sinex MASS from file', ' ', 0)

! This only needed for Glonass. We might not have Glonass satellites in the file. So this is 
! not a FATAL error
REWIND(UNIT_IN)
found = .false.
ioerr = 0
   DO WHILE(.not.found)
   READ(UNIT_IN,'(a)',IOSTAT=ioerr) record
      IF(ioerr/=0 ) THEN
         call report('WARNING', pgrm_name, 'read_sinex_file', &
                     'Failed to find FREQUENCY_CHANNEL SINEX block', ' ', ioerr)
         found = .true.
      ELSEIF(record(12:15)=='FREQ') THEN
         found = .true.
      END IF
   END DO 
   found = .false.
   IF (ioerr==0) then
   DO WHILE (.not.found)
      READ(UNIT_IN,'(a)',IOSTAT=ioerr) record
      IF(record(1:4)=='-SAT'.or.ioerr.ne.0) THEN
         If( ioerr/=0 ) call report('WARNING', pgrm_name, 'read_sinex_file', &
                 'No termination of FREQUENCY_CHANNEL SINEX block', ' ', ioerr)
         found = .true.
      ELSE IF(record(1:1)==' ' .and. record(11:11)== ':') THEN
         !IF(record(1:1)==' ') THEN
         READ(record,'(1x,a4,2(1x,i4,1x,i3,1x,f5.0),1x,i4)',IOSTAT=ioerr) satsvn,yr1,doy1,sod1,yr2,doy2,sod2,frqchn
         IF(ioerr/=0) call report('WARNING', pgrm_name, 'read_sinex_file', &
                 'Reading GLONASS frequency:'//record, ' ',ioerr)
         IF(yr2==0000) THEN
         yr2  = 2100
         doy2 = 365
         sod2 = 86400
         END IF
         time1 = yr1 + doy1/365.d0 + sod1/86400.d0/365.d0
         time2 = yr2 + doy2/365.d0 + sod2/86400.d0/365.d0
         
      do i = 1, SAT_COUNT
        if (satellites(i)%SVN /= satsvn .or. satellites(i)%TSTOP < time1 .or. satellites(i)%TSTART > time2) cycle
        j = i
        exit ! FIXME: if on correct svn but totally before the current window should we not create a new row
      enddo
      i = j
          if (satellites(i)%TSTART == time1 .and. satellites(i)%TSTOP == time2) then
            !yippee. It aligns. no splitting of record required. Read next record
            satellites(i)%FRQCHN = frqchn
            cycle
          endif
          if (satellites(i)%TSTART > time1) then
            ! need to append new row for this SVN, from time1 to TSTART with no PRN or MASS
            satellites(i)%FRQCHN = frqchn
            SAT_COUNT = SAT_COUNT + 1
            if (SAT_COUNT > MAX_SAT) call report('FATAL', pgrm_name, 'read_sinex_file', &
                   'Too many satellite rows in sinex file (increase MAX_SAT in mdl_param.f03)', &
                   ' ', 0)
            satellites(SAT_COUNT)%SVN = satsvn
            satellites(SAT_COUNT)%BLKTYP = satellites(i)%BLKTYP
            satellites(SAT_COUNT)%BLKID = satellites(i)%BLKID
            satellites(SAT_COUNT)%COSPAR = satellites(i)%COSPAR
            satellites(SAT_COUNT)%PRN = ' '
            satellites(SAT_COUNT)%TSTART = time1
            satellites(SAT_COUNT)%STARTYR = yr1
            satellites(SAT_COUNT)%STARTDOY = doy1
            satellites(SAT_COUNT)%STARTSOD = sod1
            satellites(SAT_COUNT)%TSTOP = satellites(i)%TSTART
            satellites(SAT_COUNT)%STOPYR = satellites(i)%STARTYR
            satellites(SAT_COUNT)%STOPDOY = satellites(i)%STARTDOY
            satellites(SAT_COUNT)%STOPSOD = satellites(i)%STARTSOD
            satellites(SAT_COUNT)%MASS = 0.d0
            satellites(SAT_COUNT)%POWER = satellites(i)%POWER
            satellites(SAT_COUNT)%FRQCHN = FRQCHN
            satellites(SAT_COUNT)%E_PX = satellites(i)%E_PX
            satellites(SAT_COUNT)%E_PY = satellites(i)%E_PY
            satellites(SAT_COUNT)%E_PZ = satellites(i)%E_PZ
            satellites(SAT_COUNT)%E_LX = satellites(i)%E_LX
            satellites(SAT_COUNT)%E_LY = satellites(i)%E_LY
            satellites(SAT_COUNT)%E_LZ = satellites(i)%E_LZ
          elseif (satellites(i)%TSTART < time1) then
            ! need to append new row for this SVN, from TSTART to time1 with PRN & MASS but no frqchn
            ! amend current row to start at time1
            SAT_COUNT = SAT_COUNT + 1
            if (SAT_COUNT > MAX_SAT) call report('FATAL', pgrm_name, 'read_sinex_file', &
                   'Too many satellite rows in sinex file (increase MAX_SAT in mdl_param.f03)', &
                   ' ', 0)
            satellites(SAT_COUNT)%SVN = satsvn
            satellites(SAT_COUNT)%BLKTYP = satellites(i)%BLKTYP
            satellites(SAT_COUNT)%BLKID = satellites(i)%BLKID
            satellites(SAT_COUNT)%COSPAR = satellites(i)%COSPAR
            satellites(SAT_COUNT)%PRN = satellites(i)%PRN
            satellites(SAT_COUNT)%TSTART = satellites(i)%TSTART
            satellites(SAT_COUNT)%STARTYR = satellites(i)%STARTYR
            satellites(SAT_COUNT)%STARTDOY = satellites(i)%STARTDOY
            satellites(SAT_COUNT)%STARTSOD = satellites(i)%STARTSOD
            satellites(SAT_COUNT)%TSTOP = time1
            satellites(SAT_COUNT)%STOPYR = yr1
            satellites(SAT_COUNT)%STOPDOY = doy1
            satellites(SAT_COUNT)%STOPSOD = sod1
            satellites(SAT_COUNT)%MASS = satellites(i)%MASS
            satellites(SAT_COUNT)%POWER = satellites(i)%POWER
            satellites(SAT_COUNT)%FRQCHN = 0
            satellites(SAT_COUNT)%E_PX = satellites(i)%E_PX
            satellites(SAT_COUNT)%E_PY = satellites(i)%E_PY
            satellites(SAT_COUNT)%E_PZ = satellites(i)%E_PZ
            satellites(SAT_COUNT)%E_LX = satellites(i)%E_LX
            satellites(SAT_COUNT)%E_LY = satellites(i)%E_LY
            satellites(SAT_COUNT)%E_LZ = satellites(i)%E_LZ
            satellites(i)%TSTART = time1
            satellites(i)%STARTYR = yr1
            satellites(i)%STARTDOY = doy1
            satellites(i)%STARTSOD = sod1
          endif 
          if (satellites(i)%TSTOP < time2) then
            ! need to append new row for this SVN, from TSTOP to time2 with no PRN or MASS
            do while (satellites(i)%SVN == satsvn .and. satellites(i)%TSTOP <= time2)
              satellites(i)%FRQCHN = frqchn
              i = i+1
            enddo
          endif
          if (satellites(i)%SVN ==satsvn .and. satellites(i)%TSTOP < time2) then
            !OK new row required
            SAT_COUNT = SAT_COUNT + 1
            if (SAT_COUNT > MAX_SAT) call report('FATAL', pgrm_name, 'read_sinex_file', &
                 'Too many satellite rows in sinex file (increase MAX_SAT in mdl_param.f03)', &
                 ' ', 0)
            satellites(SAT_COUNT)%SVN = satsvn
            satellites(SAT_COUNT)%BLKTYP = satellites(i)%BLKTYP
            satellites(SAT_COUNT)%BLKID = satellites(i)%BLKID
            satellites(SAT_COUNT)%COSPAR = satellites(i)%COSPAR
            satellites(SAT_COUNT)%PRN = ' '
            satellites(SAT_COUNT)%TSTART = satellites(i)%TSTOP
            satellites(SAT_COUNT)%STARTYR = satellites(i)%STOPYR
            satellites(SAT_COUNT)%STARTDOY = satellites(i)%STOPDOY
            satellites(SAT_COUNT)%STARTSOD = satellites(i)%STOPSOD
            satellites(SAT_COUNT)%TSTOP = time2
            satellites(SAT_COUNT)%STOPYR = yr2
            satellites(SAT_COUNT)%STOPDOY = doy2
            satellites(SAT_COUNT)%STOPSOD = sod2
            satellites(SAT_COUNT)%MASS = 0.d0
            satellites(SAT_COUNT)%POWER = satellites(i)%POWER
            satellites(SAT_COUNT)%FRQCHN = satellites(i)%FRQCHN
            satellites(SAT_COUNT)%E_PX = satellites(i)%E_PX
            satellites(SAT_COUNT)%E_PY = satellites(i)%E_PY
            satellites(SAT_COUNT)%E_PZ = satellites(i)%E_PZ
            satellites(SAT_COUNT)%E_LX = satellites(i)%E_LX
            satellites(SAT_COUNT)%E_LY = satellites(i)%E_LY
            satellites(SAT_COUNT)%E_LZ = satellites(i)%E_LZ
          elseif (satellites(i)%SVN == satsvn .and. satellites(i)%TSTOP > time2) then
            ! need to append new row for this SVN, from time2 to TSTOP with no FRQCHN
            satellites(i)%FRQCHN = frqchn
            SAT_COUNT = SAT_COUNT + 1
            if (SAT_COUNT > MAX_SAT) call report('FATAL', pgrm_name, 'read_sinex_file', &
                   'Too many satellite rows in sinex file (increase MAX_SAT in mdl_param.f03)', &
                   ' ', 0)
            satellites(SAT_COUNT)%SVN = satsvn
            satellites(SAT_COUNT)%BLKTYP = satellites(i)%BLKTYP
            satellites(SAT_COUNT)%BLKID = satellites(i)%BLKID
            satellites(SAT_COUNT)%COSPAR = satellites(i)%COSPAR
            satellites(SAT_COUNT)%PRN = satellites(i)%PRN
            satellites(SAT_COUNT)%TSTART = time2
            satellites(SAT_COUNT)%STARTYR = yr2
            satellites(SAT_COUNT)%STARTDOY = doy2
            satellites(SAT_COUNT)%STARTSOD = sod2
            satellites(SAT_COUNT)%TSTOP = satellites(i)%TSTOP
            satellites(SAT_COUNT)%STOPYR = satellites(i)%STOPYR
            satellites(SAT_COUNT)%STOPDOY = satellites(i)%STOPDOY
            satellites(SAT_COUNT)%STOPSOD = satellites(i)%STOPSOD
            satellites(SAT_COUNT)%MASS = satellites(i)%MASS
            satellites(SAT_COUNT)%POWER = satellites(i)%POWER
            satellites(SAT_COUNT)%FRQCHN = 0
            satellites(SAT_COUNT)%E_PX = satellites(i)%E_PX
            satellites(SAT_COUNT)%E_PY = satellites(i)%E_PY
            satellites(SAT_COUNT)%E_PZ = satellites(i)%E_PZ
            satellites(SAT_COUNT)%E_LX = satellites(i)%E_LX
            satellites(SAT_COUNT)%E_LY = satellites(i)%E_LY
            satellites(SAT_COUNT)%E_LZ = satellites(i)%E_LZ
            satellites(i)%TSTOP = time2
            satellites(i)%STOPYR=yr2
            satellites(i)%STOPDOY=doy2
            satellites(i)%STOPSOD=sod2
          endif
          cycle ! read next line
      END IF
   END DO 
   END IF

! sort array based on svn (first) then start time
call Qsort_Sinex_Svn(satellites(1:SAT_COUNT))
ioerr = 0
call report('STATUS', pgrm_name, 'read_sinex_file', &
        'Finished reading Sinex FRQCHN from file', ' ', 0)

! Get the transmitter power ( SATELLITE/TX_POWER block)
REWIND(UNIT_IN)
found = .false.
DO WHILE (.not.found)
   READ(UNIT_IN,'(a)',IOSTAT=ioerr) record
   IF(ioerr/=0 ) THEN
      call report('WARNING', pgrm_name, 'read_sinex_file', &
              'Failed to find TX_POWER SINEX block', ' ', ioerr)
      found = .true.
   ELSEIF( record(12:15)=='TX_P' ) THEN
   found = .true.
   END IF
END DO
k = 0
! if we didn't find the section skip to next optional segment
if (ioerr == 0) found = .false.
ioerr = 0
DO WHILE (.not.found)
   READ(UNIT_IN,'(a)',IOSTAT=ioerr) record
   IF(record(1:4)=='-SAT'.or.ioerr/=0) THEN
      If(ioerr/=0) call report('WARNING', pgrm_name, 'read_sinex_file', &
              'No termination of TX_POWER SINEX block', ' ',ioerr)
   found = .true.
   ELSE IF(record(1:1)==' '.and.record(11:11)== ':') THEN
      k = k+1
      READ(record,'(1x,a4,2(1x,i4,1x,i3,1x,F5.0),1x,i4)',IOSTAT=ioerr) satsvn,yr1,doy1,sod1,yr2,doy2,sod2,s_POWER
!      write (message, '("SVN :",a,", POWER:", i4)') satsvn, s_POWER
!      call report('STATUS', pgrm_name, 'read_sinex_file', &
!              message, ' ', 0)
      IF(yr2==0000) THEN
         yr2 = 2100
        doy2 = 365
        sod2 = 86400
      END IF 
      IF(ioerr/=0)  call report ('WARNING', pgrm_name, 'read_sinex_file', &
              'Error reading satellite TX power ', ' ', ioerr)
      time1 = yr1 + doy1/365.d0 + sod1/86400.d0/365.d0
      time2 = yr2 + doy2/365.d0 + sod2/86400.d0/365.d0
      ioerr = 0
      do i = 1, SAT_COUNT
        if (satellites(i)%SVN /= satsvn .or. satellites(i)%TSTOP < time1 .or. satellites(i)%TSTART > time2) cycle ! Not in this window 
        j = i
        exit
      enddo
      i = j
      if (satellites(i)%TSTART == time1 .and. satellites(i)%TSTOP == time2) then
         !yippee. It aligns. no splitting of record required. Read next record
         satellites(i)%POWER=s_POWER
         cycle
      endif
          if (satellites(i)%TSTART > time1) then
            ! need to append new row for this SVN, from time1 to TSTART with no PRN, MASS, or FRQCHN
            satellites(i)%POWER = s_POWER
            SAT_COUNT = SAT_COUNT + 1
            if (SAT_COUNT > MAX_SAT) call report('FATAL', pgrm_name, 'read_sinex_file', &
                   'Too many satellite rows in sinex file (increase MAX_SAT in mdl_param.f03)', &
                   ' ', 0)
            satellites(SAT_COUNT)%SVN = satsvn
            satellites(SAT_COUNT)%BLKTYP = satellites(i)%BLKTYP
            satellites(SAT_COUNT)%BLKID = satellites(i)%BLKID
            satellites(SAT_COUNT)%COSPAR = satellites(i)%COSPAR
            satellites(SAT_COUNT)%PRN = ' '
            satellites(SAT_COUNT)%TSTART = time1
            satellites(SAT_COUNT)%STARTYR = yr1
            satellites(SAT_COUNT)%STARTDOY = doy1
            satellites(SAT_COUNT)%STARTSOD = sod1
            satellites(SAT_COUNT)%TSTOP = satellites(i)%TSTART
            satellites(SAT_COUNT)%STOPYR = satellites(i)%STARTYR
            satellites(SAT_COUNT)%STOPDOY= satellites(i)%STARTDOY
            satellites(SAT_COUNT)%STOPSOD = satellites(i)%STARTSOD
            satellites(SAT_COUNT)%MASS = 0.d0
            satellites(SAT_COUNT)%POWER = s_POWER
            satellites(SAT_COUNT)%FRQCHN = 0
            satellites(SAT_COUNT)%E_PX = satellites(i)%E_PX
            satellites(SAT_COUNT)%E_PY = satellites(i)%E_PY
            satellites(SAT_COUNT)%E_PZ = satellites(i)%E_PZ
            satellites(SAT_COUNT)%E_LX = satellites(i)%E_LX
            satellites(SAT_COUNT)%E_LY = satellites(i)%E_LY
            satellites(SAT_COUNT)%E_LZ = satellites(i)%E_LZ
          elseif (satellites(i)%TSTART < time1) then
            ! need to append new row for this SVN, from TSTART to time1 with PRN, MASS & FRQCHN
            ! amend current row to end at time1
            SAT_COUNT = SAT_COUNT + 1
            if (SAT_COUNT > MAX_SAT) call report('FATAL', pgrm_name, 'read_sinex_file', &
                   'Too many satellite rows in sinex file (increase MAX_SAT in mdl_param.f03)', &
                   ' ', 0)
            satellites(SAT_COUNT)%SVN = satsvn
            satellites(SAT_COUNT)%BLKTYP = satellites(i)%BLKTYP
            satellites(SAT_COUNT)%BLKID = satellites(i)%BLKID
            satellites(SAT_COUNT)%COSPAR = satellites(i)%COSPAR
            satellites(SAT_COUNT)%PRN = satellites(i)%PRN
            satellites(SAT_COUNT)%TSTART = time1
            satellites(SAT_COUNT)%STARTYR = yr1
            satellites(SAT_COUNT)%STARTDOY = doy1
            satellites(SAT_COUNT)%STARTSOD = sod1
            satellites(SAT_COUNT)%TSTOP = satellites(i)%TSTOP
            satellites(SAT_COUNT)%STOPYR = satellites(i)%STARTYR
            satellites(SAT_COUNT)%STOPDOY= satellites(i)%STARTDOY
            satellites(SAT_COUNT)%STOPSOD = satellites(i)%STARTSOD
            satellites(SAT_COUNT)%MASS = satellites(i)%MASS
            satellites(SAT_COUNT)%POWER = s_POWER
            satellites(SAT_COUNT)%FRQCHN = satellites(i)%FRQCHN
            satellites(SAT_COUNT)%E_PX = satellites(i)%E_PX
            satellites(SAT_COUNT)%E_PY = satellites(i)%E_PY
            satellites(SAT_COUNT)%E_PZ = satellites(i)%E_PZ
            satellites(SAT_COUNT)%E_LX = satellites(i)%E_LX
            satellites(SAT_COUNT)%E_LY = satellites(i)%E_LY
            satellites(SAT_COUNT)%E_LZ = satellites(i)%E_LZ
            satellites(i)%TSTOP = time1
            satellites(i)%STOPYR = yr1
            satellites(i)%STOPDOY = doy1
            satellites(i)%STOPSOD = sod1
          endif 
          if (satellites(i)%TSTOP < time2) then
            ! check next row(s) and update mass again, we may also need to create a new row for any trailing period
            ! note that sort ordering means we only need to check forwards
            Do while (satellites(i)%SVN == satsvn .and. satellites(i)%TSTOP <= time2)
                satellites(i)%POWER = s_POWER
                i=i+1
            enddo
          endif
          if (satellites(i)%TSTOP < time2) then
            ! need to append new row for this SVN, from TSTOP to time2 with no PRN, MASS or FRQCHN
            satellites(i)%POWER = s_POWER
            SAT_COUNT = SAT_COUNT + 1
            if (SAT_COUNT > MAX_SAT) call report('FATAL', pgrm_name, 'read_sinex_file', &
                   'Too many satellite rows in sinex file (increase MAX_SAT in mdl_param.f03)', &
                   ' ', 0)
            satellites(SAT_COUNT)%SVN = satsvn
            satellites(SAT_COUNT)%BLKTYP = satellites(i)%BLKTYP
            satellites(SAT_COUNT)%BLKID = satellites(i)%BLKID
            satellites(SAT_COUNT)%COSPAR = satellites(i)%COSPAR
            satellites(SAT_COUNT)%PRN = ' '
            satellites(SAT_COUNT)%TSTART = satellites(i)%TSTOP
            satellites(SAT_COUNT)%STARTYR = satellites(i)%STOPYR
            satellites(SAT_COUNT)%STARTDOY= satellites(i)%STOPDOY
            satellites(SAT_COUNT)%STARTSOD = satellites(i)%STOPSOD
            satellites(SAT_COUNT)%TSTOP = time2
            satellites(SAT_COUNT)%STOPYR = yr2
            satellites(SAT_COUNT)%STOPDOY = doy2
            satellites(SAT_COUNT)%STOPSOD = sod2
            satellites(SAT_COUNT)%MASS = 0.d0
            satellites(SAT_COUNT)%POWER = s_POWER
            satellites(SAT_COUNT)%E_PX = satellites(i)%E_PX
            satellites(SAT_COUNT)%E_PY = satellites(i)%E_PY
            satellites(SAT_COUNT)%E_PZ = satellites(i)%E_PZ
            satellites(SAT_COUNT)%E_LX = satellites(i)%E_LX
            satellites(SAT_COUNT)%E_LY = satellites(i)%E_LY
            satellites(SAT_COUNT)%E_LZ = satellites(i)%E_LZ
            satellites(SAT_COUNT)%FRQCHN = 0
          elseif (satellites(i)%TSTOP > time2) then
            ! need to append new row for this SVN, from time2 to TSTOP with no POWER
            satellites(i)%POWER = s_POWER
            SAT_COUNT = SAT_COUNT + 1
            if (SAT_COUNT > MAX_SAT) call report('FATAL', pgrm_name, 'read_sinex_file', &
                   'Too many satellite rows in sinex file (increase MAX_SAT in mdl_param.f03)', &
                   ' ', 0)
            satellites(SAT_COUNT)%SVN = satsvn
            satellites(SAT_COUNT)%BLKTYP = satellites(i)%BLKTYP
            satellites(SAT_COUNT)%BLKID = satellites(i)%BLKID
            satellites(SAT_COUNT)%COSPAR = satellites(i)%COSPAR
            satellites(SAT_COUNT)%PRN = satellites(i)%PRN
            satellites(SAT_COUNT)%TSTART = time2
            satellites(SAT_COUNT)%STARTYR = yr2
            satellites(SAT_COUNT)%STARTDOY = doy2
            satellites(SAT_COUNT)%STARTSOD = sod2
            satellites(SAT_COUNT)%TSTOP = satellites(i)%TSTOP
            satellites(SAT_COUNT)%STOPYR = satellites(i)%STOPYR
            satellites(SAT_COUNT)%STOPDOY = satellites(i)%STOPDOY
            satellites(SAT_COUNT)%STOPSOD = satellites(i)%STOPSOD
            satellites(SAT_COUNT)%MASS = satellites(i)%MASS
            satellites(SAT_COUNT)%POWER = 0
            satellites(SAT_COUNT)%FRQCHN = satellites(i)%FRQCHN
            satellites(SAT_COUNT)%E_PX = satellites(i)%E_PX
            satellites(SAT_COUNT)%E_PY = satellites(i)%E_PY
            satellites(SAT_COUNT)%E_PZ = satellites(i)%E_PZ
            satellites(SAT_COUNT)%E_LX = satellites(i)%E_LX
            satellites(SAT_COUNT)%E_LY = satellites(i)%E_LY
            satellites(SAT_COUNT)%E_LZ = satellites(i)%E_LZ
            satellites(i)%TSTOP = time2
            satellites(i)%STOPYR = yr2
            satellites(i)%STOPDOY = doy2
            satellites(i)%STOPSOD = sod2
          endif
          cycle !inserted this data line now. Read next line
        endif
        END DO

! If no power for a particular satellite we blindly set it to
! 185 (BDS-IGSO). If no PRN for a satellite set it to X99 (invalid)
do i = 1, SAT_COUNT
!    if (satellites(i)%POWER == 0) satellites(i)%POWER = 185
    if (satellites(i)%PRN == "") satellites(i)%PRN="X99"
end do
call report('STATUS', pgrm_name, 'read_sinex_file', &
        'Finished reading Sinex POWER from file', ' ', k)

! sort array based on prn (first) then start time: since we always search by prn ...
call Qsort_Sinex_Prn(satellites(1:SAT_COUNT))
ioerr = 0

write(message, '(a, i4, a)') "Created satellite array with ", SAT_COUNT, " rows of data"
call report('STATUS', pgrm_name, 'read_sinex_file', &
        message, ' ', 0)

Do i = 1, SAT_COUNT
   write(message, '(i4,";",a," BLKTYP: ", a, " MASS: ", f7.2, " PRN: ", a, " POWER: ", i4,' //&
           '" Start: ", i5, " ", i5, " ", i5, " Stop: ", i5, " ", i5, " ", i5)') &
           i, satellites(i)%SVN, TRIM(satellites(i)%BLKTYP), satellites(i)%MASS, satellites(i)%PRN, &
           satellites(i)%POWER, satellites(i)%STARTYR, satellites(i)%STARTDOY, &
           satellites(i)%STARTSOD, satellites(i)%STOPYR, satellites(i)%STOPDOY, satellites(i)%STOPSOD
   call report ('STATUS', pgrm_name, 'read_sinex_file', &
                message, ' ', 0)
enddo

                   
! repeat for com block (com_x, com_y, com_z)


END SUBROUTINE

END

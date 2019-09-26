      Subroutine read_svsinex( lun,idir,iyr,iday,ihr,imin,gnss,isat
     .                       , satid,frqchn,antbody,blkid,svid,sbmass
     .                       , yawbias,yawrate,power)

c PURPOSE: Read an IGS Satellite metadata SINEX file, returning the values
c          needed for the SV at the input epoch

c PARAMETERS:
c       IN:   lun     : logical unit number for read               I*4
c             idir    : idir = 1 SVN to PRN, idir = -1 PRN to SVN  I*4
c             iyr     : 4-digit year                               I*4
c             iday    : 3-digit day of year                        I*4
c             ihr     : 2-digit hr                                 I*4
c             imin    : 2-digit minute                             I*4
c             gnss    : 1-character GNSS code (G R E C J I)        C*1
c             isat    : input SVN or PRN number                    I*4
c
c        OUT: satid   : output SFN or PRN number                   I*4 
c             frqchn  : frquency channel (GLONASS only)            I*4
c             antbody : antenna/body-type (rcvant_tab/ANTEX std)   C*20
c             svid    : SINEX name for satellite                   C*10 - not yet used
c             sbmass  : S/C mass in grams (if old-style will be 0) R*8
c             yawbias : Yaw bias                                   C*1
c             yawrate : maximum yaw rate of the input prn or svn   R*8 
c             power   : transmitted power in watts                 I*4

c CREATED 14 November 2017 by R. King

        
      implicit none

      INTEGER*4    SVNID, BLKID            
      character*1  gnss,gnss_tmp,yawbias 
      character*3  prn,          ! PRN string read from sinex file
     .             prn_in        ! PRN string for satellite we are looking for
                         ! or the correct one when svn found (depends on idir)                        

      character*4  svn,          ! SVN string read from sinex file
     .             svn_in        ! SVN string for satellite we are looking for
                         ! or the correct one when prn found.

      character*10 cospar_id   ! ID not used in GAMIT
      character*6  SatCat      ! Catalog number; not used in GAMIT.
      character*20 antbody_in  ! Body type read from metadata snx

      integer*4 idir     ! Direction of PRN->SVN or visa-vers
      integer*4 yr1,  yr2   ! Year read from sinex file
      integer*4 doy1, doy2  ! DOY read from sinex file
      real*8    sod1, sod2  ! Seconds of day read from sinex file
      real*8    secs        ! Needed for call to ds2hms: Added TAH 190701

      integer*8 itimdif     ! GAMIT function to return time difference in seconds?
                                                    
      character*10 svid 
      character*20 antbody 
      character*80 prog_name   
      character*128 record
      character*256 message

      integer*4 lun,isat,satid,frqchn,power,iyr,iday,ihr,imin,
     .          rcpar,ioerr,i
c        
      real*4 version, time_id, time1, time2
      real*8 sbmass,yawrate

      logical found  ! Used to indicate that blocks and SVNs have been found.
      logical debug / .false. /   ! Turn on ouput.
                    
                                                                    
c Put the requested time into an array for comparing with file entires
      time_id=iyr+iday/365.d0 +(ihr*3600.d0+imin*60.d0)/86400.d0/365.d0                  

c     times read from file have no seconds
      frqchn = 0   ! Initialize value for when GLONASS no used.

c Get the SVN/PRN correspondences  ( SATELLITE/PRN  block ) 
* Start read the sinex file looking for the start of SATELLITE/PRN  block
*     Once found, look for our specific prn or svn. 
      rewind(lun)
      found = .false.
      do while(.not.found) 
        read(lun,'(a)',iostat=ioerr) record
        if(ioerr.ne.0 ) then 
            print*, 'Failed to find SATELLITE/PRN block',ioerr
        else
            if( record(12:14).eq.'PRN' ) then
               found = .true.
            endif 
        endif 
      enddo  

* MOD TAH 190627: Don't check direction yet we will do that below as
*     the PRN or SVN is found.  We need to find the SVN number because
*     most the file is based on SVN not (non-unique) PRN. 
C     if( idir.eq.-1 ) then
c       PRN to SVN
*     Depending on direction either assign PRN straing (G04) or the 
*     svn string (G074) so we can search the records.
*     Assuming here isat in input PRN or SVN and satid is output PRN or SCN
      if( idir.eq.-1 ) then  ! PRN to SVN
         write(prn_in,'(a1,I2.2)') gnss, isat
      else                   ! SVN to PRN 
         write(svn_in,'(a1,I3.3)') gnss, isat
      endif

*     Now read the file to see if can find the entry 
      found = .false.
      rewind(lun)
      do while ( .not.found ) 
        read(lun,'(a)',iostat=ioerr) record
        if(record(1:3).eq.'-SAT') then !.or.ioerr.ne.0 ) then 
           print*,'Failed to find SV in PRN block of SINEX file',ioerr    
        else
*          Only try to decode record if it is not a comment
           if(record(1:1).eq.' ' ) then
              read(record,'(1x,a4,2(1x,i4,1x,i3,1x,f5.0),1x,a3)'
     .            ,iostat=ioerr) svn,yr1,doy1,sod1,yr2,doy2,sod2,prn 

c              If( ioerr.ne.0 ) 
c     .        print*,'Error decoding PRN block of SINEX file',ioerr
              if( yr2.eq.0000 ) then 
                  yr2 = 2100
                  doy2 = 365
                  sod2 = 86400
              time2 = yr2+doy2/365.d0+sod2/86400.d0/365.d0
              endif
              time1 = yr1+doy1/365.d0+sod1/86400.d0/365.d0 
              time2 = yr2+doy2/365.d0+sod2/86400.d0/365.d0
*             See if have a match to the prn_in or svn_in we are
*             looking for
              if( idir.eq.-1 ) then     ! PRN_in passed, see if match
c                PRN to SVN: Compare prn string from file with prn_in
                 if(prn.eq.prn_in) then
                    if(time_id>=time1 .and. time_id<=time2) then
                       found = .true. 
                       read(svn(2:4),'(i3)') satid
                       svn_in = svn   ! Save so that we use to find other 
                                      ! entries such as power, mass and frequencies.
                    endif
                 endif  
              elseif( idir.eq.1 ) then
c               SVN to PRN: Compare svn string from file with svn_in
                if(svn.eq.svn_in) then
                   if(time_id>=time1 .and. time_id<=time2) then
                      found = .true. 
                      read(prn(2:3),'(i2)') satid
                      prn_in = prn    ! Save in case we need.
                   endif
                endif     
              else 
                If( ioerr.ne.0 )   
     .          print*,'idir not -1 or 1 in call to read_svsinex',ioerr
              endif
           end if
        endif
      enddo
         

c Get the frequency channel for Glonass ( SATELLITE/FREQUENCY_CHANNEL block)
*     This only needed if gnss == R for Glonass.
      if( gnss(1:1).eq.'R' ) then 
         rewind(lun)
         found = .false.
         do while(.not.found) 
           read(lun,'(a)',iostat=ioerr) record              
           if( ioerr.ne.0 ) then
              print*,
     .           'Failed to find FREQUENCY_CHANNEL SINEX block',ioerr
           elseif( record(12:15).eq.'FREQ' ) then
             found = .true.
           endif          
         enddo
         found = .false.                                 
         do while(.not.found) 
            read(lun,'(a)',iostat=ioerr) record
            if( record(1:3).eq.'-SAT'.or.ioerr.ne.0 ) then 
                If( ioerr.ne.0 )   
     .          print*,'Failed to find SV in FREQUENCY_CHANNEL 
     .                  SINEX block',ioerr
            else             
*           Only try to decode record if it is not a comment
               if( record(1:1).eq.' ' ) then
                 read(record,'(1x,a4,2(1x,i4,1x,i3,1x,f5.0),1x,i4)',
     .              iostat=ioerr) svn,yr1,doy1,sod1,
     .                                yr2,doy2,sod2, frqchn
                 if(ioerr.ne.0) 
     .           print*,'Reading GLONASS frequency',ioerr                   
                 if( yr2.eq.0000 ) then 
                    yr2 = 2100
                    doy2 = 365
                    sod2 = 86400
                 time2 = yr2+doy2/365.d0+sod2/86400.d0/365.d0
                 endif 
                 time1 = yr1+doy1/365.d0+sod1/86400.d0/365.d0
                 time2 = yr2+doy2/365.d0+sod2/86400.d0/365.d0
                 if(svn.eq.svn_in) then 
                    if(time_id>=time1 .and. time_id<=time2) then
                     found = .true. 
                    endif
                  endif  
               endif
            endif
         enddo
      endif         ! gnss = R

c Get the SATELLITE/IDENTIFIER block
c     +SATELLITE/IDENTIFIER                                                           
c     *                                                                               
c     *SVN_ COSPAR ID SatCat Block__________ Comment__________________________________
c     *                                                                               
c      G001 1978-020A  10684 GPS-I           Launched 1978-02-22; NAVSTAR 1
c      G002 1978-047A  10893 GPS-I           Launched 1978-05-13; NAVSTAR 2
c      G003 1978-093A  11054 GPS-I           Launched 1978-10-07; NAVSTAR 3

      rewind(lun)
      found = .false.
      do while(.not.found) 
        read(lun,'(a)',iostat=ioerr) record   
        if(ioerr.ne.0 ) then
           print*,'Failed to find IDENTIFIER SINEX block',ioerr
        elseif( record(12:21).eq.'IDENTIFIER' ) then
          found = .true.
        endif          
      enddo
      found = .false.                                 
      do while(.not.found) 
        read(lun,'(a)',iostat=ioerr) record
        if( record(1:3).eq.'-SAT'.or.ioerr.ne.0 ) then 
           If( ioerr.ne.0 )
     .     print*,'Failed to find SV in IDENTIFIER SINEX block',ioerr
        else             
*         Only try to decode record if it is not a comment
          if( record(1:1).eq.' ' ) then
             read(record,'(1x,a4,1x,a9,1x,A6,1x,a15)'
     .           ,iostat=ioerr) svn, cospar_id,SatCat, antbody_in
             if(ioerr.ne.0) 
     .       print*,'Fail to reading satellite id ',ioerr 
*            Only one line per satellite (unique ID)  
             if( svn.eq.svn_in ) then
                 found = .true. 
                 antbody = antbody_in
             endif  
          endif
        endif 
      enddo
             

c Get the satellite mass ( SATELLITE/MASS block )

      rewind(lun)
      found = .false.
      do while(.not.found) 
        read(lun,'(a)',iostat=ioerr) record   
        if(ioerr.ne.0 ) then
           print*,'Failed to find MASS SINEX block',ioerr
        elseif( record(12:15).eq.'MASS' ) then
          found = .true.
        endif          
      enddo
      found = .false.                                 
      do while(.not.found) 
        read(lun,'(a)',iostat=ioerr) record
        if( record(1:3).eq.'-SAT'.or.ioerr.ne.0 ) then 
           If( ioerr.ne.0 )   
     .     print*,'Failed to find SV in MASS SINEX block',ioerr
        else             
*         Only try to decode record if it is not a comment
          if( record(1:1).eq.' ' ) then
             read(record,'(1x,a4,2(1x,i4,1x,i3,1x,f5.0),1x,f9.3)'
     .           ,iostat=ioerr) svn,yr1,doy1,sod1,yr2,doy2,sod2,sbmass 
             if(ioerr.ne.0) 
     .       print*,'Reading satellite mass ',ioerr                   
             if( yr2.eq.0000 ) then 
                 yr2 = 2100
                 doy2 = 365
                 sod2 = 86400
              time2 = yr2+doy2/365.d0+sod2/86400.d0/365.d0
             endif
              time1 = yr1+doy1/365.d0+sod1/86400.d0/365.d0
              time2 = yr2+doy2/365.d0+sod2/86400.d0/365.d0
             if(svn.eq.svn_in) then 
                if(time_id>=time1 .and. time_id<=time2) then
                   sbmass = sbmass 
                   found = .true.
                endif
             endif  
          endif
        endif 
      enddo

c Get the yaw bias and rate.
* MOD TAH 190701: Implemented with new SATELLITE/YAW_BIAS_RATE block added
*     but TAH to igs_metadata.snx (using script sh_svnav_yaw_to_igsmeta)

c Switch off the following block. When the yaw rate and yaw bias
c information is ready in the metadat file, then the block should be
c turn on.(09-09-2019 Tzuapng Tseng)

c      rewind(lun)
c      found = .false.
c      do while(.not.found) 
c        read(lun,'(a)',iostat=ioerr) record   
c        if( record(1:3).eq.'-SAT'.or.ioerr.ne.0 ) then 
c           if(ioerr.ne.0 ) then
c           print*,'Failed to find YAW_BIAS_RATE SINEX block',ioerr
c           endif 
c        elseif( record(12:19).eq.'YAW_BIAS' ) then
c          found = .true.
c        endif          
c      enddo
c      found = .false.                                 
c      do while(.not.found) 
c        read(lun,'(a)',iostat=ioerr) record
c        if( record(1:3).eq.'-SAT'.or.ioerr.ne.0 ) then 
c          If( ioerr.ne.0 )  
c     .    print*,'Failed to find SV in YAW block',ioerr
c        else             
*         Only try to decode record if it is not a comment
c          if( record(1:1).eq.' ' ) then
c             read(record,'(1x,a4,2(1x,i4,1x,i3,1x,F5.0),4x,a1,1x,f8.4)'
c     .           ,iostat=ioerr) svn,yr1,doy1,sod1,yr2,doy2,sod2,
c     .                          yawbias, yawrate
c             if(ioerr.ne.0)  print*,'Fail to reading yaw rate ',ioerr   
c             if( yr2.eq.0000 ) then 
c                 yr2 = 2100
c                 doy2 = 365
c                 sod2 = 86400
c              time2 = yr2+doy2/365+sod2/86400/365
c             endif
c              time1 = yr1+doy1/365+sod1/86400/365
c              time2 = yr2+doy2/365+sod2/86400/365
c             if(svn.eq.svn_in) then
c               if(time_id>=time1 .and. time_id<=time2) then
c                 found = .true. 
c               endif
c             endif  
c           endif
c         end if
c      enddo


c Get the transmitter power ( SATELLITE/TX_POWER block)

      rewind(lun)
      found = .false.
      do while(.not.found) 
        read(lun,'(a)',iostat=ioerr) record 
        if(ioerr.ne.0 ) then
        print*,'svnav.dat','Failed to find TX_POWER SINEX block',ioerr
        elseif( record(12:15).eq.'TX_P' ) then
          found = .true.
        endif 
      enddo   
      found = .false.                                 
      do while(.not.found) 
        read(lun,'(a)',iostat=ioerr) record 
        if( record(1:3).eq.'-SAT'.or.ioerr.ne.0 ) then 
           If( ioerr.ne.0 )  
     .     print*,'Failed to find SV in TX_POWER block',ioerr
        else             
*          Only try to decode record if it is not a comment
           if( record(1:1).eq.' ' ) then
              read(record,'(1x,a4,2(1x,i4,1x,i3,1x,F5.0),1x,i4)'
     .           ,iostat=ioerr) svn,yr1,doy1,sod1,yr2,doy2,sod2,power  
              if( yr2.eq.0000 ) then 
                  yr2 = 2100
                  doy2 = 365
                  sod2 = 86400
              time2 = yr2+doy2/365.d0+sod2/86400.d0/365.d0
              endif
              if(ioerr.ne.0) print*,'Check metadata file'
              time1 = yr1+doy1/365.d0+sod1/86400.d0/365.d0
              time2 = yr2+doy2/365.d0+sod2/86400.d0/365.d0
              if(svn.eq.svn_in) then
                if(time_id>=time1 .and. time_id<=time2) then
                  found = .true.
                endif
              endif  
           endif
         endif 
      enddo

      if( debug ) then
         write(*,*) svn, prn, antbody, iyr, iday, frqchn, sbmass, 
     .                power, yawbias, yawrate 
      endif
      SVNID = satid
      IF(antbody=='GPS-I')      BLKID = 1
      IF(antbody=='GPS-II')     BLKID = 2
      IF(antbody=='GPS-IIA')    BLKID = 3
      IF(antbody=='GPS-IIR')    BLKID = 4
      IF(antbody=='GPS-IIR-A')  BLKID = 5
      IF(antbody=='GPS-IIR-B')  BLKID = 6
      IF(antbody=='GPS-IIR-M')  BLKID = 7
      IF(antbody=='GPS-IIF')    BLKID = 8
      IF(antbody=='GPS-IIIA')   BLKID = 9
      IF(antbody=='GLO')        BLKID = 101
      IF(antbody=='GLO-M'  .or.antbody == 'GLO-M+')  BLKID = 102
      IF(antbody=='GLO-K1A'.or.antbody == 'GLO-K1B') BLKID = 103
      IF(antbody=='GLA-1')     BLKID = 201 ! Galileo (IOV)
      IF(antbody=='GLA-2')     BLKID = 202 ! Galileo (FOC)

      return
      end


subroutine read_cmdline

! This subroutine reads commnad line options and returns variables defined in mdl_config
! Uses getopt command line parsing functions from f90getopt.F90

USE f90getopt
USE mdl_config

! ----------------------------------------------------------------------
! Command line local variables
integer            :: len_optarg
character (LEN=80) :: pgm_name

! Set number of long command line options available
type(option_s) :: opts(15)

! Current mdl_config varaible options
! ----------------------------------------------------------------------
! POD_fname_cfg                -c
! POD_MODE_cfg                 -m
! EQM_fname_cfg                -e
! VEQ_fname_cfg                -v
! pseudobs_orbit_filename_cfg  -s
! ext_orbit_filename_cfg       -o
! orbit_determination_arc_cfg  -a
! orbit_prediction_arc_cfg     -p
! EOP_solution_cfg             -t
! EOP_fname_cfg                -s
!! ERP_fname_cfg               -- not implemented on cmdline 
!! EOP_Nint_cfg                -- not implemented on cmdline 
! iau_model_cfg                -n
! Estimator_Iterations_cfg     -i
! sp3_velocity_cfg             -u
! IC_MODE_cfg                  -q
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! Read the Command line
! ----------------------------------------------------------------------
! Define command line options
!                    long_opt    argument  short_opt
opts(1)  = option_s( "config",   .true.,      'c' )
opts(2)  = option_s( "podmode",  .true.,      'm' )
opts(3)  = option_s( "pobs",     .true.,      's' )
opts(4)  = option_s( "cobs",     .true.,      'o' )
opts(5)  = option_s( "eqm",      .true.,      'e' )
opts(6)  = option_s( "veq",      .true.,      'v' )
opts(7)  = option_s( "arclen",   .true.,      'a' )
opts(8)  = option_s( "predlen",  .true.,      'p' )
opts(9)  = option_s( "eopf",     .true.,      'r' )
opts(10) = option_s( "eopsol",   .true.,      't' )
opts(11) = option_s( "nutpre",   .true.,      'n' )
opts(12) = option_s( "estiter",  .true.,      'i' )
opts(13) = option_s( "sp3vel",   .false.,     'u' )
opts(14) = option_s( "icmodel",  .true.,      'q' )
opts(15) = option_s( "help",     .false.,     'h' )

! Get the program name
call get_command_argument( 0, pgm_name )

! Read Command line from beginning
optind=1

! If no comand line options given provide some help [By default run the POD with Default POD.in file available]
!if (command_argument_count() .eq. 0 ) then
!   print*, trim(pgm_name),' -h or ',trim(pgm_name),' --help for command line help'
!   stop
!end if

! Set master configuration file name to be the default defined name (set in main_pod.f03)
POD_fname_cfg = 'DEFAULT'

! Process options given sequentially
do
   select case(getopt("c:m:s:o:e:v:a:p:r:t:n:i:u:q:h",opts))
      case( char(0) )
         exit
      case( 'c' )
!         print *, 'option config/c=', optarg
   	      POD_fname_cfg = trim(optarg)
      case( 'm' )
!         print *, 'option podmode/m=', optarg
          len_optarg = len_trim(optarg)
          read(optarg(1:len_optarg),'(i2)')	POD_MODE_cfg
      case( 's' )
!          print *, 'option pobs/s=', optarg
	      pseudobs_orbit_filename_cfg = trim(optarg)
      case( 'o' )
!         print *, 'option cobs/o=', optarg
	      ext_orbit_filename_cfg = trim(optarg)
      case( 'e' )
!         print *, 'option eqm/e=', optarg
	      EQM_fname_cfg = trim(optarg)
      case( 'v' )
!         print *, 'option veq/v=', optarg
	      VEQ_fname_cfg = trim(optarg)
      case( 'a' )
!         print *, 'option arclen/a=', optarg
          len_optarg = len_trim(optarg)
          read(optarg(1:len_optarg),'(f14.6)') orbit_determination_arc_cfg
      case( 'p' )
!         print *, 'option predlen/p=', optarg
          len_optarg = len_trim(optarg)
          read(optarg(1:len_optarg),'(f14.6)') orbit_prediction_arc_cfg
      case( 'r' )
!         print *, 'option eopf/r=', optarg
	      EOP_fname_cfg = trim(optarg)
      case( 't' )
!         print *, 'option eopsol/t=', optarg
          len_optarg = len_trim(optarg)
          read(optarg(1:len_optarg),'(i3)') EOP_solution_cfg
      case( 'n' )
!      print *, 'option nutpre/n=', optarg
          len_optarg = len_trim(optarg)
          read(optarg(1:len_optarg),'(i4)') iau_model_cfg
      case( 'i' )
!         print *, 'option estiter/i=', optarg
          len_optarg = len_trim(optarg)
          read(optarg(1:len_optarg),'(i2)') Estimator_Iterations_cfg
      case( 'u' )
!      print *, 'option sp3vel/u=', optarg
          len_optarg = len_trim(optarg)
          read(optarg(1:len_optarg),'(i4)') sp3_velocity_cfg
      case( 'q' )
!      print *, 'option icmode/u=', optarg
          len_optarg = len_trim(optarg)
          read(optarg(1:len_optarg),'(i4)') IC_MODE_cfg
      case( 'h' )
          print*,'Default master POD config file = POD.in'
		  print*,'To run from default config file: ',trim(pgm_name),' or ',trim(pgm_name),' -c POD.in'
          print*,''
		  print*,'POD.in config file options by defaut can be overridden on the command line'
          print*,''
          print*,'Command line: ',trim(pgm_name),' -c -m -s -o -e -v -a -p -r -t -n -i -u -h '
          print*,''
          print*,'Where: '
          print*,'      -c --config  = Config file name [Default POD.config]'
          print*,'      -m --podmode = POD Mode:'
          print*,'				1 - Orbit Determination (pseudo-observations; orbit fitting)'
          print*,'				2 - Orbit Determination and Prediction'
          print*,'				3 - Orbit Integration (Equation of Motion only)'
          print*,'				4 - Orbit Integration and Partials (Equation of Motion and Variational Equations)'
          print*,'      -s --pobs    = Pseudo observations orbit .sp3 file name'
          print*,'      -o --cobs    = Comparison orbit .sp3 file name'
          print*,'      -e --eqm     = EQuations of Motion input file name  [Default: EQM.in]'
          print*,'      -v --veq     = Variatinal EQuations input file name [Default: VEQ.in]'
          print*,'      -a --arclen  = Orbit Estimation Arc length (hours)'
          print*,'      -p --predlen = Orbit Prediction Arc length (hours)'
	      print*,'      -r --eopf    = Earth Orientation Paramaeter (EOP) values file'
          print*,'      -t --eopsol  = Earth Orientation Paramaeter file type: (1,2)'
          print*,'				1 - IERS C04 EOP'
          print*,'				2 - IERS RS/PC Daily EOP'
          print*,'				3 - IGS RP + IERS RS/PC Daily (dX,dY)'
          print*,'      -n --nutpre  = IAU Precession / Nutation model'
          print*,'				2000 - IAU2000A'
          print*,'				2006 - IAU2006/2000A'
          print*,'      -i --estiter = Orbit Estimatimation Iterations (1 or greater)'
          print*,'      -u --sp3vel  = Output .sp3 file with velocities'
          print*,'      -q --icmode  = Initial condition from parameter estimation procedure'
		  print*,'				0 - Do not write Velocity vector to sp3 orbit'
		  print*,'				1 - Write Velocity vector to sp3 orbit'  
          print*,'      -h --help.   = Print program help'
          print*,''
          print*,'Examples:'
          print*,''
          print*,'       ',trim(pgm_name),' -m 1 -s igs16403.sp3 '		  
          print*,'       ',trim(pgm_name),' -m 2 -s igs16403.sp3 -e EQMx.in -v VEQx.in -p 12'
          stop
   end select
end do

return
end

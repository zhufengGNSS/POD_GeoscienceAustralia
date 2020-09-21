# Brief introduction of precise orbit determination 

The precise orbit determination (POD) provides a refernce orbit of GNSS satellite for the GNSS measurement.
The reference orbit is precisely modeled with all considered forces. 
In the following, we list the observation type, earth orientation parameters, 
orbital force models and the numerical integration method, all of which are used for POD at Geoscience Australia:

---------------------------------------------------------------------------
Observation tupe: pseudo observation from SP3 file

----------------------------------------------------------------------------
Earth Orientation: IERS C04 product
 
Precession-Nutation model: IAU 2000

----------------------------------------------------------------------------
Gravitational Forces
 
 -geopotential model: goco05s.gfc, up to degree/order 15	
 
 -Planetary effect  : JPL DE430	  
 
 -Tidal effects     : solid tides non-frequency terms
                      solid tides frequency terms				 
                      solid earth pole tide	
                      
 -Ocean tides model : FES2004 ocean tide model, up to degree/order 15
                      ocean pole tide

 -Relativistic effect 
 
---------------------------------------------------------------------------
Non-gravitational Effects
 -Solar radiation : ECOM-based model	
 
 -Earth radiation : Analytical model 		  
 
 -Antenna thrust  : Simple model
 
---------------------------------------------------------------------------
Numerical integration methods: Runge-Kutta-Nystrom 7th order 

 

# Installtion
# Analysis Centre Software - POD

`ACS POD Version 1.1.0 beta release`

### Supported Platforms
The POD is supported on the following Platforms

* Linux
* Mac OSX
  
### Dependencies
The POD has several software dependencies

* C/C++ and Fortran compiler. We use and recommend [gcc-g++ and gfortran](https://gcc.gnu.org/git.html)
* BLAS and [LAPACK](https://github.com/Reference-LAPACK/lapack) linear algebra libraries. We use and recommend [OpenBlas](https://www.openblas.net/)
* Cmake Version3 
* Python Version3 (including Numpy and Matplotlib modules)

### Download
To downlaod ACS POD source from bitbucket, you need to first install [Git](https://www.atlassian.com/git) and [Git-LFS](https://git-lfs.github.com/)

Then you can download the POD source code using git clone:

    $git clone git@bitbucket.org:geoscienceaustralia/pod.git

## Directory Structure
    pod/
    ├── README.md			! General README information
    ├── LICENSE.md		    ! Software License information
    ├── aws			        ! Amazon Web Services config
    ├── bin			        ! Binary executables directory
    ├── build			    ! Cmake build directory
    ├── CMakeLists.txt		! Cmake build file
    ├── config			    ! Default/Master configuration files directoy 
    ├── docs			    ! Documentation directory
    ├── examples ├          ! POD examples directory
                 ├── EX1    ! POD example 1
                 ├── EX2    ! POD example 2
    ├── lib			        ! Compiled objectlibrary directory 
    ├── scripts			    ! Auxillary Python and Shell scripts 
    ├── src			        ! Source code directoy
    ├── tables			    ! Metadata and auxillary tables directory
    ├── test                ! POD test directory


### Build
To build the `POD` ...

    $ cd pod
    $ mkdir build
    $ cd build
    $ cmake3 .. >cmake.out 2>cmake.err
    $ make >make.out 2>make.err
    $ less make.err (to verify everything was built correctly)

You should now have the executables in the bin directory: pod crs2trs brdc2ecef

### User Environment (Shell) setup

In you home directory create a link to the pod installation directory ...

    $ cd ~
    $ ln -s /data/software/acs/pod/install/directory pod
    
Add the pod installation executable diretory (~/pod/bin) to your execuatble search path.

Bash shell example:
    
    export PATH=~/pod/bin:$PATH

### Test 

To test your build of the  `POD` ...

    $ cd ~/pod/test
    $ sh_test_pod

At the completion of the test run, the sh_test_pod script will return any differences to the standard test resuts

---------------------------------------------------------------------------

# Testing examples
# Analysis Centre Software - POD

## POD

The `ACS` Version 1.1.0 beta release supports:

## POD Examples
In the following POD examples, the command line options are provided in controlling observation source, 
earth rotation parameter source, solar radiation pressure model, empirical model and so on.

The options in the major configuration file 'POD.in' can be overridden on the command line:

 Command line: ./pod -c -m -s -o -e -v -a -p -r -t -n -i -u -q -k -w -h

 Where:
       -c --config  = Config file name [Default POD.config]
       -m --podmode = POD Mode:
                                1 - Orbit Determination (pseudo-observations; orbit fitting)
                                2 - Orbit Determination and Prediction
                                3 - Orbit Integration (Equation of Motion only)
                                4 - Orbit Integration and Partials (Equation of Motion and Variational Equations)
       -s --pobs    = Pseudo observations orbit .sp3 file name
       -o --cobs    = Comparison orbit .sp3 file name
       -e --eqm     = EQuations of Motion input file name  [Default: EQM.in]
       -v --veq     = Variatinal EQuations input file name [Default: VEQ.in]
       -a --arclen  = Orbit Estimation Arc length (hours)
       -p --predlen = Orbit Prediction Arc length (hours)
       -r --eopf    = Earth Orientation Paramaeter (EOP) values file
       -t --eopsol  = Earth Orientation Paramaeter file type: (1,2)
                                1 - IERS C04 EOP
                                2 - IERS RS/PC Daily EOP
                                3 - IGS RP + IERS RS/PC Daily (dX,dY)
       -n --nutpre  = IAU Precession / Nutation model
                                2000 - IAU2000A
                                2006 - IAU2006/2000A
       -i --estiter = Orbit Estimatimation Iterations (1 or greater)
       -u --sp3vel  = Output .sp3 file with velocities
       -q --icmode  = Initial condition from parameter estimation procedure
                                0 - Do not write Velocity vector to sp3 orbit
                                1 - Write Velocity vector to sp3 orbit
       -k --srpmodel= 1: ECOM1, 2:ECOM2, 12: ECOM1+ECOM2, 3:SBOX
       -w --empmodel= 1: activated, 0: no estimation
       -d --verbosity = output verbosity level [Default: 0]
       -h --help.   = Print program help


Examples:
        ./pod -m 1 -q 1 -k 1 -w 0 -s igs16403.sp3 -o igs16403.sp3
        ./pod -m 2 -q 1 -k 1 -w 0 -s igs16403.sp3 -e EQMx.in -v VEQx.in -p 12

 For orbit updates using Parameter Estimation Algorithm (PEA):
        ./pod -m 4 -q 2 -k 1 -w 0 -s igs16403.sp3 -o igs16403.sp3

---------------------------------------------------------------------------

## Broadcast example

Conversion from broadcast dynamic elements to earth-center earth-fixed (ecef) coordinates
The program can work for multi-GNSS constellations, except for GLONASS that will be handled soon.
 To run main program : ./brdc2ecef

 Command line: ./brdc2ecef -f -o -c -h

 Where:
       -f --inputfile  = input file name
       -o --outputfile = output file name
       -c --constellation = GNSS TYPE, e.g., A: All GNSS; G: GPS; R: GLONASS
                                             E: GALILEO;  C: BDS; J: QZSS
       -l --leap second file = leap.second
       -r --EOP file
       -h --help.   = Print program help

 Examples:

        ./brdc2ecef -f BRDC00IGS_R_20163210000_01D_MN.rnx -o brdc16321.sp3 -c A -l leap.second -r finals2000A.daily

-----------------------------------------------------------------------------

In the pod/examples directory, there are six POD user exmples: 

## Example 1 - (pod/examples/ex1):

* Simple GPS IGS final SP3 file orbit fitting with ECOM1 (same as pod/test example)

* Command line: ./pod -m 1 -q 1 -k 1 -w 0 -s igs16403.sp3 -o igs16403.sp3

## Example 2 - (pod/examples/ex2):

* Multi GNSS CODE MGEX SP3 file orbit fitting with ECOM1

* Command line: ./pod -m 1 -q 1 -k 1 -w 0 -s COD0MGXFIN_20191990000_01D_05M_ORB.SP3 -o COD0MGXFIN_20191990000_01D_05M_ORB.SP3
## Example 3 - (pod/examples/ex3):

* GPS IGS SP3 file orbit fitting, orbit prediction and comparison to next IGS SP3 file

* Command line: /pod -m 2 -q 1 -k 1 -w 0 -s igs20000.sp3 -o igs20001.sp3

## Example 4 - (pod/examples/ex4):

* Integration of POD initial conditions file generated by the PEA

* Please make sure whether the a priori SRP option is consistent.
* Command line: ./pod -m 4 -q 2 -k 1 -w 0 -o igs20624.sp3

## Example 5 - (pod/examples/ex5):

* ECOM1+ECOM2 hybrid model in POD

* Command line: /pod -m 1 -q 1 -k 12 -w 1 -s igs20000.sp3 -o igs20000.sp3

## Example 6 - (pod/examples/ex6):

* SP3 orbital file generated from broadcast ephemeris 

* Command line: ulimit -s unlimited
* Command line: ./brdc2ecef -f BRDC00IGS_R_20163210000_01D_MN.rnx -o brdc16321.sp3 -c A -l leap.second -r finals2000A.daily

In each example directory (ex1/ex2/ex3/ex4) there is a sh_ex? script that when exectuted 
will run the example and compare the output with the expected solution. 

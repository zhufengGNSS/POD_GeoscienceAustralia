# Analysis Centre Software - POD

The `ACS` Version 1.1.0 beta release supports:

1. The `POD` 

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

### Supported Platforms
The POD is supported on the following Platforms:
  Linux
  Mac OSX
  
### Dependencies
The POD has several software dependencies:

  C/C++ and Fortran compiler. We use and recommend [gcc-g++ and gfortran](https://gcc.gnu.org/git.html)
  BLAS and [LAPACK](https://github.com/Reference-LAPACK/lapack) linear algebra libraries. We use and recommend [OpenBlas](https://www.openblas.net/)
  Cmake version3 
  Python version3 (including Numpy and Matplotlib modules) 

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

In you home directory create a link to the pod installation directory
    $ cd ~
    $ ln -s /data/software/acs/pod/install/directory pod
    
And/Or
    Add the pod installation executable diretory (~/pod/bin) to your execuatble search path.

    Bash shell example:
    
      export PATH=~/pod/bin:$PATH

### Test 

To test your build of the  `POD` ...

    $ cd ~/pod/test
    $ sh_test_pod

At the completion of the test run, the sh_test_pod script will return any differences to the standard test resuts

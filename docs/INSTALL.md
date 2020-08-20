# Analysis Centre Software - POD

`ACS POD Version 1.1.0 beta release`

### Supported Platforms

The POD is supported on the following Platforms

* Linux
* Mac OSX
  
### Dependencies

The POD has several software dependencies:

* C/C++ and Fortran compiler. We use and recommend [gcc-g++ and gfortran](https://gcc.gnu.org/git.html)
* BLAS and [LAPACK](https://github.com/Reference-LAPACK/lapack) linear algebra libraries. We use and recommend [OpenBlas](https://www.openblas.net/)
* Cmake Version3 
* Python Version3 (including Numpy and Matplotlib modules)

### Installing dependencies with Ubuntu

    $ sudo apt -y install git gobjc gobjc++ gfortran openssl curl net-tools openssh-server cmake make gzip vim libssl1.0-dev
    $ sudo apt -y install libopenblas-dev git-lfs
    $ sudo apt -y install python3-matplotlib 

### Download

To downlaod ACS POD source from bitbucket, you need to first install [Git](https://www.atlassian.com/git) and [Git-LFS](https://git-lfs.github.com/)

Then you can download the POD source code using git clone:

    $ git clone git@bitbucket.org:geoscienceaustralia/pod.git
    
Then download all of the example data using git-lfs:

    $ git-lfs pull 

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
    $ cmake .. >cmake.out 2>cmake.err
    $ make >make.out 2>make.err
    $ less make.err (to verify everything was built correctly)

You should now have the executables in the bin directory: 
    ├── bin			        

    * pod 
    * crs2trs 
    * brdc2ecef

### User Environment (Shell) setup

In you home directory create a link to the pod installation directory ...

    $ cd ~
    $ ln -s /data/software/acs/pod/install/directory pod
    
Add the pod installation executable diretory (~/pod/bin) to your execuatble search path.

Bash shell example:
    
    export PATH=~/pod/bin:$PATH

### Test 

To test your build of the  `POD` ...

    $ cd ~/pod
    $python ./scripts/download_tables.py
    $python ./scripts/download_examples.py -p ./examples
    $cd examples/ex1
    $chmod oug+x sh_ex1
    $sh_ex1

At the completion of the test run, the sh_ex1 script will return any differences to the standard test resuts


62 directories, 1135 files

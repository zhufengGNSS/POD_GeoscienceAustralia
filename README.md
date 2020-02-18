# Analysis Centre Software - POD

## Overview

The *Analysis Centre Software (ACS)* is a processing package being developed to processes GNSS observations for geodetic 
applications.  

We currently support the processing of:

* the American Global Positioning System (`GPS`)
* the Russian GLONASS system ('GLONASS')
* The European Gallileo system ('Gallileo')
* the Chinese Navigation Satellite System ('Beidou`); and
* The Japanese QZSS develop system ('QZSS')

We are actively developing the ACS to have the following capabilities and features:

* Precise Orbit & Clock determination of GNSS satellites (GNSS POD).
* Precise Point Positioning (PPP) of GNSS stations in network and individual mode.
* Real-Time corrections for PPP users.
* Analyse full, single and multi-frequency, multi-GNSS data.
* Delivering atmospheric products such as ionosphere and troposphere models.
* Servicing a wide range of users and receiver types.
* Delivering outputs usable and accessible by non-experts.
* Providing both a real-time and off-line processing capability.
* Delivering both position and integrity information.
* Routinely produce IGS final, rapid, ultra-rapid and real-time (RT) products. 

The software is broken into two main components:

* the Network Parameter Estimation Algorithm (`PEA-N`); and 
* the Precise Orbit Determination (`POD`).

## POD

The `ACS` Version 0.0.1 beta release supports:

1. The `POD` 

## Directory Structure

    pod/
    ├── README.md

### Dependencies

1. The lapack numerical linear algebra library (lapack.x86_64) (You may need to run the command ln -s /usr/lib64/liblapack.so.3 /usr/lib64/liblapack.so)
2. The basic linear algebra library (blas.x86_64,liblas-libs.x86_64) (You may need to run the command ln -s /usr/lib64/libblas.so.3 /usr/lib64/libblas.so)
3. A working C compiler (gcc will do), a working C++ compiler (gcc-g++ will do) and a fortran compiler (we have used gfortran)
4. If the flags set in CMakeLists.txt do not work with your compiler please remove the ones that don't

### Build

To build the `POD` ...

    $ cd pod
    $ mkdir build
    $ cd build
    $ cmake3 .. >cmake.out 2>cmake.err
    $ make >make.out 2>make.err
    $ less make.err (to verify everything was built correctly)

You should now have the executables in the bin directory: pod crs2trs brdc2ecef


### Configuration File

The `POD` Precise Orbit Determination (`./bin/pod`) uses the configuration file:
    ├── EQM.in (Full force model equation of motion)
    ├── VEQ.in (For variational equations)
    ├── POD.in (For all other config)


### Processing Example #1

In this example the pod will perform a dynamic orbit determination for PRN04 over a 6 hour arc. The full gravitational force models are applied, with a cannonball model SRP model.

  
To run the `POD` ...

    $ bin/pod

This should output the following to `stdout`...

    Orbit Determination
    Orbit residuals in ICRF : RMS(XYZ)   1.6754034501980351E-002   5.2908718335411935E-002   1.5676115599034774E-002
    Orbit Determination: Completed
    CPU Time (sec)   298.48134399999998
    External Orbit comparison
    Orbit comparison: ICRF
    RMS RTN   2.8094479714173427E-002   2.4358145601708528E-002   4.4097979280889953E-002
    RMS XYZ   1.6754034501980351E-002   5.2908718335411935E-002   1.5676115599034774E-002
    Orbit comparison: ITRF
    RMS XYZ   3.9069978513805753E-002   3.9343671258381237E-002   1.5660654272651970E-002
    Write orbit matrices to output files
    CPU Time (sec)   349.19307899999995

The results above show that our orbits arcs, over 6 hours, are currently within 2-5 cm of the final combined IGS orbit. 

The prcessing also produces the following output files...

    ├── DE.430            planetary ephemris intermediate file
    ├── Amatrix.out       design matrix
    ├── Wmatrix.out       reduced observation matrix
    ├── orbext_ICRF.out   intermediary file for the IGS orbit solution in ICRF for comparison purposes
    ├── orbext_ITRF.out   intermediary file for the IGS orbit solution in ITRFfor comparison purposes
    ├── dorb_icrf.out     differences in solutions in ICRF
    ├── dorb_RTN.out      differences in solutions in orbital frame components radial, tangential and normal (RTN)
    ├── dorb_Kepler.out   differences in solutions in keperian elements 
    ├── dorb_itrf.out     differences in solutions in ITRF 
    ├── orb_icrf.out      the final estimated orbit in ICRF
    ├── orb_itrf.out      the final estimated orbit in ITRF
    ├── VEQ_Smatrix.out   State transition matrix from the variational equations solution
    ├── VEQ_Pmatrix.out   Sensitivity matrix from the variational equations solution


### Processing Example #2 - ECOM2 SRP

In this example we will change the SRP model to use the ECOM2 model. 

Edit the EQM.in file so that the Solar Radiation Pressure configuration section now looks:

! Solar Radiation Pressure model:
! 1. Cannonball model
! 2. Box-wing model
! 3. ECOM (D2B1) model
SRP_model 3

Then edit VEQ.in, so that the Non-gravitational forces now looks like:

%% Non-gravitational Effects
Solar_radiation           0
Earth_radiation           0
Antenna_thrust            0

! Solar Radiation Pressure model:
! 1. Cannonball model
! 2. Box-wing model
! 3. ECOM (D2B1) model
SRP_model                 3

run the `POD` ...

    $ bin/pod

This should output the following to `stdout`...

    Orbit Determination
    Orbit residuals in ICRF : RMS(XYZ)   2.0336204859568077E-002   8.4715644601919167E-003   3.9687932322714677E-002
    Orbit Determination: Completed
    CPU Time (sec)   299.68054799999999
    External Orbit comparison
    Orbit comparison: ICRF
    RMS RTN   2.8182836396022540E-002   2.4598832384842121E-002   2.5879201921952168E-002
    RMS XYZ   2.0336204859568077E-002   8.4715644601919167E-003   3.9687932322714677E-002
    Orbit comparison: ITRF
    RMS XYZ   1.8757217704973204E-002   1.1635302426688266E-002   3.9702619816620370E-002
    Write orbit matrices to output files
    CPU Time (sec)   350.88653299999999
    




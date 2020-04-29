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
    ├── LICENSE.md
    ├── INSTALL.md
    ├── README.md
    ├── src/
    ├── bin/  (created)
    ├── lib/  (created)
    ├── config/
    ├── tables/
    ├── scripts/

### Dependencies

1. The open basic linear algebra library (Openblas.x86_64,liblas-libs.x86_64) (You may need to run the command ln -s /usr/lib64/libopenblas.so.3 /usr/lib64/libopenblas.so)
2. A working C compiler (gcc will do), a working C++ compiler (gcc-g++ will do) and a fortran compiler (we have used gfortran)
3. Cmake (from cmake.org) at least version 2.8
4. If the flags set in CMakeLists.txt do not work with your compiler please remove/replace the ones that don't

### Build

To build the `POD` ...

    $ cd pod
    $ mkdir build
    $ cd build
    $ cmake3 .. 
    $ make >make.out 2>make.err
    $ less make.err (to verify everything was built correctly)

You should now have the executables in the bin directory: pod crs2trs brdc2ecef

### Test 

To test your build of the  `POD` ... - You may not need the ulimit command but we found it necessary

    $ cd ../pod/test
    $ ulimit -s unlimited
    $ ./sh_test_pod

At the completion of the test run, the sh_test_pod script will return any differences to the standard test resuts

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
    

## Acknowledgements

In this section we wish to acknowledge the use of and heritage of some of the source code that we have used to help develop the POD.

### Eclipse Routine

The routines to calculate the eclipsing times for GPS satellites were based off the original routines written by Jan Kouba, they have since been heavily modified.

### JPL Planetary Ephemerides

We are using the Jet Propulison (JPL) Planetary and Lunar Ephemerides processing program ((ftp://ssd.jpl.nasa.gov/pub/eph/planets/fortran/ ), in particular the routines:

- CONST.f
- FSIZER3.f
- INTERP.f
- PLEPH.f
- SPLIT.f

We have modified the following subroutines:

- asc2eph.f90
- STATE.f90

so thate there is no longer a dependency on a binary file produced in the original JPL form.

### Standards of Fundamental Astronomy (SOFA) routines

We have used a number of routines obtained from SOFA, http://www.iausofa.org/ :

- anp.for
- bi00.for
- bpn2xy.for
- bpn2xy.for
- c2ixys.for
- c2tcio.for
- cal2jd.for
- cp.for
- cr.for
- era00.for
- fad03.for
- fae03.for
- faf03.for
- faju03.for
- fal03.for
- falp03.for
- fama03.for
- fame03.for
- fane03.for
- faom03.for
- fapa03.for
- fasa03.for
- faur03.for
- fave03.for
- gmst00.for
- gmst06.for
- gmst_iers.f03
- ir.for
- jd2cal.for
- jdcalf.for
- numat.for
- nut00a.for
- obl80.for
- pn00a.for
- pn00.for
- pnm00a.for
- pnm06a.for
- pom00.for
- pr00.for
- rx.for
- rxr.for
- ry.for
- rz.for
- s00.for
- s06.for
- sp00.for
- taiutc.for
- tide_pole_oc.f90
- tide_pole_se.f90
- time_GPS.f90
- time_TAI.f90
- time_TT.f90
- time_TT_sec.f90
- time_UTC.f90
- tr.for
- xy06.for
- xys00a.for
- xys06a.for

### International Earth Rotation Service (IERS) routines 

The following routines we originally sourced fromthe IERS:

- interp_iers.f
- CNMTX.F
- FUNDARG.F
- LAGINT.f
- ORTHO_EOP.F
- PMSDNUT2.F
- RG_ZONT2.F
- UTLIBR.F
- IERS_CMP_2015.F








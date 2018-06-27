# Analysis Centre Software - POD

## Overview

The *Analysis Centre Software (ACS)* is a processing package being developed to processes GNSS observations for geodetic 
applications.  

We currently support the processing of:

* the American Global Positioning System (`GPS`); and
* the Chinese Navigation Satellite System (`Beidou`)

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

The `ACS` Version 1.0 beta release supports:

1. The `POD` 

## Directory Structure

    pod/
    ├── README.md

### Dependencies

1. The lapack numerical linear algebra library (lapack.x86_64) (You may need to run the command ln -s /usr/lib64/liblapack.so.3 /usr/lib64/liblapack.so)
1. The basic linear algebra library (blas.x86_64,liblas-libs.x86_64) (You may need to run the command ln -s /usr/lib64/libblas.so.3 /usr/lib64/libblas.so)

### Build

To build the `POD` ...

    $ cd pod

You should now have...

    pod/
    ├── README.md

### Configuration File

The `POD` Precise Orbit Determination (`./main_orb.e`) uses the configuration file:
    ├── EQM.in




### Run
  
To run the `POD` ...

    $ ./main_orb.e

This should output the following to `stdout`...

 Orbit Determination
 Orbit residuals in ICRF : RMS(XYZ)   5.0362114805437154E-002  0.37740266043432752        1.0641007845262487    
 Orbit Determination: Completed
 External Orbit comparison
 Orbit comparison: ICRF
 RMS RTN  0.59203966644540662       0.91129368973505365       0.32532340451976832
 RMS XYZ   5.0362114805437154E-002  0.37740266043432752        1.0641007845262487
 Orbit comparison: ITRF
 RMS XYZ  0.37973061226403487        2.6051541426323321E-002   1.0641454646428574
 Write orbit matrices to output files
 CPU Time (sec)   266.74518600000005

and produce the following output files...

    ├── DE.430
    ├── Amatrix.out
    ├── Wmatrix.out
    ├── orbext_ICRF.out
    ├── orbext_ITRF.out
    ├── dorb_icrf.out
    ├── dorb_RTN.out
    ├── dorb_Kepler.out
    ├── dorb_itrf.out
    ├── orb_icrf.out
    ├── orb_itrf.out
    ├── VEQ_Smatrix.out
    ├── VEQ_Pmatrix.out

    $ tail -n 31 output/
    
Performance on a t2.medium (2 virtual CPUS and 4GB of RAM) amazon server: 

    $ time ./main_orb.e 
 Orbit Determination
 Orbit residuals in ICRF : RMS(XYZ)   5.0362114805437154E-002  0.37740266043432752        1.0641007845262487    
 Orbit Determination: Completed
 External Orbit comparison
 Orbit comparison: ICRF
 RMS RTN  0.59203966644540662       0.91129368973505365       0.32532340451976832
 RMS XYZ   5.0362114805437154E-002  0.37740266043432752        1.0641007845262487
 Orbit comparison: ITRF
 RMS XYZ  0.37973061226403487        2.6051541426323321E-002   1.0641454646428574
 Write orbit matrices to output files
 CPU Time (sec)   266.74518600000005

real    4m26.779s
user    4m24.284s
sys     0m2.465s


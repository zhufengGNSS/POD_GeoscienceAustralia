# Analysis Centre Software - POD

## Precise Orbit Determination (POD)

The POD configuration refers to the overall orbit modelling applied and is performed through three configuration files as listed:

1. `POD.in` is the master configuration file for the overall parameterisation. It provides the major parameters to be written to the following two files. 
2. `EQM.in` Equation of Motion parameterisation (full force model) 
3. `VEQ.in` Variational Equations parameterisation (minimum force model) 

It should be clarified that the EQM.in and VEQ.in configuration files are identical but they are applied to different sets of the orbital equations i.e. the Equation of Motion and Variational Equations. 

Therefore, the parameters are set to different values.


The description of the three configuration files is given by the three corresponding files 'POD.md', 'EQM.md' and 'VEQ.md'.

The description is provided per each configurable parameter that is defined within the three aforementioned configuration files of the POD.


# VEQ.in configuration file's parameters description

The VEQ.in configuration file is used for parameterisation of the Variational Equations solution


** Parameter name:	`Satellite_PRN`  **

Description: 	GNSS satellite PRN number as composed by one letter denoting the GNSS constellation and a two digits number. 
				The parameter is being automatically written here based on the configuration of the master configuration file POD.in
				
Input Options:	-
	
Input example:	G23


# 
** Parameter name:	`Reference_frame` **

Description: 	Reference frame of the satellite's initial state vector.  
				The parameter is being automatically written here based on the configuration of the master configuration file POD.in
				
Input Options:	

				1. ITRF   : International Terrestrial Reference Frame 
				2. ICRF   : International Celestial Reference Frame
				3. Kepler : Kepler elements
				
Input example:	ITRF


#
** Parameter name:	`Time_scale` ** 

Description: 	Time scale of the satellite's initial epoch 

Input Options:	

				1. GPS : GPS time system 
				2. TT  : Terrestrial Time
				3. UTC : Universal Coordinated Time
				
Input example:	GPS


# 
** Parameter name:	`Orbit_arc_length`  ** 

Description: 	Orbit arc length in seconds. 
				The parameter is being automatically written here based on the configuration of the master configuration file POD.in 
				
Input Options:	Value in seconds greater than the numerical integrator step 

Input example:	86400


 
# Initial Conditions: Epoch and State Vector

## Initial Epoch:

** Parameter name:	`Year` ** 

Description: 	Year of the initial conditions' epoch. The parameter is being written here based on the configuration of the master configuration file POD.in

Input Options:	- 

Input example:	2015


# 
** Parameter name:	`Month` **

Description: 	Month of the initial conditions' epoch. The parameter is being written here based on the configuration of the master configuration file POD.in

Input Options:	1 to 12 

Input example:	4


# 
** Parameter name:	`Day` ** 

Description: 	Calendar Day of the initial conditions' epoch. The parameter is being written here based on the configuration of the master configuration file POD.in 

Input Options:	1 to 31 

Input example:	12


# 
** Parameter name:	`Seconds`  **

Description: 	Seconds since the start of the day (00 hour) of the initial conditions' epoch. The parameter is being written here based on the configuration of the master configuration file POD.in 

Input Options:	0 to any positive value 

Input example:	0.000000000


 
## Initial State Vector:

** Parameter name:	`state_vector`  **

Description: 	Initial state vector' six coodinates values.  The parameter is being written here based on the configuration of the master configuration file POD.in

Input Options:	

				1. Position vector X,Y,Z in meters and Velcoity vector (Vx, Vy, Vz) in meter/seconds
				2. Kepler elements values: 
				
Input example:	6804571.8859999999       -15232535.475000001       -20431678.779000003        1697.1661907155067        2026.6151679363102       -982.28787769004703


  
# Earth Orientation modelling

** Parameter name:	`EOP_data_sol`  **

Description: 	Earth Orientation Parameters (EOP) data solution.  The parameter is being written here based on the configuration of the master configuration file POD.in 

Input Options:	

				1. 1 : IERS C04 combined solution
				2. 2 : IERS Rapid Service/Prediciton (RS/PC) daily solution 
				3. 3 : Earth Rotation Parameters (ERP) solution by IGS ultra-rapid products  

Input example:	1 or 2 or 3


# 
** Parameter name:	`EOP_filename`   **

Description: 	EOP data file name.  The parameter is being written here based on the configuration of the master configuration file POD.in  

Input Options:	  

Input example:	eopc04_14_IAU2000.62-now


#  
** Parameter name:	`ERP_filename`  **

Description: 	ERP data file name. The parameter is being written here based on the configuration of the master configuration file POD.in  

Input Options:	  

Input example:	igu18543_12.erp


# 
** Parameter name:	`EOP_interpolation_points`  ** 

Description: 	Number of EOP data epochs (days) used for numerical interpolation based on Lagrange polynomials. The parameter is being written here based on the configuration of the master configuration file POD.in 

Input Options:	  

Input example:	4


# 
** Parameter name:	`iau_pn_model`  **

Description: 	Precession-Nutation model by International Astronomical Union (IAU). The parameter is being written here based on the configuration of the master configuration file POD.in

Input Options:	

				1. 2000 : IAU2000A model
				2. 2006 : IAU2006/2000A model   
				
Input example:	2000


 
# Numerical integration methods

** Parameter name:	`integrator_meth`  **

Description: 	Numerical Integration method

Input Options:	

				1. RKN7 : Runge-Kutta-Nystrom 7th order RKN7(6)8 
				2. RK8  : Runge-Kutta 8th order RK8(7)13 
				3. RK4  : Runge-Kutta 4th order	 
				
Input example:	RKN7


# 
** Parameter name:	`integrator_step` **

Description: 	Numerical Integrator' stepsize in seconds

Input Options:	

Input example:	900 sec

 

# Gravitational Forces

** Parameter name:	`Gravity_field`  ** 

Description: 	Earth gravity field effect 

Input Options:	

				1. 1 : Effect is considered
				2. 0 : Effect is not considered

Input example:	1

 
# 
** Parameter name:	`Planets_perturbations`  **

Description: 	Planetary and Lunar perturbations effect 

Input Options:	

				1. 1 : Effect is considered
				2. 0 : Effect is not considered

Input example:	1


# 
** Parameter name:	`Tidal_effects` ** 

Description: 	Tides effect (overall including Solid Earth tides, Ocean tides, pole tide) 

Input Options:	

				1. 1 : Effect is considered
				2. 0 : Effect is not considered

Input example:	0


# 
** Parameter name:	`Relativistic_effects` **

Description: 	General relativistic corrections to the equation of motion 

Input Options:	

				1. 1 : Effect is considered
				2. 0 : Effect is not considered

Input example:	0



## Earth Gravity Field

** Parameter name:	`gravity_model` **

Description: 	Earth gravity field model 

Input Options:	

				1. 1 : Static global gravity model     	  
				2. 2 : Time-variable global gravity model 
				3. 3 : IERS conventional geopotential model 	 

Input example:	2


# 
** Parameter name:	`gravity_model_filename` **

Description: 	Gravity model file name that provides the spherical harmonic coefficients 

Input Options:	File in gfc format (gravity field coefficients)

Input example:	goco05s.gfc


# 
** Parameter name:	`degree_max`  **

Description: 	Gravity model maximum degree/order (d/o) of spherical harmonics expansion 

Input Options:	Values range varies from 1 up to Model's maximum degree

Input example:	5


# 
** Parameter name:	`degree_max_timevar`  **

Description: 	Time-variable gravity coefficients' maximum degree/order (d/o); This parameters is applied in the case of a time-variable gravity model (gravity_model=2)

Input Options:	1 to maximum degree of the time-variable part only

Input example:	5


  
## Planetary/Lunar ephemeris

** Parameter name:	`DE_fname_header`  **

Description: 	The ephemeris is provided by the JPL DE (Development Ephemeris) through two files. This parameter refers to the ephemeris' header file name

Input Options:	 

Input example:	header.430_229


#  
** Parameter name:	`DE_fname_data`  **

Description: 	The ephemeris is provided by the JPL DE (Development Ephemeris) through two files. This parameter refers to the ephemeris' data file name

Input Options:	 

Input example:	ascp1950.430


 
## Tidal Effects:

** Parameter name:	`solid_tides_nonfreq`  **

Description: 	Solid Earth Tides frequency-independent terms 

Input Options:	

				1. 1 : Effect is considered
				2. 0 : Effect is not considered

Input example:	1


# 
** Parameter name:	`solid_tides_freq`  **

Description: 	Solid Earth Tides frequency-dependent terms 

Input Options:	

				1. 1 : Effect is considered
				2. 0 : Effect is not considered

Input example:	1


# 
** Parameter name:	`ocean_tides`  **

Description: 	Ocean Tides 

Input Options:	

				1. 1 : Effect is considered
				2. 0 : Effect is not considered

Input example:	1


# 
** Parameter name:	`solid_earth_pole_tide` **

Description: 	Solid Earth Pole tide 

Input Options:	

				1. 1 : Effect is considered
				2. 0 : Effect is not considered

Input example:	1


# 
** Parameter name:	`ocean_pole_tide` **

Description: 	Ocean Pole tide 

Input Options:	

				1. 1 : Effect is considered
				2. 0 : Effect is not considered

Input example:	1


 
## Ocean Tides:

** Parameter name:	`ocean_tides_model_file`  **

Description: 	Ocean Tides model file name 

Input Options:	

Input example:	fes2004_Cnm-Snm.dat   (This refers to the FES2004 ocean tide model that is applied as standard here) 


# 
** Parameter name:	`ocean_tides_model_deg` ** 

Description: 	Maximum degree/order of the spherical harmonics expansion of the ccean tides effect 

Input Options:	Values range varies from 0 up to Model's maximum degree

Input example:	15 (It is suggested to set the same value with the gravity field' maximum degree expansion)


 
# Non-gravitational Effects

** Parameter name:	`Solar_radiation` ** 

Description: 	Solar radiation pressure effect 

Input Options:	

				1. 1 : Effect is considered
				2. 0 : Effect is not considered

Input example:	1


# 
** Parameter name:	`Earth_radiation`  **

Description: 	Earth radiation pressure effect 

Input Options:	

				1. 1 : Effect is considered
				2. 0 : Effect is not considered

Input example:	0


# 
** Parameter name:	`Antenna_thrust`  **

Description: 	Antenna thrust effect 

Input Options:	

				1. 1 : Effect is considered
				2. 0 : Effect is not considered

Input example:	0



## Solar Radiation Pressure (SRP) models

** Parameter name:	`ECOM_param`  **

Description: 	The configurable parameter defines the option of the SRP models and the relevant parameters to be estimated

Input Options:	

				1. 1 : ECOM1 model (parameters estimated on D,Y,B directions)
				2. 2 : ECOM2 model (parameters estimated on D,Y,B directions)
				3. 3 : Simple Box-wing (SBOXW) model (paramters estimated on D,Y,B,X,Z directions)
				4. 0 : Cannonball model (No parameters are estimated)

Input example:	1 or 2 or 3 or 0


 
# SRP ECOM models parameters' directions to be estimated


## Bias accelerations in D, Y and B directions

** Parameter name:	`bias_D`  **

Description: 	Bias acceleration estimated in D direction 

Input Options:	

				1. 1 : Parameter is estimated
				2. 0 : Parameter is not estimated

Input example:	1


# 
** Parameter name:	`bias_Y`  **

Description: 	Bias acceleration estimated in Y direction 

Input Options:	

				1. 1 : Parameter is estimated
				2. 0 : Parameter is not estimated

Input example:	1


# 
** Parameter name:	`bias_B`  **

Description: 	Bias acceleration estimated in B direction 

Input Options:	

				1. 1 : Parameter is estimated
				2. 0 : Parameter is not estimated

Input example:	1


 
## Cycle-per-revolution accelerations in D, Y and B directions

** Parameter name:	`cpr_D`  **

Description: 	Cycle per revolution term estimated in D direction 

Input Options:	

				1. 1 : Parameter is estimated
				2. 0 : Parameter is not estimated

Input example:	1


# 
** Parameter name:	`cpr_Y`  **

Description: 	Cycle per revolution term estimated in Y direction 

Input Options:	

				1. 1 : Parameter is estimated
				2. 0 : Parameter is not estimated

Input example:	1


# 
** Parameter name:	`cpr_B`  **

Description: 	Cycle per revolution term estimated in B direction 

Input Options:	

				1. 1 : Parameter is estimated
				2. 0 : Parameter is not estimated

Input example:	1

 
 
# Empirical forces modelling (bias accelerations and one cycle per revolution terms)

** Parameter name:	`EMP_param`  **

Description: 	Empirical forces modelling in orbital frame. Empirical forces include bias accelerations and one cycle per revolution 

Input Options:	

				1. 1 : Effect is considered
				2. 0 : Effect is not considered

Input example:	0


## Bias accelerations per radial, along-track and cross-track directions

** Parameter name:	`bias_r`  **

Description: 	Bias acceleration estimated in radial direction 

Input Options:	

				1. 1 : Parameter is estimated
				2. 0 : Parameter is not estimated

Input example:	1


# 
** Parameter name:	`bias_t`  **

Description: 	Bias acceleration estimated in along-track direction 

Input Options:	

				1. 1 : Parameter is estimated
				2. 0 : Parameter is not estimated

Input example:	1


# 
** Parameter name:	`bias_n`  **

Description: 	Bias acceleration estimated in cross-track direction 

Input Options:	

				1. 1 : Parameter is estimated
				2. 0 : Parameter is not estimated

Input example:	1


## Cycle-per-revolution accelerations per radial, along-track and cross-track directions

** Parameter name:	`cpr_r`  **

Description: 	Cycle per revolution coefficient estimated in radial direction 

Input Options:	

				1. 1 : Parameter is estimated
				2. 0 : Parameter is not estimated

Input example:	1


# 
** Parameter name:	`cpr_t`  **

Description: 	Cycle per revolution coefficient estimated in along-track direction 

Input Options:	

				1. 1 : Parameter is estimated
				2. 0 : Parameter is not estimated

Input example:	1


# 
** Parameter name:	`cpr_n`  **

Description: 	Cycle per revolution coefficient estimated in cross-track direction 

Input Options:	

				1. 1 : Parameter is estimated
				2. 0 : Parameter is not estimated

Input example:	1



# Observation Model (external orbit used to form pseudo-observations)

** Parameter name:	`pseudobs_filename`  **

Description: 	GNSS orbit' file name in sp3 format. The parameter is being written here based on the configuration of the master configuration file POD.in 

Input Options:	sp3 orbits as provided by International GNSS Service (IGS) products

Input example:	igs18400.sp3


## Numerical interpolation of the sp3 orbit based on Lagrange polynomials  

** Parameter name:	`pseudobs_interp_points`  ** 

Description: 	Number of data points used in the Lagrange interpolation 

Input Options:	Values greater than zero

Input example:	12 (12 is suggested for the GNSS orbits)


# 
** Parameter name:	`pseudobs_interp_step`  ** 

Description: 	Interval of the interpolated orbit in seconds 

Input Options:	Values should be equal or smaller than the rate of the sp3 orbit

Input example:	900

 
# Variational Equations

** Parameter name:	`VEQ_integration`  **

Description: 	Variational Equations numerical integration 

Input Options:	The parameter is being written here based on the configuration of the master configuration file POD.in

Input example:	0

 
# Parameter Estimation

** Parameter name:	`Estimator_procedure`  **

Description: 	Orbit parameter estimation based on least-squares method.  The parameter is being written here based on the configuration of the master configuration file POD.in 

Input Options:	

				1. 1 : Orbit parameter estimator is applied
				2. 0 : Orbit parameter estimator is not applied

Input example:	1 or 0


# 
** Parameter name:	`Estimator_Iterations`  **

Description: 	Number of iterations of the orbit parameter estimator (least-squares method).  The parameter is being written here based on the configuration of the master configuration file POD.in 

Input Options:	Values range 0 (no iterations) to any positive integer. 1 or 2 iterations are usually sufficient in terms of orbit accuracy.

Input example:	2

 
# External Orbit Comparison

** Parameter name:	`orbit_external_opt` ** 

Description: 	Type of external orbit to be applied for orbit comparison 

Input Options:	1. 2 : Interpolated orbit based on Lagrange interpolation of sp3 orbit

Input example:	2


# 
** Parameter name:	`orbit_ext_frame`  **

Description: 	Reference frame of external orbit

Input Options:	ITRF 

Input example:	ITRF


# 
** Parameter name:	`orbit_filename`  ** 

Description: 	External orbit' file name (sp3 format). The parameter is being written here based on the configuration of the master configuration file POD.in

Input Options:	Orbit file in sp3 format 

Input example:	igs18400.sp3


# 
** Parameter name:	`orbit_interp_step`  ** 

Description: 	Interval of the interpolated orbit in seconds 

Input Options:	Values should be equal or smaller than the rate of the sp3 orbit

Input example:	900


# 
** Parameter name:	`orbit_interp_points`  **

Description: 	Number of data points used in the Lagrange interpolation 

Input Options:	Values greater than zero

Input example:	12 (12 is suggested for the GNSS orbits)

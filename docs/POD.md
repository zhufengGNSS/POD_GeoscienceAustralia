# Precise Orbit Determination (POD) Tool configuration 

The 'POD.in' is the master configuration file of the POD
This document provides the description of the configuration file parameters 


## POD mode:
** Parameter name:	`POD_MODE_cfg` **
 
Description: 	POD basic modes options

Input Options:	

  				1. 1 : Orbit Determination (orbit fitting and parameter estimation using pseudo-observations)
				2. 2 : Orbit Determination and Prediction
				3. 3 : Orbit Integration (Equation of Motion numerical integration solution)				
				4. 4 : Orbit Integration and Partials (Equation of Motion and Variational Equations numerical inegration solution)	
				
Input example:	2


## Initial Conditions (IC):

** Parameter name:	`IC_input` **

Description: 	Initial Conditions input mode

Input Options:	

				1. 1 : Input a-priori orbit in sp3 format (applied as pseudo-observations)
				2. 2 : Input file with Initial Conditions (State Vector and Parameters at initial epoch per satellite) based on internal POD format 
				
Input example:	1
		
  
#    
** Parameter name:	`IC_refsys` **

Description: 	Initial Conditions reference frame

Input Options:	
				
				1. ICRF : International Celestial Reference Frame
				2. ITRF : International Terrestrial Reference Frame 
				
Input example:	ICRF
  
  
#    
** Parameter name:	`IC_filename_cfg` **

Description: 	Initial Conditions file name based on POD format

Input Options:	File name charachters

Input example:	orb_pea.out



## Configuration files of Orbit modelling (2 basic initial files):

** Parameter name:	`EQM_fname_cfg`  **

Description: 	Configuration file for the Equation of Motion based on POD format

Input Options:	File name charachters

Input example:	EQM.in 
  
  
#    
** Parameter name:	`VEQ_fname_cfg`  **

Description: 	Configuration file for the Variational Equations based on POD format 

Input Options:	File name charachters

Input example:	VEQ.in 



## Orbit arc length

** Parameter name:	`orbit_determination_arc_cfg`  **

Description: 	Orbit Estimation arc length in hours

Input Options:	

Input example:	24 
   
   
#      
** Parameter name:	`orbit_prediction_arc_cfg`  **

Description: 	Orbit Prediction arc length in hours (in the case that POD mode is set to 2, POD_MODE_cfg=2 ) 

Input Options:
	
Input example:	12 
   
   
#    
** Parameter name:	`orbit_backwards_arc_cfg`  **

Description: 	Arc length (in hours) of the backwards orbit numerical integration 

Input Options:	

Input example:	2


## Earth Orientation modelling

** Parameter name:	`EOP_solution_cfg`  **

Description: 	Earth Orientation Parameters (EOP) data solution.  

Input Options:	

				1. 1 : IERS C04 combined solution
				2. 2 : IERS Rapid Service/Prediciton (RS/PC) daily solution 				
				3. 3 : Earth Rotation Parameters (ERP) solution by IGS ultra-rapid products  
				
Input example:	1
  
   
#     
** Parameter name:	`EOP_fname_cfg`  **

Description: 	EOP data file name 

Input Options:	  

Input example:	eopc04_14_IAU2000.62-now
   
   
#    
** Parameter name:	`ERP_fname_cfg`  **

Description: 	ERP data file name 

Input Options:	  

Input example:	igu18543_12.erp
   
   
#     
** Parameter name:	`EOP_Nint_cfg`  **

Description: 	Number of EOP data epochs (days) used for numerical interpolation based on Lagrange polynomials. 

Input Options:	  

Input example:	4
   
   
#      
** Parameter name:	`iau_model_cfg`  **

Description: 	Precession-Nutation model by International Astronomical Union (IAU). 

Input Options:	
	
				1. 2000 : IAU2000A model
				2. 2006 : IAU2006/2000A model   
				
Input example:	2000



## Observation Model (external orbit used to form pseudo-observations)
** Parameter name:	`pseudobs_orbit_filename_cfg`  **

Description: 	GNSS orbit' file name in sp3 format 

Input Options:	sp3 orbits as provided by International GNSS Service (IGS) products

Input example:	igs18400.sp3



## External Orbit Comparison
** Parameter name:	`ext_orbit_filename_cfg` **

Description: 	File name of the external orbit (sp3 format) to be applied for orbit comparison

Input Options:	Orbit file in sp3 format 

Input example:	igs18400.sp3



## Orbit parameter estimation based on least-squares method (in case that POD_MODE_cfg = 1 or 2)
** Parameter name:	`Estimator_Iterations_cfg`  **

Description: 	Number of iterations of the orbit parameter estimator (least-squares method).  

Input Options:	Values range 0 (no iterations) to any positive integer. 1 or 2 iterations are usually sufficient in terms of orbit accuracy.

Input example:	2
   
   
#    
** Parameter name:	`VEQ_REFSYS_cfg`  **

Description: 	Reference System of the partial derivatives of the Variational Equations' solution 

Input Options:	
				
				1. ICRS : Celestial Reference System 
				2. ITRS : Terrestrial Reference System
				
Input example:	ITRS
   
   
#     
** Parameter name:	`sp3_velocity_cfg`  **

Description: 	Option for write or not the satellite Velocity vector to the orbit sp3 format 

Input Options:	

				1. 0 : Do not write Velocity vector to sp3 orbit
				2. any value > 0 : Write Velocity vector to sp3 orbit

Input example:	0


#    
** Parameter name:	`partials_velocity_cfg`  **

Description: 	Option for write or not the partials derivatives of the velocity vector w.r.t. parameters into the orbits_partials output file in POD format

Input Options:	

				1. 0 : partials_velocity_cfg = 0 :: Do not write Velocity vector's partials elements
				2. any value > 0 : Write Velocity vector's partials elements
				
Input example:	0
   
   
#    
** Parameter name:	`leapsec_filename_cfg`  **

Description: 	Leap seconds file name (internal format) as provided by POD package.

Input Options:	 

Input example:	leap.second
   
   
#    
** Parameter name:	`satsinex_filename_cfg` **

Description: 	File name of the satellite metadata in SINEX format  

Input Options:	 

Input example:	igs_metadata_2063.snx
   
   
#   
** Parameter name:	`SRP_MOD_arp` **

Description: 	A-priori model for the Solar Radiation Pressure effect 

Input Options:	

				1. 1 : use a cannonball f0 model 
				2. 2 : use simple box-wing model 
				3. 3 : use box-wing model from repro3 routines 
				4. 0 : no a priori model  
				
Input example:	0
 

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

 




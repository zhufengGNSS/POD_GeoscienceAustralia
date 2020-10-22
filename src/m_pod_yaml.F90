! -----------------------------------------------------------------------------
! This file is part of Fortran-YAML: a lightweight YAML parser written in
! object-oriented Fortran.
!
! Official repository: https://github.com/BoldingBruggeman/fortran-yaml
!
! Copyright 2013-2016 Bolding & Bruggeman ApS.
!
! This is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation (https://www.gnu.org/licenses/gpl.html). It is distributed in the
! hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
! implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! A copy of the license is provided in the COPYING file.
! -----------------------------------------------------------------------------

module pod_yaml

   use yaml_version, only: yaml_commit_id=>git_commit_id, &
                           yaml_branch_name=>git_branch_name
   use yaml
   use yaml_types
   use, intrinsic :: iso_fortran_env

   implicit none

   character(error_length) :: error
   character(256) :: path
   class (type_node), pointer :: root
   type (type_error), target :: my_error
   type (type_error), pointer :: my_error_p
   type (type_dictionary), pointer :: root_dict, pod_data_dict, pod_options_dict, eqm_options_dict, veq_options_dict
   type (type_dictionary), pointer :: srp_dict, srp_parameters_dict, integ_dict
   type (type_dictionary), pointer :: gravity_dict, gravity_model_dict, planets_dict, tides_dict, rel_dict, non_grav_dict
   type (type_dictionary), pointer :: overrides_dict
   logical yml_pod_data_enabled, yml_ext_orbit_enabled, yml_estimate_params, yml_write_sp3_velocities, yml_write_partial_velocities
   logical yml_veq_integration
   integer*2 yml_pod_mode
   integer*2 yml_ic_input_format
   integer*2 yml_ic_input_refsys
   integer*2 yml_ext_orbit_opt
   integer*2 yml_ext_orbit_frame
   integer*2 yml_iau_model
   integer*2 yml_eqm_integrate_method, yml_veq_integrate_method
   integer*2 yml_ECOM_estimator
   integer*2 yml_apriori_srp
   integer*2 yml_eqm_gravity_model, yml_veq_gravity_model
   ! srp_parameters are a bitfield type for all variables, 6 (DYB,RTN) ECOM biases, 8 (DYB,RTN,D2,D4) ECOM cprs, 3 (RTN) EMP biases, 
   ! 3 (RTN) EMP cprs.
   ! NB EMP bits are 17-22
   integer*4 yml_default_srp_parameters, yml_srp_parameters_defined
   ! tidal effects are a bitfield type, 1 = solid_nonfreq, 2 = solid_freq, 3 = ocean, 4 = solid_pole, 5 = ocean_pole
   integer*2 yml_eqm_tidal_effects, yml_veq_tidal_effects
   integer*2 yml_eqm_gravity_max_degree, yml_eqm_gravity_time_max_degree, yml_eqm_tides_max_degree
   integer*2 yml_veq_gravity_max_degree, yml_veq_gravity_time_max_degree, yml_veq_tides_max_degree
   ! non grav effects are a bitfield type, 1 = thrust, 2 = earth, 3 = solar
   integer*2 yml_eqm_non_grav_effects, yml_veq_non_grav_effects
   logical   yml_eqm_planetary_perturbations_enabled, yml_eqm_tidal_effects_enabled, yml_eqm_rel_effects_enabled
   logical   yml_eqm_non_grav_effects_enabled, yml_veq_non_grav_effects_enabled
   logical   yml_veq_planetary_perturbations_enabled, yml_veq_tidal_effects_enabled, yml_veq_rel_effects_enabled
   logical   yml_eqm_gravity_enabled, yml_veq_gravity_enabled

   character(4) yml_pod_data_prn
   character(8) yml_pod_data_ref_frame
   character(8) yml_pod_data_time_scale
   integer*8    yml_pod_data_orbit_arc_seconds
   character(20) yml_pod_data_initial_epoch
   character(512) yml_pod_data_state_vector
   character(128) yml_orbit_filename, yml_ext_orbit_filename, yml_satsinex_filename, yml_leapsecond_filename, yml_eop_filename
   character(128) yml_gravity_filename, yml_ephemeris_header, yml_ephemeris_data_file, yml_ocean_tides_file

   integer*4 yml_orbit_step, yml_orbit_points, yml_orbit_arc_determination, yml_orbit_arc_prediction, yml_orbit_arc_backwards
   integer*4 yml_ext_orbit_steps, yml_ext_orbit_points, yml_eop_option, yml_eop_int_points, yml_eqm_integ_stepsize
   integer*4 yml_estimator_iterations, yml_estimator_procedure, yml_veq_integ_stepsize

   type :: srp_override
      integer*4 parameters, parameters_defined
      integer*2 apriori_model
   end type

   type :: integration_override
      real*8 veq_stepsize, veq_iterations, eqm_stepsize, eqm_iterations
      logical veq_enabled, eqm_enabled
   end type

   type :: override
      character(20) name ! system or block or prn
      type (srp_override) ::  srp
      type (integration_override) ::  integ
   end type

   integer*4 max_sys_overrides, max_block_overrides, max_prn_overrides

   parameter (max_sys_overrides = 8)
   parameter (max_block_overrides = 256)
   parameter (max_prn_overrides = 2048)

   integer*4 sys_override_count, block_override_count, prn_override_count

   type (override) :: yml_sys_overrides (max_sys_overrides), yml_block_overrides (max_block_overrides)
   type (override) :: yml_prn_overrides (max_prn_overrides)

   contains

subroutine get_yaml()

   nullify(root_dict)
   nullify(pod_data_dict)
   nullify(pod_options_dict)
   nullify(eqm_options_dict)
   nullify(veq_options_dict)
   nullify(srp_dict)
   nullify(srp_parameters_dict)
   nullify(gravity_dict)
   nullify(gravity_model_dict)
   nullify(planets_dict)
   nullify(tides_dict)
   nullify(rel_dict)
   nullify(non_grav_dict)
   nullify(overrides_dict)
   nullify(integ_dict)

   my_error_p => my_error
   yml_pod_mode = -1
   yml_ic_input_format = -1
   yml_ic_input_refsys = -1

   sys_override_count = 0
   block_override_count = 0
   prn_override_count = 0

   write(*,*)
   write(*,*) 'YAML version:   ',yaml_commit_id,' (',yaml_branch_name,' branch)'
   write(*,*)

   call get_command_argument(1, path)
   if (path=='') then
      write (*,*) 'ERROR: path to YAML file not provided.'
      stop 2
   end if
   root => parse(path,unit=100,error=error)
   if (error/='') then
      write (*,*) 'PARSE ERROR: '//trim(error)
      stop 1
   end if
   !call root%dump(unit=output_unit,indent=0)

   if (associated(root)) then
      select type (root)
      type is (type_dictionary)
         root_dict => root
      class default
         ! do nothing
      end select 
   end if

   !TODO: check for errors each step of the way
   if (.not.associated(root_dict)) then
        write(*,*) "could not find any data in YAML config"
        STOP
   end if
   pod_data_dict = root_dict%get_dictionary("pod_data", .true., my_error_p)
   if (.not.associated(pod_data_dict)) then
      write(*,*) "could not find pod_data label in YAML config"
      STOP
   else
      call get_pod_data(pod_data_dict, my_error, yml_pod_data_enabled, yml_pod_data_prn, yml_pod_data_ref_frame,&
              yml_pod_data_time_scale, yml_pod_data_orbit_arc_seconds, yml_pod_data_initial_epoch,&
              yml_pod_data_state_vector)  
   end if

   pod_options_dict = root_dict%get_dictionary("pod_options", .true., my_error_p);
   if (.not.associated(pod_options_dict)) then
      write(*,*) "could not find pod_options label in YAML config"
      STOP
   else
      yml_pod_mode = get_pod_mode(pod_options_dict, my_error)
      yml_ic_input_format = get_input_format(pod_options_dict, my_error)
      yml_ic_input_refsys = get_input_refsys(pod_options_dict, my_error)
      call get_pseudoobs(pod_options_dict, my_error, yml_orbit_filename, yml_orbit_step, yml_orbit_points);
      call get_orbitarcs(pod_options_dict, my_error, yml_orbit_arc_determination, yml_orbit_arc_prediction,&
              yml_orbit_arc_backwards)
      yml_ext_orbit_enabled = pod_options_dict%get_logical("ext_orbit_enabled", .false., my_error_p)
      if (yml_ext_orbit_enabled) then
         yml_ext_orbit_opt = get_ext_orbit_opt(pod_options_dict, my_error)
         yml_ext_orbit_frame = get_ext_orbit_frame(pod_options_dict, my_error)
         call get_extorbit_comp(pod_options_dict, my_error, yml_ext_orbit_filename, yml_ext_orbit_steps,&
                 yml_ext_orbit_points)
      end if
      yml_iau_model = get_iau_model(pod_options_dict, my_error)
      call get_earth_orientation_params(pod_options_dict, my_error, yml_eop_option, yml_eop_filename, yml_eop_int_points)
      yml_satsinex_filename = pod_options_dict%get_string("satsinex_filename", "", my_error_p)
      yml_leapsecond_filename = pod_options_dict%get_string("leapsecond_filename", "", my_error_p)
      yml_gravity_filename = pod_options_dict%get_string("gravity_model_file", "", my_error_p)
      yml_ephemeris_header = pod_options_dict%get_string("DE_fname_header", "", my_error_p)
      yml_ephemeris_data_file = pod_options_dict%get_string("DE_fname_data", "", my_error_p)
      yml_ocean_tides_file = pod_options_dict%get_string("ocean_tides_model_file", "", my_error_p)
      yml_write_sp3_velocities = pod_options_dict%get_logical("sp3_velocity", .false., my_error_p)
      yml_write_partial_velocities = pod_options_dict%get_logical("partials_velocity", .false., my_error_p)
      yml_estimator_iterations = pod_options_dict%get_integer("estimator_iterations", -1, my_error_p)
      srp_dict = pod_options_dict%get_dictionary("srp_apriori_model", .true., my_error_p)
      if (.not.associated(srp_dict)) then
         write (*,*) "cannot find srp_apriori_model label in pod_options"
         STOP
      endif
      yml_apriori_srp = get_apriori_srp(srp_dict, my_error)
      srp_parameters_dict = pod_options_dict%get_dictionary("srp_parameters", .true., my_error_p)
      if (.not.associated(srp_parameters_dict)) then
         write (*,*) "cannot find srp_parameters label in pod_options"
         STOP
      end if
      yml_default_srp_parameters = get_srp_parameters(srp_parameters_dict, .true., yml_srp_parameters_defined, my_error)
      ! do we actually need this?
      yml_ECOM_estimator = get_ECOM_estimator(pod_options_dict, my_error)
      if (yml_ECOM_estimator .ge. 3) then
         yml_estimate_params = .true.
      else
         yml_estimate_params = .false.
      endif
      ! what is estimator_procedure???
      yml_estimator_procedure = pod_options_dict%get_integer("estimatr_procedure", -1, my_error_p)
      yml_veq_integration = pod_options_dict%get_logical("veq_integration", .false., my_error_p)
   end if

   eqm_options_dict = root_dict%get_dictionary("eqm_options", .true., my_error_p)
   if (.not.associated(pod_options_dict)) then
      write(*,*) "could not find eqm_options label in YAML config"
      STOP
   else
      gravity_dict = eqm_options_dict%get_dictionary("gravity_field", .true., my_error_p)
      if (.not.associated(gravity_dict)) then
         write (*,*) "could not find gravity_field label in eqm_options"
         STOP
      end if
      gravity_model_dict = gravity_dict%get_dictionary("gravity_model", .true., my_error_p)
      if (.not.associated(gravity_model_dict)) then
         write (*,*) "could not find gravity_model label in eqm_options"
         STOP
      end if
      planets_dict = eqm_options_dict%get_dictionary("planetary_perturbations", .true., my_error_p)
      if (.not.associated(planets_dict)) then
         write (*,*) "could not find planetary_perturbations label in eqm_options"
         STOP
      end if
      tides_dict = eqm_options_dict%get_dictionary("tidal_effects", .true., my_error_p)
      if (.not.associated(tides_dict)) then
         write (*,*) "could not find tidal_effects label in eqm_options"
         STOP
      end if
      rel_dict = eqm_options_dict%get_dictionary("relativistic_effects", .true., my_error_p)
      if (.not.associated(rel_dict)) then
         write (*,*) "could not find relativistic_effects label in eqm_options"
         STOP
      end if
      non_grav_dict = eqm_options_dict%get_dictionary("non_gravitational_effects", .true., my_error_p)
      if (.not.associated(non_grav_dict)) then
         write (*,*) "could not find non_gravitational_effects label in eqm_options"
         STOP
      end if
      integ_dict = eqm_options_dict%get_dictionary("integration_options", .true., my_error_p)
      if (.not.associated(integ_dict)) then
         write (*,*) "could not find integration_options label in eqm_options"
         STOP
      end if

      yml_eqm_gravity_enabled = gravity_dict%get_logical("enabled", .false., my_error_p)
      yml_eqm_gravity_model = get_gravity_model(gravity_model_dict, my_error, "eqm")
      yml_eqm_gravity_max_degree = gravity_dict%get_integer("gravity_degree_max", -1, my_error_p)
      yml_eqm_gravity_time_max_degree = gravity_dict%get_integer("timevar_degree_max", -1, my_error_p)
      yml_eqm_planetary_perturbations_enabled = planets_dict%get_logical("perturbations_enabled", .false., my_error_p)
      yml_eqm_tidal_effects_enabled = tides_dict%get_logical("enabled", .false., my_error_p)
      yml_eqm_tidal_effects = get_tidal_effects(tides_dict, my_error, "eqm")
      yml_eqm_tides_max_degree = tides_dict%get_integer("ocean_tides_degree_max", -1, my_error_p)
      yml_eqm_rel_effects_enabled = rel_dict%get_logical("enabled", .false., my_error_p)
      yml_eqm_non_grav_effects_enabled = non_grav_dict%get_logical("enabled", .false., my_error_p)
      yml_eqm_non_grav_effects = get_non_grav_effects(non_grav_dict, my_error, "eqm")
      yml_eqm_integrate_method = get_integrator_method(integ_dict, my_error, yml_eqm_integ_stepsize)
   end if

   veq_options_dict = root_dict%get_dictionary("veq_options", .true., my_error_p)
   if (.not.associated(veq_options_dict)) then
      write(*,*) "could not find veq_options label in YAML config"
      STOP
   else
      gravity_dict = veq_options_dict%get_dictionary("gravity_field", .true., my_error_p)
      if (.not.associated(gravity_dict)) then
         write (*,*) "could not find gravity_field label in veq_options"
         STOP
      end if
      gravity_model_dict = gravity_dict%get_dictionary("gravity_model", .true., my_error_p)
      if (.not.associated(gravity_model_dict)) then
         write (*,*) "could not find gravity_model label in veq_options"
         STOP
      end if
      planets_dict = veq_options_dict%get_dictionary("planetary_perturbations", .true., my_error_p)
      if (.not.associated(planets_dict)) then
         write (*,*) "could not find planetary_perturbations label in veq_options"
         STOP
      end if
      tides_dict = veq_options_dict%get_dictionary("tidal_effects", .true., my_error_p)
      if (.not.associated(tides_dict)) then
         write (*,*) "could not find tidal_effects label in veq_options"
         STOP
      end if
      rel_dict = veq_options_dict%get_dictionary("relativistic_effects", .true., my_error_p)
      if (.not.associated(rel_dict)) then
         write (*,*) "could not find relativistic_effects label in veq_options"
         STOP
      end if
      non_grav_dict = veq_options_dict%get_dictionary("non_gravitational_effects", .true., my_error_p)
      if (.not.associated(non_grav_dict)) then
         write (*,*) "could not find non_gravitational_effects label in veq_options"
         STOP
      end if
      integ_dict = veq_options_dict%get_dictionary("integration_options", .true., my_error_p)
      if (.not.associated(integ_dict)) then
         write (*,*) "could not find integration_options label in veq_options"
         STOP
      end if
      yml_veq_gravity_enabled = gravity_dict%get_logical("enabled", .false., my_error_p)
      yml_veq_gravity_model = get_gravity_model(gravity_model_dict, my_error, "veq")
      yml_veq_gravity_max_degree = gravity_dict%get_integer("gravity_degree_max", -1, my_error_p)
      yml_veq_gravity_time_max_degree = gravity_dict%get_integer("timevar_degree_max", -1, my_error_p)
      yml_veq_planetary_perturbations_enabled = planets_dict%get_logical("perturbations_enabled", .false., my_error_p)
      yml_veq_tidal_effects_enabled = tides_dict%get_logical("enabled", .false., my_error_p)
      yml_veq_tidal_effects = get_tidal_effects(tides_dict, my_error, "veq")
      yml_veq_tides_max_degree = tides_dict%get_integer("ocean_tides_degree_max", -1, my_error_p)
      yml_veq_rel_effects_enabled = rel_dict%get_logical("enabled", .false., my_error_p)
      yml_veq_non_grav_effects_enabled = non_grav_dict%get_logical("enabled", .false., my_error_p)
      yml_veq_non_grav_effects = get_non_grav_effects(non_grav_dict, my_error, "veq")
      yml_veq_integrate_method = get_integrator_method(integ_dict, my_error, yml_veq_integ_stepsize)
   end if

end subroutine get_yaml

function get_non_grav_effects(dict, error, label)
   type (type_dictionary) :: dict
   type (type_error) :: error
   type (type_error), pointer :: e
   logical solar, earth, thrust
   character (*) label
   integer*2 get_non_grav_effects

   nullify (e)
   get_non_grav_effects = 0;

   solar = dict%get_logical("solar_radiation", .false., e)
   earth = dict%get_logical("earth_radiation", .false., e)
   thrust = dict%get_logical("antenna_thrust", .false., e)

   if (associated(e)) then
      error = e
      error%message = label // ':' // error%message
   end if

   if (solar) then
      get_non_grav_effects = 1
   end if
   get_non_grav_effects = 2 * get_non_grav_effects
   if (earth) then
      get_non_grav_effects = get_non_grav_effects + 1
   end if
   get_non_grav_effects = 2 * get_non_grav_effects
   if (thrust) then
      get_non_grav_effects = get_non_grav_effects + 1
   end if

   return
end function get_non_grav_effects

function get_tidal_effects(dict, error, label)
   type (type_dictionary) :: dict
   type (type_error) :: error
   type (type_error), pointer :: e
   logical solid_nonfreq, solid_freq, ocean, solid_pole, ocean_pole
   character (*) label
   integer*2 get_tidal_effects

   nullify (e)
   get_tidal_effects = 0;

   solid_nonfreq = dict%get_logical("solid_tides_nonfreq", .false., e)
   solid_freq = dict%get_logical("solid_tides_freq", .false., e)
   ocean = dict%get_logical("ocean_tides", .false., e)
   solid_pole = dict%get_logical("solid_earth_pole_tides", .false., e)
   ocean_pole = dict%get_logical("ocean_pole_tide", .false., e)

   if (associated(e)) then
      error = e
      error%message = label // ':' // error%message
   end if

   if (solid_nonfreq) then
      get_tidal_effects = 1
   end if
   get_tidal_effects = get_tidal_effects * 2
   if (solid_freq) then
      get_tidal_effects = get_tidal_effects + 1
   end if
   get_tidal_effects = get_tidal_effects * 2
   if (ocean) then
      get_tidal_effects = get_tidal_effects + 1
   end if
   get_tidal_effects = get_tidal_effects * 2
   if (solid_pole) then
      get_tidal_effects = get_tidal_effects + 1
   end if
   get_tidal_effects = get_tidal_effects * 2
   if (ocean_pole) then
      get_tidal_effects = get_tidal_effects + 1
   end if

   return
end function get_tidal_effects

function get_srp_parameters(dict, listall, defined, error)
   type (type_dictionary) :: dict
   type (type_error) :: error
   logical listall
   integer*4 defined
   integer*4 get_srp_parameters

   defined = 0
   get_srp_parameters = 0

   call get_logical_parameter(dict, error, "ECOM_D_bias", listall, get_srp_parameters, 1, defined)
   call get_logical_parameter(dict, error, "ECOM_Y_bias", listall, get_srp_parameters, 2, defined)
   call get_logical_parameter(dict, error, "ECOM_B_bias", listall, get_srp_parameters, 3, defined)
   call get_logical_parameter(dict, error, "ECOM_R_bias", listall, get_srp_parameters, 4, defined)
   call get_logical_parameter(dict, error, "ECOM_T_bias", listall, get_srp_parameters, 5, defined)
   call get_logical_parameter(dict, error, "ECOM_N_bias", listall, get_srp_parameters, 6, defined)
   call get_logical_parameter(dict, error, "ECOM_D_cpr", listall, get_srp_parameters, 7, defined)
   call get_logical_parameter(dict, error, "ECOM_Y_cpr", listall, get_srp_parameters, 8, defined)
   call get_logical_parameter(dict, error, "ECOM_B_cpr", listall, get_srp_parameters, 9, defined)
   call get_logical_parameter(dict, error, "ECOM_R_cpr", listall, get_srp_parameters, 10, defined)
   call get_logical_parameter(dict, error, "ECOM_T_cpr", listall, get_srp_parameters, 11, defined)
   call get_logical_parameter(dict, error, "ECOM_N_cpr", listall, get_srp_parameters, 12, defined)
   call get_logical_parameter(dict, error, "ECOM_D_2_cpr", listall, get_srp_parameters, 13, defined)
   call get_logical_parameter(dict, error, "ECOM_D_4_cpr", listall, get_srp_parameters, 14, defined)
   call get_logical_parameter(dict, error, "EMP_R_bias", listall, get_srp_parameters, 17, defined)
   call get_logical_parameter(dict, error, "EMP_T_bias", listall, get_srp_parameters, 18, defined)
   call get_logical_parameter(dict, error, "EMP_N_bias", listall, get_srp_parameters, 19, defined)
   call get_logical_parameter(dict, error, "EMP_R_cpr", listall, get_srp_parameters, 20, defined)
   call get_logical_parameter(dict, error, "EMP_T_cpr", listall, get_srp_parameters, 21, defined)
   call get_logical_parameter(dict, error, "EMP_N_cpr", listall, get_srp_parameters, 22, defined)

   return
end function get_srp_parameters

! temp: use the math library version instead
function pow(base, idx)
   integer*4 base
   integer*4 idx

   integer*4 pow
   integer*4 i

   pow = 1

   if (idx .eq. 0) then
      return
   end if

   do i=1, idx
      pow = pow * base
   end do

   return
end function pow

subroutine get_logical_parameter(dict, error, label, listall, parameters, idx, defined)
   type (type_dictionary) :: dict
   type (type_error) :: error
   type (type_error), pointer :: e
   type (type_scalar), pointer :: s
   logical listall
   integer*4 parameters, defined
   integer*4 idx
   character(*) label
   logical value

   nullify(e)
   nullify(s)

   s = dict%get_scalar(label, listall, e)
   if (associated(e)) then
      if (listall) then
         write (*,*) "Must supply a default " // label // " srp parameter"
         error = e
         STOP
      else
         return
      end if
   end if

   value = dict%get_logical(label, .false., e)
   if (associated(e)) then
      error = e
      write (*,*) "label "//label//" found but cannot interpret as logical value"
      STOP
   endif
   defined = defined + pow(2, idx-1)

   if (value) then
      parameters = parameters + pow(2, idx-1)
   end if

   return
end subroutine get_logical_parameter

function get_gravity_model(dict, error, label)
   type (type_dictionary) :: dict
   type (type_error) :: error
   type (type_error), pointer :: e
   logical central, static_model, time_model, iers
   integer*2 get_gravity_model
   character(*) label

   nullify (e)

   get_gravity_model = -1

   central = dict%get_logical("central_force", .false., e)
   static_model = dict%get_logical("static_gravity_model", .false., e)
   time_model = dict%get_logical("time_variable_model", .false., e)
   iers = dict%get_logical("iers_geopotential_model", .false., e)

   if (central) then
      if (static_model .or. time_model .or. iers) then
         write (*,*) "Can only have one " // label // " gravity model selected"
         STOP
      end if
      get_gravity_model = 0
   else if (static_model) then
      if (time_model .or. iers) then
         write (*,*) "Can only have one " // label // " gravity model selected"
         STOP
      end if
      get_gravity_model = 1
   else if (time_model) then
      if (iers) then
         write (*,*) "Can only have one " // label // " gravity model selected"
         STOP
      end if
      get_gravity_model = 2
   else if (iers) then
      get_gravity_model = 3
   end if

   if (get_gravity_model < 0) then
      write (*,*) "Must choose a gravity model for " // label
      STOP
   end if

   if (associated(e)) then
      error = e
   end if

   return
end function get_gravity_model

function get_apriori_srp(dict, error)
   type (type_dictionary) :: dict
   type (type_error) :: error
   type (type_error), pointer :: e
   logical cannonball, simple_bw, full_bw, no_model
   integer*2 get_apriori_srp

   nullify (e)

   get_apriori_srp = -1

   cannonball = dict%get_logical("cannon_ball_model", .false., e)
   simple_bw = dict%get_logical("simple_boxwing_model", .false., e)
   full_bw = dict%get_logical("full_boxwing_model", .false., e)
   no_model = dict%get_logical("no_model", .false., e)

   if (no_model) then
      if (cannonball .or. simple_bw .or. full_bw) then
         write (*,*) "Can only have one srp_apriori_model selected"
         STOP
      end if
      get_apriori_srp = 0
   else if (cannonball) then
      if (simple_bw .or. full_bw) then
         write (*,*) "Can only have one srp_apriori_model selected"
         STOP
      end if
      get_apriori_srp = 1
   else if (simple_bw) then
      if (full_bw) then
         write (*,*) "Can only have one srp_apriori_model selected"
         STOP
      end if
      get_apriori_srp = 2 
   else if (full_bw) then
      get_apriori_srp = 3 
   end if

   if (get_apriori_srp < 0) then
      write (*,*) "Must choose a default srp apriori model (even if its 'no_model')"
      STOP
   end if

   if (associated(e)) then
      error = e
   end if

   return
end function get_apriori_srp

function get_ECOM_estimator(dict, error)
   type (type_dictionary) :: dict
   type (type_error) :: error
   type (type_error), pointer :: e
   logical ecom1, ecom2, hybrid, sboxw, none
   integer*2 get_ECOM_estimator

   nullify (e)

   ecom1 = dict%get_logical("ECOM1", .false., e)
   ecom2 = dict%get_logical("ECOM2", .false., e)
   hybrid = dict%get_logical("hybrid", .false., e)
   sboxw = dict%get_logical("SBOXW", .false., e)

   get_ECOM_estimator = 0
   if (ecom1) then
      if (ecom2 .or. hybrid .or. sboxw) then
         write (*,*) "Can choose at most 1 ECOM estimation model"
         STOP
      end if
      get_ECOM_estimator = 1
   else if (ecom2) then
      if (hybrid .or. sboxw) then
         write (*,*) "Can choose at most 1 ECOM estimation model"
         STOP
      end if
      get_ECOM_estimator = 2
   else if (hybrid) then
      if (sboxw) then
         write (*,*) "Can choose at most 1 ECOM estimation model"
         STOP
      end if
      get_ECOM_estimator = 12
   else if (sboxw) then
      get_ECOM_estimator = 3     
   end if

   if (associated(e)) then
      error = e
   end if

   return
end function get_ECOM_estimator

subroutine get_pod_data(pod_data_dict, error, pod_data_enabled, prn, ref_frame, time_scale, orbit_arc_secs,&
                init_epoch, state_vector)
   type (type_dictionary) :: pod_data_dict
   type (type_error) :: error
   type (type_error), pointer :: e
   logical pod_data_enabled
   character(*) prn, ref_frame, time_scale, init_epoch, state_vector
   character(6) var_string
   integer*8 orbit_arc_secs

   nullify (e)

   pod_data_enabled = pod_data_dict%get_logical("enable", .false., e)

   if (pod_data_enabled) then
      var_string = "TRUE"
      prn = pod_data_dict%get_string("satellite_PRN", "", e)
      ref_frame = pod_data_dict%get_string("reference_frame", "", e)
      time_scale = pod_data_dict%get_string("time_scale", "", e)
      orbit_arc_secs = pod_data_dict%get_integer("orbit_arc_length", -1, e)
      init_epoch = pod_data_dict%get_string("initial_epoch", "", e)
      state_vector = pod_data_dict%get_string("state_vector", "", e)

      write(*,*) "pod_data(enable) is "//var_string
      if (prn /= "") then
         write(*,*) "pod_data(prn) is "//prn
      else
         write(*,*) "pod_data(prn) is missing"
      end if
      if (ref_frame /= "") then
         write(*,*) "pod_data(reference frame) is "//ref_frame
      else
         write(*,*) "pod+data(reference frame) is missing"
      end if
      if (time_scale /= "") then
         write(*,*) "pod_data(time scale) is "//time_scale
      else
         write(*,*) "pod_data(time scale) is missing"
      end if
      if (init_epoch /= "") then
         write(*,*) "pod_data(initial epoch) is "//init_epoch
      else
         write(*,*) "pod_data(initial epoch) is missing"
      end if
      if (state_vector /= "") then
         write(*,*) "pod_data(state vector) is " //trim(state_vector)
      else
         write(*,*) "pod_data(state vector) is missing"
      end if
   else
      ! not specified is as good as false here
      var_string = "FALSE"
      write(*,*) "pod_data(enable) is "//var_string
   end if

   if (associated(e)) then
      error = e
   end if

   return
end subroutine get_pod_data

subroutine get_earth_orientation_params(dict, error, eop_option, filename, npoints)
   type (type_dictionary) :: dict
   type (type_error) :: error
   type (type_error), pointer :: e
   character(*) filename
   integer*4 npoints, eop_option
   logical EOP_c04, EOP_rapid, EOP_ultra_rapid

   nullify(e)

   EOP_c04 = dict%get_logical("EOP_soln_c04", .false., e)
   EOP_rapid = dict%get_logical("EOP_soln_rapid", .false., e)
   EOP_ultra_rapid = dict%get_logical("EOP_soln_igs", .false., e)

   npoints = dict%get_integer("EOP_soln_interp_points", -1, e)

   eop_option = 0
   if (EOP_c04) then
      eop_option = 1
      filename = dict%get_string("EOP_soln_c04_file", "", e)
      if (EOP_rapid .or. EOP_ultra_rapid) then
         write (*,*) "Can only choose one EOP option"
         STOP
      end if
   end if
   if (EOP_rapid) then
      eop_option = 2
      filename = dict%get_string("EOP_soln_rapid_filename", "", e)
      if (EOP_ultra_rapid) then
         write (*,*) "Can only choose one EOP option"
         STOP
      end if
   end if
   if (EOP_ultra_rapid) then
      eop_option = 3
      filename = dict%get_string("ERP_soln_igs_filename", "", e)
   end if
   if (eop_option == 0) then
      write (*,*) "Must choose an EOP option"
      STOP
   end if
   if (filename == "") then
      write (*,*) "Must specify EOP file relevant to your choice of option"
      STOP
   end if

   if (associated(e)) then
      error = e
   end if

end subroutine get_earth_orientation_params

subroutine get_extorbit_comp(dict, error, filename, steps, points)
   type (type_dictionary) :: dict
   type (type_error) :: error
   type (type_error), pointer :: e
   integer*4 :: steps, points
   character(*) filename

   nullify(e)

   filename = dict%get_string("ext_orbit_filename", "", e);
   steps = dict%get_integer("ext_orbit_interp_step", -1, e)
   points = dict%get_integer("ext_orbit_interp_points", -1, e)

   if (associated(e)) then
      error = e
   end if

   return
end subroutine get_extorbit_comp

integer*2 function get_iau_model(dict, error)
   use yaml_types
   type (type_dictionary) :: dict
   type (type_error) :: error
   type (type_error), pointer :: e
   logical iau_model_2000, iau_model_2006
   
   nullify(e)
   get_iau_model = -1

   iau_model_2000 = dict%get_logical("iau_model_2000", .false., e);
   iau_model_2006 = dict%get_logical("iau_model_2006", .false., e)

   if (associated(e)) then
      error = e
   endif

   if (iau_model_2000) then
      get_iau_model = 2000
      if (iau_model_2006) then
         write(*,*) "Can only select 1 nutation model"
         STOP
      end if
   else if (iau_model_2006) then
      get_iau_model = 2006
   end if

   if (get_iau_model == -1) then
      write (*,*) "Must select a nutation model"
      STOP
   end if

   return
end function get_iau_model

integer*2 function get_integrator_method(dict, error, stepsize)
   type (type_dictionary) :: dict
   type (type_error) :: error
   type (type_error), pointer :: e
   logical rk4, rkn7, rk8
   integer*4 stepsize

   nullify(e)
   get_integrator_method = 0

   rk4 = dict%get_logical("RK4_integrator_method", .false., e)
   rkn7 = dict%get_logical("RKN7_integrator_method", .false., e)
   rk8 = dict%get_logical("RK8_integrator_method", .false., e)

   stepsize = dict%get_integer("integrator_step", -1, e)

   if (associated(e)) then
      error = e
   end if

   if (rk4) then
      get_integrator_method = 4
      if (rkn7 .or. rk8) then
         write (*,*) "Can only choose 1 integrator method"
         STOP
      endif
   else if (rkn7) then
      get_integrator_method = 7
      if (rk8) then
         write (*,*) "Can only choose 1 integrator method"
         STOP
      end if
   else if (rk8) then
      get_integrator_method = 8
   end if

   if (get_integrator_method == 0) then
      write (*,*) "Must choose an integration method"
      STOP
   end if

   if (stepsize == -1) then
      write (*,*) "Must set the integrator stepsize"
      STOP
   end if

   return
end function get_integrator_method

integer*2 function get_ext_orbit_frame(dict, error)
   type (type_dictionary) :: dict
   type (type_error) :: error
   type (type_error), pointer :: e
   logical ext_orbit_refsys_itrf, ext_orbit_refsys_icrf

   nullify(e)
   get_ext_orbit_frame = 0

   ext_orbit_refsys_itrf = dict%get_logical("ext_orbit_frame_itrf", .false., e)
   ext_orbit_refsys_icrf = dict%get_logical("ext_orbit_frame_icrf", .false., e)

   if (associated(e)) then
      error = e
   end if

   if (ext_orbit_refsys_itrf) then
      get_ext_orbit_frame = 1
   end if
   if (ext_orbit_refsys_icrf) then
      if (get_ext_orbit_frame > 0) then
         write (*,*) "can only select one external orbit frame"
         STOP
      end if
      get_ext_orbit_frame = 2
   end if

   if (get_ext_orbit_frame == 0) then
      write (*,*) "must select an external orbit frame"
      STOP
   end if

   return
end function get_ext_orbit_frame

integer*2 function get_ext_orbit_opt(dict, error)
   type (type_dictionary) :: dict
   type (type_error) :: error
   type (type_error), pointer :: e
   logical ext_orbit_sp3, ext_orbit_interp, ext_orbit_kepler

   nullify(e)
   get_ext_orbit_opt = 0

   ext_orbit_sp3 = dict%get_logical("ext_orbit_type_sp3", .false., e)
   ext_orbit_interp = dict%get_logical("ext_orbit_type_iterp", .false., e)
   ext_orbit_kepler = dict%get_logical("ext_orbit_type_kepler", .false., e)

   if (associated(e)) then
      error = e
   end if

   if (ext_orbit_sp3) then
      get_ext_orbit_opt = 1
   end if
   if (ext_orbit_interp) then
      if (get_ext_orbit_opt > 0) then
         write (*,*) "Cannot select more than one external orbit type"
         STOP
      end if
      get_ext_orbit_opt = 2
   end if
   if (ext_orbit_kepler) then
      if (get_ext_orbit_opt > 0) then
         write (*,*) "Cannot select more than one external orbit type"
         STOP
      end if
      get_ext_orbit_opt = 3
   end if
   if (get_ext_orbit_opt == 0) then
      write (*,*) "Must select an external orbit type"
      STOP
   end if

   ! Options 1 and 3 no longer valid ...
   if (get_ext_orbit_opt /= 2) then
      if (get_ext_orbit_opt == 1) then
         write (*,*) "ext_orbit_type_sp3 :"
      else if (get_ext_orbit_opt == 3) then
         write (*,*) "ext_orbit_type_kepler :"
      endif
      write (*,*) "Option discontinued. Please enable 'ext_orbit_type_iterp'"
      STOP
   end if
   return
end function get_ext_orbit_opt

integer*2 function get_pod_mode(dict, error)
   type(type_dictionary) :: dict
   type (type_error) :: error
   type (type_error), pointer :: e
   logical pod_mode_fit, pod_mode_predict, pod_mode_eqm_int, pod_mode_ic_int

   nullify(e)
   get_pod_mode = -1

   pod_mode_fit = dict%get_logical("pod_mode_fit", .false., e)
   !TODO check for error
   pod_mode_predict = dict%get_logical("pod_mode_predict", .false., e)
   !TODO check for error
   pod_mode_eqm_int = dict%get_logical("pod_mode_eqm_int", .false., e)
   !TODO check for error
   pod_mode_ic_int = dict%get_logical("pod_mode_ic_int", .false., e)
   !TODO check for error

   if (associated(e)) then
      error = e
   end if

   if (pod_mode_fit) then
      get_pod_mode = 0
   end if
   if (pod_mode_predict) then
      if (get_pod_mode >= 0) then
         write (*,*) "Cannot select more than one pod mode"
         STOP
      end if
      get_pod_mode = 1
   end if
   if (pod_mode_eqm_int) then
      if (get_pod_mode >= 0) then
         write (*,*) "Cannot select more than one pod mode"
         STOP
      end if
      get_pod_mode = 2
   end if
   if (pod_mode_ic_int) then
      if (get_pod_mode >= 0) then
         write (*,*) "Cannot select more than one pod mode"
         STOP
      end if
      get_pod_mode = 3
   end if
   if (get_pod_mode < 0) then
      write (*,*) "Must select a pod mode"
      STOP
   end if

   return
end function get_pod_mode

integer*2 function get_input_format(dict, error)
   type(type_dictionary) :: dict
   type (type_error) :: error
   type (type_error), pointer :: e
   logical ic_input_sp3, ic_input_icf

   nullify(e)
   get_input_format = -1

   ic_input_sp3 = dict%get_logical("ic_input_format_sp3", .false., e)
   !TODO check for error
   ic_input_icf = dict%get_logical("ic_input_format_icf", .false., e)
   !TODO check for error

   if (associated(e)) then
      error = e
   end if

   if (ic_input_sp3) then
      get_input_format = 0
   end if
   if (ic_input_icf) then
      if (get_input_format >= 0) then
         write (*,*) "Cannot select more than one input format type"
         STOP
      else
         get_input_format = 1
      end if
   end if
   if (get_input_format < 0) then
      write (*,*) "Must select an input format type"
      STOP
   end if

   return
end function get_input_format

integer*2 function get_input_refsys(dict, error)
   type(type_dictionary) :: dict
   type (type_error) :: error
   type (type_error), pointer :: e
   logical ic_input_itrf, ic_input_icrf

   get_input_refsys = -1

   ic_input_itrf = dict%get_logical("ic_input_refsys_itrf", .false., e)
   !TODO check for error
   ic_input_icrf = dict%get_logical("ic_input_refsys_icrf", .false., e)
   !TODO check for error

   if (associated(e)) then
      error = e
   end if

   if (ic_input_itrf) then
      get_input_refsys = 0
   end if
   if (ic_input_icrf) then
      if (get_input_refsys >= 0) then
         write (*,*) "Cannot select more than one input refsys type"
         STOP
      else
         get_input_refsys = 1
      end if
   end if
   if (get_input_refsys < 0) then
      write (*,*) "Must select an input format type"
      STOP
   end if

   return
end function get_input_refsys

subroutine get_pseudoobs(dict, error, filename, step, points)
   type(type_dictionary) :: dict
   type (type_error) :: error
   character(*) filename
   integer*4 step
   integer*4 points
   type (type_error), pointer :: e

   nullify(e)

   filename = dict%get_string("pseudoobs_orbit_filename", "", e);
   step = dict%get_integer("pseudoobs_interp_step", -1, e)
   points = dict%get_integer("pseudoobs_interp_points", -1, e)

   if (associated(e)) then
      error = e
   end if

   return
end subroutine get_pseudoobs

subroutine get_orbitarcs(dict, error, determination, prediction, backwards)
   type (type_dictionary) :: dict
   type (type_error) :: error
   integer*4 determination, prediction, backwards
   type (type_error), pointer :: e

   nullify(e)

   determination=dict%get_integer("orbit_arc_determination", -1, e)
   prediction=dict%get_integer("orbit_arc_prediction", -1, e)
   backwards=dict%get_integer("orbit_arc_backwards", -1, e)

   if (associated(e)) then
      error = e
   end if

   return
end subroutine get_orbitarcs

end module pod_yaml

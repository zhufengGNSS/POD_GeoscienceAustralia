cmake_minimum_required(VERSION 2.8)

project(pod)

#project (add_vectors Fortran)
enable_language(Fortran)

#set (CMAKE_Fortran_FLAGS "-g -O3 -frecursive -fall-intrinsics -Wmaybe-uninitialized -fcheck=all -fbacktrace")
set (CMAKE_Fortran_FLAGS "-O3 -frecursive -fall-intrinsics -Winteger-division -Wmaybe-uninitialized -Wunused-dummy-argument -fcheck=all -finit-local-zero -fbacktrace")
find_package(LAPACK REQUIRED)
find_package(BLAS REQUIRED)
#message("Fortran default extensions = ${CMake_Fortran_SOURCE_FILE_EXTENSIONS}")
set (CMAKE_Fortran_SOURCE_FILE_EXTENSIONS "f03;F90;f90;for;f;F")
file(GLOB_RECURSE sources src/*.f90 src/*.F90 src/*.f03 src/*.f src/*.for src/*.F)
foreach (file ${sources})
	execute_process(COMMAND /usr/bin/basename -z ${file} OUTPUT_VARIABLE srcfile)
	set_source_files_properties(src/${srcfile} PROPERTIES LANGUAGE Fortran)
	#set_source_files_properties(src/${srcfile} PROPERTIES COMPILE_FLAGS "-O3 -frecursive")
	#message("set src/${srcfile} to fortran compile")
endforeach (file ${sources})

add_library(pod_lib
src/mdl_precision.f90 
src/mdl_num.f90 
src/mdl_config.f03 
src/mdl_param.f03 
src/mdl_eop.f90 
src/mdl_planets.f90 
src/mdl_tides.f90 
src/mdl_arr.f90 
src/m_integrEQM.f03 
src/m_matrixRxR.f03 
src/m_legendre.f03 
src/m_legendre1.f03 
src/m_force_tides.f03 
src/m_tides_ocean.f03 
src/m_legendre2.f03 
src/m_pd_geopotential.f03 
src/m_force_gfm.f03 
src/pd_forceZ.f03 
src/m_get_lambda.f 
src/m_shadow.f90 
src/m_pd_empirical.f03 
src/m_pd_ECOM.f90 
src/m_satinfo.f90 
src/m_pd_force.f03 
src/m_veq_rkn768.f03 
src/m_veqC2T.f03
src/m_integrVEQ.f03 
src/m_betainfo.f90 
src/m_orbinteg.f03 
src/m_matrixinv.f03 
src/m_matrixreverse.f03 
src/m_matrixmerge.f03 
src/m_orb_estimator.f03 
src/m_orbC2T.f03 
src/m_orbT2C.f03 
src/m_statist.f03 
src/m_statdelta.f03 
src/m_orbinfo.f90 
src/m_statorbit.f03 
src/m_statorbit2.f03 
src/m_writearray.f03 
src/m_writearray2.f03 
src/m_orbresize.f03 
src/m_orbext.f03 
src/m_orbext2.f03 
src/m_gfc.f03 
src/m_gfc3.f03 
src/m_obsorbT2C.f03 
src/m_write_prmfile_init.f03 
src/m_writedata.f03 
src/m_writeorbit.f03 
src/m_read_svsinex.f03
src/m_read_satsnx.f03
src/m_writeorbit_multi.f03 
src/m_write_orb2sp3.f03 
src/m_rso.f03 
src/m_sp3.f03 
src/m_sp3_PRN.f03 
src/m_lagrange.f03 
src/m_interporb.f03 
src/m_eop_cor.f03 
src/m_eop_igu.f03 
src/m_eop_data.f03 
src/m_keplerorb.f03 
src/m_orbdet.f03 
src/m_orbitmain.f03 
src/m_orbitIC.f03
src/m_pod_gnss.f03 
src/gblinit.f03
)

add_executable(pod
src/m_read_leapsec.f03
src/m_clock_read.f03
src/write_prmfile_init0.f03 
src/writeparam.f03 
src/writeparam1.f03 
src/write_prmfile.f03 
src/write_prmfile2.f03 
src/readparam.f03  
src/arctan.f90 
src/productdot.f90 
src/productcross.f90 
src/coord_r2sph.f90 
src/matrix_inv3.f90 
src/matrix_Rr.f90 
src/matrix_RxR.f90 
src/jd2cal.for 
src/cal2jd.for 
#src/dat.for 
src/fad03.for 
src/faf03.for 
src/fal03.for 
src/fama03.for 
src/fane03.for  
src/fapa03.for 
src/faur03.for 
src/fae03.for  
src/faju03.for 
src/falp03.for 
src/fame03.for 
src/faom03.for 
src/fasa03.for 
src/fave03.for 
src/xy06.for 
src/s06.for 
src/c2ixys.for 
src/xys00a.for 
src/pnm00a.for 
src/bpn2xy.for 
src/s00.for 
src/pn00a.for 
src/nut00a.for 
src/pn00.for 
src/pr00.for 
src/obl80.for 
src/bp00.for 
src/numat.for 
src/bi00.for 
src/ir.for 
src/rz.for 
src/ry.for 
src/rx.for 
src/anp.for 
src/era00.for 
src/gmst00.for 
src/gmst06.for 
src/sp00.for 
src/pom00.for 
src/taiutc.for
src/ORTHO_EOP.F 
src/CNMTX.F 
src/UTLIBR.F 
src/PMSDNUT2.F 
src/FUNDARG.F 
src/RG_ZONT2.F 
src/interp_iers.f 
src/rxr.for 
src/tr.for 
src/cr.for 
src/cp.for 
src/time_TT.f90 
src/time_GPS.f90 
src/time_UTC.f90 
src/time_TAI.f90 
src/time_GPSweek.f90 
src/time_GPSweek3.f90
src/time_TT_sec.f90 
src/crs_trs.f90 
src/eop_rd.f90 
src/eop_c04.f90 
src/eop_finals2000A.f90 
src/erp_igu.f90 
src/interp_lin.f90 
src/EOP.f90 
src/era_matrix.f90 
src/orb_frame.f90 
src/crf_bff.f90 
src/kepler_eq.f90 
src/kepler_k2z.f90 
src/kepler_z2k.f90 
src/force_gm.f90 
src/GM_de.f90 
src/CATfile.f90 
src/asc2eph.f90 
src/STATE.f90 
src/PLEPH.f 
src/CONST.f 
src/SPLIT.f 
src/INTERP.f 
src/ERPFBOXW.f
src/SURFBOXW.f
src/PROPBOXW.f
src/SRPFBOXW.f
src/FSIZER3.f 
src/force_gm3rd.f90 
src/indirectJ2.f90 
src/delaunay.f90
src/gmst_iers.f03 
src/tides_solid1.f90 
src/tides_solid2.f90 
src/tide_perm.f90 
src/tides_fes2004.f90 
src/IERS_CMP_2015.F 
src/tide_pole_se.f90 
src/tide_pole_oc.f90 
src/rel_schwarzschild.f90 
src/rel_LenseThirring.f90 
src/rel_deSitter.f90 
src/force_srp.f90 
src/prn_shift.f03 
src/surfprop.f90 
src/cross_product.f90 
src/R3.for 
src/R1.for 
src/scan0orb.f03
src/force_erp.f90 
src/force_ant.f90 
src/empirical_init_file.f03 
src/empirical_init.f03 
src/empirical_cor.f03 
src/beta_angle.f90 
src/yaw_nom.f90
src/yaw_angle.f90 
src/yaw_attitude.f90 
src/eclips.f 
src/eclips2017.f 
src/eclips201707.f 
src/beta_pred.f90 
src/yaw_bds.f90 
src/yaw_gal.f90 
src/yawdyn.f90 
src/attitude.f03 
src/force_sum.f03 
src/integr_rkn768.f03 
src/integr_rk87.f03 
src/integr_rk4.f03 
src/pd_gm.f03 
src/pd_forceZ.f03 
src/eclipse_integstep.f03 
src/write_prmfile.f03 
src/writeparam.f03 
src/writeparam1.f03 
src/time_TT_sec.f90 
src/prm_main.f03 
src/prm_emp.f03 
src/empirical_cor.f03 
src/empirical_init.f03 
src/ecom_init.f03 
src/prm_srp.f03 
src/prm_grav.f03 
src/prm_gravity.f03 
src/prm_nongrav.f03 
src/prm_ocean.f03 
src/prm_orbext.f03 
src/prm_planets.f03 
src/prm_pseudobs.f03 
src/f90getopt.F90 
src/read_cmdline.f03 
src/main_pod.f03
)

link_directories(/usr/lib64)

target_include_directories(pod PUBLIC
        src/
        ${LAPACK_INCLUDE_DIRS}
        ${BLAS_INCLUDE_DIRS})

target_link_libraries(pod pod_lib ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})

set_target_properties(pod PROPERTIES LINKER_LANGUAGE Fortran)

target_compile_definitions(pod PRIVATE LAPACK BLAS)

set(CMAKE_VERBOSE_MAKEFILE ON)

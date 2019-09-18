BIN = ../bin
# Fortan compiler
#FC = /opt/intel/bin/ifort
FC = gfortran

# Compiler flags

FCFLAGS = -O3 -llapack -lblas
#-O3 -Wuninitialized -fno-f2c -ffast-math -fno-automatic -fno-backslash -m64 -mcmodel=large 
# Linker flags
FLFLAGS = -L/usr/lib64

SRCS = mdl_precision.f90 mdl_num.f90 mdl_config.f03 mdl_param.f03 mdl_eop.f90 mdl_planets.f90 mdl_tides.f90 mdl_arr.f90 \
m_integrEQM.f03 m_matrixRxR.f03 m_legendre.f03 m_legendre1.f03 m_force_tides.f03 m_tides_ocean.f03 m_legendre2.f03 \
m_pd_geopotential.f03 m_force_gfm.f03 pd_forceZ.f03 m_get_lambda.f m_shadow.f90 m_pd_empirical.f03 m_pd_ECOM.f90 \
m_satinfo.f90 m_pd_force.f03 m_veq_rkn768.f03 m_veqC2T.f03 m_integrVEQ.f03 m_betainfo.f90 m_orbinteg.f03 m_matrixinv.f03 m_orb_estimator.f03 \
m_orbC2T.f03 m_orbT2C.f03 m_statist.f03 m_statdelta.f03 m_orbinfo.f90 m_statorbit.f03 m_statorbit2.f03 m_writearray.f03 m_writearray2.f03 m_orbresize.f03 \
m_orbext.f03 m_orbext2.f03 m_gfc.f03 m_gfc3.f03 m_obsorbT2C.f03 m_write_prmfile_init.f03 m_writedata.f03 m_writeorbit.f03 m_writeorbit_multi.f03 \
m_write_orb2sp3.f03 m_rso.f03 m_sp3.f03 m_sp3_PRN.f03 m_lagrange.f03 m_interporb.f03 m_eop_cor.f03 m_eop_igu.f03 m_eop_data.f03 m_keplerorb.f03 \
m_orbdet.f03 m_orbitmain.f03 m_pod_gnss.f03 m_read_leapsec.f03\
write_prmfile_init0.f03 writeparam.f03 writeparam1.f03 write_prmfile.f03 write_prmfile2.f03 readparam.f03  \
arctan.f90 productdot.f90 productcross.f90 coord_r2sph.f90 matrix_inv3.f90 matrix_Rr.f90 matrix_RxR.f90 \
jd2cal.for cal2jd.for taiutc.for fad03.for faf03.for fal03.for fama03.for fane03.for  fapa03.for faur03.for \
fae03.for  faju03.for falp03.for fame03.for faom03.for fasa03.for fave03.for xy06.for s06.for c2ixys.for \
xys00a.for pnm00a.for bpn2xy.for s00.for pn00a.for nut00a.for pn00.for pr00.for obl80.for bp00.for numat.for bi00.for \
ir.for rz.for ry.for rx.for anp.for era00.for gmst00.for gmst06.for sp00.for pom00.for \
ORTHO_EOP.F CNMTX.F UTLIBR.F PMSDNUT2.F FUNDARG.F RG_ZONT2.F interp_iers.f rxr.for tr.for cr.for cp.for time_TT.f90 time_GPS.f90 \
time_UTC.f90 time_TAI.f90 time_GPSweek.f90 time_TT_sec.f90 crs_trs.f90 eop_rd.f90 eop_c04.f90 eop_finals2000A.f90 erp_igu.f90 interp_lin.f90 \
EOP.f90 era_matrix.f90 orb_frame.f90 crf_bff.f90 kepler_eq.f90 kepler_k2z.f90 kepler_z2k.f90 \
force_gm.f90 GM_de.f90 CATfile.f90 asc2eph.f90 STATE.f90 PLEPH.f CONST.f SPLIT.f INTERP.f FSIZER3.f \
force_gm3rd.f90 indirectJ2.f90 delaunay.f90 gmst_iers.f03 tides_solid1.f90 tides_solid2.f90 tide_perm.f90 \
tides_fes2004.f90 IERS_CMP_2015.F tide_pole_se.f90 tide_pole_oc.f90 rel_schwarzschild.f90 rel_LenseThirring.f90 rel_deSitter.f90 \
force_srp.f90 prn_shift.f03 surfprop.f90 cross_product.f90 R3.for R1.for force_erp.f90 force_ant.f90 \
empirical_init_file.f03 empirical_init.f03 empirical_cor.f03 beta_angle.f90 yaw_nom.f90 yaw_angle.f90 \
yaw_attitude.f90 eclips.f eclips2017.f eclips201707.f beta_pred.f90 yaw_bds.f90 yaw_gal.f90 yawdyn.f90 \
attitude.f03 force_sum.f03 integr_rkn768.f03 integr_rk87.f03 integr_rk4.f03 pd_gm.f03 pd_forceZ.f03 eclipse_integstep.f03 \
write_prmfile.f03 writeparam.f03 writeparam1.f03 time_TT_sec.f90 prm_main.f03 prm_emp.f03 empirical_cor.f03 empirical_init.f03 \
ecom_init.f03 prm_srp.f03 prm_grav.f03 prm_gravity.f03 prm_nongrav.f03 prm_ocean.f03 scan0orbext.f03 prm_orbext.f03 prm_planets.f03 scan0orb.f03 prm_pseudobs.f03 \
f90getopt.F90 read_cmdline.f03 \
main_pod.f03

CRS2TRSSRCS = mdl_precision.f90 mdl_num.f90 mdl_param.f03 mdl_eop.f90 mdl_arr.f90 m_writearray.f03 \
arctan.f90 productdot.f90 productcross.f90 coord_r2sph.f90 \
m_matrixinv.f03 matrix_Rr.f90 matrix_RxR.f90 m_matrixRxR.f03 \
jd2cal.for cal2jd.for taiutc.for \
fad03.for  faf03.for   fal03.for   fama03.for  fane03.for  fapa03.for  faur03.for \
fae03.for  faju03.for  falp03.for  fame03.for  faom03.for  fasa03.for  fave03.for \
xy06.for s06.for c2ixys.for \
xys00a.for pnm00a.for bpn2xy.for s00.for pn00a.for nut00a.for pn00.for pr00.for obl80.for bp00.for numat.for bi00.for \
ir.for rz.for ry.for rx.for \
anp.for era00.for gmst00.for gmst06.for \
sp00.for pom00.for \
ORTHO_EOP.F CNMTX.F UTLIBR.F PMSDNUT2.F FUNDARG.F RG_ZONT2.F \
interp_iers.f \
rxr.for tr.for cr.for cp.for \
time_TT.f90 time_GPS.f90 time_UTC.f90 time_TAI.f90 time_GPSweek.f90 \
crs_trs.f90 m_eop_data.f03 m_eop_cor.f03 m_eop_igu.f03 eop_rd.f90 eop_c04.f90 eop_finals2000A.f90 erp_igu.f90 interp_lin.f90 \
era_matrix.f90 \
main_crs2trs.f03

BRDC2ECEF = mdl_precision.f90 mdl_num.f90 mdl_param.f03 mdl_eop.f90 mdl_arr.f90 m_writearray.f03 \
brdc_cmdline.f03 mdl_brdconfig.f03 f90getopt.F90 \
arctan.f90 productdot.f90 productcross.f90 coord_r2sph.f90 \
m_matrixinv.f03 matrix_Rr.f90 matrix_RxR.f90 m_matrixRxR.f03 \
m_write_brd2sp3.f03 prn2str.f03 prn_shift_brdc.f03 chkbrdc.f \
readbrdcheader.f03 readbrdcmessg.f90 reformbrdc.f brdc2ecef.f \
m_read_leapsec.f03 taiutc.for time_UTC.f90 time_GPSweek.f90 jd2cal.for cal2jd.for time_GPSweek2.f90 \
brdc_qc_gps.f90 brdc_qc_gal.f90 m_meanstd.f90 \
m_antoffset.f90 glnorbint.f90 glnacc.f \
\
EOP.f90 \
fad03.for  faf03.for   fal03.for   fama03.for  fane03.for  fapa03.for  faur03.for \
fae03.for  faju03.for  falp03.for  fame03.for  faom03.for  fasa03.for  fave03.for \
xy06.for s06.for c2ixys.for \
xys00a.for pnm00a.for bpn2xy.for s00.for pn00a.for nut00a.for pn00.for pr00.for obl80.for bp00.for numat.for bi00.for \
ir.for rz.for ry.for rx.for \
anp.for era00.for gmst00.for gmst06.for \
sp00.for pom00.for \
ORTHO_EOP.F CNMTX.F UTLIBR.F PMSDNUT2.F FUNDARG.F RG_ZONT2.F \
interp_iers.f \
rxr.for tr.for cr.for cp.for \
time_TT.f90 time_GPS.f90 time_TAI.f90 \
crs_trs.f90 m_eop_data.f03 m_eop_cor.f03 m_eop_igu.f03 eop_rd.f90 eop_c04.f90 eop_finals2000A.f90 erp_igu.f90 interp_lin.f90 \
era_matrix.f90 \
\
CATfile.f90 asc2eph.f90 mdl_planets.f90 STATE.f90 PLEPH.f CONST.f SPLIT.f INTERP.f FSIZER3.f \
\
main_brdcorbit.f90


PROGRAM = pod

all: $(BIN) $(PROGRAM) crs2trs brd2ecef

directories: $(BIN)

$(BIN): 
	mkdir -p $(BIN)

$(PROGRAM): $(SRCS)  
	$(FC) $(FCFLAGS) -o $(BIN)/$@ $^

%.o: %.F
	$(FC) $(FLFLAGS) -o $@ $<

%.mod: %.h
	$(FC) $(FLFLAGS) -o $@ $<

crs2trs: $(CRS2TRSSRCS) $(BIN)
	$(FC) $(FCFLAGS) $(FLFLAGS) $(CRS2TRSSRCS) -o $(BIN)/crs2trs
brd2ecef: $(BRDC2ECEF) $(BIN)
	$(FC) $(FCFLAGS) $(FLFLAGS) $(BRDC2ECEF) -o $(BIN)/brd2ecef

clean:
	rm -f *.o *.mod $(BIN)/$(PROGRAM) $(BIN)/crs2trs


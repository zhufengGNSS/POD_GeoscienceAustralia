BIN = ./bin
# Fortan compiler
#FC = /opt/intel/bin/ifort
FC = gfortran

# Compiler flags
FCFLAGS = -llapack -lblas

# Linker flags
FLFLAGS = -L/usr/lib64

SRCS = mdl_precision.f90 mdl_num.f90 mdl_param.f03 mdl_eop.f90 mdl_planets.f90 mdl_tides.f90 mdl_arr.f90 m_writearray.f03 \
arctan.f90 productdot.f90 productcross.f90 coord_r2sph.f90 \
m_matrixinv.f03 matrix_Rr.f90 matrix_RxR.f90 m_matrixRxR.f03 \
jd2cal.for cal2jd.for dat.for \
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
m_eop_cor.f03 m_eop_igu.f03 m_eop_data.f03\
crs_trs.f90 eop_rd.f90 eop_c04.f90 eop_finals2000A.f90 erp_igu.f90 interp_lin.f90 EOP.f90 \
era_matrix.f90 \
orb_frame.f90 \
kepler_eq.f90 kepler_k2z.f90 kepler_z2k.f90 m_keplerorb.f03 \
m_rso.f03 m_sp3.f03 m_gnssp3.f03 m_sat_ini_vet.f03 \
m_lagrange.f03 m_interporb.f03 \
force_gm.f90 \
m_legendre.f03 m_legendre1.f03 m_force_gfm.f03 \
GM_de.f90 CATfile.f90 asc2eph.f90 STATE.f90 PLEPH.f CONST.f SPLIT.f INTERP.f FSIZER3.f \
force_gm3rd.f90 indirectJ2.f90 \
delaunay.f90 gmst_iers.f03 \
tides_solid1.f90 tides_solid2.f90 tide_perm.f90 \
tides_fes2004.f90 m_tides_ocean.f03 \
IERS_CMP_2015.F  tide_pole_se.f90 tide_pole_oc.f90 \
m_force_tides.f03 \
rel_schwarzschild.f90 rel_LenseThirring.f90 rel_deSitter.f90 \
force_srp.f90 prn_shift.f03 surfprop.f90 cross_product.f90 R3.for R1.for \
force_sum.f03 \
integr_rkn768.f03 integr_rk87.f03 integr_rk4.f03 \
pd_gm.f03 m_legendre2.f03 m_pd_geopotential.f03 \
pd_forceZ.f03 \
m_veq_rkn768.f03 \
m_integrEQM.f03 m_integrVEQ.f03 m_orbinteg.f03 \
m_orb_estimator.f03 \
m_gfc.f03 m_gfc3.f03 \
m_orbC2T.f03 m_orbT2C.f03 m_obsorbT2C.f03 \
prm_read.f03 prm_gravity.f03 prm_planets.f03 prm_ocean.f03 \
prm_orbext.f03 prm_pseudobs.f03 \
prm_main.f03 prm_grav.f03 prm_nongrav.f03 \
m_statist.f03 m_statdelta.f03 m_statorbit.f03 \
m_orbdet.f03 m_orbext.f03 \
main_orb.f03

CRS2TRSSRCS = mdl_precision.f90 mdl_num.f90 mdl_param.f03 mdl_eop.f90 mdl_arr.f90 m_writearray.f03 \
arctan.f90 productdot.f90 productcross.f90 coord_r2sph.f90 \
m_matrixinv.f03 matrix_Rr.f90 matrix_RxR.f90 m_matrixRxR.f03 \
jd2cal.for cal2jd.for dat.for \
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

PROGRAM = main_orb.e

all: $(PROGRAM) crs2trs

directories: $(BIN)

$(BIN): 
	mkdir -p $(BIN)

$(PROGRAM): $(SRCS)
	$(FC) $(FCFLAGS) -o $@ $^

%.o: %.F
	$(FC) $(FLFLAGS) -o $@ $<

%.mod: %.h
	$(FC) $(FLFLAGS) -o $@ $<

crs2trs: $(CRS2TRSSRCS) $(BIN)
	$(FC) $(FCFLAGS) $(FLFLAGS) $(CRS2TRSSRCS) -o $(BIN)/crs2trs

clean:
	rm -f *.o *.mod $(PROGRAM) $(BIN)/crs2trs


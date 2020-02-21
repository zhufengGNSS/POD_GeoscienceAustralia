#!/bin/bash

# Test script for POD install

# Remove old test files
\rm -r gag19424.sp3 pod.out >& /dev/null

# Create necessary links
ln -s ../tables/ascp1950.430 .             >& /dev/null
ln -s ../tables/fes2004_Cnm-Snm.dat .      >& /dev/null
ln -s ../tables/goco05s.gfc .              >& /dev/null
ln -s ../tables/header.430_229 .           >& /dev/null
ln -s ../tables/igs_metadata_2063.snx .    >& /dev/null
ln -s ../tables/leap.second .              >& /dev/null
ln -s ../tables/eopc04_14_IAU2000.62-now . >& /dev/null

# Run the POD
../bin/pod | tee pod.out

echo 'Plotting: 1] orbit fit residual time series'
python3 ../scripts/res_plot.py     -i gag19424_igs19424_orbdiff_rtn.out -d . -c G >& /dev/null
echo '          2] orbit fit statistics'
python3 ../scripts/rms_bar_plot.py -i gag19424_igs19424_orbdiff_rtn.out -d . -c G >& pod.rms

# Diff output
diff pod.out pod.out.good
diff pod.rms pod.rms.good
diff gag19424.sp3 gag19424.sp3.good

# Cleanup 
\rm -r EQM0* VEQ0* ECOM1_*.in emp*.in DE.430 

exit 

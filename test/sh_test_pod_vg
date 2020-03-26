#!/bin/bash

# Test script for POD install

# Remove old test output files
\rm -r gag19424.sp3 pod.out >& /dev/null

# Create necessary links
ln -s ../tables/ascp1950.430 .        >& /dev/null
ln -s ../tables/fes2004_Cnm-Snm.dat . >& /dev/null
ln -s ../tables/goco05s.gfc .         >& /dev/null
ln -s ../tables/header.430_229 .      >& /dev/null
ln -s ../tables/igs_metadata_2063.snx . >& /dev/null
ln -s ../tables/leap.second .           >& /dev/null
ln -s ../tables/eopc04_14_IAU2000.62-now . >& /dev/null

# Run the POD
valgrind --tool=memcheck -s --leak-check=full --show-leak-kinds=all --track-origins=yes --log-file=valgrind.errors ../bin/pod > pod.out

# Diff output
diff pod.out pod.out.good
diff gag19424.sp3 gag19424.sp3.good

# Cleanup 
\rm -r EQM0* VEQ0* ECOM1_*.in emp*.in DE.430 

exit 

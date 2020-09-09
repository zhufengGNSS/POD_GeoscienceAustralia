import os
import argparse
import shutil
import urllib.request

import re
#===============================================================================
# Command line argument
parser = argparse.ArgumentParser(description="Download POD examples")

parser.add_argument('-d',dest='DestDir',type=str, required=False, default='/data/acs/pod/examples', help='Destination directory path (default: /data/acs/pod/examples/)')
args = parser.parse_args()

# This will remove a trailing slash, which we can add on later
# when we define the download and detsination paths
destDir = os.path.dirname(args.DestDir)
#===============================================================================
# 
#===============================================================================
def downloadfiles(srcdir,destdir,srcfiles):

    # Get all of the data subdirectory
    for filename in srcfiles:
        url = srcdir + filename
        destfile = destdir + filename
    
        if not os.path.exists(os.path.dirname(destfile)):
            os.makedirs(os.path.dirname(destfile))
            print("Making the directories needed for: ",destfile)
        
        # Check to see if the file exists already before trying to download            
        if(not os.path.isfile(destfile) ):
            print("Downloading to:<",destfile,">")
        
            # Download the file from `url` and save it locally under `file_name`:
            with urllib.request.urlopen(url) as response, open(destfile, 'wb') as out_file:
                shutil.copyfileobj(response, out_file)
        else:
            print("Already have file, not downloading: ",destfile)

#===============================================================================
# Define source paths
#===============================================================================
S3bucket = 'https://peanpod.s3-ap-southeast-2.amazonaws.com/pod_test/'

tablesSrcDir = S3bucket+'tables/'

ex1SrcDir = S3bucket + 'examples/ex1/'
ex1SrcSolDir = ex1SrcDir + 'solution/'

ex2SrcDir = S3bucket + 'examples/ex2/'
ex2SrcSolDir = ex2SrcDir + 'solution/'

ex3SrcDir = S3bucket + 'examples/ex3/'
ex3SrcSolDir = ex3SrcDir + 'solution/'

ex4SrcDir = S3bucket + 'examples/ex4/'
ex4SrcSolDir = ex4SrcDir + 'solution/'

ex5SrcDir = S3bucket + 'examples/ex5/'
ex5SrcSolDir = ex5SrcDir + 'solution/'

#===============================================================================
# Form destination paths
#===============================================================================
tablesDestDir = destDir + '/tables/' 
ex1DestDir = destDir + '/ex1/'
ex1DestSolDir = ex1DestDir + 'solution/' 

ex2DestDir = destDir + '/ex2/'
ex2DestSolDir = ex2DestDir + 'solution/' 

ex3DestDir = destDir + '/ex3/'
ex3DestSolDir = ex3DestDir + 'solution/' 

ex4DestDir = destDir + '/ex4/'
ex4DestSolDir = ex4DestDir + 'solution/' 

ex5DestDir = destDir + '/ex5/'
ex5DestSolDir = ex5DestDir + 'solution/' 
#===============================================================================
# Define the files to be downloaded
#===============================================================================
#ex03Datafile - Here are some hints on how to recreate the list if the data is updated :
#ls -1 /data/acs/pea/examples/EX03/data/ | xargs -L1 -I% echo "'%'," > data.listingt
#  perl -i -pe 's%\n%%g' data.listing
# Then edit data.listing to remove the last comma and put in the brackets at the start and en dof the list, then copy and paste below 
#
# table files
#
tablesFiles = ('ascp1950.430','fes2004_Cnm-Snm.dat','goco05s.gfc','header.430_229','igs_metadata_2063.snx','leap.second', 'eopc04_IAU2000.62-now')
#
# EX1 files
#
ex1Files = ('EQM.in','POD.in','VEQ.in','igs19424.sp3')
ex1SolFiles = ('pod.out','pod.rms')
#
# EX2 files
#
ex2Files = ('EQM.in','POD.in','VEQ.in','igs19424.sp3','COD0MGXFIN_20191990000_01D_05M_ORB.SP3')
ex2SolFiles = ('pod.out','pod_C.rms','pod_E.rms','pod_G.rms','pod_R.rms')
#
# EX3 files
#
ex3Files = ('EQM.in','POD.in','VEQ.in','igs20010.sp3', 'igs20011.sp3')
ex3SolFiles = ('pod.out','pod.rms')
#
# EX4 files
#
ex4Files = ('EQM.in','POD.in','VEQ.in','igs20624.sp3')
ex4SolFiles = ('pod.out','pod.rms')
#
# EX5 files
#
ex5Files = ('EQM.in','POD.in','VEQ.in','igs19424.sp3')
ex5SolFiles = ('pod.out','pod.rms')

#       
#===============================================================================
# Syantax for downlfiles:
# downloadfiles(src_directory,dest_directory,list_of_files_to_download)
downloadfiles(tablesSrcDir,tablesDestDir,tablesFiles)
#
downloadfiles(ex1SrcDir,ex1DestDir,ex1Files)
downloadfiles(ex1SrcSolDir,ex1DestSolDir,ex1SolFiles)
#
downloadfiles(ex2SrcDir,ex2DestDir,ex2Files)
downloadfiles(ex2SrcSolDir,ex2DestSolDir,ex2SolFiles)
#
downloadfiles(ex3SrcDir,ex3DestDir,ex3Files)
downloadfiles(ex3SrcSolDir,ex3DestSolDir,ex3SolFiles)
#
downloadfiles(ex4SrcDir,ex4DestDir,ex4Files)
downloadfiles(ex4SrcSolDir,ex4DestSolDir,ex4SolFiles)
#
downloadfiles(ex5SrcDir,ex5DestDir,ex5Files)
downloadfiles(ex5SrcSolDir,ex5DestSolDir,ex5SolFiles)

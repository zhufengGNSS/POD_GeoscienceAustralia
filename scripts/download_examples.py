import os
import argparse
import shutil
import urllib.request

import re
#===============================================================================
# Command line argument
parser = argparse.ArgumentParser(description="Download POD examples")

parser.add_argument('-d',dest='DestDir',type=str, required=False, default='/data/acs/pod/examples/', help='Destination directory path (default: /data/acs/pod/examples/)')
args = parser.parse_args()

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
# Download example EX03
#===============================================================================
S3bucket = 'https://peanpod.s3-ap-southeast-2.amazonaws.com/pod/examples/'

ex1SrcDir = S3bucket + 'ex1/'
ex1SrcSolDir = ex1SrcDir + 'solution/'

ex1DestDir = args.DestDir + 'ex1/'
ex1DestSolDir = ex1DestDir + 'solution/' 

#ex03Datafile - Here are some hints on how to recreate the list if the data is updated :
#ls -1 /data/acs/pea/examples/EX03/data/ | xargs -L1 -I% echo "'%'," > data.listingt
#  perl -i -pe 's%\n%%g' data.listing
# Then edit data.listing to remove the last comma and put in the brackets at the start and en dof the list, then copy and paste below 
ex1Files = ('EQM.in','POD.in','VEQ.in','ex1.clean','igs19424.sp3','sh_ex1')
ex1SolFiles = ('gag19424.sp3','gag19424_igs19424_orbdiff_rtn.out','gag19424_igs19424_orbitstat_N.out','gag19424_igs19424_orbitstat_R.out','gag19424_igs19424_orbitstat_T.out','gag19424_orbits_partials.out','orb_icrf.out','orb_itrf.out','orbres_gag19424_igs19424.out_G.png','orbrms_gag19424_igs19424G.png','pod.out','pod.rms','pseudobs_ICRF.out','pseudobs_ITRF.out')

# Syantax for downlfiles:
# downloadfiles(src_directory,dest_directory,list_of_files_to_download)
downloadfiles(ex1SrcDir,ex1DestDir,ex1Files)
downloadfiles(ex1SrcSolDir,ex1DestSolDir,ex1SolFiles)


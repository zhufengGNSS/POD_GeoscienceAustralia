import os
import argparse
import shutil
import urllib.request

import re
#===============================================================================
# Command line argument
parser = argparse.ArgumentParser(description="Download POD tables")

parser.add_argument('-d',dest='DestDir',type=str, required=False, default='tables', help='Destination directory path (default: tables)')
args = parser.parse_args()

#===============================================================================
# 
#===============================================================================
def downloadfiles(srcdir,destdir,srcfiles):

    # Get all of the data subdirectory
    for filename in srcfiles:
        url = srcdir + filename
        destfile = destdir + '/' + filename
    
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
# Download pod/tables files
#===============================================================================
S3bucket = 'https://peanpod.s3-ap-southeast-2.amazonaws.com/pod/tables/'

tables_SrcDir = S3bucket

tables_DestDir = args.DestDir

tables_Files = ('ascp1950.430','fes2004_Cnm-Snm.dat','goco05s.gfc','header.430_229','igs_metadata_2063.snx','leap.second','eopc04_14_IAU2000.62-now')

# Syantax for downlfiles:
# downloadfiles(src_directory,dest_directory,list_of_files_to_download)
downloadfiles(tables_SrcDir,tables_DestDir,tables_Files)


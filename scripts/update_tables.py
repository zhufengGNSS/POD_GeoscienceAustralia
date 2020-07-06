import os
import argparse
from ftplib import FTP
#import urllib3, shutil

#===============================================================================
# Command line argument
parser = argparse.ArgumentParser(description="Download POD examples")

parser.add_argument('-d',dest='outdir',type=str, required=False, default='/data/acs/pod/tables', help='Output directory path')
args = parser.parse_args()

#===============================================================================
filename = 'eopc04_14_IAU2000.62-now'
destfilename = args.outdir + "/" + filename

# Create the directory path if needed
if not os.path.exists(os.path.dirname(destfilename)):
    print("Making the directory"+destfilename)
    os.makedirs(os.path.dirname(destfilename))

# Connect to the IERS server
ftp = FTP('ftp.iers.org')
ftp.login()
ftp.cwd('products/eop/long-term/c04_14/iau2000')

with open(destfilename, 'wb') as fp:
    ftp.retrbinary('RETR '+filename, fp.write)

# The below only works for http addresses    
#url = "ftp://ftp.iers.org/products/eop/long-term/c04_14/iau2000/eopc04_14_IAU2000.62-now"
#c = urllib3.PoolManager()
# HTTP example
#with c.request('GET', url, preload_content=False) as res, open(filename, 'wb') as out_file:
#    shutil.copyfileobj(res, out_file)
    

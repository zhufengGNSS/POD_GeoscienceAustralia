from boto3.session import Session
import boto3
import os
import argparse

# Key is for PEANPODPUBLIC user
# Will only allow download from the peanpod  S3 bucket
ACCESS_KEY = 'AKIAYZV65R4UQHCAKYWT'
SECRET_KEY = 'lyFk1VgveMsR4uGZ08sUUkPjPWm9VNTrfa7BubcV'


# Command line argument
parser = argparse.ArgumentParser(description="Download POD examples")

parser.add_argument('-d',dest='outdir',type=str, required=False, default='/data/acs/', help='Output directory path')

#==============================================================================
def download_dir(client, resource, dist, local='/tmp', bucket='peanpod'):
    paginator = client.get_paginator('list_objects')
    for result in paginator.paginate(Bucket=bucket, Delimiter='/', Prefix=dist):
        if result.get('CommonPrefixes') is not None:
            for subdir in result.get('CommonPrefixes'):
                download_dir(client, resource, subdir.get('Prefix'), local, bucket)

        for file in result.get('Contents', []):
            dest_pathname = os.path.join(local, file.get('Key'))
            if not os.path.exists(os.path.dirname(dest_pathname)):
                os.makedirs(os.path.dirname(dest_pathname))
            
           
            if(file.get('Key')[-1] != '/' and not os.path.isfile(dest_pathname)  ):
                print("Downloading to:",dest_pathname)
                resource.meta.client.download_file(bucket, file.get('Key'), dest_pathname)
#==============================================================================
args = parser.parse_args()

session = Session(aws_access_key_id=ACCESS_KEY,aws_secret_access_key=SECRET_KEY)
s3      = session.resource('s3')
client  = session.client('s3')

print("Outdir:",args.outdir)
download_dir(client, s3, 'pod/examples', args.outdir, bucket='peanpod')

download_dir(client, s3, 'pod/tables', args.outdir, bucket='peanpod')

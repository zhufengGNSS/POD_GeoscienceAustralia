import boto3
import os
import argparse


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

resource = boto3.resource('s3')
client  = boto3.client('s3', region_name='ap-southeast-2')

print("Outdir:",args.outdir)
download_dir(client, resource, 'pod/examples', args.outdir, bucket='peanpod')

download_dir(client, resource, 'pod/tables', args.outdir, bucket='peanpod')

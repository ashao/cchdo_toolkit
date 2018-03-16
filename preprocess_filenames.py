#!/usr/bin/env python
import argparse
from netCDF4 import Dataset
from glob import glob
import os

parser = argparse.ArgumentParser(description="Checks all bottle files in the given directory and prepends the expocode found to the file name if not already included")
parser.add_argument('filepath', type=str, help="Path to where the bottle files are stored")
parser.add_argument('-t', action='store_true', help="Only write the output, but do not actually rename the files")
args = parser.parse_args()
print(args)

filenames = glob(os.path.join(args.filepath+'*.nc'))

for file in filenames:
    data = Dataset(os.path.join(args.filepath,file))
    expocode = getattr(data,'EXPOCODE')
    data.close()

    if not expocode in file:
        # If only testing, write out the old and new names
        name = os.path.basename(file)
        if args.t:
            print(os.path.join(args.filepath,name), os.path.join(args.filepath,expocode+'_'+name))
        else:
            os.rename(os.path.join(args.filepath,name), os.path.join(args.filepath,expocode+'_'+name))


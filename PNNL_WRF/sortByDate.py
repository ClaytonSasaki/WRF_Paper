#!/usr/bin/python3

import os
import shutil

indir = '/home/disk/monsoon/relampago/raw/wrf'

for file in os.listdir(indir):
    if file.startswith('wrfout'):
        print(file)
        [prefix,domain,date,time] = file.split('_')
        newDate = date.replace('-','')
        if not os.path.isdir(indir+'/'+newDate):
            os.mkdir(indir+'/'+newDate)
        shutil.move(indir+'/'+file,
                    indir+'/'+newDate+'/'+file)
        
        

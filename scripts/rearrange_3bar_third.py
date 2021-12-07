import os, errno
import glob
import h5py
import argparse
import xml.etree.ElementTree as ET
import numpy as np
from read_column import *
import bootstrap
parser = argparse.ArgumentParser()
parser.add_argument('--input_dir', help="input directory name", type=str, required=True)
parser.add_argument('--output', '-o', help="output directory name", type=str, required=True)
parser.add_argument('--mom', '-p', help="momentum", type=str, required=True)
parser.add_argument('--NT', '-nt', help="NT", type=int, required=True)

args = parser.parse_args()
input_dir=args.input_dir
output=args.output
mom=args.mom #'0 0 0' , '1 0 0'

MOM=mom.split(' ')
MOM2=[]
for m in MOM:
    MOM2.append(str(-1*int(m)))
mom2=' '.join(MOM2) #'0 0 0','-1 0 0' #opposit momentum

NT=args.NT

confs=[]  
with open('conf.list','r') as c:
    for line in c:
        confs.append(line.strip())

for y in range(NT): #loop for y
    DATA3=[]
    for iconf, conf in enumerate(confs):
        print(iconf+1, conf)
        #Third Term
        filelist_third=sorted(glob.glob(input_dir+'/THIRD/{}/{}_y{}_X*.xml'.format(conf,conf,y)))
        Data3=[]
        for f3 in filelist_third:
            data3=[]
            try:
                tree=ET.parse(f3)
                root=tree.getroot()
                for meson in root.findall("./Wilson_hadron_measurements/elem/Shell_Shell_Wilson_Mesons/elem/[gamma_value='15']/momenta/elem/[sink_mom='{}']/mesprop/elem/re".format(mom2)):
                    data3.append(float(meson.text)) #data_TrSC[y']
                    Data3.append(data3[y]) #pick y=y' Data3[Nsrc]
            except:
                print(f3)
        DATA3.append(np.mean(Data3)) #DATA3[nconf]
    hf=h5py.File(output+'.h5','a')
    group='THIRD/MOM_{}'.format('_'.join(mom.split(' ')))
    if group in hf.keys():
        g=hf[group]
    else:
        g=hf.create_group(group)
    if 'y_{}'.format(y) in g:
        del g['y_{}'.format(y)]
    g.create_dataset('y_{}'.format(y),data=DATA3,compression="gzip", compression_opts=9)
    hf.close()

    hf=h5py.File('boot_'+output+'.h5','a')
    group='THIRD/MOM_{}'.format('_'.join(mom.split(' ')))
    if group in hf.keys():
        g=hf[group]
    else:
        g=hf.create_group(group)
    if 'y_{}'.format(y) in g:
        del g['y_{}'.format(y)]
    g.create_dataset('y_{}'.format(y),data=bootstrap.resampling(DATA3,1000),compression="gzip", compression_opts=9)
    hf.close()

              







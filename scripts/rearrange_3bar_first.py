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

DATA1=[]
for iconf, conf in enumerate(confs):
    print(iconf+1, conf)

    #First Term TrSS*TrSC
    filelist_TrSS=sorted(glob.glob(input_dir+'/FIRST/{}/TrSS/{}_*.xml'.format(conf,conf)))
    print('len(filelist_TrSS) : ',len(filelist_TrSS))
    Data_First=[]
    for f_TrSS in filelist_TrSS:
        f_TrSC=f_TrSS.replace('TrSS','TrSC').replace('shell_sink','c_shell_sink')
#        print(f_TrSS)
#        print(f_TrSC)
        data_TrSS,data_TrSC=[],[]
        try:
            tree=ET.parse(f_TrSS)
            root=tree.getroot()
            for meson in root.findall("./Wilson_hadron_measurements/elem/Shell_Shell_Wilson_Mesons/elem/[gamma_value='15']/momenta/elem/[sink_mom='{}']/mesprop/elem/re".format(mom)):
                data_TrSS.append(float(meson.text)) #data_TrSC[NT]
        except:
            print('f_TrSS : ',f_TrSS)
#        print('len(data_TrSS) : ',len(data_TrSS))
        try:
            tree=ET.parse(f_TrSC)
            root=tree.getroot()
            for meson in root.findall("./Wilson_hadron_measurements/elem/Shell_Shell_Wilson_Mesons/elem/[gamma_value='15']/momenta/elem/[sink_mom='{}']/mesprop/elem/re".format(mom2)):
                data_TrSC.append(float(meson.text)) #data_TrSC[NT]
        except:
            print('f_TrSC : ',f_TrSC)
#        print('len(data_TrSC) : ',len(data_TrSC))
        Data_First.append(np.array(data_TrSS) * np.array(data_TrSC))#Data_First[Nsrc][NT]
    print('len(Data_First) : ',len(Data_First))
    Data_First=np.array(Data_First).transpose() #Data_First[NT][Nsrc]
    Data=[]
    for y in range(NT):
        Data.append(np.mean(Data_First[y])) #Data[NT]
    DATA1.append(np.array(Data)) #DATA1[Nconf][NT]
print('len(DATA1) nconf : ',len(DATA1))
DATA1=np.array(DATA1).transpose() #DATA1[NT][Nconf]
print('len(DATA1) NT : ',len(DATA1))


hf=h5py.File(output+'.h5','a')
group='FIRST/MOM_{}'.format('_'.join(mom.split(' ')))
if group in hf.keys():
    g=hf[group]
else:
    g=hf.create_group(group)

for y in range(NT):
    if 'y_{}'.format(y) in g:
        del g['y_{}'.format(y)]
    g.create_dataset('y_{}'.format(y),data=DATA1[y],compression="gzip", compression_opts=9)
hf.close()

hf=h5py.File('boot_'+output+'.h5','a')
group='FIRST/MOM_{}'.format('_'.join(mom.split(' ')))
if group in hf.keys():
    g=hf[group]
else:
    g=hf.create_group(group)

for y in range(NT):
    if 'y_{}'.format(y) in g:
        del g['y_{}'.format(y)]
    g.create_dataset('y_{}'.format(y),data=bootstrap.resampling(DATA1[y],1000),compression="gzip", compression_opts=9)
hf.close()


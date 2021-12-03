import sys
import os
from sobol import *

src_start = 0
src_stop = 128
ndim = 4

nx = 32
nt = 64

doQQBAR = False
doQQQNUC = True

ms = sys.argv[1]  #-0.013
mc = sys.argv[2]  #0.25
beta = sys.argv[3]  #3.6
cfgNum = int(sys.argv[4])
direc= sys.argv[5]  # directory where things are stored

# here is an example of the full path to some configuration
#cfg = '/p/fastdata/slnpp/GREGORY/CONFIGS/NF3P1/B3.6/B3.6_M-0.013M0.25_L32T64/NS8_LS1_G2/test_b3.6_m-0.013m0.25_l32t64_nf3p1_cfg_810.lime'

cfgDir = '/p/fastdata/slnpp/GREGORY/CONFIGS/NF3P1/B'+str(beta)+'/B'+str(beta)+'_M'+str(ms)+'M'+str(mc)+'_L'+str(nx)+'T'+str(nt)+'/NS8_LS1_G2/'  # base directory of configurations

cfg = cfgDir + 'test_b'+str(beta)+'_m'+str(ms)+'m'+str(mc)+'_l'+str(nx)+'t'+str(nt)+'_nf3p1_cfg_'+str(cfgNum)+'.lime'

src=[]

#gamma = [1,2,4,15]
gamma = [15]

seed = src_start
for _ in range(src_start, src_stop):
    [srn, seed] = i4_sobol(ndim, seed)
    src.append([int(nx*srn[0]+.5)%nx,int(nx*srn[1]+.5)%nx,int(nx*srn[2]+.5)%nx,int(nt*srn[3]+.5)%nt])

src_start = seed

def printPreamble(cfg):
    print(
        '<?xml version="1.0"?>\n'+
        '<chroma>\n'+
        '  <Cfg>\n'+
        '    <cfg_type>SZINQIO</cfg_type>\n'+
        '    <cfg_file>'+cfg+'</cfg_file>\n'+
        '    <parallel_io>true</parallel_io>\n'+
        '  </Cfg>\n'+
        '  <Param>\n'+
        '    <nrow>'+str(nx)+' '+str(nx)+' '+str(nx)+' '+ str(nt)+'</nrow>\n'+
        '    <InlineMeasurements>\n'
    )

def printFinale():
    print(
        '    </InlineMeasurements>\n'+
        '  </Param>\n'+
        '</chroma>\n'
    )
        
def printMakeSource(tsrc):
    print(
        '      <elem>\n'+
        '       <Name>MAKE_SOURCE</Name>\n'+
        '       <Frequency>1</Frequency>\n'+
        '       <Param>\n'+
        '        <version>6</version>\n'+
        '         <Source>\n'+
        '          <version>2</version>\n'+
        '          <SourceType>SHELL_SOURCE</SourceType>\n'+
        '          <j_decay>3</j_decay>\n'+
        '          <t_srce>'+str(tsrc[0])+' '+str(tsrc[1])+' '+str(tsrc[2])+' '+str(tsrc[3])+'</t_srce>\n'+
        '          <SmearingParam>\n'+
        '            <wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind>\n'+
        '            <wvf_param>2.0</wvf_param>\n'+
        '            <wvfIntPar>5</wvfIntPar>\n'+
        '            <no_smear_dir>3</no_smear_dir>\n'+
        '          </SmearingParam>\n'+
        '         </Source>\n'+
        '       </Param>\n'+
        '       <NamedObject>\n'+
        '         <gauge_id>default_gauge_field</gauge_id>\n'+
        '         <source_id>sm_source.'+str(tsrc[0])+'_'+str(tsrc[1])+'_'+str(tsrc[2])+'_'+str(tsrc[3])+'</source_id>\n'+
        '       </NamedObject>\n'+
        '      </elem>\n'
        )

def printPropagator(tsrc,quark,mass):
    print(
        '      <elem>\n'+
        '       <Name>PROPAGATOR</Name>\n'+
        '       <Frequency>1</Frequency>\n'+
        '       <Param>\n'+
        '         <version>10</version>\n'+
        '         <quarkSpinType>FULL</quarkSpinType>\n'+
        '         <obsvP>false</obsvP>\n'+
        '         <FermionAction>\n'+
        '           <FermAct>CLOVER</FermAct>\n'+
        '           <Mass>'+str(mass)+'</Mass>\n'+
        '           <clovCoeff>1.0</clovCoeff>\n'+
        '           <FermState>\n'+
        '             <Name>STOUT_FERM_STATE</Name>\n'+
        '             <rho>0.11</rho>\n'+
        '             <orthog_dir>-1</orthog_dir>\n'+
        '             <n_smear>6</n_smear>\n'+
        '             <FermionBC>\n'+
        '               <FermBC>SIMPLE_FERMBC</FermBC>\n'+
        '               <boundary>1 1 1 -1</boundary>\n'+
        '             </FermionBC>\n'+
        '           </FermState>\n'+
        '         </FermionAction>\n'+
        '         <InvertParam>\n'+
        '         <invType>QPHIX_CLOVER_INVERTER</invType>\n'+
        '           <SolverType>BICGSTAB</SolverType>\n'+
        '           <MaxIter>10000</MaxIter>\n'+
        '           <RsdTarget>1e-9</RsdTarget>\n'+
        '            <CloverParams>\n'+
        '              <Mass>'+str(mass)+'</Mass>\n'+
        '              <clovCoeff>1.0</clovCoeff>\n'+
        '            </CloverParams>\n'+
        '            <AntiPeriodicT>true</AntiPeriodicT>\n'+
        '            <RsdToleranceFactor>50</RsdToleranceFactor>\n'+
        '          </InvertParam>\n'+
        '       </Param>\n'+
        '       <NamedObject>\n'+
        '         <gauge_id>default_gauge_field</gauge_id>\n'+
        '         <source_id>sm_source.'+str(tsrc[0])+'_'+str(tsrc[1])+'_'+str(tsrc[2])+'_'+str(tsrc[3])+'</source_id>\n'+
        '         <prop_id>'+quark+'_sm_src_prop.'+str(tsrc[0])+'_'+str(tsrc[1])+'_'+str(tsrc[2])+'_'+str(tsrc[3])+'</prop_id>\n'+
        '       </NamedObject>\n'+
        '      </elem>\n'+
        '\n'+
        '      <elem>\n'+
        '        <Name>SINK_SMEAR</Name>\n'+
        '        <Frequency>1</Frequency>\n'+
        '        <Param>\n'+
        '          <version>5</version>\n'+
        '          <Sink>\n'+
        '            <version>2</version>\n'+
        '            <SinkType>POINT_SINK</SinkType>\n'+
        '            <j_decay>3</j_decay>\n'+
        '          </Sink>\n'+
        '        </Param>\n'+
        '        <NamedObject>\n'+
        '          <gauge_id>default_gauge_field</gauge_id>\n'+
        '          <prop_id>'+quark+'_sm_src_prop.'+str(tsrc[0])+'_'+str(tsrc[1])+'_'+str(tsrc[2])+'_'+str(tsrc[3])+'</prop_id>\n'+
        '          <smeared_prop_id>'+quark+'_sm_src_pt_sink_prop.'+str(tsrc[0])+'_'+str(tsrc[1])+'_'+str(tsrc[2])+'_'+str(tsrc[3])+'</smeared_prop_id>\n'+
        '        </NamedObject>\n'+
        '      </elem>\n'
    )
    return quark+'_sm_src_pt_sink_prop.'+str(tsrc[0])+'_'+str(tsrc[1])+'_'+str(tsrc[2])+'_'+str(tsrc[3])

def printSU3(sprop,cprop,directory,gamma1,gamma2):
    print(
        '      <elem>\n'+
        '	<annotation> NOW DO SOME SU3 MEASUREMENTS </annotation>\n'+
        '	<Name>SU3_SYSTEM</Name>\n'+
        '	<Frequency>1</Frequency>\n'+
        '	<Param>\n'+
        '         <file_dir>'+directory+'/</file_dir>\n'+
        '         <ConsistencyTests>false</ConsistencyTests>\n'+
        '         <Rep15>true</Rep15>\n'+
	'         <Rep6>true</Rep6>\n'+
	'         <Rep3bar>false</Rep3bar>\n'+
        '         <Gamma1>'+str(gamma1)+'</Gamma1>\n'+
        '         <Gamma2>'+str(gamma2)+'</Gamma2>\n'+
        '	</Param>\n'+
        '	<NamedObject>\n'+
        '          <gauge_id>default_gauge_field</gauge_id>\n'+
        '          <light_prop_id>'+sprop+'</light_prop_id>\n'+
        '          <heavy_prop_id>'+cprop+'</heavy_prop_id>\n'+
        '	</NamedObject>\n'+
        '	<xml_file></xml_file>\n'+
        '      </elem>\n'
        )

def printEraseObject(prop):
    print(
        '      <elem>\n'+
        '        <Name>ERASE_NAMED_OBJECT</Name>\n'+
        '        <Frequency>1</Frequency>\n'+
        '        <NamedObject>\n'+
        '          <object_id>'+prop+'</object_id>\n'+
        '        </NamedObject>\n'+
        '      </elem>\n'
    )

def printApeSmear(tsrc,quark):
    print(
        '      <elem>\n'+
        '	<Name>SINK_SMEAR</Name>\n'+
        '	<Frequency>1</Frequency>\n'+
        '	<Param>\n'+
        '          <version>5</version>\n'+
        '          <Sink>\n'+
        '            <version>2</version>\n'+
        '            <SinkType>SHELL_SINK</SinkType>\n'+
        '            <j_decay>3</j_decay>\n'+
        '            <Displacement>\n'+
        '              <version>1</version>\n'+
        '              <DisplacementType>NONE</DisplacementType>\n'+
        '            </Displacement>\n'+
        '            <SmearingParam>\n'+
        '              <wvf_kind>GAUGE_INV_GAUSSIAN</wvf_kind>\n'+
        '              <wvf_param>2.0</wvf_param>\n'+
        '              <wvfIntPar>5</wvfIntPar>\n'+
        '              <no_smear_dir>3</no_smear_dir>\n'+
        '            </SmearingParam>\n'+
        '            <LinkSmearing>\n'+
        '              <LinkSmearingType>APE_SMEAR</LinkSmearingType>\n'+
        '              <link_smear_fact>2.5</link_smear_fact>\n'+
        '              <link_smear_num>1</link_smear_num>\n'+
        '              <no_smear_dir>3</no_smear_dir>\n'+
        '            </LinkSmearing>\n'+
        '          </Sink>\n'+
        '	</Param>\n'+
        '	<NamedObject>\n'+
        '          <gauge_id>default_gauge_field</gauge_id>\n'+
        '          <prop_id>'+quark+'_sm_src_prop.'+str(tsrc[0])+'_'+str(tsrc[1])+'_'+str(tsrc[2])+'_'+str(tsrc[3])+'</prop_id>\n'+
        '          <smeared_prop_id>'+quark+'_sm_src_sm_sink_prop.'+str(tsrc[0])+'_'+str(tsrc[1])+'_'+str(tsrc[2])+'_'+str(tsrc[3])+'</smeared_prop_id>\n'+
        '	</NamedObject>\n'+
        '      </elem>\n'
    )
    return quark+'_sm_src_sm_sink_prop.'+str(tsrc[0])+'_'+str(tsrc[1])+'_'+str(tsrc[2])+'_'+str(tsrc[3])

def printQQBAR(prop1,prop2,directory):
    print(
        '      <elem>\n'+
        '       <Name>QQBAR</Name>\n'+
        '       <Frequency>1</Frequency>\n'+
        '       <Param>\n'+
        '        <version>4</version>\n'+
        '        <Dirac_basis>true</Dirac_basis>\n'+
        '        <nrow>'+str(nx)+' '+str(nx)+' '+str(nx)+' '+ str(nt)+'</nrow>\n'+
        '       </Param>\n'+
        '       <NamedObject>\n'+
        '        <gauge_id>default_gauge_field</gauge_id>\n'+
        '        <prop_ids>\n'+
        '          <elem>'+prop1+'</elem>\n'+
        '          <elem>'+prop2+'</elem>\n'+
        '        </prop_ids>\n'+
        '        <qqbar_file>'+directory+'/qqbar_'+prop1+'_'+prop2+'</qqbar_file>\n'+
        '       </NamedObject>\n'+
        '      </elem>\n'
        )

def printQQQNUC(prop1,prop2,directory):
    print(
        '      <elem>\n'+
        '       <Name>QQQ_NUCNUC</Name>\n'+
        '       <Frequency>1</Frequency>\n'+
        '       <Param>\n'+
        '        <version>2</version>\n'+
        '        <max_p2>9</max_p2>\n'+
        '       </Param>\n'+
        '       <NamedObject>\n'+
        '        <gauge_id>default_gauge_field</gauge_id>\n'+
        '        <prop_ids>\n'+
        '          <elem>'+prop1+'</elem>\n'+
        '          <elem>'+prop2+'</elem>\n'+
        '        </prop_ids>\n'+
        '       </NamedObject>\n'+
        '       <qqbar_file>'+directory+'/qqbar_'+prop1+'_'+prop2+'</qqbar_file>\n'+
        '      </elem>\n'
        )

    
def main():
    
    printPreamble(cfg)
    for n in range(len(src)): # loop through the number of sources
        printMakeSource(src[n])
        cprop = printPropagator(src[n],'charm',mc) # charm quark propagator and point smear
        sprop = printPropagator(src[n],'strange',ms) # strange quark propagator and point smear
        if doQQBAR:
            printQQBAR(sprop,sprop,direc)
            printQQBAR(sprop,cprop,direc)
            printQQBAR(cprop,cprop,direc)
        if doQQQNUC:
            printQQQNUC(sprop,cprop,direc)
        for i in range(len(gamma)):
            for	j in range(len(gamma)):
                printSU3(sprop,cprop,direc,gamma[i],gamma[j]) # this does the SU3 measurements
        printEraseObject(cprop) # now erase point propagators since they are no longer needed
        printEraseObject(sprop)
        cprop = printApeSmear(src[n],'charm')
        sprop = printApeSmear(src[n],'strange')
        if doQQBAR:
            printQQBAR(sprop,sprop,direc)
            printQQBAR(sprop,cprop,direc)
            printQQBAR(cprop,cprop,direc)
        if doQQQNUC:
            printQQQNUC(sprop,cprop,direc)
        for i in range(len(gamma)):
            for	j in range(len(gamma)):
                printSU3(sprop,cprop,direc,gamma[i],gamma[j]) # this does the SU3 measurements
        printEraseObject(cprop) # now erase point propagators since they are no longer needed
        printEraseObject(sprop)
        printEraseObject('charm_sm_src_prop.'+str(src[n][0])+'_'+str(src[n][1])+'_'+str(src[n][2])+'_'+str(src[n][3])) # and erase the source propagators
        printEraseObject('strange_sm_src_prop.'+str(src[n][0])+'_'+str(src[n][1])+'_'+str(src[n][2])+'_'+str(src[n][3]))
    printFinale()

if __name__ == "__main__":
    main()

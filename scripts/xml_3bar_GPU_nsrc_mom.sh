#!/bin/bash
. parameter.sh
BASE=`pwd -P`
templete_dir=${BASE}/xml.templete
conf_dir=${BASE}/conf

DATA=${BASE}/scratch
DATA_3bar=${DATA}/3bar

XML_DIR=${BASE}/xml.input/3bar

mkdir -p ${XML_DIR}
mkdir -p ${DATA_3bar}

iconf=$1
tag=`head -n $iconf ./conf.list | tail -1`
conf=${conf_dir}/${tag}
XML=${XML_DIR}/${tag}_$1_$2.xml

NSRC=$3
SRC_STR=$2
SRC_END=$((SRC_STR+NSRC-1))

gamma_index=15

#MAKE XML INPUT
echo "<?xml version=\"1.0\"?>" > ${XML}
echo "<chroma>" >> ${XML}
echo "   <Param>" >> ${XML}
echo "      <InlineMeasurements>" >> ${XML}
for nsrc in $(seq $SRC_STR $SRC_END)
do
  SRC=`python3 ${BASE}/sobol.py ${NS} ${NS} ${NS} ${NT} ${nsrc} | tail -n 1`
  X=`echo $SRC | awk '{print $1}'`
  Y=`echo $SRC | awk '{print $2}'`
  Z=`echo $SRC | awk '{print $3}'`
  T=`echo $SRC | awk '{print $4}'`

#point source at (X,Y,Z,T)
sed -e s%X1%${X}% -e s%X2%${Y}% -e s%X3%${Z}% -e s%X4%${T}% -e s%SMEAR_PARAM%${SMEAR_PARAM}% -e s%N_SMEAR%${N_SMEAR}% -e s%SOURCE_ID%sh_source% ${templete_dir}/shell_source.xml >> ${XML}

#S(y,x)
sed -e s%MASS%${mass}% -e s%CLOVCOEFF%${clov}% -e s%STABILIZED%${STABILIZED}% -e s%ACTION%${ACTION}% -e s%SOURCE_ID%sh_source% -e s%PROP_ID%sh_prop% ${templete_dir}/invert_mass_gpu.xml >> ${XML}
  
#C(y,x)
sed -e s%MASS%${mass_c}% -e s%CLOVCOEFF%${clov}% -e s%STABILIZED%${STABILIZED}% -e s%ACTION%${ACTION}% -e s%SOURCE_ID%sh_source% -e s%PROP_ID%sh_prop_c% ${templete_dir}/invert_mass_gpu.xml >> ${XML}

sed -e s%OBJECT_ID%sh_source% ${templete_dir}/erase.xml >> ${XML}


sed -e s%SINK_TYPE%SHELL_SINK% -e s%SMEAR_PARAM%${SMEAR_PARAM}% -e s%N_SMEAR%${N_SMEAR}% -e s%GAUGE_ID%default_gauge_field% -e s%PROP_ID%sh_prop% -e s%SMEAR_PROP%smear_sh_prop% ${templete_dir}/sink_smear_APE_SMEAR.xml >> ${XML}

sed -e s%SINK_TYPE%SHELL_SINK% -e s%SMEAR_PARAM%${SMEAR_PARAM}% -e s%N_SMEAR%${N_SMEAR}% -e s%GAUGE_ID%default_gauge_field% -e s%PROP_ID%sh_prop_c% -e s%SMEAR_PROP%smear_sh_prop_c% ${templete_dir}/sink_smear_APE_SMEAR.xml >> ${XML}

#compute the first term Tr[Gamma S(x,y) Gamma S(y,x)], Tr[Gamma S(x,y) Gamma C(y,x)]
mkdir -p ${DATA_3bar}/FIRST/${tag}/TrSS
sed -e s%GAUGE_ID%default_gauge_field% -e s%PROP1_ID%smear_sh_prop% -e s%PROP2_ID%smear_sh_prop% -e s%XML_FILE%${DATA_3bar}/FIRST/${tag}/TrSS/${tag}_shell_sink_X${X}_Y${Y}_Z${Z}_T${T}.xml% -e s%MOM2_MAX%5% -e s%AVG_EQUIV_MOM%false% ${templete_dir}/hadron_spectrum.xml >> ${XML}
mkdir -p ${DATA_3bar}/FIRST/${tag}/TrSC
sed -e s%GAUGE_ID%default_gauge_field% -e s%PROP1_ID%smear_sh_prop_c% -e s%PROP2_ID%smear_sh_prop% -e s%XML_FILE%${DATA_3bar}/FIRST/${tag}/TrSC/${tag}_c_shell_sink_X${X}_Y${Y}_Z${Z}_T${T}.xml% -e s%MOM2_MAX%5% -e s%AVG_EQUIV_MOM%false% ${templete_dir}/hadron_spectrum.xml >> ${XML}

#compute the second term 1/3Tr[Gamma S(x,y) Gamma S(x,y') Gamma S(y',x) Gamma C(y,x)]
mkdir -p ${DATA_3bar}/SECOND/${tag}
sed -e s%FILE%${DATA_3bar}/SECOND/${tag}/% -e s%REP15%false% -e s%REP6%false% -e s%REP3BAR%true% -e s%GAMMA1%15% -e s%GAMMA2%15% -e s%GAUGE_ID%default_gauge_field% -e s%LIGHT_PROP%smear_sh_prop% -e s%HEAVY_PROP%smear_sh_prop_c% ${templete_dir}/SU3_system.xml >> ${XML}

#S(x,x)
sed -e s%PROP_IN%smear_sh_prop% -e s%PROP_OUT%Sxx% -e s%X1%${X}% -e s%X2%${Y}% -e s%X3%${Z}% -e s%X4%${T}% ${templete_dir}/PeekSpinColor.xml >> ${XML}
#Gamma
sed -e s%GAUGE_ID%default_gauge_field% -e s%PROP_ID%sh_prop% -e s%GAMMA_ID%GAMMA% -e s%T_OP%1000% -e s%INDEX%${gamma_index}% ${templete_dir}/gamma_Op.xml >> ${XML}

#Gamma S(x,x)
sed -e s%FACTOR%1.0% -e s%PROP_A%GAMMA% -e s%PROP_B%Sxx% -e s%PROP_C%GAMMA_Sxx% ${templete_dir}/qmul.xml >> ${XML}

MOM[0]="0 0 0"
MOM[1]="1 0 0"
MOM[2]="1 1 0"
MOM[3]="1 1 1"
MOM[4]="2 0 0"
MOM[5]="2 1 0"

MOM_tag[0]=000
MOM_tag[1]=100
MOM_tag[2]=110
MOM_tag[3]=111
MOM_tag[4]=200
MOM_tag[5]=210

MOM_LIST="1" #"0 1 2 3 4 5"

for imom in ${MOM_LIST}
do

  for y in $(seq 0 63)
  do
    #sequential source Gamma C(y,x)  
    sed -e s%T_OPER%${y}% -e s%SINK_TYPE%SHELL_SINK% -e s%SMEAR_PARAM%${SMEAR_PARAM}% -e s%N_SMEAR%${N_SMEAR}% -e s%GAMMA_INDEX%${gamma_index}% -e s%SEQ_QUARKS%false% -e s%PROP_ID%sh_prop_c% -e s%SEQSOURCE_ID%gamma_prop_c_source% -e s%MOM%"${MOM[${imom}]}"% ${templete_dir}/gamma_seqsource.xml >> ${XML}
    # S(y',y) Gamma C(y,x)
    sed -e s%MASS%${mass}% -e s%CLOVCOEFF%${clov}% -e s%STABILIZED%${STABILIZED}% -e s%ACTION%${ACTION}% -e s%SOURCE_ID%gamma_prop_c_source% -e s%PROP_ID%prop_gamma_prop_c% ${templete_dir}/invert_mass_gpu.xml >> ${XML}

    sed -e s%OBJECT_ID%gamma_prop_c_source% ${templete_dir}/erase.xml >> ${XML}

    sed -e s%FACTOR_A%0.0% -e s%FACTOR_B%1.0% -e s%PROP_A%smear_sh_prop% -e s%PROP_B%prop_gamma_prop_c% -e s%PROP_C%smear_prop_gamma_prop_c% ${templete_dir}/qadd.xml >> ${XML}

    sed -e s%OBJECT_ID%prop_gamma_prop_c% ${templete_dir}/erase.xml >> ${XML}

    #S(y',y) Gamma C(y,x) * Gamma S(x,x)
    sed -e s%FACTOR%1.0% -e s%PROP_A%smear_prop_gamma_prop_c% -e s%PROP_B%GAMMA_Sxx% -e s%PROP_C%smear_prop_gamma_prop_c_gamma_Sxx% ${templete_dir}/qmul.xml >> ${XML}

    sed -e s%OBJECT_ID%smear_prop_gamma_prop_c% ${templete_dir}/erase.xml >> ${XML}
    #contraction S(y',y) Gamma C(y,x) * Gamma S(x,x) and S(y',x)
    mkdir -p ${DATA_3bar}/THIRD/mom_${MOM_tag[${imom}]}/${tag}
    sed -e s%GAUGE_ID%default_gauge_field% -e s%PROP1_ID%smear_prop_gamma_prop_c_gamma_Sxx% -e s%PROP2_ID%smear_sh_prop% -e s%XML_FILE%${DATA_3bar}/THIRD/mom_${MOM_tag[${imom}]}/${tag}/${tag}_y${y}_X${X}_Y${Y}_Z${Z}_T${T}.xml% -e s%MOM2_MAX%5% -e s%AVG_EQUIV_MOM%false% ${templete_dir}/hadron_spectrum.xml >> ${XML}
    sed -e s%OBJECT_ID%smear_prop_gamma_prop_c_gamma_Sxx% ${templete_dir}/erase.xml >> ${XML}
  done
done

sed -e s%OBJECT_ID%sh_prop% ${templete_dir}/erase.xml >> ${XML}
sed -e s%OBJECT_ID%sh_prop_c% ${templete_dir}/erase.xml >> ${XML}
sed -e s%OBJECT_ID%smear_sh_prop% ${templete_dir}/erase.xml >> ${XML}
sed -e s%OBJECT_ID%smear_sh_prop_c% ${templete_dir}/erase.xml >> ${XML}
sed -e s%OBJECT_ID%Sxx% ${templete_dir}/erase.xml >> ${XML}
sed -e s%OBJECT_ID%GAMMA% ${templete_dir}/erase.xml >> ${XML}
sed -e s%OBJECT_ID%GAMMA_Sxx% ${templete_dir}/erase.xml >> ${XML}



done #nsrc

echo "     </InlineMeasurements>" >> ${XML}
echo "     <nrow>${NS} ${NS} ${NS} ${NT}</nrow>" >> ${XML}
echo "   </Param>" >> ${XML}

sed -e s%CFG_TYPE%${CFG_TYPE}% -e s%CFG_FILE%${conf}% ${templete_dir}/configuration.xml >> ${XML}

echo "</chroma>" >> ${XML}

echo ${XML}

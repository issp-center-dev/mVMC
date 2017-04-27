#!/bin/sh
#QSUB -queue F18acc
#QSUB -node 1
#QSUB -mpi 24
#QSUB -omp 1
#QSUB -place pack
#QSUB -over false
#PBS -l walltime=23:30:00  
#PBS -N hybrid
cd ${PBS_O_WORKDIR}

MPI="mpijob"
VMC="./vmc.out -F 5"
DEF1="xnamelist.def  zinitial.def"
DEF2="xnamelist_pro.def  zqp1_opt.dat"
DEFaft="xnamelist_aft.def  zqp_opt.dat"
DEFene="xnamelist_ene.def  zqp_opt.dat"
DEFSC1="xnamelist_aft_SC1.def  zqp_opt.dat"
DEFtest="xnamelist_test.def  zinitial.def"

perl -w R_CisAjsCktAltDC.pl 
perl -w New_R_SC1.pl 


cd ./V_Def
 cp ../input.txt .
 sh def_VMC.sh
 cp *def ..
cd ..

date
echo "opt"
${MPI} ${VMC} ${DEF1}
date

date
echo "pro"
${MPI} ${VMC} ${DEF2}
date

date
echo "aft"
${MPI} ${VMC} ${DEFaft}
sh Aft.sh
date


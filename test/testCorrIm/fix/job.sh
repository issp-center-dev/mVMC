#!/bin/sh

#$ -S /bin/sh
#$ -cwd
#$ -V 
#$ -q gr3.q
#$ -pe openmpi   16
#$ -w n

export OMP_NUM_THREADS=1

#L=12
#Sx=2
#Sy=2
#N=58
#tp=-0.3
#U=7.0
#U=7.0
#J=0.2
AF=0.0
SC=1

#./makeDefFile.py Lx Ly Sx Sy Ne U
#./makeDef_Hubbard2dAPSCAft.py input_makedef.txt 
#./makeDef_Hubbard2dAPSCAft.py  ${L} ${L} ${Sx} ${Sy} ${N} ${U} 
#./makeDef_Hubbard2dAPSCAft.py  ${L} ${L} ${Sx} ${Sy} ${N} ${U} ${J}
#            Lx Ly Sx Sy Ne APFlag CDWFlag SymSc GapSc GapAF
#./initHFB.py 12 12 12  2 62   1       0      2d    0     0
#./normal3.py  Lx   Ly    Nup   tp  AF APFlag
#./normal3.py  ${L} ${L}  ${N} ${tp} 0.2   1
#./scafForAP_ver3.0.py  Lx   Ly  Nup    tp   AF    SC
#./scafForAP_ver3.0.py  input_makedef.txt    ${AF} ${SC}
#./normal3.py  Lx   Ly  Nup      AF 
#./normal3.py input_makedef.txt ${AF}

#rm moto -rf
#mkdir moto
#cp * moto
#cd moto
#./makeDef_StripeSCAft.py zinitpara.def 

#date
#mpirun -np $NSLOTS ./vmc.out xnamelist.def zqp_init.dat
#date

#date
#mpirun -np $NSLOTS ./vmc.out xnamelist_aft.def zqp_opt.dat
#mpirun -np $NSLOTS ./vmc.out xnamelist_aft.def zqp_init.dat
#date
#./average.py xnamelist_aft.def 0 0 1 0

#cp zqp_opt.dat ../zqp_init.dat
#cd ../

date
mpirun -np $NSLOTS ./vmc.out namelist.def
date

#date
#mpirun -np $NSLOTS ./vmc.out xnamelist_aft.def zqp_opt.dat
#date
#./average.py xnamelist_aft.def 0 0 1 0

#date

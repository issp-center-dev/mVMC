#[s] definitions of executions
MPI=" "
VMC="path2vmc.out"       
VMCDRY="path2vmcdry.out" 
UHF="path2UHF" 
#[e] definitions of executions

mkdir random
mv opt random
mv aft random
mv Ene.dat MaxSq.dat  Real.dat  SqNq.dat occ.dat  MaxNq.dat Nk.dat  Sij.dat  random

#[s] opt
  ${VMCDRY} ./stan_opt.in
  ${UHF}  namelist.def
  echo " InOrbital zqp_APOrbital_opt.dat " >> namelist.def
  ${MPI} ${VMC} namelist.def 
  cp ./output/zqp_opt.dat . 
  mv output opt
#[e] opt

#[s] aft
  ${VMCDRY} ./stan_aft.in
  cp green1 greenone.def 
  cp green2 greentwo.def
  ${MPI} ${VMC} namelist.def ./zqp_opt.dat
  mv output aft
#[e] aft

#[s] post process
  python3 VMClocal.py  input.toml
  python3 VMCcor.py    input.toml
#[e] post process

mkdir uhf
mv opt uhf
mv aft uhf
mv Ene.dat MaxSq.dat  Real.dat  SqNq.dat occ.dat  MaxNq.dat Nk.dat  Sij.dat  uhf

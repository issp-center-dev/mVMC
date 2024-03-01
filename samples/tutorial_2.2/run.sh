#[s] definitions of executions
MPI=" "
VMC="path2vmc.out"       
VMCDRY="path2vmcdry.out" 
#[e] definitions of executions

#python3 MakeInput.py input.toml
#[s] opt
  ${VMCDRY} ./stan_opt.in
  ${MPI} ${VMC} namelist.def 
  cp ./output/zqp_opt.dat . 
  mv output opt
#[e] opt

#[s] aft
  python3 SCGreen.py input.toml
  ${VMCDRY}     ./stan_aft.in
  cp green1     greenone.def 
  cp SC_1swave  greentwo.def
  ${MPI} ${VMC} namelist.def ./zqp_opt.dat
  mv output aft
#[e] aft

#[s] post process
  python3 VMClocal.py  input.toml
  python3 CalcSC.py    input.toml
#[e] post process

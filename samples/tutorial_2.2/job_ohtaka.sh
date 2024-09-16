#!/bin/bash

#SBATCH -J TMTTF
#SBATCH -p i8cpu
#SBATCH --time=00:30:00
#SBATCH -N 4
#SBATCH -n 64
#SBATCH -c 8
#SBATCH -o log.%j
#SBATCH -e log.%j

source /home/issp/materiapps/intel/mvmc/mvmcvars.sh
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}

#[s] definitions of executions
MPI="srun "
VMC="vmc.out"       #Pre-installed
VMCDRY="vmcdry.out" #Pre-installed
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

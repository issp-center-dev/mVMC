#!/bin/bash
#PJM -L "node=4"
#PJM --mpi "proc=8"
#PJM -L "elapse=10:00"
#PJM -L "rscgrp=debug"

#------- Program execution -------# 
export OMP_NUM_THREADS=1
fipp -C -I call,hwm -d test mpiexec ./vmc.out xnamelist.def zqp_init.dat

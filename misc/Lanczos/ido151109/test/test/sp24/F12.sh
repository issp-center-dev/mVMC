#!/bin/bash
#------ pjsub option --------#
#PJM -L "rscgrp=F12"
#PJM -L "node=2x5"
#PJM --mpi "proc=40"
#PJM -L "elapse=1440:00"
#PJM -j
#------- Program execution -------# 
export OMP_NUM_THREADS=4
fipp -C -I call,hwm -d test mpiexec ./vmc.out xnamelist.def zqp_opt.dat

#!/bin/bash
#------ pjsub option --------#
#PJM -L "rscgrp=F96"
#PJM -L "node=4x2x10"
#PJM --mpi "proc=320"
#PJM -L "elapse=1440:00"
#PJM -j
#------- Program execution -------# 
export OMP_NUM_THREADS=4
fipp -C -I call,hwm -d test mpiexec ./vmc.out xnamelist.def zqp_opt.dat

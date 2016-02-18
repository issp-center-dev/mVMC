#!/bin/bash -x
#
#PJM --bulk
#PJM --sparam "0-5"
#PJM --rsc-list "node=8x8x8"
#PJM --mpi "rank-map-bynode=XYZ"
#PJM --rsc-list "elapse=8:00:00"
#PJM --rsc-list "rscgrp=large"
#PJM --mpi assign-online-node
#PJM --stg-transfiles all
#PJM -s
#PJM --mpi "use-rankdir"
#
#PJM --stgin "rank=* ./vmc.out %r:./"
#PJM --stgin "rank=* ./%b/zinitial.def %r:./"
#PJM --stgin "rank=0 ./%b/*.def %r:./"
#PJM --stgout "rank=0  %r:./*.dat ./%b/"


. /work/system/Env_base
#. /work/system/Env_base_1.2.0-16-3

export PARALLEL=8
export OMP_NUM_THREADS=$PARALLEL

MPI="mpiexec"
VMC="./vmc.out -F 5"
DEF1="xnamelist.def  zinitial.def"
DEF2="xnamelist_pro.def  zqp_1_opt.dat"
DEFaft="xnamelist_aft.def  zqp_opt.dat"

date
echo "opt"
${MPI} ${VMC} ${DEF1}
date

date
echo "opt_pro"
${MPI} ${VMC} ${DEF2}
date

date
echo "aft"
${MPI} ${VMC} ${DEFaft}
date


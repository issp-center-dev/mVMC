#[s] definitions of executions
MPI=" "
VMC="path2vmc.out"       
#[e] definitions of executions

target=$1

mkdir ${target}
${MPI} ${VMC} -s stan_${target}_S0.in 
mv output $target/output_s0

${MPI} ${VMC} -s stan_${target}_S1.in 
mv output ${target}/output_s1

E0=$(awk ' {print +$1}' ${target}/output_s0/zqp_opt.dat)
E1=$(awk ' {print +$1}' ${target}/output_s1/zqp_opt.dat)
Delta_s=`echo "${E1} - ${E0}" | bc -l`

echo "#E0_from_zqp_opt.dat   E1_from_zqp_opt.dat   Delta_s=E1-E0" > spin_gap.dat
echo $E0 $E1 $Delta_s >> spin_gap.dat
mv spin_gap.dat ${target}

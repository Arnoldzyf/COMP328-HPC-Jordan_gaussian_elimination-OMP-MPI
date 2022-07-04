#!/bin/bash -l
#SBATCH -p course -t 5
#SBATCH -D ./
#SBATCH --export=ALL
#SBATCH --exclusive
#SBATCH -c 40

## sbatch --exclusive -t 3 run_omp.sh gaussian_OMP.c gaussianLib.c

## !! need to change the value of #define N_cons and #define N_var in gaussianLib.h according to the size !!
size_row=240
size_col=600
opt=2

loop=5 ## the number of loops to check time

## cd COMP328/Assignment/compare
module load compilers/intel 
## module load compilers/intel/2019u5
module load mpi/intel-mpi/2019u5/bin

echo
echo "test case 1 of size ${size_row} * ${size_col}"
echo "using the maximum of ${SLURM_CPUS_PER_TASK} cores"
echo
## sbatch --exclusive -t 3 run_omp.sh gaussian_OMP.c gaussianLib.c
SRC_OMP=$1
EXE_OMP=${SRC_OMP%%.c}.exe
SRC_lib=$2


export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1} # if '-c' not used then default to 1 # making #threads=#cores

## echo
## echo
## echo "******************** OMP ***********************"
## rm -f ${EXE_OMP}
## icc -qopenmp -O0 $SRC_OMP $SRC_lib -o $EXE_OMP
## ## echo '============= print test case ================='
## ./${EXE_OMP} -2 ## debug=-2, print nothing (too lengthy to print the result three times...change to 2 if you want to see the result)
## for i in $(seq 1 $loop)
## do
##   ./${EXE_OMP} -1
## done

echo
echo
echo "******************** OMP ***********************"
rm -f ${EXE_OMP}
icc -qopenmp -O${opt} $SRC_OMP $SRC_lib -o $EXE_OMP
export OMP_NUM_THREADS=2
## echo '============= print test case ================='
./${EXE_OMP} -2 ## debug=-2, print nothing (...change to 2 if you want to see the result)
for k in $(seq 1 ${SLURM_CPUS_PER_TASK}); 
do
  export OMP_NUM_THREADS=${k}
  echo "============= OMP with ${OMP_NUM_THREADS} threads ================="
  for i in $(seq 1 $loop)
  do
     ./${EXE_OMP} -1
  done
done

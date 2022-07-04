#!/bin/bash -l
#SBATCH -p course -t 5
#SBATCH -D ./
#SBATCH --export=ALL
#SBATCH --exclusive
#SBATCH -n 1

## sbatch --exclusive -t 3 opt_MPI.sh gaussian_MPI.c gaussianLib.c

## !! need to change the value of #define N_cons and #define N_var in gaussianLib.h according to the size !!
size_row=240
size_col=240

loop=5 ## the number of loops to check time

## cd COMP328/Assignment/compare
module load compilers/intel 
## module load compilers/intel/2019u5
module load mpi/intel-mpi/2019u5/bin

echo
echo "test case of size ${size_row} * ${size_col}"
echo "using ${SLURM_NTASKS} cores"
echo
## sbatch --exclusive -t 3 opt_MPI.sh gaussian_MPI.c gaussianLib.c
SRC_MPI=$1
EXE_MPI=${SRC_MPI%%.c}.exe
SRC_lib=$2

echo
echo
echo "******************** MPI ***********************"
rm -f ${EXE_MPI}
mpiicc -O0 $SRC_MPI $SRC_lib -o $EXE_MPI 
## echo '============= print test case ================='
export numMPI=${SLURM_NTASKS}
mpirun -np ${numMPI} ./${EXE_MPI} -2 ## debug=-2, print nothing (too lengthy to print the result three times...change to 2 if you want to see the result)
for opt in $(seq 0 3); 
do
  echo "============= opt=${opt} MPI with ${SLURM_NTASKS} processors ================="
  if ((${size_row} % ${SLURM_NTASKS} != 0))
  then 
    echo "!! N_cons ${size_row} is indivisible by comm_size ${SLURM_NTASKS} , please change the number of processors !!"
  else
    rm -f ${EXE_MPI}
    mpiicc -O${opt} $SRC_MPI $SRC_lib -o $EXE_MPI 
    for i in $(seq 1 $loop)
    do
       mpirun -np ${numMPI} ./${EXE_MPI} -1
    done
  fi
done
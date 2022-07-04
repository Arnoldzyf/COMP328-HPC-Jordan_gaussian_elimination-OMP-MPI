#!/bin/bash -l
#SBATCH -p course -t 5
#SBATCH -D ./
#SBATCH --export=ALL
#SBATCH --exclusive
#SBATCH -n 40

## sbatch --exclusive -t 3 run_MPI.sh gaussian_MPI.c gaussianLib.c

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
echo "using the maximum of ${SLURM_NTASKS} cores"
echo
## sbatch --exclusive -t 3 run_MPI.sh gaussian_MPI.c gaussianLib.c
SRC_MPI=$1
EXE_MPI=${SRC_MPI%%.c}.exe
SRC_lib=$2

echo
echo
echo "******************** MPI ***********************"
rm -f ${EXE_MPI}
mpiicc -O${opt} $SRC_MPI $SRC_lib -o $EXE_MPI
export numMPI=2 
## echo '============= print test case ================='
mpirun -np ${numMPI} ./${EXE_MPI} -2 ## debug=-2, print nothing (...change to 2 if you want to see the result)
for k in $(seq 1 ${SLURM_NTASKS})
do
  echo "============= MPI with ${k} processors ================="
  if ((${size_row} % ${k} != 0))
  then 
    echo "!! N_cons ${size_row} is indivisible by comm_size ${k} , please change the number of processors !!"
  else
    export numMPI=${k}
    for i in $(seq 1 $loop)
    do
       mpirun -np ${numMPI} ./${EXE_MPI} -1
    done
  fi
done
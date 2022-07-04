#!/bin/bash -l
#SBATCH -p course -t 5
#SBATCH -D ./
#SBATCH --export=ALL
#SBATCH --exclusive

## sbatch --exclusive -t 3 opt_serial.sh gaussian_serial.c gaussianLib.c

## !! need to change the value of #define N_cons and #define N_var in gaussianLib.h according to the size !!
size_row=240
size_col=240

loop=5 ## the number of loops to check time

## cd COMP328/Assignment/compare
module load compilers/intel 
## module load compilers/intel/2019u5
module load mpi/intel-mpi/2019u5/bin

echo
echo "test optimization"
echo "test case of size ${size_row} * ${size_col}"
echo
## sbatch --exclusive -t 3 opt_serial.sh gaussian_serial.c gaussianLib.c
SRC_serial=$1
EXE_serial=${SRC_serial%%.c}.exe
SRC_lib=$2

## echo "******************** Serial *********************"
rm -f ${EXE_serial}
icc -O0 $SRC_serial $SRC_lib -o $EXE_serial
## echo '============= print test case ================='
./${EXE_serial} -2 ## if set to 2(debug value) will print complete result
## echo "Each test case is printed to a single file"
echo
echo "******************** Serial *********************"
for opt in $(seq 0 3)
do
  rm -f ${EXE_serial}
  icc -O${opt} $SRC_serial $SRC_lib -o $EXE_serial 
  echo "============= opt=${opt} run time =============="
  for i in $(seq 1 $loop)
  do ./${EXE_serial} -1 ## debug=-1, only print time
  done
done

#! /bin/bash
#PBS -l nodes=1:ppn=24
#PBS -N rtv3.00u1.19
#PBS -q default
#PBS -j oe
export MKL_NUM_THREADS=24
export OMP_NUM_THREADS=24
cd $PBS_O_WORKDIR
cat $PBS_NODEFILE > pbsnodes
echo $PWD >>/home/rabh/qsub-$(pwd| head -c 4| tail -c 3)info.txt
time ./ov0.50u.out

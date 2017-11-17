#!/bin/sh

PBS_O_WORKDIR=/home/sarun002/Proj3

cd $PBS_O_WORKDIR

module purge
module load gcc-4.6.2
module load mvapich2-1.9a2/gnu-4.6.2

mpicc seivingprime.c -o seiving -lm
mpicc spart2.c -o spart2 -lm
mpicc spart3.c -o spart3 -lm

first=($qsub jobseive)
second=($qsub -W depend=afterany:($first) jobseive2) 
third=($qsub -W depend=afterany:($second) jobseive4)
fourth=($qsub -W depend=afterany:($third) jobseive8)
fifth=($qsub -W depend=afterany:($fourth) jobspart)
sixth=($qsub -W depend=afterany:($fifth) jobspart2)
seventh=($qsub -W depend=afterany:($sixth) jobspart4)
eighth=($qsub -W depend=afterany:($seventh) jobspart8)
ninth=($qsub -W depend=afterany:($eighth) job3part)
tenth=($qsub -W depend=afterany:($ninth) job3part2)
eleth=($qsub -W depend=afterany:($tenth) job3part4)
$qsub -W depend=afterany:($eleth) job3part8




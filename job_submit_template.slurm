#!/bin/bash
#SBATCH --job-name=adm_test      # create a short name for your job
#SBATCH --nodes=1               # node count
# SBATCH --exclusive              # fully dedicate single node 
#SBATCH -p physics
#SBATCH --ntasks-per-node=25
#SBATCH --cpus-per-task=1        # cpu-cores per task (>1 if multi-threaded tasks)
#SBATCH --mem=5G               # memory for entire job
#SBATCH --time=18:00:00          # total run time limit (HH:MM:SS)
#SBATCH --mail-type=all          # send email on job start, end, and fail                                                            
#SBATCH --mail-user=sandiproy@princeton.edu                                                                                                            

module purge
module load openmpi/gcc/1.10.2/64 gsl/2.4 fftw/gcc/openmpi-1.10.2/3.3.4 hdf5/gcc/1.8.16

MYPATH=/home/sr39/gizmo-adm
#srun $MYPATH/GIZMO /tigress/sr39/gas_disk_test_Feb2020/gas_spin_analhalo_w1/gas_spin_w1.params 1
srun $MYPATH/GIZMO /tigress/sr39/gizmo_adm_testing/all_adm_v2/all_adm_v2.params

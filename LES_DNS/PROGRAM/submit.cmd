# job name (default is name of pbs script file)
#PBS -N ns3d
# resource limits: number of CPUs to be used
#PBS -l ncpus=8
# resource limits: amount of memory to be used
#PBS -l mem=4000mb, vmem=4000mb
# resource limits: max. wall clock time during which job can be running
#PBS -l walltime=4:00:00

cd $PBS_O_WORKDIR

ulimit -s 10000000
export OMP_STACKSIZE=100M

export OMP_NUM_THREADS=8
#export OMP_SCHEDULE="STATIC,8"

export OMP_PROC_BIND=true


./ns3d > ns3d_test.out
 

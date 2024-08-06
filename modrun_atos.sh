#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128
#SBATCH --qos=np
#SBATCH --account=noemep
#SBATCH --time=12:00:0
#SBATCH --job-name=emepctm
#SBATCH --output=%x.%j.out --error=%x.%j.out

YEAR=(2018)
year=${YEAR[0]}

sim="rowcle50_bs15ch4"

DIRLOCA=${year}_$sim
DIRPATH="/ec/res4/hpcperm/fawc/Output_EMEP/acp_methane/specials/spinup/"

# git directory where to run emep model from
RUNDIR="/home/fawc/EMEP_codes/ACP_Methane/emep-mscw/"

mkdir -p $DIRPATH$DIRLOCA

# store logs and run model
cd $DIRPATH$DIRLOCA

sed "s:RUNYEAR:$year:g;s:SPNYEAR:$spin:g;s:SIMTYPE:$sim:g" $RUNDIR/config_emep_spinup.nml > config_emep.nml

mpirun $RUNDIR/emepctm

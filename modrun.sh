#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=128
#SBATCH --account=noemep
#SBATCH --qos=np
#SBATCH --time=24:00:0
#SBATCH --job-name=emepctm
#SBATCH --output=%x.%j.out --error=%x.%j.out


### Minimalistic script for run the Unified EMEP model

# user=sm_gunla
# user_id=GL
# run=rv4_51

DIRLOCA="EMEP_out/PBAP_test"
DIRPATH="/ec/res4/scratch/nor4796/"

# git directory where to run emep model from
RUNDIR="/home/nor4796/emep-mscw"

#if [ ! -d "$DIRPATH$DIRLOCA" ]; then
mkdir -p $DIRPATH$DIRLOCA
# fi

# store git information in a small git log
pwd                              >>  $DIRPATH$DIRLOCA/GitLog.txt
date                             >>  $DIRPATH$DIRLOCA/GitLog.txt
git describe                     >>  $DIRPATH$DIRLOCA/GitLog.txt
echo "git branch: "              >>  $DIRPATH$DIRLOCA/GitLog.txt
git status --branch --short -uno >>  $DIRPATH$DIRLOCA/GitLog.txt
echo "git hash: "                >>  $DIRPATH$DIRLOCA/GitLog.txt
git rev-parse HEAD               >>  $DIRPATH$DIRLOCA/GitLog.txt
echo " "                         >>  $DIRPATH$DIRLOCA/GitLog.txt
echo " ======= DIFFS =========:" >>  $DIRPATH$DIRLOCA/GitLog.txt
git diff ':!config_emep.nml' ':!modrun.sh' >>  $DIRPATH$DIRLOCA/GitLog.txt

# store logs and run model
cd $DIRPATH$DIRLOCA

cp $RUNDIR/emepctm .
cp $RUNDIR/config_emep.nml config_emep.nml

# run the model
mpirun emepctm


#!/bin/bash

### job options for Slurm/sbatch
#SBATCH -A emep
#SBATCH --job-name=emepctm
#SBATCH --output=%x.%j.out --error=%x.%j.out
#SBATCH --nodes=8 --ntasks-per-node=32 --time=24:00:00

### Minimalistic script for run the Unified EMEP model

# user=sm_gunla
# user_id=GL
# run=rv4_51

DIRLOCA="EMEP_output/PBAP_month_new_dep/"
DIRPATH="/home/sm_gunla/"

# git directory where to run emep model from
RUNDIR="/home/sm_gunla/emep-mscw"

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
mpprun emepctm


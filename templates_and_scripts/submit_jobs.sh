#!/bin/bash

for runscript in `ls run*.sh`
do
  echo $subfolder
  #script=$subfolder + "/run_gs2.sh"
  echo $runscript
  sbatch $runscript
  sleep 0.1
done

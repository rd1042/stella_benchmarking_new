#!/bin/bash

for subfolder in `ls beta_*/ -d`
do
  echo $subfolder
  cd $subfolder
  #script=$subfolder + "/run_gs2.sh"
  sbatch run_stella.sh
  cd ../
  sleep 0.1
done

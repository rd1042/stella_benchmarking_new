#!/bin/bash

for subfolder in `ls fprim*/ -d`
do
  echo $subfolder
  #script=$subfolder + "/run_gs2.sh"
  cd $subfolder
  sbatch run_stella.sh
  cd ..
  sleep 0.1
done

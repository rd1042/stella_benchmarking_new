#!/bin/bash

for subfolder in `ls */ -d`
do
  echo $subfolder
  cd $subfolder
  runscript="run_stella.sh"
  #echo $runscript
  sbatch $runscript 
  cd ../
  sleep 0.1
done

#!/bin/bash

for runscript in `ls run*.sh`
do
  echo $runscript
  sbatch $runscript
  sleep 0.1
done

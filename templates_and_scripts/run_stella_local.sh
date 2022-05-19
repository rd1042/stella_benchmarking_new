#!/bin/bash

for inputfile in `ls *.in`
do
  echo $inputfile
  mpirun -n 4 ~/stella/stella $inputfile
  sleep 0.1
done

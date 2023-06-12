#!/bin/bash

program_name="./h2o_omp"

argument_pairs=(
  "-s 0"
  "-s 1"
  "-s 2"
  "-s 3"
  "-s 4"
  "-d 0"
  "-d 1"
  "-d 2"
  "-d 3"
  "-d 4"
)

export PYTHONPATH=$PYTHONPATH:$HOME/.xacc
export LD_LIBRARY_PATH=$HOME/.xacc/lib:$LD_LIBRARY_PATH
for arguments in "${argument_pairs[@]}"; do
  log_file="log/log_${arguments// /_}.txt"
  
  # Execute the program with the argument pair and redirect output to log file using nohup
  nohup $program_name $arguments > $log_file 2>&1 &
  
  echo "Started: $program_name $arguments and saved output to $log_file"
done

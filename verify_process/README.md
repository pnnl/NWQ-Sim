This version can produce the correct results with only C1 Gate. 

## Run the code to verify the correctness
#### Compile the project
```
cd NWQ-sim
source environment/setup_frontier.sh
mkdir build
cd build
cmake ..
make
```

#### Run small example
Ask for 2 Node 
```
salloc -N 2 -A <Account> -t <time> 
srun -N 2  --ntasks-per-node=1 --gpus-per-task=1 --gpu-bind=none ./qasm/nwq_qasm -backend AMDGPU_MPI -sim sv -q ./verify_process/test.qasm > ./verify_process/sv_hip_c1.txt
diff sv_hip_c1.txt sv_cpu_c1.txt
```


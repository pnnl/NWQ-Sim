This version can produce the correct results with only C1 and C2 Gate. 

## Run the code to verify the correctness

#### Compile the project
```bash
cd NWQ-sim
source environment/setup_frontier.sh
mkdir build
cd build
cmake ..
make
```

#### Run small example
Ask for 2 Node 
```bash
salloc -N 2 -A <Account> -t <time> 
srun -N 2  --ntasks-per-node=1 --gpus-per-task=1 --gpu-bind=none ./qasm/nwq_qasm -backend AMDGPU_MPI -sim sv -q ./verify_process/test.qasm > ./verify_process/sv_hip_c2.txt
diff sv_hip_c2.txt sv_cpu_c2.txt
```

## Troubleshooting
This code is not reliable
Any changes in the code will result in different answer
For example, in line 1504 - 1516 adding any `grid.sync()` will result in incorrect answer even if no threads will enter in the branch. 

We need to investigate how the AMD implement the `grid.sync()` and `roc_shmem_ctx_wg_barrier(ctx)` to make the code run correctly.
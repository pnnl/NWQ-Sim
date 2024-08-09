## Before starting
To ensure the reproducibility, change the 63rd line `rng.seed(time(0));` in `dm_cpu.hpp` as `rng.seed(12345);`. 

## Generate the superoperator for basis gates
1. Change `DEVICE_CONFIG_PATH` in `default_config.json` as `../../data/device/`
2. Generate superoperators
```bash
module load gcc/<versions>
g++ -std=c++17 -I../include test_interface.cpp -o test
srun -N 1 <confirguration paras> ./test
```

This will generate a file which will store the superoperators for basis gates define in the device file. 
## Use the customized superoperator for simulation
1. Comment `add_subdirectory(vqe)` in `CMakeList.txt`
```bash
# add_subdirectory(vqe)
```
2. Build NWQ-sim:
```bash
cd ~/NWQ-Sim
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```
3. Change `DEVICE_CONFIG_PATH` in `default_config.json` as `../data/device/`
4. Run simulator without using the customized superoperators
- Change `CUSTOMIZED_SP` as `false`
- Run
```bash
srun -N 1 ./qasm/nwq_qasm -backend CPU -shots 2048 -sim dm -q ../data/openqasm_basis/adder_n10.qasm
```
5. Run simulator without using the customized superoperators
- Change `CUSTOMIZED_SP` as `true`
```
srun -N 1 ./qasm/nwq_qasm -backend CPU -shots 2048 -sim dm -q ../data/openqasm_basis/adder_n10.qasm
```

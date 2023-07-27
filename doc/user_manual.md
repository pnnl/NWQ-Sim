# User Manual

NWQ-Sim in general only requires a C++ compiler. However, in order to build for GPUs or scaling (up and out) or using other APIs (python, qir, qiskit), we need the following libraries:

## Dependencies

NWQ-Sim has certain dependencies based on the features you want to utilize. The core requirements are:

* C++ Compiler (GCC 7.0 or later recommended)
* CMake 3.20 or later (Project Build Tool)

Additional dependencies depend on the specific features you require:

* OpenMP (for local single-node execution)
* MPI (for multi-node execution)
* CUDA 11.0 or later (for NVIDIA GPU Backend)
* NVSHMEM 2.9.0 or later (for NVIDIA GPU Cluster Scale-out)
* ROCM 3.1.0 (for AMD GPU Backend)
* XACC (for Frontend VQE Algorithm)
* Python 3.4 (for Python API)
* Pybind11 2.5.0 (for Python API)
* mpi4py 3.0.3 (for Python API on a cluster)
* Qiskit 0.20.0 (for Qiskit interface)
* Q# runtime (for Q#/QIR interface)

## <a id="build_base"></a>Build from Source

To build NWQ-Sim from source, follow these steps:

1. Clone the NWQ-Sim repository:
```bash
git clone https://github.com/pnnl/NWQ-Sim.git
cd NWQ-Sim
```
2. Create a build directory and navigate into it:
```bash
mkdir build
cd build
```
3. Configure the project using CMake and build
```bash
cmake ..
make
```

## Build on NERSC Perlmutter Supercomputer
NVSHMEM is required to enable multi-gpu execution on Perlmutter HPC. The current NVSHMEM library has an issue where it prevents multi-GPU execution from utilizing more than 2GB of GPU memory on each GPU. We have provided a fix for this issue in the source code.

Follow these steps to build NWQ-Sim on the NERSC Perlmutter Supercomputer:

1. Initialize the environment with provided script

```bash
source ~/NWQ-Sim/environment/setup_perlmutter.sh
```

2. Build NVSHMEM
* Download and extract the NVSHMEM txz archive from [here](https://developer.download.nvidia.com/compute/redist/nvshmem/). For example, to download and extract NVSHMEM 2.9.0:
```bash
wget https://developer.download.nvidia.com/compute/redist/nvshmem/2.9.0/source/nvshmem_src_2.9.0-2.tar.xz

tar -xf nvshmem_src_2.9.0-2.tar.xz
```
* Replace the mem.cpp file in nvshmem_src

```bash
cp ~/NWQ-Sim/nvshmem_util/mem.cpp ~/nvshmem_src_2.9.0-2/src/mem/mem.cpp
```

* Copy the provided NVSHMEM build script to nvshmem_src folder and then build it
```bash
cp ~/NWQ-Sim/nvshmem_util/scripts/build_nvshmem_perlmutter.sh ~/nvshmem_src_2.9.0-2/
cd ~/nvshmem_src_2.9.0-2
./build_nvshmem_perlmutter.sh
```

Finally, build NWQ-Sim using the steps in [Build from Source](#build_base)


## Configure and run on ORNL Frontier Supercomputer

TO BE ADDED.

## Configure and run on ORNL Summit Supercomputer

Follow these steps to build NWQ-Sim on the OLCF Summit Supercomputer:

1. Initialize the environment with provided script

```bash
source ~/NWQ-Sim/environment/setup_summit.sh
```

2. Build NVSHMEM
* Download and extract the NVSHMEM txz archive from [here](https://developer.download.nvidia.com/compute/redist/nvshmem/). For example, to download and extract NVSHMEM 2.9.0:
```bash
wget https://developer.download.nvidia.com/compute/redist/nvshmem/2.9.0/source/nvshmem_src_2.9.0-2.tar.xz

tar -xf nvshmem_src_2.9.0-2.tar.xz
```
* Replace the mem.cpp file in nvshmem_src

```bash
cp ~/NWQ-Sim/nvshmem_util/mem.cpp ~/nvshmem_src_2.9.0-2/src/mem/mem.cpp
```

* Copy the provided NVSHMEM build script to nvshmem_src folder and then build it
```bash
cp ~/NWQ-Sim/nvshmem_util/scripts/build_nvshmem_summit.sh ~/nvshmem_src_2.9.0-2/
cd ~/nvshmem_src_2.9.0-2
./build_nvshmem_summit.sh
```

Finally, build NWQ-Sim using the steps in [Build from Source](#build_base)


## Program Runtime Configuration Options

This guide provides detailed instructions on how to execute the compiled program along with the available command-line arguments to configure the program runtime.

**Location:** Navigate to the `build` directory in your local project workspace.

**Execution:** Run the executable program with the desired command-line arguments to adjust program behaviors as needed. Here is a comprehensive list of the command-line arguments:

- `-q`: Executes a simulation with the given QASM file.

- `-t <index>`: Runs the testing benchmarks for the specific index provided.

- `-a`: Runs all testing benchmarks. 

- `-backend list`: Lists all the available backends. The list of available backends are:
  - CPU
  - OpenMP
  - MPI
  - NVGPU
  - NVGPU_MPI


- `-backend <name>`: Sets the backend for your program to the specified one. The backend name string is case-insensitive.

- `-shots <value>`: Configures the total number of shots.

- `-basis`: Activates the program to run benchmark circuits using only basis gates.

- `-sim <method>`: Sets the simulation method. There are two available options:
  - `sv`: Stochastic Vector simulation.
  - `dm`: Density Matrix simulation. Please note, when running with `dm`, the given circuit can only contain IBM basis gates and 2-qubit gates that are included in the device configuration file specified in the default_configuration.json file.

**Example Usage:** To run the qasm frontend from the `build` directory with a specific backend, a total number of shots, and a simulation method, use the following command: 

```
./qasm/nwq_sim -backend <name> -shots <value> -sim <method> -q <path/to/qasm>
```

Replace `<name>`, `<value>`, `<method>`, and `<path/to/qasm>` with your desired backend name, number of shots, and simulation method respectively.

Please ensure that you replace `/qasm/nwq_sim` with the actual name of your compiled executable file if not using the qasm frontend.

### Running on Perlmutter Supercomputer
To run NWQ-Sim on the Perlmutter Supercomputer, initilize the environment first
```bash
source ~/NWQ-Sim/environment/setup_perlmutter.sh
```

Launch multi-GPU execution for regular or interactive jobs:
```bash
srun -C gpu -N <NODES> -n <GPUS> -c 1 --gpus-per-task=1 --gpu-bind=single:1 <NWQ-Sim Command> -backend NVGPU_MPI
```

Replace `<NODES>` with the number of compute nodes, `<GPUS>` with the total number of GPUs, and `<NWQ-Sim Command>` with the NWQ-Sim execution command.

## XACC Frontend
To use NWQ-Sim as an execution backend for XACC, follow these steps:
1. Install XACC by following the steps outlined in the [XACC repository](https://github.com/eclipse/xacc#build-from-source).

2. Include the NWQ-Sim backend implementation in your code:
```c++
#include "nwq_accelerator.hpp"
```
Create an NWQAccelerator object:
```c++
auto nwq_acc = std::make_shared<xacc::quantum::NWQAccelerator>();
``` 

Utilize the NWQAccelerator with XACC. For example, you can run XACC-VQE:
```c++
 xacc::Initialize(argc, argv);

// Get reference to the Accelerator
auto nwq_acc = std::make_shared<xacc::quantum::NWQAccelerator>();

nwq_acc->updateConfiguration(
  { std::make_pair("shots", 4096),
    std::make_pair("backend", "cpu"),
    std::make_pair("sim-method", "sv"),
  });

// Create the N=2 deuteron Hamiltonian
auto H_N_2 = xacc::quantum::getObservable(
    "pauli", std::string("5.907 - 2.1433 X0X1 "
                          "- 2.1433 Y0Y1"
                          "+ .21829 Z0 - 6.125 Z1"));

auto optimizer = xacc::getOptimizer("mlpack");

// JIT map Quil QASM Ansatz to IR
xacc::qasm(R"(
.compiler xasm
.circuit deuteron_ansatz
.parameters theta
.qbit q
X(q[0]);
Ry(q[1], theta);
CNOT(q[1],q[0]);
)");

auto ansatz = xacc::getCompiled("deuteron_ansatz");

// Get the VQE Algorithm and initialize it
auto vqe = xacc::getAlgorithm("vqe");
vqe->initialize({std::make_pair("ansatz", ansatz),
                  std::make_pair("observable", H_N_2),
                  std::make_pair("accelerator", accelerator),
                  std::make_pair("optimizer", optimizer)});

// Allocate some qubits and execute
auto buffer = xacc::qalloc(2);
vqe->execute(buffer);
xacc::Finalize();
```

Replace the target source file in `NWQ-Sim/xacc/CMakeList.txt` and build the project. The executable will be located at `NWQ-Sim/build/xacc/nwq_xacc`.

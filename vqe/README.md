# NWQ-VQE
VQE simulation based on the NWQ-Sim platform
## Dependencies
- [NLOpt](https://github.com/stevengj/nlopt): To simplify installation, NLOpt is a git submodule of NWQ-Sim. After cloning `NWQ-Sim`, we need to sync and build `vqe/nlopt` prior to building NWQ-Sim. The installation steps are:

```shell
  git submodule init
  git submodule update
  cd vqe/nlopt
  mkdir build;cd build
  cmake ..
  make
```

Some users may need to define C and C++ compilers for CMake
```shell
  cmake .. -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DCMAKE_C_COMPILER=/usr/bin/gcc
```

The specific CXX/C compiler paths can be altered as needed. The main emphasis is to ensure that the optimization code has access to C++ STD libraries.
On OLCF Frontier, you may need to load the cray-python module before using cmake for NLOpt.

On OLCF Summit, you may use gcc/9.3 for the compilation of NLOpt:
```shell
  module load gcc/9.3.0-compiler_only
  cmake .. -DCMAKE_CXX_COMPILER=/sw/summit/gcc/9.3.0-2/bin/g++ -DCMAKE_C_COMPILER=/sw/summit/gcc/9.3.0-2/bin/gcc
```


## Installation/Configuration Directions
After installing NLOpt, build NWQ-Sim as normal (see the [User Manual](doc/user_manual.md)).

Note that debug messages showing the Fermionic operator indices and Hamiltonian Pauli strings will be printed *unless* the project is built in `Release` mode (add `-DCMAKE_BUILD_TYPE=Release`). 

Running `make` from the `build` directory will compile the binary `build/vqe/nwq_vqe`, a shared library `build/vqe/libvqe.dylib`, and the example binaries under `/build/vqe/examples`. 


## Command Line Usage
NWQ-VQE can be used from the command line via the `nwq_vqe` executable. The executable supports a number of arguments to configure runtime parameters. 

The arguments below are only for displaying information.
```shell
INFORMATIONAL
-h, --help            Show help menu.
-l, --list-backends   List available backends and exit.
```

To run VQE or QFlow on a provided Hamiltonian file, only two arguments are required: the path to the XACC/DUCC-formatted Hamiltonian (see [example_hamiltonians](example_hamiltonians) for examples) and the number of electrons:
```shell
REQUIRED
-f, --hamiltonian     Path to the input Hamiltonian file (formatted as a sum of Fermionic operators, see examples).
-p, --nparticles      Number of electrons in molecule.
```
The default indexing scheme is XACC/Qiskit. This can be altered with flags below.
```shell
OPTIONAL
--ducc                Use DUCC indexing scheme. Defaults to false (defaults to use XACC/Qiskit scheme).
```

The default ansatz is UCCSD and comes with 3 levels of symmetries to adjust the number of parameters. All three levels are supposed to provide the same accuracy without the discrepancy from the trotterization and optimization. There is also an optional singlet and triplet ansatz. It has a much larger number of parameters and is mainly designed for ADAPT-VQE (but also compatible with VQE if necessary).
``` shell
OPTIONAL
--sym, --symm         UCCSD Symmetry level (0->none, 1->spin symmetry, 2->also orbital symmetry). Defaults to 0.
--gsd                 Use singlet GSD ansatz for ADAPT-VQE. Default to false.
```

The first class of optional arguments are used for configuring the simulator: selecting the backend, listing initial amplitude, setting random seeds, and passing NWQ-Sim configuration JSON files
```shell
OPTIONAL
-b, --backend,        Simulation backend. Defaults to CPU.
--seed                Random seed for initial point and empirical gradient estimation. Defaults to time(NULL).
--config              Path to NWQ-Sim config file. Defaults to "../default_config.json".
```

NWQ-VQE uses the [NLOpt](https://nlopt.readthedocs.io/en/latest/) library for optimization algorithm implementations. To configure the NLOpt object, we support setting termination criteria, upper/lower bounds, and optimizer algorithms via command line. Algorithm-specific options can be set via a JSON file, see [mma_config.json](mma_config.json) for an example. To get a full list of valid NLOpt optimizer strings, refer to the [NLOpt documentation](https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/). An example of a correctly formatted string is `--optimizer LN_COBYLA`.
```shell
-v, --verbose         Print optimization callback on each iteration. Defaults to false.
--opt-config          Path to config file for NLOpt optimizer parameters
-o, --optimizer       NLOpt optimizer name. Defaults to LN_COBYLA. Other examples are LN_NEWUOA and LD_LBFGS
-lb, --lbound         Optimizer lower bound. Defaults to -2Pi
-ub, --ubound         Optimizer upper bound. Defaults to 2Pi
--reltol              Relative tolerance termination criterion. Defaults to -1 (off)
--abstol              Relative tolerance termination criterion. Defaults to -1 (off)
--maxeval             Maximum number of function evaluations for optimizer (only for VQE). Defaults to 100
--maxtime             Maximum optimizer time (seconds). Defaults to -1.0 (off)
--stopval             Cutoff function value for optimizer. Defaults to -MAXFLOAT (off)
```

Following are the options to set circuit simulation backend related parameters:
```shell
--num_threads         Specify the number of OMP threads. Defaults to use all hardware threads.
--disable_fusion      Disable gate fusion. Defaults to enabled.
```

To run the $\mathrm{H_4}$ example using command line, run:
```shell
./vqe/nwq_vqe -f ../vqe/example_hamiltonians/H4_4_0.9_xacc.hamil -p 4 --maxeval 1000
```
where we have now increased the evaluation limit to 1000.


### ADAPT-VQE Usage
NWQ-Sim also supports ADAPT-VQE simulations, using both Fermionic operators ([Ref](https://www.nature.com/articles/s41467-019-10988-2)) and single Pauli strings ([Ref](https://journals.aps.org/prxquantum/abstract/10.1103/PRXQuantum.2.020310)). The executable supports the following ADAPT-specific options. `--adapt-maxeval`, `--adapt-gradtol`, and `--adapt-fvaltol` set termination criteria for the ADAPT loop, while the `--adapt` and `--qubit` flags enable ADAPT and Qubit-ADAPT respectively. To enable Qubit-ADAPT, you must pass both the `--adapt` and `--qubit` flags. Qubit-ADAPT also supports randomized operator pool subsampling, with the size set by the `--adapt-pool` flag (-1 indicates the full pool). 
```shell
ADAPT-VQE OPTIONS
--adapt               Use ADAPT-VQE for dynamic ansatz construction. Defaults to false
-ag, --adapt-gradtol  Cutoff absolute tolerance for operator gradient norm. Defaults to 1e-3
-af, --adapt-fvaltol  Cutoff absolute tolerance for function value. Defaults to 1e-6
-am, --adapt-maxeval  Set a maximum iteration count for ADAPT-VQE. Defaults to 100
```


To run the $\mathrm{H_4}$ problem using Fermionic ADAPT, we add the ``--adapt'' flag:
```shell
./vqe/nwq_vqe -f ../vqe/example_hamiltonians/H4_4_0.9_xacc.hamil -p 4 --maxeval 200 --adapt
```
where `--maxeval` now indicates the number of iterations for the VQE subsolver, not the number of ADAPT iterations (Default 100).

To increase the number of ADAPT iterations, we can add the flag `--adapt-maxeval` flag:
```shell
./vqe/nwq_vqe -f ../vqe/example_hamiltonians/H4_4_0.9_xacc.hamil -p 4 --maxeval 200 --adapt --adapt-maxeval 500
```
to increase the iteration limit to 500.

ADAPT-VQE also supports termination criteria based on gradient norm and function value convergence. The gradient norm is $\sqrt{\sum_{i=1}\langle \psi|[H,O_i]|\psi\rangle^2}$ for each $O_i$ in the operator pool.

To set a gradient norm tolerance cutoff of $0.001$, we can add the `--adapt-gradtol` flag:
```shell
./vqe/nwq_vqe -f ../vqe/example_hamiltonians/H4_4_0.9_xacc.hamil -p 4 --maxeval 200 --adapt --qubit --adapt-maxeval 400 --adapt-gradtol 1e-3
```
Conversely, we can add a tolerance cutoff for the function value difference between iterations with the `--adapt-fvaltol` flag:
```shell
./vqe/nwq_vqe -f ../vqe/example_hamiltonians/H4_4_0.9_xacc.hamil -p 4 --maxeval 200 --adapt --qubit --adapt-maxeval 400 --adapt-fvaltol 1e-3
```

#### Qubit-ADAPT VQE
Qubit-ADAPT also supports custom configuration options:
```
--qubit               Uses Qubit instead of Fermionic operators for ADAPT-VQE. Defaults to false
--adapt-pool          Sets the pool size for Qubit operators. Defaults to -1
``` 

The `--qubit` flag enables a Pauli string operator pool (based on UCCSD ansatz) rather than using Fermionic operators e.g.:
```shell
./vqe/nwq_vqe -f ../vqe/example_hamiltonians/H4_4_0.9_xacc.hamil -p 4 --maxeval 200 --adapt --qubit --adapt-maxeval 400 --adapt-fvaltol 1e-3 --qubit
```
However, this constructs an operator pool with all Pauli strings from the JW-mapped Fermionic pool, which may be excessive. To randomly subsample from the operator pool, use the `--adapt-pool` to set the pool size:
```shell
./vqe/nwq_vqe -f ../vqe/example_hamiltonians/H4_4_0.9_xacc.hamil -p 4 --maxeval 200 --adapt --qubit --adapt-maxeval 400 --adapt-fvaltol 1e-3 --qubit --adapt-pool 50
```
Note that `--adapt-pool` has no effect on the Fermionic ADAPT solver.




## Running the Example(s) and General Workflow
To add NWQ-VQE functionality:
1. Add the include paths:
  - Makefile: Add the `-I<Path_To_NWQSIM>/vqe/include` flag 
  - CMake: Add `include_directories(BEFORE PRIVATE <path_to_NWQSIM>/vqe/include)` to your CMakeLists.txt file.
2. Link to the shared library:
  - Makefile: Add flags for the library directory `-L<Path_To_NWQSIM>/build/vqe` and library `-lvqe`
  - CMake: Add the following lines to specify the directory and executable dependencies:
  ```cmake
  link_directories(BEFORE PRIVATE <Path_To_NWQSIM>/build/vqe)
  link_libraries(vqe) 
  ```

The file [examples/basic_example.cpp](src/example.cpp) contains basic example code, explained below. The corresponding executable is compiled to `build/vqe/examples/basic_example`. Note that the example code

The NWQ-VQE API can be broken into 3 main steps: problem definition, ansatz construction, and solver configuration.
1. **Problem Definition**: We start by specifying the target Hamiltonian. The program expects a sum of Fermionic operators with complex coefficients. [example_hamiltonians/h2O.hamil](examples/h2O.hamil) provides an example generated by `xacc` for Quantum Chemistry ([source](https://github.com/npbauman/DUCC-Hamiltonians/blob/main/H2O/cc-pv5z/single/6/out-xacc)). 
```c++
NWQSim::IdxType n_particles = 4; // Set the number of particles
  // Note: path relative to presumed build directory
  std::string hamiltonian_path = "../vqe/example_hamiltonians/H4_4_0.9_xacc.hamil"; //  Effective Hamiltonian file path

  std::shared_ptr<Hamiltonian> hamil = std::make_shared<Hamiltonian>(hamiltonian_path, n_particles, true); // Build the Hamiltonian object (used for energy calculation)
```
2. **Ansatz Construction**: For VQE we need a parameterized wavefunction in the form of a quantum ansatz. The NWQ-VQE package supports the UCCSD ansatz with a Jordan-Wigner (JW) Fermion-Qubit mapping. We specify first specify the JW tranformation function, then call the UCCSD ansatz constructor (subclassed from the NWQ-Sim circuit class). Note that NWQ-Sim backends expect circuits to be wrapped by a `std::shared_ptr`. The `Ansatz` class contains the general functionality needed for optimization, while `UCCSD` provides the specific circuit constructor. After declaring the circuit object, invoke the `buildAnsatz` function to initialize the circuit
```c++
  Transformer jw_transform = getJordanWignerTransform; // Choose a transformation function

  // Build the ansatz circuit using the Hamiltonian Molecular environment and JW mapping
  //      (shared_ptr used to match baseline NWQ-Sim functionality)
  std::shared_ptr<Ansatz> ansatz = std::make_shared<UCCSD>(hamil->getEnv(), jw_transform, 1);
  ansatz->buildAnsatz();
```

3. **Solver Configuration**: To solve the problem, we pass the ansatz pointer and the Hamiltonian object to a `SV_CPU_VQE` object to run the main optimization loop. Here we use the COBYLA algorithm implemented by NLOpt, see (NLOpt Algorithms)[https://nlopt.readthedocs.io/en/latest/NLopt_Algorithms/] for a full list of algorithms. If the algorithm requires a gradient, the SPSA scheme is used to compute the emprical gradient.
```c++
  std::string config_path = "../default_config.json";
  // Build the Quantum State object
  NWQSim::VQE::SV_CPU_VQE state(ansatz, // reference to ansatz
                                hamil,  // reference to Hamiltonian
                                nlopt::algorithm::LN_COBYLA, // NLOpt algorithm for optimization
                                callback_function, // Callback function for each energy evaluation
                                config_path, // path to config file
                                0 // Random seed (passed to the SPSA gradient estimator for random perturbations)
                                );
```
The `callback_function` is specified above, and it permits per-iteration logging. Note that this is called on every function evalution, which may occur multiple times per iteration:
```c++
// Callback function, requires signature (void*) (const std::vector<NWQSim::ValType>&, NWQSim::ValType, NWQSim::IdxType)
double last = 0.0;
void callback_function(const std::vector<NWQSim::ValType>& x, NWQSim::ValType fval, NWQSim::IdxType iteration) {
  double relchange = abs(last) > 1e-6 ? 0.0 : (fval - last) / last;
  std::cout << "Iteration " << iteration << ", fval = " << fval << ", relchange " << relchange << std::endl;
  last = fval;
}
```

With those 3 steps done, all we need to do is specify a starting point (all zero here)
```c++
  // Set initial parameters to 0
  std::vector<double> parameters(ansatz->numParams(), 0.0);

  // Return destination for the function value
  double fval;

  // Start the VQE optimization
  state.optimize(parameters, fval);
  std::cout << "Final Parameters: " << parameters << std::endl;
  std::cout << "Final Energy: " << fval << std::endl;
```
With the given parameters, the example output is:
```shell
Final Parameters: [-0.017729, 0.00190356, -0.00687456, -0.00595677, -0.0362287, 0.0676943, 0.00566315, -0.00576158, -0.00107353, 0.00854818, 0.0418309, -0.00767838, 0.0104276, 0.168057, 0.0427336, 0.00160509, -0.00336476, 0.0288875]
Final Energy: -2.16938
```
The final energy (-2.16938) is about 0.11 mHartree away from the true ground state (-2.18031661). 

The example file [examples/config_example.cpp](examples/config_example.cpp) is similar, but demonstrates how to specify `NLOpt` algorithm parameters and cutoff criteria within NWQ-VQE.


### Backends
Note that [examples/basic_example_cuda_mpi.cu](examples/basic_example_cuda_mpi.cu) separates the `SV_CUDA_MPI_VQE` constructor call and the call to `MPI_Finalize()` into two different functions. This is necessary to ensure that the `SV_CUDA_MPI` destructor (which invokes `nvshmem_finalize()`) is called before `MPI_Finalize()`.

To locally run the example CUDA, MPI, and examples respectively, run (from `NWQ-Sim/build`):
```shell
./vqe/examples/vqe/examples/basic_example_cuda
mpirun -n <nproc> ./vqe/examples/vqe/examples/basic_example_mpi
```
To test on NERSC Perlmutter for the CPU MPI version using 4 nodes, 128 cores:
```shell
salloc --nodes 4 --qos interactive -t 60 --constraint cpu --account=m4243
source ../environment/setup_perlmutter.sh
srun -N4 -n128 ./vqe/examples/basic_example_mpi
```

To test on NERSC Perlmutter for the CUDA MPI version using 4 nodes with 1 GPU each:
```shell
salloc --nodes 4 --qos interactive -t 60 --constraint gpu --account=m4243
source ../environment/setup_perlmutter.sh
srun -C gpu -N 4 -n 4 -c 1 --gpus-per-task=1 --gpu-bind=single:1 ./vqe/examples/basic_example_cuda_mpi
```

To test on OLCF Frontier for the CPU MPI version using 4 nodes, 128 cores:
```shell
salloc -N 4 -A CSC528 -t 30 -q debug
source ../environment/setup_frontier.sh
srun -N4 -n128 ./vqe/examples/basic_example_mpi
```
To test on OLCF Summit for the CPU MPI version using 4 nodes, 16 cores per node:
```shell
bsub -W 60 -nnodes 4 -P CHP125 -Is /bin/bash
source ../environment/setup_summit.sh
jsrun -n 64 -a1 -g0 -c1 -r16 ./vqe/examples/basic_example_mpi
```
To run the command-line solver with the aforementioned backends, replace the `./vqe/examples/...` executable with the command line tool and add the appropriate `--backend` flag. For example, to run multi-GPU on Perlmutter:
```shell
salloc --nodes 4 --qos interactive -t 60 --constraint gpu --account=m4243
source ../environment/setup_perlmutter.sh
srun -C gpu -N 4 -n 4 -c 1 --gpus-per-task=1 --gpu-bind=single:1 ./vqe/nwq_vqe -f ../vqe/example_hamiltonians/H4_4_0.9_xacc.hamil -p 4 --maxeval 1000 --backend NVGPU_MPI
```
where we have added the `NVGPU_MPI` backend option to utilize NVSHMEM-linked NVIDIA GPU accelerators.


# QFlow
In the QFlow algorithm, we repeatedly perform single-direction gradient descent over downfolded Hamiltonians. To perform a QFlow gradient descent for an effective Hamiltonian, NWQ-VQE provides a command line interface.

## Command Line
The `build/vqe/nwq_qflow` binary provides a command line interface, the same as `nwq_vqe` (before ADAPT-VQE portion). The specific options differ for local graident follower:
```
NWQ-VQE QFlow Options
--localfd             Use local gradient follower pipeline (Finite Differences).
--localspsa           Use local gradient follower pipeline (SPSA).
-g, --grad-samples    SPSA gradient samples.
--delta               Perturbation magnitude for SPSA and Central Finite Difference, [f(x+delta)-f(x-delta)]/2delta. Defaults to 1e-2."
--eta                 Gradient descent step size. Defaults to 1.
```
`--delta` and `--eta` are QFlow-specific parameters which control the gradient descent procedure (TODO: implement adaptive stepsize line search). `--delta` is used to perturb parameter vectors to compute the empirical gradient using SPSA, whereas `--eta` controls the descent stepsize. The algorithm descends the gradient from the (random) initial point until it finds a minimum, then returns.

# Result Summary in the Output

After successful execution of `nwq_vqe`, the output will contain a summary of the results at the end of the printout. For example, the command
```shell
./vqe/nwq_vqe -f ./H415.hamil -p 4 -v --abstol 1e-8 -lb -3.1 -ub 3.1 --maxeval 1000 -o LN_COBYLA --sym 2
```
gives
```
--------- Result Summary ---------
Method                 : VQE, Symmetry Level = 2
Ansatz                 : UCCSD Minimal
# Ham. Pauli Strings   : 197
Operator Stats         : 26 operators, 15 parameters, and 3428 Gates
Circuit Stats          : 1934 depth, 1988 1q gates, 1440 2q gates, 0.315 gate density
Optimization terminated: Function tolerance reached
# function eval.       : 673
Evaluation Time        : 0 hrs 0 mins 1.3581 secs
Final objective value  : -1.9947639715805670
Final parameters:
  2^ 0 :: -0.0045417680471343
  6^ 4 :: -0.0045417680471343
  3^ 0 :: -0.0000072708605015
  7^ 4 :: -0.0000072708605015
  2^ 1 :: -0.0000086624570434
  6^ 5 :: -0.0000086624570434
  3^ 1 :: -0.0072011151790247
  7^ 5 :: -0.0072011151790247
  3^ 2^ 1 0 :: 0.0559262141883256
  7^ 6^ 5 4 :: 0.0559262141883256
  6^ 2^ 4 0 :: -0.0953014244963014
  7^ 2^ 4 0 :: -0.0000130689432080
  6^ 3^ 4 0 :: -0.0000130689432080
  6^ 2^ 5 0 :: 0.0000125812063835
  6^ 2^ 4 1 :: 0.0000125812063835
  7^ 2^ 5 0 :: 0.1753959843323752
  6^ 3^ 4 1 :: 0.1753959843323752
  7^ 3^ 4 0 :: -0.1781937728268054
  6^ 3^ 5 0 :: 0.1262203587657208
  7^ 2^ 4 1 :: 0.1262203587657208
  7^ 3^ 5 0 :: 0.0000102348263700
  7^ 3^ 4 1 :: 0.0000102348263700
  6^ 2^ 5 1 :: -0.3725823551270946
  7^ 2^ 5 1 :: 0.0000003964811323
  6^ 3^ 5 1 :: 0.0000003964811323
  7^ 3^ 5 1 :: -0.0810860032080360
```
The default orbital ordering is XACC/Qiskit. For H4 example, the spin orbital ordering is like
```
 3 ----- 7     <= Spatial 3
 2 ----- 6     <= Spatial 2
 1 ----- 5     <= Spatial 1
 0 ----- 4     <= Spatial 0
 α       β
```
That is, `2^ 0` represents $t^{2}_0$ or a single excitation from spin-orbital 0 (spatial 0, alpha) to spin-orbital 2 (spatial 2, alpha). And similarly, `6^ 3^ 4 1` represents $t^{6 3}_{4 1}$ or a double excitation from spin-orbitals 4 (spatial 0, beta) and 1 (spatial 1, alpha) to 6 (spatial 2, beta) and 3 (spatial 3, alpha).



<!-- 
## Python
A pybind11-enabled API allows for a direct Python interface. Two functions are provided: `optimize_effective_hamiltonian` and `get_param_count`. The latter returns the number of unique parameters, accounting for symmetry reparameterizations.

A basic example using `qiskit-nature` utilities is shown below (assuming in root `NWQ-Sim` directory):
```python
%load_ext autoreload
%autoreload 2
import sys
sys.path.insert(0, './build/vqe')
import nwqflow
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer
from qiskit_nature.second_q.mappers import JordanWignerMapper
import numpy as np

# Example Hamiltonian from PySCF
driver = PySCFDriver(
    atom='Li 0.0, 0.0, 0.0; H 0.0, 0.0 1.2',
    basis='sto-3g',
    charge=0,
    spin=0,
    unit=DistanceUnit.ANGSTROM
)
problem = driver.run()
problem = ActiveSpaceTransformer(problem.num_particles, 5, range(5)).transform(problem)

# Convert the PySCF format to XACC (e.g. "+_3 +_2 -_0 -_1" -> 3^ 2^ 0 1)
oplist = []
for op, coeff in problem.second_q_ops()[0].items():
    new_op = []
    for o in op.split():
        ostr = o[2:]
        if o[0:2] == '+_':
            ostr+="^"
        new_op.append(ostr)
    oplist.append((' '.join(new_op), coeff))

n_params = nwqflow.get_param_count(num_particles = 4, num_spatial_orbitals=5)  # number of optimization parameters
initial_point = np.random.rand(n_params)  # initial parameter point
# Perform single-direction gradient descent
result = nwqflow.optimize_effective_hamiltonian(operators=oplist,
                                                num_particles=4,
                                                x0=initial_point,
                                                backend="CPU",
                                                xacc=True)
print(result)
```


Documentation for `optimize_effective_hamiltonian` is shown below:
```python
    optimize_effective_hamiltonian(
      operators: list[tuple[str, complex]], 
      num_particles: int, 
      x0: list[float], 
      backend: str = 'CPU',
      xacc: bool = False, 
      seed: int = -1, 
      delta: float = 0.001, 
      eta: float = 0.001,
      num_trotter: int = 1,
      tnum_samples: int = 1) -> list[tuple[str, float]]

    Perform single-direction gradient descentreturn the locally-optimal parameters
            Arguments:
                    operators (Iterable[(str, complex)]): List of xacc-formatted operator strings with coefficients
                    num_particles (int): Number of electrons (assumed to be equal number of alpha/beta)
                    x0 (Iterable[float]): Initial parameter values
                    backend (str): NWQ-Sim backend for simulation. Defaults to "CPU"
                    xacc (bool): Use XACC operator indexing, otherwise use DUCC. Defaults to False
                    seed (int): Random seed for optimizer and SPSA perturbation. Defaults to time(NULL)
                    delta (float): Magnitude of SPSA perturbation. Defaults to 1e-3
                    eta (float): Gradient descent stepsize. Defaults to 1e-3
                    num_trotter (int): Number of Trotter steps (linearly increases number of parameters). Defaults to 1
                    num_samples (int): Number of gradient samples for SPSA average. Defaults to 1

``` -->
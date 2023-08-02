User Manual
===========

NWQ-Sim in general only requires a C++ compiler. However, in order to
build for GPUs or scaling (up and out) or using other APIs (python, qir,
qiskit), we need the following libraries:

Dependencies
------------

NWQ-Sim has certain dependencies based on the features you want to
utilize. The core requirements are:

-  C++ Compiler (GCC 7.0 or later recommended)
-  CMake 3.20 or later (Project Build Tool)

Additional dependencies depend on the specific features you require:

-  OpenMP (for local single-node execution)
-  MPI (for multi-node execution)
-  CUDA 11.0 or later (for NVIDIA GPU Backend)
-  NVSHMEM 2.9.0 or later (for NVIDIA GPU Cluster Scale-out)
-  ROCM 3.1.0 (for AMD GPU Backend)
-  XACC (for Frontend VQE Algorithm)
-  Python 3.4 (for Python API)
-  Pybind11 2.5.0 (for Python API)
-  mpi4py 3.0.3 (for Python API on a cluster)
-  Qiskit 0.20.0 (for Qiskit interface)
-  Q# runtime (for Q#/QIR interface)

Build from Source
-----------------

NWQ-Sim uses CMake for building, which automatically detects the
execution environment, determines which backends to build, and includes
the appropriate libraries accordingly. To build NWQ-Sim from source,
follow these steps:

1. Clone the NWQ-Sim repository:

.. code:: bash

   git clone https://github.com/pnnl/NWQ-Sim.git
   cd NWQ-Sim

2. Create a build directory and navigate into it:

.. code:: bash

   mkdir build
   cd build

3. Configure the project using CMake and build

.. code:: bash

   cmake ..
   make

Build on HPC Systems
--------------------

For multi-GPU execution on HPCs with NVIDIA GPUs, NWQ-Sim requires the
NVSHMEM library. However, the current version of NVSHMEM library has a
known issue that restricts each GPU to utilize no more than 2GB of GPU
memory. To overcome this limitation, we have incorporated a fix in our
source code.

ORNL Frontier HPC
~~~~~~~~~~~~~~~~~

Currently, for AMD CPU, single/OpenMP/MPI work, for AMD GPU, only single
AMD MI250X GPU works.

Follow these steps to build NWQ-Sim on the OLCF Frontier HPC:

1. Initialize the environment with provided script

.. code:: bash

   source ~/NWQ-Sim/environment/setup_frontier.sh

Then, build NWQ-Sim using the steps in `Build from
Source <#build_base>`__

ORNL Summit HPC
~~~~~~~~~~~~~~~

Follow these steps to build NWQ-Sim on the OLCF Summit HPC:

1. Initialize the environment with provided script

.. code:: bash

   source ~/NWQ-Sim/environment/setup_summit.sh

2. Build NVSHMEM

-  Download and extract the NVSHMEM txz archive from
   `here <https://developer.download.nvidia.com/compute/redist/nvshmem/>`__.
   For example, to download and extract NVSHMEM 2.9.0:

.. code:: bash

   wget https://developer.download.nvidia.com/compute/redist/nvshmem/2.9.0/source/nvshmem_src_2.9.0-2.tar.xz

   tar -xf nvshmem_src_2.9.0-2.tar.xz

-  Replace the mem.cpp file in nvshmem_src

.. code:: bash

   cp ~/NWQ-Sim/nvshmem_util/mem.cpp ~/nvshmem_src_2.9.0-2/src/mem/mem.cpp

-  Copy the provided NVSHMEM build script to nvshmem_src folder and then
   build it

.. code:: bash

   cp ~/NWQ-Sim/nvshmem_util/scripts/build_nvshmem_summit.sh ~/nvshmem_src_2.9.0-2/
   cd ~/nvshmem_src_2.9.0-2
   ./build_nvshmem_summit.sh

Finally, build NWQ-Sim using the steps in `Build from
Source <#build_base>`__

NERSC Perlmutter HPC
~~~~~~~~~~~~~~~~~~~~

Follow these steps to build NWQ-Sim on the NERSC Perlmutter HPC:

1. Initialize the environment with provided script

.. code:: bash

   source ~/NWQ-Sim/environment/setup_perlmutter.sh

2. Build NVSHMEM

-  Download and extract the NVSHMEM txz archive from
   `here <https://developer.download.nvidia.com/compute/redist/nvshmem/>`__.
   For example, to download and extract NVSHMEM 2.9.0:

.. code:: bash

   wget https://developer.download.nvidia.com/compute/redist/nvshmem/2.9.0/source/nvshmem_src_2.9.0-2.tar.xz

   tar -xf nvshmem_src_2.9.0-2.tar.xz

-  Replace the mem.cpp file in nvshmem_src

.. code:: bash

   cp ~/NWQ-Sim/nvshmem_util/mem.cpp ~/nvshmem_src_2.9.0-2/src/mem/mem.cpp

-  Copy the provided NVSHMEM build script to nvshmem_src folder and then
   build it

.. code:: bash

   cp ~/NWQ-Sim/nvshmem_util/scripts/build_nvshmem_perlmutter.sh ~/nvshmem_src_2.9.0-2/
   cd ~/nvshmem_src_2.9.0-2
   ./build_nvshmem_perlmutter.sh

Finally, build NWQ-Sim using the steps in `Build from
Source <#build_base>`__

Program Runtime Configuration Options
-------------------------------------

This guide provides detailed instructions on how to execute the compiled
program along with the available command-line arguments to configure the
program runtime.

**Location:** Navigate to the ``build`` directory in your local project
workspace.

**Execution:** Run the executable program with the desired command-line
arguments to adjust program behaviors as needed. Here is a comprehensive
list of the command-line arguments:

-  ``-q``: Executes a simulation with the given QASM file.

-  ``-t <index>``: Runs the testing benchmarks for the specific index
   provided.

-  ``-a``: Runs all testing benchmarks.

-  ``-backend list``: Lists all the available backends. The list of
   available backends are:

   -  CPU
   -  OpenMP
   -  MPI
   -  NVGPU
   -  NVGPU_MPI
   -  AMDGPU

-  ``-backend <name>``: Sets the backend for your program to the
   specified one. The backend name string is case-insensitive.

-  ``-shots <value>``: Configures the total number of shots.

-  ``-basis``: Activates the program to run benchmark circuits using
   only basis gates.

-  ``-sim <method>``: Sets the simulation method. There are two
   available options:

   -  ``sv``: Stochastic Vector simulation.
   -  ``dm``: Density Matrix simulation. Please note, when running with
      ``dm``, the given circuit can only contain IBM basis gates and
      2-qubit gates that are included in the device configuration file
      specified in the default_configuration.json file.

**Example Usage:** To run the qasm frontend from the ``build`` directory
with a specific backend, a total number of shots, and a simulation
method, use the following command:

::

   ./qasm/nwq_sim -backend <name> -shots <value> -sim <method> -q <path/to/qasm>

Replace ``<name>``, ``<value>``, ``<method>``, and ``<path/to/qasm>``
with your desired backend name, number of shots, and simulation method
respectively.

Please ensure that you replace ``/qasm/nwq_sim`` with the actual name of
your compiled executable file if not using the qasm frontend.

Running on Frontier HPC
~~~~~~~~~~~~~~~~~~~~~~~

To run NWQ-Sim on the Frontier or Crusher Supercomputer, initilize the
environment first

.. code:: bash

   source ~/NWQ-Sim/environment/setup_frontier.sh

Launch multi-CPU execution for regular or interactive jobs:

.. code:: bash

   srun -N<nodes> -n<CPUS> ./qasm/nwq_sim <NWQ-Sim Command> -backend MPI

Running on Summit HPC
~~~~~~~~~~~~~~~~~~~~~

To run NWQ-Sim on the Summit Supercomputer, initilize the environment
first

.. code:: bash

   source ~/NWQ-Sim/environment/setup_summit.sh

Launch multi-GPU execution for regular or interactive jobs:

.. code:: bash

   jsrun -n<GPUS> -a1 -g1 -c1 -brs <NWQ-Sim Command> -backend NVGPU_MPI

Replace with the total number of GPUs, and with the NWQ-Sim execution
command.

Running on Perlmutter HPC
~~~~~~~~~~~~~~~~~~~~~~~~~

To run NWQ-Sim on the Perlmutter Supercomputer, initilize the
environment first

.. code:: bash

   source ~/NWQ-Sim/environment/setup_perlmutter.sh

Launch multi-GPU execution for regular or interactive jobs:

.. code:: bash

   srun -C gpu -N <NODES> -n <GPUS> -c 1 --gpus-per-task=1 --gpu-bind=single:1 <NWQ-Sim Command> -backend NVGPU_MPI

Replace ``<NODES>`` with the number of compute nodes, ``<GPUS>`` with
the total number of GPUs, and ``<NWQ-Sim Command>`` with the NWQ-Sim
execution command.

NWQ-Sim for Chemistry Simulations
---------------------------------

NWQ-Sim is also capable of conducting chemistry simulations using the
XACC frontend, such as Variational Quantum Eigensolver (VQE)
simulations. This allows for a range of complex quantum chemical
computations using NWQ-Sim.

Below is an example of how to use NWQ-Sim with the XACC frontend for a
VQE simulation:

1. Install XACC by following the steps outlined in the `XACC
   repository <https://github.com/eclipse/xacc#build-from-source>`__.

Note, to successfully install and run XACC on Summit, you need to:

.. code:: bash

   module load openblas

Also, do not use all threads to build (make -j$(nproc –all) install)
which draws error, use:

.. code:: bash

   make -j8 install

To successfully install XACC on Frontier, you need to load the two
modules (the default cray-python/3.9 won’t work)

.. code:: bash

   module load openblas/0.3.17
   module load cray-python/3.10.10

2. Navigate to /NWQSim/xacc folder and create a source file.
3. Include the NWQ-Sim backend implementation in your code:

.. code:: cpp

   #include "nwq_accelerator.hpp"

4. Create an NWQAccelerator object:

.. code:: cpp

   auto nwq_acc = std::make_shared<xacc::quantum::NWQAccelerator>();

5. Utilize the NWQAccelerator with XACC. For example, you can run
   XACC-VQE:

.. code:: cpp

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

Replace the target source file in ``NWQ-Sim/xacc/CMakeList.txt`` and
build the project. The executable will be located at
``NWQ-Sim/build/xacc/nwq_xacc``.

Example Execution
~~~~~~~~~~~~~~~~~

Here, we illustrate an execution of the Adapt VQE simulation on a water
molecule using NWQ-Sim. The chart below depicts the variation in delta
energy per iteration of the algorithm. As observed, the desired chemical
accuracy is achieved around the 14th iteration, demonstrating the
effectiveness of the approach.

.. figure:: adapt_vqe.png
   :alt: Adapt VQE Delta Energy Chart

   Adapt VQE Delta Energy Chart

Please note, this is an example; actual results may vary based on the
specific quantum chemistry problem and the precision of your
Hamiltonian.

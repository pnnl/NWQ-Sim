# StabSim How To

This directory contains an example file for using STAB_CPU and STAB_GPU w/ 1D and 2D parallelization on NWQSim

## Building

# Executing the example file is streamlined by cmake. First ensure a compatible environment with gcc, cmake, and cuda.

Ex.
cuda/12.1
gcc/10.3.0
cmake/3.26.3

# Then create a build folder in the top level under NWQ-SIM.

mkdir build
cd build
cmake ..

# Make all projects or navigate to a project, stabilizer simulation in this case.

cd stabilizer/STABSim_exmaple
make

# Run the generated executable.

./example

## Creating your own file to simulate

# STABSim uses the same structure as other simulators on NWQSim, such as SVSim and DMSim. As such, it supports some of the same operations, and many others specific to stabilizers. Some examples are given in example.cpp, but the basic structure follows.

1. Create a circuit
    auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);

2. Create a state
    std::string backend = "nvgpu";
    std::string sim_method = "stab";
    auto state = BackendManager::create_state(backend, n_qubits, sim_method);

3. Simulate a circuit on the state. The state object retains state information, so different circuits can continue to evolve the same state.

    double timer = 0;
    state->sim(circuit, timer);

    # or, for a 2D stabilizer simulation

    double timer = 0;
    std::vector<int> gate_chunks = {16, 16, 16};
    state->sim2D(circuit, gate_chunks, timer);

    # where gate_chunks is a vector that defines how many sequential gates in 'circuit' can be parallelized at a time

4. View the completed state in the terminal

    state->print_res_state();

5. Store measurement results (measure_all is a peak rather than a strong measurement, which can be done with M gates)

    int shots = 1000;
    long long int* results = state->measure_all(shots);

    # or, if n_qubits is too long to fit in a long long use measure_all_long. The first 64 qubit results will be at the first column, next 64 at the second column, etc.

    int shots = 1000;
    long long int** results = state->measure_all_long(shots);


## Future Work

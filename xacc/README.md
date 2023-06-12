# NWQSim and XACC Integration 

This directory contains the components necessary to make our NWQSim compatible with the XACC framework. As of now, it is not ready to be directly compiled and installed as a custom accelerator in the XACC root folder. However, we have provided examples demonstrating its functionality with XACC. 

## Directory Contents

- **backend_runners/**: Runner files for different backends.
- **examples/**: A collection of example programs showing how to use NWQSim as a customer accelerator for XACC.

## Build Examples

Executing examples is made straightforward with the provided Makefile in the examples folder.

```bash
# Navigate to the examples directory.
cd examples

# Make example vqe and qpe testing prgrams, running with NWQSim-CPU backend. This generates examples including: xacc_vqe xacc_qpe xacc_vqe_adapt xacc_qpe_tester.
make

# The molecule used in Adapt VQE is a BE_5 molecule. It's downfolded Hamiltonian is stored in examples/be-5.txt and the program automatically loads it and executes it. To test Adapt VQE with different backends, run one of the following commands:
make xacc_adapt_cpu 
make xacc_adapt_omp 
make xacc_adapt_nvgpu 
make xacc_adapt_nvgpu_mpi
```
This will compile the necles and run the examples on your local machine. This does not install the simulator as a plugin in the XACC framework.essary fi


## Future Work
We are actively working on enhancing NWQSim's compatibility with XACC. Our current focus is on the integration of the DMSim backend, which necessitates the inclusion of a custom gate Instruction Representation (IR) for the native SX gate required by NWQSim's DMSim.

In addition to this, we aim to facilitate direct compilation and installation into XACC's root folder. This enhancement will simplify the process for users aiming to utilize our custom simulator within the XACC framework.

Future updates are directed towards simplifying this integration process and expanding NWQSim's versatility. We are also planning to expand our 'examples' repository to offer a broader range of use-cases demonstrating the NWQSim's functionality within the XACC environment. 
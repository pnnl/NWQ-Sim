# NWQ-Sim: Northwest Quantum System Simulation Environment

A Density Matrix Quantum Simulation Environment for Single-GPU/CPU, Single-Node-Multi-GPUs/CPUs and Multi-Nodes GPU/CPU Cluster. It supports Intel/AMD/IBM CPUs, NVIDIA/AMD GPUs. 

NWQ-Sim is a quantum system simulation environment supporting Qiskit and Q# as the frontend and CPU/GPU/Phi. It supports single-device, scale-up (single-node-multi-devices) and scale-out (multi-node cluster). It supports NVIDIA GPUs and AMD GPUs. It supports AVX512 and AVX2.

This released version supports Qiskit with CPU backends.

![alt text](pic.png)

## Current version
Latest version: **2.6**

## Supported Gate

Refer to Qiskit Circuit Library [Standard Gates](https://qiskit.org/documentation/apidoc/circuit_library.html), [OpenQASM-2](https://arxiv.org/pdf/1707.03429.pdf), [OpenQASM-3](https://arxiv.org/pdf/2104.14722.pdf), [Q# Intrinsics](https://docs.microsoft.com/en-us/qsharp/api/qsharp/microsoft.quantum.intrinsic), and [QIR](https://github.com/microsoft/qsharp-runtime/blob/main/src/Qir/Runtime/public/QSharpSimApi_I.hpp).

The following gates are natively implemented in NWQSim:

|Gates | Meaning | Qiskit | OpenQASM | Q#/QIR |
|:---: | ------- | ------ | -------- | ------ |
|X | Pauli-X bit flip | X | x | X | 
|Y | Pauli-Y bit and phase flip  | Y | y | Y |
|Z | Pauli-Z phase flip | Z | z | Z |
|H | Hadamard  | H | h | H |
|S | sqrt(Z) phase  | S | s | S |
|SDG | conjugate of sqrt(Z) | Sdg | sdg | AdjointS |
|T | sqrt(S) phase | T | t | T |
|TDG | conjugate of sqrt(S) | Tdg | tdg | AdjointT |
|RI  | global phase gate | - | - | R(PauliI) |
|RX | X-axis rotation | RX | rx | R(PauliX) | 
|RY | Y-axis rotation | RY | ry | R(PauliY) |
|RZ | Z-axis rotation | RZ | rz | R(PauliZ) |
|SX | sqrt(X) gate | SX | sx | - |
|P  | phase gate | Phase | - | R1 | 
|U  | unitary gate | U | u | - |
|CX | controlled X | CX | cx | ControlledX |
|CY | controlled Y | CY | cy | ControlledY |
|CZ | controlled Z | CZ | cz | ControlledZ |
|CH | controlled H | CH | ch | ControlledH |
|CS | controlled S | - | -  | ControlledS |
|CSDG | controlled SDG | - | - | ControlledAdjointS |
|CT | controlled T | - | - | ControlledT |
|CTDG | controlled TDG | - | - | ControlledAdjointT |
|CRI  | controlled global phase gate | - | - | ControlledR(PauliI) |
|CRX | controlled X-axis rotation | CRX | crx | ControlledR(PauliX) |
|CRY | controlled Y-axis rotation | CRY | cry | ControlledR(PauliY) | 
|CRZ | controlled Z-axis rotation | CRZ | crz | ControlledR(PauliZ) | 
|CSX | controlled SX | CSX | - | - |
|CP  | controlled phase gate | CPhase | - | - | 
|CU  | controlled unitary gate | CU | - | - | 
|ID  | identity gate | I | id | I |
|SWAP| swap gate | SWAP | swap | SWAP |
|M  | measure a qubit | Measure | measure | M |
|MA | measure all qubits | measure_all | measure | Measure |
|RESET | reset a qubit | Reset | reset | Reset |


The following gates are also supported in the wrapper:

|Gates | Meaning | Qiskit | OpenQASM | Q#/QIR |
|:---: | ------- | ------ | -------- | ------ |
|EXP | exponential of a multi-qubit Pauli operator | - | - | Exp | 
|ControlledExp | multi-controlled EXP | - | - | ControlledExp | 
|CCX | Toffoli gate | CCX | CCX | CCNOT |
|CSWAP | Fredkin gate | CSWAP | CSWAP | - |
|U1 | 1 parameter 0 pulse 1-qubit  | U1 | u1 | - | 
|U2 | 2 parameter 1 pulse 1-qubit  | U2 | u2 | - |
|U3 | 3 parameter 2 pulse 1-qubit  | U3 | u3 | - |
|RXX| rotation about XX | RXX | rxx | - |
|RYY| rotation about YY | RYY | - | - |
|RZZ| rotation about ZZ | RZZ | rzz | - |
|C1 | Arbitrary 1-qubit gate | - | - | - |
|C2 | Arbitrary 2-qubit gate | - | - | - |

## Configuration

Update CMakeList.txt.  

If you would like multi-threading support (by default enabled), update CMakeLists.txt and enable the option of "USE_OPENMP". Delete the build folder and recompile. 

If your CPU supports AVX512, you can enable the option "USE_AVX512" (by default not enabled) in CMakeLists.txt and recompile. It may bring around 2X on the performance.

## Prerequisite
NWQ-Sim requires the following packages.

|  Dependency  | Version | Comments |
|:-----------: | ------- | -------- |
|  GCC (or XL) | 5.2 or later (16.01 for xlc)  | |
|    OpenMP    | 4.0     | For single-node scale-up |
|  Python      | 3.6     | For Python API |
|  Pybind11    | 2.5.0   | For Python API |
| Qiskit-Terra | 0.8.0   | For Qiskit API |
|    Qiskit    | 0.24.0  | For Qiskit API |
|Qsharp-runtime| 0.15.2101125897 | For Q# QIR backend |
|   CMake      | - | For compile |


## Build

For Qiskit API, 
```text
python setup.py build
python setup.py install
```
You may need to install pybind11 (need to ensure pybind/pybind11.h can be found) and cmake. If you install pybind11 using pip, you may need to upgrade to a new version of pip if there are issues with finding pybind11.h.
```text
conda install cmake
conda install pybind11 #python -m pip install pybind11
```

## Programming
For Qiskit, after install the qiskit_nwqsim_provider, simply adjust the default Aer to nwqsim provider with dmsim backend:
```text
from math import pi
from qiskit import QuantumCircuit, transpile, execute
from qiskit_nwqsim_provider import NWQSimProvider

qc = QuantumCircuit(2,2)
qc.cu(pi/4,pi/5,pi/6,pi/7,0,1)
qc.x(0)
qc.swap(0,1)
qc.measure([0,1], [0,1])

nwqsim = NWQSimProvider('DMSimSimulator')
backend = nwqsim.backends['svsim_cpu']

trans_qc = transpile(qc, backend)
job = backend.run(trans_qc)
print(job.get_counts())

## You can also use execute:
#result = execute(qc, backend, seed_transpiler=111).result()
#print(result.get_counts(qc))
```

For Python, after having the compiled library libdmsim.so,
```text
import libdmsim as dmsim
n_qubits = 4
sim = dmsim.Simulation(n_qubits)
sim.H(0)
sim.X(1)
sim.H(2)
sim.CX(0,1)
res = sim.measure_all(10) #10 shots
## Print measurement results
for i in range(len(res)):
    print ("Test-"+str(i)+": " + "{0:b}".format(res[i]).zfill(n_qubits))
```

## Noisy Simulation
We currently support 1-qubit and 2-qubit noisy simulation through the noisy model of Qiskit. The noisy is assumed to be the same for all qubits and for all 1-qubit or 2-qubit gates. See the folowing example:
```text
import numpy as np
from math import pi
from qiskit import QuantumCircuit, transpile, execute
from qiskit import ClassicalRegister, QuantumRegister
from qiskit_nwqsim_provider import NWQSimProvider
from qiskit.quantum_info import SuperOp
from qiskit.providers.aer.noise import NoiseModel
from qiskit.providers.aer.noise import depolarizing_error
from qiskit.providers.aer.noise import amplitude_damping_error
from qiskit.providers.aer.noise import phase_damping_error

nwqsim = NWQSimProvider('DMSimSimulator')
print (nwqsim.backends)
backend = nwqsim.backends['dmsim']

#Bell State
circuit = QuantumCircuit(2,2)
circuit.h(0)
circuit.cx(0,1)
circuit.measure([0,1], [0,1])
noise_model = NoiseModel()

#Add amplitude and phase damping error for 1-qubit
param_amplitude_damping = 0.05 
param_phase_damping = 0.03
amplitude_error = amplitude_damping_error(param_amplitude_damping)
phase_error = phase_damping_error(param_phase_damping)
qerror_q1 = phase_error.compose(amplitude_error)

#Add depolarizing error for 2-qubit
param_q2 = 0.08
qerror_q2 = depolarizing_error(param_q2, 2)
q1_superop = SuperOp(qerror_q1)
q2_superop = SuperOp(qerror_q2)
backend.set_noise_1q_sop(q1_superop)
backend.set_noise_2q_sop(q2_superop)

result = execute(circuit, backend, seed_transpiler=111).result()
print(result.get_counts(circuit))
print(result.get_statevector())
```

## Citation format

If you find NWQ-Sim useful, please cite our SC-20 paper:
 - Ang Li, Omer Subasi, Xiu Yang, and Sriram Krishnamoorthy. "Density Matrix Quantum Circuit Simulation via the BSP Machine on Modern GPU Clusters." In Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis, 2020.
 - Ang Li and Sriram Krishnamoorthy. "SV-Sim: Scalable PGAS-based State Vector Simulation of Quantum Circuits" In Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis, 2021.


Bibtex:
```text
@inproceedings{li2020density,
    title={Density Matrix Quantum Circuit Simulation via the BSP Machine on Modern GPU Clusters},
    author={Li, Ang and Subasi, Omer and Yang, Xiu and Krishnamoorthy, Sriram},
    booktitle={Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis},
    year={2020}
}
@inproceedings{li2021svsim,
    title={SV-Sim: Scalable PGAS-based State Vector Simulation of Quantum Circuits},
    author={Li, Ang and Krishnamoorthy, Sriram},
    booktitle={Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis},
    year={2021}
}

``` 

## License

This project is licensed under the MIT License, see [LICENSE](LICENSE) file for details.

## Acknowledgments

**PNNL-IPID: 31919-E, ECCN: EAR99, IR: PNNL-SA-143160**

**PNNL-IPID: 32166-E, ECCN: EAR99, IR: PNNL-SA-161181**

This project is fully supported by the [Quantum Science Center (QSC)](https://qscience.org/).The Pacific Northwest National Laboratory (PNNL) is operated by Battelle for the U.S. Department of Energy (DOE) under contract DE-AC05-76RL01830. 

## Contributing

Please contact us If you'd like to contribute to NWQ-Sim. See the contact in our paper or my [webpage](http://www.angliphd.com). Have fun!

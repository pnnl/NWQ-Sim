// 19-Qubit Random Circuit Sampling
OPENQASM 2.0;
include "qelib1.inc";

// Define quantum and classical registers
qreg q[19];
creg c[19];

// Layer 1: Random single-qubit gates
h q[0];
t q[1];
s q[2];
rx(pi/3) q[3];
h q[4];
t q[5];
s q[6];
rx(pi/3) q[7];
h q[8];
t q[9];
s q[10];
rx(pi/3) q[11];
h q[12];
t q[13];
s q[14];
rx(pi/3) q[15];
h q[16];
t q[17];
s q[18];

// Layer 1: Random two-qubit gates
cx q[0], q[1];
cx q[2], q[3];
cx q[4], q[5];
cx q[6], q[7];
cx q[8], q[9];
cx q[10], q[11];
cx q[12], q[13];
cx q[14], q[15];
cx q[16], q[17];

// Layer 2: Random single-qubit gates
t q[0];
s q[1];
rx(pi/3) q[2];
h q[3];
t q[4];
s q[5];
rx(pi/3) q[6];
h q[7];
t q[8];
s q[9];
rx(pi/3) q[10];
h q[11];
t q[12];
s q[13];
rx(pi/3) q[14];
h q[15];
t q[16];
s q[17];
rx(pi/3) q[18];

// Layer 2: Random two-qubit gates
cx q[1], q[2];
cx q[3], q[4];
cx q[5], q[6];
cx q[7], q[8];
cx q[9], q[10];
cx q[11], q[12];
cx q[13], q[14];
cx q[15], q[16];
cx q[17], q[18];

// Layer 3: Random single-qubit gates
s q[0];
rx(pi/3) q[1];
h q[2];
t q[3];
s q[4];
rx(pi/3) q[5];
h q[6];
t q[7];
s q[8];
rx(pi/3) q[9];
h q[10];
t q[11];
s q[12];
rx(pi/3) q[13];
h q[14];
t q[15];
s q[16];
rx(pi/3) q[17];
h q[18];

// Layer 3: Random two-qubit gates
cx q[0], q[2];
cx q[3], q[5];
cx q[6], q[8];
cx q[9], q[11];
cx q[12], q[14];
cx q[15], q[17];

// Layer 4: Random single-qubit gates
rx(pi/3) q[0];
h q[1];
t q[2];
s q[3];
rx(pi/3) q[4];
h q[5];
t q[6];
s q[7];
rx(pi/3) q[8];
h q[9];
t q[10];
s q[11];
rx(pi/3) q[12];
h q[13];
t q[14];
s q[15];
rx(pi/3) q[16];
h q[17];
t q[18];

// Layer 4: Random two-qubit gates
cx q[1], q[3];
cx q[4], q[6];
cx q[7], q[9];
cx q[10], q[12];
cx q[13], q[15];
cx q[16], q[18];

// Layer 5: Random single-qubit gates
h q[0];
t q[1];
s q[2];
rx(pi/3) q[3];
h q[4];
t q[5];
s q[6];
rx(pi/3) q[7];
h q[8];
t q[9];
s q[10];
rx(pi/3) q[11];
h q[12];
t q[13];
s q[14];
rx(pi/3) q[15];
h q[16];
t q[17];
s q[18];

// Layer 5: Random two-qubit gates
cx q[0], q[1];
cx q[2], q[3];
cx q[4], q[5];
cx q[6], q[7];
cx q[8], q[9];
cx q[10], q[11];
cx q[12], q[13];
cx q[14], q[15];
cx q[16], q[17];

// Measure all qubits
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
measure q[4] -> c[4];
measure q[5] -> c[5];
measure q[6] -> c[6];
measure q[7] -> c[7];
measure q[8] -> c[8];
measure q[9] -> c[9];
measure q[10] -> c[10];
measure q[11] -> c[11];
measure q[12] -> c[12];
measure q[13] -> c[13];
measure q[14] -> c[14];
measure q[15] -> c[15];
measure q[16] -> c[16];
measure q[17] -> c[17];
measure q[18] -> c[18];

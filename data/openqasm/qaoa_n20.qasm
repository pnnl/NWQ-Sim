// 20-Qubit, 1-Layer QAOA Circuit for Max-Cut on a Ring Graph
OPENQASM 2.0;
include "qelib1.inc";

// Define quantum and classical registers
qreg q[20];
creg c[20];

// Initialize all qubits to |+⟩ state
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
h q[7];
h q[8];
h q[9];
h q[10];
h q[11];
h q[12];
h q[13];
h q[14];
h q[15];
h q[16];
h q[17];
h q[18];
h q[19];

// Set variational parameters gamma and beta
// For this example, we use fixed values
// gamma = pi/4
// beta = pi/8

// Cost Hamiltonian U_C(gamma)
// Apply ZZ interactions between neighboring qubits
// Implemented using CNOT and RZ gates

// Edge between q[0] and q[1]
cx q[0], q[1];
rz(-pi/2) q[1];
cx q[0], q[1];

// Edge between q[1] and q[2]
cx q[1], q[2];
rz(-pi/2) q[2];
cx q[1], q[2];

// Edge between q[2] and q[3]
cx q[2], q[3];
rz(-pi/2) q[3];
cx q[2], q[3];

// Edge between q[3] and q[4]
cx q[3], q[4];
rz(-pi/2) q[4];
cx q[3], q[4];

// Edge between q[4] and q[5]
cx q[4], q[5];
rz(-pi/2) q[5];
cx q[4], q[5];

// Edge between q[5] and q[6]
cx q[5], q[6];
rz(-pi/2) q[6];
cx q[5], q[6];

// Edge between q[6] and q[7]
cx q[6], q[7];
rz(-pi/2) q[7];
cx q[6], q[7];

// Edge between q[7] and q[8]
cx q[7], q[8];
rz(-pi/2) q[8];
cx q[7], q[8];

// Edge between q[8] and q[9]
cx q[8], q[9];
rz(-pi/2) q[9];
cx q[8], q[9];

// Edge between q[9] and q[10]
cx q[9], q[10];
rz(-pi/2) q[10];
cx q[9], q[10];

// Edge between q[10] and q[11]
cx q[10], q[11];
rz(-pi/2) q[11];
cx q[10], q[11];

// Edge between q[11] and q[12]
cx q[11], q[12];
rz(-pi/2) q[12];
cx q[11], q[12];

// Edge between q[12] and q[13]
cx q[12], q[13];
rz(-pi/2) q[13];
cx q[12], q[13];

// Edge between q[13] and q[14]
cx q[13], q[14];
rz(-pi/2) q[14];
cx q[13], q[14];

// Edge between q[14] and q[15]
cx q[14], q[15];
rz(-pi/2) q[15];
cx q[14], q[15];

// Edge between q[15] and q[16]
cx q[15], q[16];
rz(-pi/2) q[16];
cx q[15], q[16];

// Edge between q[16] and q[17]
cx q[16], q[17];
rz(-pi/2) q[17];
cx q[16], q[17];

// Edge between q[17] and q[18]
cx q[17], q[18];
rz(-pi/2) q[18];
cx q[17], q[18];

// Edge between q[18] and q[19]
cx q[18], q[19];
rz(-pi/2) q[19];
cx q[18], q[19];

// Edge between q[19] and q[0] to complete the ring
cx q[19], q[0];
rz(-pi/2) q[0];
cx q[19], q[0];

// Mixer Hamiltonian U_B(beta)
// Apply RX rotations to all qubits
// beta = pi/8
rx(pi/4) q[0];
rx(pi/4) q[1];
rx(pi/4) q[2];
rx(pi/4) q[3];
rx(pi/4) q[4];
rx(pi/4) q[5];
rx(pi/4) q[6];
rx(pi/4) q[7];
rx(pi/4) q[8];
rx(pi/4) q[9];
rx(pi/4) q[10];
rx(pi/4) q[11];
rx(pi/4) q[12];
rx(pi/4) q[13];
rx(pi/4) q[14];
rx(pi/4) q[15];
rx(pi/4) q[16];
rx(pi/4) q[17];
rx(pi/4) q[18];
rx(pi/4) q[19];

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
measure q[19] -> c[19];

// 6 qubit QASM program to convert to X, Y, and Z bases
OPENQASM 2.0;
include "qelib1.inc";

// Declare 6 qubits
qreg q[6];

h q[0];

// q[0] and q[1] remain in the Z basis (no gates applied)
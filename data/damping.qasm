//2-qubit GHZ state in OpenQASM 2.0
OPENQASM 2.0;
include "qelib1.inc";

qreg q[2];
// creg c[2];

h q[0];    
cx q[0], q[1];

// measure q[0] -> c[0];

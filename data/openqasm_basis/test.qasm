OPENQASM 2.0;
include "qelib1.inc";
qreg q[2];
creg c[2];

x q[0];
cx q[0],q[1];

barrier q[1],q[0];
measure q[0] -> c[0];
measure q[1] -> c[1];

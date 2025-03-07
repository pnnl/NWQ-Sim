OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c[3];
x q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
cx q[1],q[2];
cx q[0],q[1];
cx q[1],q[0];
cx q[0],q[1];
rz(-pi/4) q[2];
cx q[1],q[2];
cx q[0],q[1];
cx q[1],q[0];
cx q[0],q[1];
rz(pi/4) q[2];
cx q[1],q[2];
rz(-pi/4) q[1];
cx q[0],q[1];
cx q[1],q[0];
cx q[0],q[1];
rz(-pi/4) q[2];
cx q[1],q[2];
cx q[1],q[0];
rz(-pi/4) q[0];
cx q[1],q[0];
rz(pi/2) q[0];
rz(pi/4) q[1];
rz(pi/4) q[2];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[1],q[0],q[2],q[3],q[4],q[5],q[6],q[7],q[8],q[9],q[10],q[11],q[12],q[13],q[14],q[15];
measure q[1] -> c[0];
measure q[0] -> c[1];
measure q[2] -> c[2];

OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
creg c0[1];
creg c1[1];
creg c2[1];
creg c3[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[1],q[2],q[3];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
measure q[0] -> c0[0];
if(c0==1) rz(pi/2) q[1];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
if(c0==1) rz(pi/4) q[2];
if(c0==1) rz(pi/8) q[3];
measure q[1] -> c1[0];
if(c1==1) rz(pi/2) q[2];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
if(c1==1) rz(pi/4) q[3];
measure q[2] -> c2[0];
if(c2==1) rz(pi/2) q[3];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
measure q[3] -> c3[0];

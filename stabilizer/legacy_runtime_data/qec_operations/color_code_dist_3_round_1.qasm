OPENQASM 2.0;
include "qelib1.inc";
gate cxyz q0 { U(pi/2, 0, pi/2) q0; }

qreg q[10];
creg rec[13];

reset q[0];
reset q[1];
reset q[2];
reset q[3];
reset q[4];
reset q[5];
reset q[6];
reset q[7];
reset q[8];
reset q[9];
barrier q;

cxyz q[0];
cxyz q[1];
cxyz q[3];
cxyz q[5];
cxyz q[6];
cxyz q[7];
cxyz q[9];
barrier q;

cx q[5], q[4];
cx q[3], q[2];
barrier q;

cx q[7], q[4];
cx q[6], q[2];
barrier q;

cx q[1], q[4];
cx q[6], q[8];
barrier q;

cx q[1], q[2];
cx q[7], q[8];
barrier q;

cx q[5], q[2];
cx q[9], q[8];
barrier q;

cx q[0], q[4];
cx q[5], q[8];
barrier q;

measure q[2] -> rec[0]; reset q[2]; // decomposed MR
measure q[4] -> rec[1]; reset q[4]; // decomposed MR
measure q[8] -> rec[2]; reset q[8]; // decomposed MR
barrier q;

cxyz q[0];
cxyz q[1];
cxyz q[3];
cxyz q[5];
cxyz q[6];
cxyz q[7];
cxyz q[9];
barrier q;

cx q[5], q[4];
cx q[3], q[2];
barrier q;

cx q[7], q[4];
cx q[6], q[2];
barrier q;

cx q[1], q[4];
cx q[6], q[8];
barrier q;

cx q[1], q[2];
cx q[7], q[8];
barrier q;

cx q[5], q[2];
cx q[9], q[8];
barrier q;

cx q[0], q[4];
cx q[5], q[8];
barrier q;

measure q[2] -> rec[3]; reset q[2]; // decomposed MR
measure q[4] -> rec[4]; reset q[4]; // decomposed MR
measure q[8] -> rec[5]; reset q[8]; // decomposed MR
s q[0]; s q[0]; s q[0]; h q[0]; measure q[0] -> rec[6]; h q[0]; s q[0]; // decomposed MY
s q[1]; s q[1]; s q[1]; h q[1]; measure q[1] -> rec[7]; h q[1]; s q[1]; // decomposed MY
s q[3]; s q[3]; s q[3]; h q[3]; measure q[3] -> rec[8]; h q[3]; s q[3]; // decomposed MY
s q[5]; s q[5]; s q[5]; h q[5]; measure q[5] -> rec[9]; h q[5]; s q[5]; // decomposed MY
s q[6]; s q[6]; s q[6]; h q[6]; measure q[6] -> rec[10]; h q[6]; s q[6]; // decomposed MY
s q[7]; s q[7]; s q[7]; h q[7]; measure q[7] -> rec[11]; h q[7]; s q[7]; // decomposed MY
s q[9]; s q[9]; s q[9]; h q[9]; measure q[9] -> rec[12]; h q[9]; s q[9]; // decomposed MY

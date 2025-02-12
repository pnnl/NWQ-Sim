OPENQASM 2.0;
include "qelib1.inc";

qreg q[26];
creg rec[17];

reset q[1];
reset q[3];
reset q[5];
reset q[8];
reset q[10];
reset q[12];
reset q[15];
reset q[17];
reset q[19];
reset q[2];
reset q[9];
reset q[11];
reset q[13];
reset q[14];
reset q[16];
reset q[18];
reset q[25];
barrier q;

h q[2];
h q[11];
h q[16];
h q[25];
barrier q;

cx q[2], q[3];
cx q[16], q[17];
cx q[11], q[12];
cx q[15], q[14];
cx q[10], q[9];
cx q[19], q[18];
barrier q;

cx q[2], q[1];
cx q[16], q[15];
cx q[11], q[10];
cx q[8], q[14];
cx q[3], q[9];
cx q[12], q[18];
barrier q;

cx q[16], q[10];
cx q[11], q[5];
cx q[25], q[19];
cx q[8], q[9];
cx q[17], q[18];
cx q[12], q[13];
barrier q;

cx q[16], q[8];
cx q[11], q[3];
cx q[25], q[17];
cx q[1], q[9];
cx q[10], q[18];
cx q[5], q[13];
barrier q;

h q[2];
h q[11];
h q[16];
h q[25];
barrier q;

measure q[2] -> rec[0]; reset q[2]; // decomposed MR
measure q[9] -> rec[1]; reset q[9]; // decomposed MR
measure q[11] -> rec[2]; reset q[11]; // decomposed MR
measure q[13] -> rec[3]; reset q[13]; // decomposed MR
measure q[14] -> rec[4]; reset q[14]; // decomposed MR
measure q[16] -> rec[5]; reset q[16]; // decomposed MR
measure q[18] -> rec[6]; reset q[18]; // decomposed MR
measure q[25] -> rec[7]; reset q[25]; // decomposed MR
measure q[1] -> rec[8];
measure q[3] -> rec[9];
measure q[5] -> rec[10];
measure q[8] -> rec[11];
measure q[10] -> rec[12];
measure q[12] -> rec[13];
measure q[15] -> rec[14];
measure q[17] -> rec[15];
measure q[19] -> rec[16];
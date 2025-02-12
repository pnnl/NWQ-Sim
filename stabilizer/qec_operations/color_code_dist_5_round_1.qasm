OPENQASM 2.0;
include "qelib1.inc";
gate cxyz q0 { U(pi/2, 0, pi/2) q0; }

qreg q[28];
creg rec[37];

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
reset q[10];
reset q[11];
reset q[12];
reset q[13];
reset q[14];
reset q[15];
reset q[16];
reset q[17];
reset q[18];
reset q[19];
reset q[20];
reset q[21];
reset q[22];
reset q[23];
reset q[24];
reset q[25];
reset q[26];
reset q[27];
barrier q;

cxyz q[0];
cxyz q[1];
cxyz q[3];
cxyz q[4];
cxyz q[6];
cxyz q[8];
cxyz q[9];
cxyz q[11];
cxyz q[12];
cxyz q[13];
cxyz q[15];
cxyz q[16];
cxyz q[18];
cxyz q[19];
cxyz q[21];
cxyz q[23];
cxyz q[24];
cxyz q[25];
cxyz q[27];
barrier q;

cx q[8], q[7];
cx q[3], q[2];
cx q[15], q[14];
cx q[23], q[22];
cx q[11], q[10];
cx q[21], q[20];
cx q[6], q[5];
barrier q;

cx q[13], q[7];
cx q[9], q[2];
cx q[19], q[14];
cx q[25], q[22];
cx q[16], q[10];
cx q[24], q[20];
cx q[12], q[5];
barrier q;

cx q[1], q[7];
cx q[9], q[14];
cx q[19], q[22];
cx q[4], q[10];
cx q[16], q[20];
cx q[24], q[26];
cx q[12], q[17];
barrier q;

cx q[1], q[2];
cx q[13], q[14];
cx q[9], q[10];
cx q[19], q[20];
cx q[25], q[26];
cx q[4], q[5];
cx q[16], q[17];
barrier q;

cx q[8], q[2];
cx q[18], q[14];
cx q[15], q[10];
cx q[23], q[20];
cx q[27], q[26];
cx q[11], q[5];
cx q[21], q[17];
barrier q;

cx q[0], q[7];
cx q[8], q[14];
cx q[18], q[22];
cx q[3], q[10];
cx q[15], q[20];
cx q[23], q[26];
cx q[11], q[17];
barrier q;

measure q[2] -> rec[0]; reset q[2]; // decomposed MR
measure q[5] -> rec[1]; reset q[5]; // decomposed MR
measure q[7] -> rec[2]; reset q[7]; // decomposed MR
measure q[10] -> rec[3]; reset q[10]; // decomposed MR
measure q[14] -> rec[4]; reset q[14]; // decomposed MR
measure q[17] -> rec[5]; reset q[17]; // decomposed MR
measure q[20] -> rec[6]; reset q[20]; // decomposed MR
measure q[22] -> rec[7]; reset q[22]; // decomposed MR
measure q[26] -> rec[8]; reset q[26]; // decomposed MR
barrier q;

cxyz q[0];
cxyz q[1];
cxyz q[3];
cxyz q[4];
cxyz q[6];
cxyz q[8];
cxyz q[9];
cxyz q[11];
cxyz q[12];
cxyz q[13];
cxyz q[15];
cxyz q[16];
cxyz q[18];
cxyz q[19];
cxyz q[21];
cxyz q[23];
cxyz q[24];
cxyz q[25];
cxyz q[27];
barrier q;

cx q[8], q[7];
cx q[3], q[2];
cx q[15], q[14];
cx q[23], q[22];
cx q[11], q[10];
cx q[21], q[20];
cx q[6], q[5];
barrier q;

cx q[13], q[7];
cx q[9], q[2];
cx q[19], q[14];
cx q[25], q[22];
cx q[16], q[10];
cx q[24], q[20];
cx q[12], q[5];
barrier q;

cx q[1], q[7];
cx q[9], q[14];
cx q[19], q[22];
cx q[4], q[10];
cx q[16], q[20];
cx q[24], q[26];
cx q[12], q[17];
barrier q;

cx q[1], q[2];
cx q[13], q[14];
cx q[9], q[10];
cx q[19], q[20];
cx q[25], q[26];
cx q[4], q[5];
cx q[16], q[17];
barrier q;

cx q[8], q[2];
cx q[18], q[14];
cx q[15], q[10];
cx q[23], q[20];
cx q[27], q[26];
cx q[11], q[5];
cx q[21], q[17];
barrier q;

cx q[0], q[7];
cx q[8], q[14];
cx q[18], q[22];
cx q[3], q[10];
cx q[15], q[20];
cx q[23], q[26];
cx q[11], q[17];
barrier q;

measure q[2] -> rec[9]; reset q[2]; // decomposed MR
measure q[5] -> rec[10]; reset q[5]; // decomposed MR
measure q[7] -> rec[11]; reset q[7]; // decomposed MR
measure q[10] -> rec[12]; reset q[10]; // decomposed MR
measure q[14] -> rec[13]; reset q[14]; // decomposed MR
measure q[17] -> rec[14]; reset q[17]; // decomposed MR
measure q[20] -> rec[15]; reset q[20]; // decomposed MR
measure q[22] -> rec[16]; reset q[22]; // decomposed MR
measure q[26] -> rec[17]; reset q[26]; // decomposed MR
s q[0]; s q[0]; s q[0]; h q[0]; measure q[0] -> rec[18]; h q[0]; s q[0]; // decomposed MY
s q[1]; s q[1]; s q[1]; h q[1]; measure q[1] -> rec[19]; h q[1]; s q[1]; // decomposed MY
s q[3]; s q[3]; s q[3]; h q[3]; measure q[3] -> rec[20]; h q[3]; s q[3]; // decomposed MY
s q[4]; s q[4]; s q[4]; h q[4]; measure q[4] -> rec[21]; h q[4]; s q[4]; // decomposed MY
s q[6]; s q[6]; s q[6]; h q[6]; measure q[6] -> rec[22]; h q[6]; s q[6]; // decomposed MY
s q[8]; s q[8]; s q[8]; h q[8]; measure q[8] -> rec[23]; h q[8]; s q[8]; // decomposed MY
s q[9]; s q[9]; s q[9]; h q[9]; measure q[9] -> rec[24]; h q[9]; s q[9]; // decomposed MY
s q[11]; s q[11]; s q[11]; h q[11]; measure q[11] -> rec[25]; h q[11]; s q[11]; // decomposed MY
s q[12]; s q[12]; s q[12]; h q[12]; measure q[12] -> rec[26]; h q[12]; s q[12]; // decomposed MY
s q[13]; s q[13]; s q[13]; h q[13]; measure q[13] -> rec[27]; h q[13]; s q[13]; // decomposed MY
s q[15]; s q[15]; s q[15]; h q[15]; measure q[15] -> rec[28]; h q[15]; s q[15]; // decomposed MY
s q[16]; s q[16]; s q[16]; h q[16]; measure q[16] -> rec[29]; h q[16]; s q[16]; // decomposed MY
s q[18]; s q[18]; s q[18]; h q[18]; measure q[18] -> rec[30]; h q[18]; s q[18]; // decomposed MY
s q[19]; s q[19]; s q[19]; h q[19]; measure q[19] -> rec[31]; h q[19]; s q[19]; // decomposed MY
s q[21]; s q[21]; s q[21]; h q[21]; measure q[21] -> rec[32]; h q[21]; s q[21]; // decomposed MY
s q[23]; s q[23]; s q[23]; h q[23]; measure q[23] -> rec[33]; h q[23]; s q[23]; // decomposed MY
s q[24]; s q[24]; s q[24]; h q[24]; measure q[24] -> rec[34]; h q[24]; s q[24]; // decomposed MY
s q[25]; s q[25]; s q[25]; h q[25]; measure q[25] -> rec[35]; h q[25]; s q[25]; // decomposed MY
s q[27]; s q[27]; s q[27]; h q[27]; measure q[27] -> rec[36]; h q[27]; s q[27]; // decomposed MY

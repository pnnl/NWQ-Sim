OPENQASM 2.0;
include "qelib1.inc";
gate cxyz q0 { U(pi/2, 0, pi/2) q0; }

qreg q[55];
creg rec[73];

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
reset q[28];
reset q[29];
reset q[30];
reset q[31];
reset q[32];
reset q[33];
reset q[34];
reset q[35];
reset q[36];
reset q[37];
reset q[38];
reset q[39];
reset q[40];
reset q[41];
reset q[42];
reset q[43];
reset q[44];
reset q[45];
reset q[46];
reset q[47];
reset q[48];
reset q[49];
reset q[50];
reset q[51];
reset q[52];
reset q[53];
reset q[54];
barrier q;

cxyz q[0];
cxyz q[1];
cxyz q[3];
cxyz q[4];
cxyz q[6];
cxyz q[7];
cxyz q[9];
cxyz q[11];
cxyz q[12];
cxyz q[14];
cxyz q[15];
cxyz q[17];
cxyz q[18];
cxyz q[19];
cxyz q[21];
cxyz q[22];
cxyz q[24];
cxyz q[25];
cxyz q[27];
cxyz q[28];
cxyz q[30];
cxyz q[31];
cxyz q[33];
cxyz q[35];
cxyz q[36];
cxyz q[38];
cxyz q[39];
cxyz q[40];
cxyz q[42];
cxyz q[43];
cxyz q[45];
cxyz q[46];
cxyz q[48];
cxyz q[50];
cxyz q[51];
cxyz q[52];
cxyz q[54];
barrier q;

cx q[11], q[10];
cx q[3], q[2];
cx q[21], q[20];
cx q[35], q[34];
cx q[14], q[13];
cx q[30], q[29];
cx q[42], q[41];
cx q[50], q[49];
cx q[6], q[5];
cx q[24], q[23];
cx q[38], q[37];
cx q[48], q[47];
cx q[17], q[16];
cx q[33], q[32];
cx q[9], q[8];
barrier q;

cx q[19], q[10];
cx q[12], q[2];
cx q[28], q[20];
cx q[40], q[34];
cx q[22], q[13];
cx q[36], q[29];
cx q[46], q[41];
cx q[52], q[49];
cx q[15], q[5];
cx q[31], q[23];
cx q[43], q[37];
cx q[51], q[47];
cx q[25], q[16];
cx q[39], q[32];
cx q[18], q[8];
barrier q;

cx q[1], q[10];
cx q[12], q[20];
cx q[28], q[34];
cx q[4], q[13];
cx q[22], q[29];
cx q[36], q[41];
cx q[46], q[49];
cx q[15], q[23];
cx q[31], q[37];
cx q[43], q[47];
cx q[51], q[53];
cx q[7], q[16];
cx q[25], q[32];
cx q[39], q[44];
cx q[18], q[26];
barrier q;

cx q[1], q[2];
cx q[19], q[20];
cx q[12], q[13];
cx q[28], q[29];
cx q[40], q[41];
cx q[4], q[5];
cx q[22], q[23];
cx q[36], q[37];
cx q[46], q[47];
cx q[52], q[53];
cx q[15], q[16];
cx q[31], q[32];
cx q[43], q[44];
cx q[7], q[8];
cx q[25], q[26];
barrier q;

cx q[11], q[2];
cx q[27], q[20];
cx q[21], q[13];
cx q[35], q[29];
cx q[45], q[41];
cx q[14], q[5];
cx q[30], q[23];
cx q[42], q[37];
cx q[50], q[47];
cx q[54], q[53];
cx q[24], q[16];
cx q[38], q[32];
cx q[48], q[44];
cx q[17], q[8];
cx q[33], q[26];
barrier q;

cx q[0], q[10];
cx q[11], q[20];
cx q[27], q[34];
cx q[3], q[13];
cx q[21], q[29];
cx q[35], q[41];
cx q[45], q[49];
cx q[14], q[23];
cx q[30], q[37];
cx q[42], q[47];
cx q[50], q[53];
cx q[6], q[16];
cx q[24], q[32];
cx q[38], q[44];
cx q[17], q[26];
barrier q;

measure q[2] -> rec[0]; reset q[2]; // decomposed MR
measure q[5] -> rec[1]; reset q[5]; // decomposed MR
measure q[8] -> rec[2]; reset q[8]; // decomposed MR
measure q[10] -> rec[3]; reset q[10]; // decomposed MR
measure q[13] -> rec[4]; reset q[13]; // decomposed MR
measure q[16] -> rec[5]; reset q[16]; // decomposed MR
measure q[20] -> rec[6]; reset q[20]; // decomposed MR
measure q[23] -> rec[7]; reset q[23]; // decomposed MR
measure q[26] -> rec[8]; reset q[26]; // decomposed MR
measure q[29] -> rec[9]; reset q[29]; // decomposed MR
measure q[32] -> rec[10]; reset q[32]; // decomposed MR
measure q[34] -> rec[11]; reset q[34]; // decomposed MR
measure q[37] -> rec[12]; reset q[37]; // decomposed MR
measure q[41] -> rec[13]; reset q[41]; // decomposed MR
measure q[44] -> rec[14]; reset q[44]; // decomposed MR
measure q[47] -> rec[15]; reset q[47]; // decomposed MR
measure q[49] -> rec[16]; reset q[49]; // decomposed MR
measure q[53] -> rec[17]; reset q[53]; // decomposed MR
barrier q;

cxyz q[0];
cxyz q[1];
cxyz q[3];
cxyz q[4];
cxyz q[6];
cxyz q[7];
cxyz q[9];
cxyz q[11];
cxyz q[12];
cxyz q[14];
cxyz q[15];
cxyz q[17];
cxyz q[18];
cxyz q[19];
cxyz q[21];
cxyz q[22];
cxyz q[24];
cxyz q[25];
cxyz q[27];
cxyz q[28];
cxyz q[30];
cxyz q[31];
cxyz q[33];
cxyz q[35];
cxyz q[36];
cxyz q[38];
cxyz q[39];
cxyz q[40];
cxyz q[42];
cxyz q[43];
cxyz q[45];
cxyz q[46];
cxyz q[48];
cxyz q[50];
cxyz q[51];
cxyz q[52];
cxyz q[54];
barrier q;

cx q[11], q[10];
cx q[3], q[2];
cx q[21], q[20];
cx q[35], q[34];
cx q[14], q[13];
cx q[30], q[29];
cx q[42], q[41];
cx q[50], q[49];
cx q[6], q[5];
cx q[24], q[23];
cx q[38], q[37];
cx q[48], q[47];
cx q[17], q[16];
cx q[33], q[32];
cx q[9], q[8];
barrier q;

cx q[19], q[10];
cx q[12], q[2];
cx q[28], q[20];
cx q[40], q[34];
cx q[22], q[13];
cx q[36], q[29];
cx q[46], q[41];
cx q[52], q[49];
cx q[15], q[5];
cx q[31], q[23];
cx q[43], q[37];
cx q[51], q[47];
cx q[25], q[16];
cx q[39], q[32];
cx q[18], q[8];
barrier q;

cx q[1], q[10];
cx q[12], q[20];
cx q[28], q[34];
cx q[4], q[13];
cx q[22], q[29];
cx q[36], q[41];
cx q[46], q[49];
cx q[15], q[23];
cx q[31], q[37];
cx q[43], q[47];
cx q[51], q[53];
cx q[7], q[16];
cx q[25], q[32];
cx q[39], q[44];
cx q[18], q[26];
barrier q;

cx q[1], q[2];
cx q[19], q[20];
cx q[12], q[13];
cx q[28], q[29];
cx q[40], q[41];
cx q[4], q[5];
cx q[22], q[23];
cx q[36], q[37];
cx q[46], q[47];
cx q[52], q[53];
cx q[15], q[16];
cx q[31], q[32];
cx q[43], q[44];
cx q[7], q[8];
cx q[25], q[26];
barrier q;

cx q[11], q[2];
cx q[27], q[20];
cx q[21], q[13];
cx q[35], q[29];
cx q[45], q[41];
cx q[14], q[5];
cx q[30], q[23];
cx q[42], q[37];
cx q[50], q[47];
cx q[54], q[53];
cx q[24], q[16];
cx q[38], q[32];
cx q[48], q[44];
cx q[17], q[8];
cx q[33], q[26];
barrier q;

cx q[0], q[10];
cx q[11], q[20];
cx q[27], q[34];
cx q[3], q[13];
cx q[21], q[29];
cx q[35], q[41];
cx q[45], q[49];
cx q[14], q[23];
cx q[30], q[37];
cx q[42], q[47];
cx q[50], q[53];
cx q[6], q[16];
cx q[24], q[32];
cx q[38], q[44];
cx q[17], q[26];
barrier q;

measure q[2] -> rec[18]; reset q[2]; // decomposed MR
measure q[5] -> rec[19]; reset q[5]; // decomposed MR
measure q[8] -> rec[20]; reset q[8]; // decomposed MR
measure q[10] -> rec[21]; reset q[10]; // decomposed MR
measure q[13] -> rec[22]; reset q[13]; // decomposed MR
measure q[16] -> rec[23]; reset q[16]; // decomposed MR
measure q[20] -> rec[24]; reset q[20]; // decomposed MR
measure q[23] -> rec[25]; reset q[23]; // decomposed MR
measure q[26] -> rec[26]; reset q[26]; // decomposed MR
measure q[29] -> rec[27]; reset q[29]; // decomposed MR
measure q[32] -> rec[28]; reset q[32]; // decomposed MR
measure q[34] -> rec[29]; reset q[34]; // decomposed MR
measure q[37] -> rec[30]; reset q[37]; // decomposed MR
measure q[41] -> rec[31]; reset q[41]; // decomposed MR
measure q[44] -> rec[32]; reset q[44]; // decomposed MR
measure q[47] -> rec[33]; reset q[47]; // decomposed MR
measure q[49] -> rec[34]; reset q[49]; // decomposed MR
measure q[53] -> rec[35]; reset q[53]; // decomposed MR
s q[0]; s q[0]; s q[0]; h q[0]; measure q[0] -> rec[36]; h q[0]; s q[0]; // decomposed MY
s q[1]; s q[1]; s q[1]; h q[1]; measure q[1] -> rec[37]; h q[1]; s q[1]; // decomposed MY
s q[3]; s q[3]; s q[3]; h q[3]; measure q[3] -> rec[38]; h q[3]; s q[3]; // decomposed MY
s q[4]; s q[4]; s q[4]; h q[4]; measure q[4] -> rec[39]; h q[4]; s q[4]; // decomposed MY
s q[6]; s q[6]; s q[6]; h q[6]; measure q[6] -> rec[40]; h q[6]; s q[6]; // decomposed MY
s q[7]; s q[7]; s q[7]; h q[7]; measure q[7] -> rec[41]; h q[7]; s q[7]; // decomposed MY
s q[9]; s q[9]; s q[9]; h q[9]; measure q[9] -> rec[42]; h q[9]; s q[9]; // decomposed MY
s q[11]; s q[11]; s q[11]; h q[11]; measure q[11] -> rec[43]; h q[11]; s q[11]; // decomposed MY
s q[12]; s q[12]; s q[12]; h q[12]; measure q[12] -> rec[44]; h q[12]; s q[12]; // decomposed MY
s q[14]; s q[14]; s q[14]; h q[14]; measure q[14] -> rec[45]; h q[14]; s q[14]; // decomposed MY
s q[15]; s q[15]; s q[15]; h q[15]; measure q[15] -> rec[46]; h q[15]; s q[15]; // decomposed MY
s q[17]; s q[17]; s q[17]; h q[17]; measure q[17] -> rec[47]; h q[17]; s q[17]; // decomposed MY
s q[18]; s q[18]; s q[18]; h q[18]; measure q[18] -> rec[48]; h q[18]; s q[18]; // decomposed MY
s q[19]; s q[19]; s q[19]; h q[19]; measure q[19] -> rec[49]; h q[19]; s q[19]; // decomposed MY
s q[21]; s q[21]; s q[21]; h q[21]; measure q[21] -> rec[50]; h q[21]; s q[21]; // decomposed MY
s q[22]; s q[22]; s q[22]; h q[22]; measure q[22] -> rec[51]; h q[22]; s q[22]; // decomposed MY
s q[24]; s q[24]; s q[24]; h q[24]; measure q[24] -> rec[52]; h q[24]; s q[24]; // decomposed MY
s q[25]; s q[25]; s q[25]; h q[25]; measure q[25] -> rec[53]; h q[25]; s q[25]; // decomposed MY
s q[27]; s q[27]; s q[27]; h q[27]; measure q[27] -> rec[54]; h q[27]; s q[27]; // decomposed MY
s q[28]; s q[28]; s q[28]; h q[28]; measure q[28] -> rec[55]; h q[28]; s q[28]; // decomposed MY
s q[30]; s q[30]; s q[30]; h q[30]; measure q[30] -> rec[56]; h q[30]; s q[30]; // decomposed MY
s q[31]; s q[31]; s q[31]; h q[31]; measure q[31] -> rec[57]; h q[31]; s q[31]; // decomposed MY
s q[33]; s q[33]; s q[33]; h q[33]; measure q[33] -> rec[58]; h q[33]; s q[33]; // decomposed MY
s q[35]; s q[35]; s q[35]; h q[35]; measure q[35] -> rec[59]; h q[35]; s q[35]; // decomposed MY
s q[36]; s q[36]; s q[36]; h q[36]; measure q[36] -> rec[60]; h q[36]; s q[36]; // decomposed MY
s q[38]; s q[38]; s q[38]; h q[38]; measure q[38] -> rec[61]; h q[38]; s q[38]; // decomposed MY
s q[39]; s q[39]; s q[39]; h q[39]; measure q[39] -> rec[62]; h q[39]; s q[39]; // decomposed MY
s q[40]; s q[40]; s q[40]; h q[40]; measure q[40] -> rec[63]; h q[40]; s q[40]; // decomposed MY
s q[42]; s q[42]; s q[42]; h q[42]; measure q[42] -> rec[64]; h q[42]; s q[42]; // decomposed MY
s q[43]; s q[43]; s q[43]; h q[43]; measure q[43] -> rec[65]; h q[43]; s q[43]; // decomposed MY
s q[45]; s q[45]; s q[45]; h q[45]; measure q[45] -> rec[66]; h q[45]; s q[45]; // decomposed MY
s q[46]; s q[46]; s q[46]; h q[46]; measure q[46] -> rec[67]; h q[46]; s q[46]; // decomposed MY
s q[48]; s q[48]; s q[48]; h q[48]; measure q[48] -> rec[68]; h q[48]; s q[48]; // decomposed MY
s q[50]; s q[50]; s q[50]; h q[50]; measure q[50] -> rec[69]; h q[50]; s q[50]; // decomposed MY
s q[51]; s q[51]; s q[51]; h q[51]; measure q[51] -> rec[70]; h q[51]; s q[51]; // decomposed MY
s q[52]; s q[52]; s q[52]; h q[52]; measure q[52] -> rec[71]; h q[52]; s q[52]; // decomposed MY
s q[54]; s q[54]; s q[54]; h q[54]; measure q[54] -> rec[72]; h q[54]; s q[54]; // decomposed MY

OPENQASM 2.0;
include "qelib1.inc";

qreg q[64];
creg rec[49];

reset q[1];
reset q[3];
reset q[5];
reset q[7];
reset q[9];
reset q[12];
reset q[14];
reset q[16];
reset q[18];
reset q[20];
reset q[23];
reset q[25];
reset q[27];
reset q[29];
reset q[31];
reset q[34];
reset q[36];
reset q[38];
reset q[40];
reset q[42];
reset q[45];
reset q[47];
reset q[49];
reset q[51];
reset q[53];
reset q[2];
reset q[6];
reset q[13];
reset q[15];
reset q[17];
reset q[19];
reset q[21];
reset q[22];
reset q[24];
reset q[26];
reset q[28];
reset q[30];
reset q[35];
reset q[37];
reset q[39];
reset q[41];
reset q[43];
reset q[44];
reset q[46];
reset q[48];
reset q[50];
reset q[52];
reset q[59];
reset q[63];
barrier q;

h q[2];
h q[6];
h q[15];
h q[19];
h q[24];
h q[28];
h q[37];
h q[41];
h q[46];
h q[50];
h q[59];
h q[63];
barrier q;

cx q[2], q[3];
cx q[24], q[25];
cx q[46], q[47];
cx q[15], q[16];
cx q[37], q[38];
cx q[6], q[7];
cx q[28], q[29];
cx q[50], q[51];
cx q[19], q[20];
cx q[41], q[42];
cx q[23], q[22];
cx q[45], q[44];
cx q[14], q[13];
cx q[36], q[35];
cx q[27], q[26];
cx q[49], q[48];
cx q[18], q[17];
cx q[40], q[39];
cx q[31], q[30];
cx q[53], q[52];
barrier q;

cx q[2], q[1];
cx q[24], q[23];
cx q[46], q[45];
cx q[15], q[14];
cx q[37], q[36];
cx q[6], q[5];
cx q[28], q[27];
cx q[50], q[49];
cx q[19], q[18];
cx q[41], q[40];
cx q[12], q[22];
cx q[34], q[44];
cx q[3], q[13];
cx q[25], q[35];
cx q[16], q[26];
cx q[38], q[48];
cx q[7], q[17];
cx q[29], q[39];
cx q[20], q[30];
cx q[42], q[52];
barrier q;

cx q[24], q[14];
cx q[46], q[36];
cx q[15], q[5];
cx q[37], q[27];
cx q[59], q[49];
cx q[28], q[18];
cx q[50], q[40];
cx q[19], q[9];
cx q[41], q[31];
cx q[63], q[53];
cx q[12], q[13];
cx q[34], q[35];
cx q[25], q[26];
cx q[47], q[48];
cx q[16], q[17];
cx q[38], q[39];
cx q[29], q[30];
cx q[51], q[52];
cx q[20], q[21];
cx q[42], q[43];
barrier q;

cx q[24], q[12];
cx q[46], q[34];
cx q[15], q[3];
cx q[37], q[25];
cx q[59], q[47];
cx q[28], q[16];
cx q[50], q[38];
cx q[19], q[7];
cx q[41], q[29];
cx q[63], q[51];
cx q[1], q[13];
cx q[23], q[35];
cx q[14], q[26];
cx q[36], q[48];
cx q[5], q[17];
cx q[27], q[39];
cx q[18], q[30];
cx q[40], q[52];
cx q[9], q[21];
cx q[31], q[43];
barrier q;

h q[2];
h q[6];
h q[15];
h q[19];
h q[24];
h q[28];
h q[37];
h q[41];
h q[46];
h q[50];
h q[59];
h q[63];
barrier q;

measure q[2] -> rec[0]; reset q[2]; // decomposed MR
measure q[6] -> rec[1]; reset q[6]; // decomposed MR
measure q[13] -> rec[2]; reset q[13]; // decomposed MR
measure q[15] -> rec[3]; reset q[15]; // decomposed MR
measure q[17] -> rec[4]; reset q[17]; // decomposed MR
measure q[19] -> rec[5]; reset q[19]; // decomposed MR
measure q[21] -> rec[6]; reset q[21]; // decomposed MR
measure q[22] -> rec[7]; reset q[22]; // decomposed MR
measure q[24] -> rec[8]; reset q[24]; // decomposed MR
measure q[26] -> rec[9]; reset q[26]; // decomposed MR
measure q[28] -> rec[10]; reset q[28]; // decomposed MR
measure q[30] -> rec[11]; reset q[30]; // decomposed MR
measure q[35] -> rec[12]; reset q[35]; // decomposed MR
measure q[37] -> rec[13]; reset q[37]; // decomposed MR
measure q[39] -> rec[14]; reset q[39]; // decomposed MR
measure q[41] -> rec[15]; reset q[41]; // decomposed MR
measure q[43] -> rec[16]; reset q[43]; // decomposed MR
measure q[44] -> rec[17]; reset q[44]; // decomposed MR
measure q[46] -> rec[18]; reset q[46]; // decomposed MR
measure q[48] -> rec[19]; reset q[48]; // decomposed MR
measure q[50] -> rec[20]; reset q[50]; // decomposed MR
measure q[52] -> rec[21]; reset q[52]; // decomposed MR
measure q[59] -> rec[22]; reset q[59]; // decomposed MR
measure q[63] -> rec[23]; reset q[63]; // decomposed MR
measure q[1] -> rec[24];
measure q[3] -> rec[25];
measure q[5] -> rec[26];
measure q[7] -> rec[27];
measure q[9] -> rec[28];
measure q[12] -> rec[29];
measure q[14] -> rec[30];
measure q[16] -> rec[31];
measure q[18] -> rec[32];
measure q[20] -> rec[33];
measure q[23] -> rec[34];
measure q[25] -> rec[35];
measure q[27] -> rec[36];
measure q[29] -> rec[37];
measure q[31] -> rec[38];
measure q[34] -> rec[39];
measure q[36] -> rec[40];
measure q[38] -> rec[41];
measure q[40] -> rec[42];
measure q[42] -> rec[43];
measure q[45] -> rec[44];
measure q[47] -> rec[45];
measure q[49] -> rec[46];
measure q[51] -> rec[47];
measure q[53] -> rec[48];

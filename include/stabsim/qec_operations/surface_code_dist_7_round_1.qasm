OPENQASM 2.0;
include "qelib1.inc";

qreg q[118];
creg rec[97];

reset q[1];
reset q[3];
reset q[5];
reset q[7];
reset q[9];
reset q[11];
reset q[13];
reset q[16];
reset q[18];
reset q[20];
reset q[22];
reset q[24];
reset q[26];
reset q[28];
reset q[31];
reset q[33];
reset q[35];
reset q[37];
reset q[39];
reset q[41];
reset q[43];
reset q[46];
reset q[48];
reset q[50];
reset q[52];
reset q[54];
reset q[56];
reset q[58];
reset q[61];
reset q[63];
reset q[65];
reset q[67];
reset q[69];
reset q[71];
reset q[73];
reset q[76];
reset q[78];
reset q[80];
reset q[82];
reset q[84];
reset q[86];
reset q[88];
reset q[91];
reset q[93];
reset q[95];
reset q[97];
reset q[99];
reset q[101];
reset q[103];
reset q[2];
reset q[6];
reset q[10];
reset q[17];
reset q[19];
reset q[21];
reset q[23];
reset q[25];
reset q[27];
reset q[29];
reset q[30];
reset q[32];
reset q[34];
reset q[36];
reset q[38];
reset q[40];
reset q[42];
reset q[47];
reset q[49];
reset q[51];
reset q[53];
reset q[55];
reset q[57];
reset q[59];
reset q[60];
reset q[62];
reset q[64];
reset q[66];
reset q[68];
reset q[70];
reset q[72];
reset q[77];
reset q[79];
reset q[81];
reset q[83];
reset q[85];
reset q[87];
reset q[89];
reset q[90];
reset q[92];
reset q[94];
reset q[96];
reset q[98];
reset q[100];
reset q[102];
reset q[109];
reset q[113];
reset q[117];
barrier q;

h q[2];
h q[6];
h q[10];
h q[19];
h q[23];
h q[27];
h q[32];
h q[36];
h q[40];
h q[49];
h q[53];
h q[57];
h q[62];
h q[66];
h q[70];
h q[79];
h q[83];
h q[87];
h q[92];
h q[96];
h q[100];
h q[109];
h q[113];
h q[117];
barrier q;

cx q[2], q[3];
cx q[32], q[33];
cx q[62], q[63];
cx q[92], q[93];
cx q[19], q[20];
cx q[49], q[50];
cx q[79], q[80];
cx q[6], q[7];
cx q[36], q[37];
cx q[66], q[67];
cx q[96], q[97];
cx q[23], q[24];
cx q[53], q[54];
cx q[83], q[84];
cx q[10], q[11];
cx q[40], q[41];
cx q[70], q[71];
cx q[100], q[101];
cx q[27], q[28];
cx q[57], q[58];
cx q[87], q[88];
cx q[31], q[30];
cx q[61], q[60];
cx q[91], q[90];
cx q[18], q[17];
cx q[48], q[47];
cx q[78], q[77];
cx q[35], q[34];
cx q[65], q[64];
cx q[95], q[94];
cx q[22], q[21];
cx q[52], q[51];
cx q[82], q[81];
cx q[39], q[38];
cx q[69], q[68];
cx q[99], q[98];
cx q[26], q[25];
cx q[56], q[55];
cx q[86], q[85];
cx q[43], q[42];
cx q[73], q[72];
cx q[103], q[102];
barrier q;

cx q[2], q[1];
cx q[32], q[31];
cx q[62], q[61];
cx q[92], q[91];
cx q[19], q[18];
cx q[49], q[48];
cx q[79], q[78];
cx q[6], q[5];
cx q[36], q[35];
cx q[66], q[65];
cx q[96], q[95];
cx q[23], q[22];
cx q[53], q[52];
cx q[83], q[82];
cx q[10], q[9];
cx q[40], q[39];
cx q[70], q[69];
cx q[100], q[99];
cx q[27], q[26];
cx q[57], q[56];
cx q[87], q[86];
cx q[16], q[30];
cx q[46], q[60];
cx q[76], q[90];
cx q[3], q[17];
cx q[33], q[47];
cx q[63], q[77];
cx q[20], q[34];
cx q[50], q[64];
cx q[80], q[94];
cx q[7], q[21];
cx q[37], q[51];
cx q[67], q[81];
cx q[24], q[38];
cx q[54], q[68];
cx q[84], q[98];
cx q[11], q[25];
cx q[41], q[55];
cx q[71], q[85];
cx q[28], q[42];
cx q[58], q[72];
cx q[88], q[102];
barrier q;

cx q[32], q[18];
cx q[62], q[48];
cx q[92], q[78];
cx q[19], q[5];
cx q[49], q[35];
cx q[79], q[65];
cx q[109], q[95];
cx q[36], q[22];
cx q[66], q[52];
cx q[96], q[82];
cx q[23], q[9];
cx q[53], q[39];
cx q[83], q[69];
cx q[113], q[99];
cx q[40], q[26];
cx q[70], q[56];
cx q[100], q[86];
cx q[27], q[13];
cx q[57], q[43];
cx q[87], q[73];
cx q[117], q[103];
cx q[16], q[17];
cx q[46], q[47];
cx q[76], q[77];
cx q[33], q[34];
cx q[63], q[64];
cx q[93], q[94];
cx q[20], q[21];
cx q[50], q[51];
cx q[80], q[81];
cx q[37], q[38];
cx q[67], q[68];
cx q[97], q[98];
cx q[24], q[25];
cx q[54], q[55];
cx q[84], q[85];
cx q[41], q[42];
cx q[71], q[72];
cx q[101], q[102];
cx q[28], q[29];
cx q[58], q[59];
cx q[88], q[89];
barrier q;

cx q[32], q[16];
cx q[62], q[46];
cx q[92], q[76];
cx q[19], q[3];
cx q[49], q[33];
cx q[79], q[63];
cx q[109], q[93];
cx q[36], q[20];
cx q[66], q[50];
cx q[96], q[80];
cx q[23], q[7];
cx q[53], q[37];
cx q[83], q[67];
cx q[113], q[97];
cx q[40], q[24];
cx q[70], q[54];
cx q[100], q[84];
cx q[27], q[11];
cx q[57], q[41];
cx q[87], q[71];
cx q[117], q[101];
cx q[1], q[17];
cx q[31], q[47];
cx q[61], q[77];
cx q[18], q[34];
cx q[48], q[64];
cx q[78], q[94];
cx q[5], q[21];
cx q[35], q[51];
cx q[65], q[81];
cx q[22], q[38];
cx q[52], q[68];
cx q[82], q[98];
cx q[9], q[25];
cx q[39], q[55];
cx q[69], q[85];
cx q[26], q[42];
cx q[56], q[72];
cx q[86], q[102];
cx q[13], q[29];
cx q[43], q[59];
cx q[73], q[89];
barrier q;

h q[2];
h q[6];
h q[10];
h q[19];
h q[23];
h q[27];
h q[32];
h q[36];
h q[40];
h q[49];
h q[53];
h q[57];
h q[62];
h q[66];
h q[70];
h q[79];
h q[83];
h q[87];
h q[92];
h q[96];
h q[100];
h q[109];
h q[113];
h q[117];
barrier q;

measure q[2] -> rec[0]; reset q[2]; // decomposed MR
measure q[6] -> rec[1]; reset q[6]; // decomposed MR
measure q[10] -> rec[2]; reset q[10]; // decomposed MR
measure q[17] -> rec[3]; reset q[17]; // decomposed MR
measure q[19] -> rec[4]; reset q[19]; // decomposed MR
measure q[21] -> rec[5]; reset q[21]; // decomposed MR
measure q[23] -> rec[6]; reset q[23]; // decomposed MR
measure q[25] -> rec[7]; reset q[25]; // decomposed MR
measure q[27] -> rec[8]; reset q[27]; // decomposed MR
measure q[29] -> rec[9]; reset q[29]; // decomposed MR
measure q[30] -> rec[10]; reset q[30]; // decomposed MR
measure q[32] -> rec[11]; reset q[32]; // decomposed MR
measure q[34] -> rec[12]; reset q[34]; // decomposed MR
measure q[36] -> rec[13]; reset q[36]; // decomposed MR
measure q[38] -> rec[14]; reset q[38]; // decomposed MR
measure q[40] -> rec[15]; reset q[40]; // decomposed MR
measure q[42] -> rec[16]; reset q[42]; // decomposed MR
measure q[47] -> rec[17]; reset q[47]; // decomposed MR
measure q[49] -> rec[18]; reset q[49]; // decomposed MR
measure q[51] -> rec[19]; reset q[51]; // decomposed MR
measure q[53] -> rec[20]; reset q[53]; // decomposed MR
measure q[55] -> rec[21]; reset q[55]; // decomposed MR
measure q[57] -> rec[22]; reset q[57]; // decomposed MR
measure q[59] -> rec[23]; reset q[59]; // decomposed MR
measure q[60] -> rec[24]; reset q[60]; // decomposed MR
measure q[62] -> rec[25]; reset q[62]; // decomposed MR
measure q[64] -> rec[26]; reset q[64]; // decomposed MR
measure q[66] -> rec[27]; reset q[66]; // decomposed MR
measure q[68] -> rec[28]; reset q[68]; // decomposed MR
measure q[70] -> rec[29]; reset q[70]; // decomposed MR
measure q[72] -> rec[30]; reset q[72]; // decomposed MR
measure q[77] -> rec[31]; reset q[77]; // decomposed MR
measure q[79] -> rec[32]; reset q[79]; // decomposed MR
measure q[81] -> rec[33]; reset q[81]; // decomposed MR
measure q[83] -> rec[34]; reset q[83]; // decomposed MR
measure q[85] -> rec[35]; reset q[85]; // decomposed MR
measure q[87] -> rec[36]; reset q[87]; // decomposed MR
measure q[89] -> rec[37]; reset q[89]; // decomposed MR
measure q[90] -> rec[38]; reset q[90]; // decomposed MR
measure q[92] -> rec[39]; reset q[92]; // decomposed MR
measure q[94] -> rec[40]; reset q[94]; // decomposed MR
measure q[96] -> rec[41]; reset q[96]; // decomposed MR
measure q[98] -> rec[42]; reset q[98]; // decomposed MR
measure q[100] -> rec[43]; reset q[100]; // decomposed MR
measure q[102] -> rec[44]; reset q[102]; // decomposed MR
measure q[109] -> rec[45]; reset q[109]; // decomposed MR
measure q[113] -> rec[46]; reset q[113]; // decomposed MR
measure q[117] -> rec[47]; reset q[117]; // decomposed MR
measure q[1] -> rec[48];
measure q[3] -> rec[49];
measure q[5] -> rec[50];
measure q[7] -> rec[51];
measure q[9] -> rec[52];
measure q[11] -> rec[53];
measure q[13] -> rec[54];
measure q[16] -> rec[55];
measure q[18] -> rec[56];
measure q[20] -> rec[57];
measure q[22] -> rec[58];
measure q[24] -> rec[59];
measure q[26] -> rec[60];
measure q[28] -> rec[61];
measure q[31] -> rec[62];
measure q[33] -> rec[63];
measure q[35] -> rec[64];
measure q[37] -> rec[65];
measure q[39] -> rec[66];
measure q[41] -> rec[67];
measure q[43] -> rec[68];
measure q[46] -> rec[69];
measure q[48] -> rec[70];
measure q[50] -> rec[71];
measure q[52] -> rec[72];
measure q[54] -> rec[73];
measure q[56] -> rec[74];
measure q[58] -> rec[75];
measure q[61] -> rec[76];
measure q[63] -> rec[77];
measure q[65] -> rec[78];
measure q[67] -> rec[79];
measure q[69] -> rec[80];
measure q[71] -> rec[81];
measure q[73] -> rec[82];
measure q[76] -> rec[83];
measure q[78] -> rec[84];
measure q[80] -> rec[85];
measure q[82] -> rec[86];
measure q[84] -> rec[87];
measure q[86] -> rec[88];
measure q[88] -> rec[89];
measure q[91] -> rec[90];
measure q[93] -> rec[91];
measure q[95] -> rec[92];
measure q[97] -> rec[93];
measure q[99] -> rec[94];
measure q[101] -> rec[95];
measure q[103] -> rec[96];

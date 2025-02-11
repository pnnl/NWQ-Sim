OPENQASM 2.0;
include "qelib1.inc";
gate cxyz q0 { U(pi/2, 0, pi/2) q0; }

qreg q[91];
creg rec[121];

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
reset q[55];
reset q[56];
reset q[57];
reset q[58];
reset q[59];
reset q[60];
reset q[61];
reset q[62];
reset q[63];
reset q[64];
reset q[65];
reset q[66];
reset q[67];
reset q[68];
reset q[69];
reset q[70];
reset q[71];
reset q[72];
reset q[73];
reset q[74];
reset q[75];
reset q[76];
reset q[77];
reset q[78];
reset q[79];
reset q[80];
reset q[81];
reset q[82];
reset q[83];
reset q[84];
reset q[85];
reset q[86];
reset q[87];
reset q[88];
reset q[89];
reset q[90];
barrier q;

cxyz q[0];
cxyz q[1];
cxyz q[3];
cxyz q[4];
cxyz q[6];
cxyz q[7];
cxyz q[9];
cxyz q[10];
cxyz q[12];
cxyz q[14];
cxyz q[15];
cxyz q[17];
cxyz q[18];
cxyz q[20];
cxyz q[21];
cxyz q[23];
cxyz q[24];
cxyz q[25];
cxyz q[27];
cxyz q[28];
cxyz q[30];
cxyz q[31];
cxyz q[33];
cxyz q[34];
cxyz q[36];
cxyz q[37];
cxyz q[39];
cxyz q[40];
cxyz q[42];
cxyz q[43];
cxyz q[45];
cxyz q[47];
cxyz q[48];
cxyz q[50];
cxyz q[51];
cxyz q[53];
cxyz q[54];
cxyz q[55];
cxyz q[57];
cxyz q[58];
cxyz q[60];
cxyz q[61];
cxyz q[63];
cxyz q[64];
cxyz q[66];
cxyz q[67];
cxyz q[69];
cxyz q[71];
cxyz q[72];
cxyz q[74];
cxyz q[75];
cxyz q[76];
cxyz q[78];
cxyz q[79];
cxyz q[81];
cxyz q[82];
cxyz q[84];
cxyz q[86];
cxyz q[87];
cxyz q[88];
cxyz q[90];
barrier q;

cx q[14], q[13];
cx q[3], q[2];
cx q[27], q[26];
cx q[47], q[46];
cx q[17], q[16];
cx q[39], q[38];
cx q[57], q[56];
cx q[71], q[70];
cx q[6], q[5];
cx q[30], q[29];
cx q[50], q[49];
cx q[66], q[65];
cx q[78], q[77];
cx q[86], q[85];
cx q[20], q[19];
cx q[42], q[41];
cx q[60], q[59];
cx q[74], q[73];
cx q[84], q[83];
cx q[9], q[8];
cx q[33], q[32];
cx q[53], q[52];
cx q[69], q[68];
cx q[23], q[22];
cx q[45], q[44];
cx q[12], q[11];
barrier q;

cx q[25], q[13];
cx q[15], q[2];
cx q[37], q[26];
cx q[55], q[46];
cx q[28], q[16];
cx q[48], q[38];
cx q[64], q[56];
cx q[76], q[70];
cx q[18], q[5];
cx q[40], q[29];
cx q[58], q[49];
cx q[72], q[65];
cx q[82], q[77];
cx q[88], q[85];
cx q[31], q[19];
cx q[51], q[41];
cx q[67], q[59];
cx q[79], q[73];
cx q[87], q[83];
cx q[21], q[8];
cx q[43], q[32];
cx q[61], q[52];
cx q[75], q[68];
cx q[34], q[22];
cx q[54], q[44];
cx q[24], q[11];
barrier q;

cx q[1], q[13];
cx q[15], q[26];
cx q[37], q[46];
cx q[4], q[16];
cx q[28], q[38];
cx q[48], q[56];
cx q[64], q[70];
cx q[18], q[29];
cx q[40], q[49];
cx q[58], q[65];
cx q[72], q[77];
cx q[82], q[85];
cx q[7], q[19];
cx q[31], q[41];
cx q[51], q[59];
cx q[67], q[73];
cx q[79], q[83];
cx q[87], q[89];
cx q[21], q[32];
cx q[43], q[52];
cx q[61], q[68];
cx q[75], q[80];
cx q[10], q[22];
cx q[34], q[44];
cx q[54], q[62];
cx q[24], q[35];
barrier q;

cx q[1], q[2];
cx q[25], q[26];
cx q[15], q[16];
cx q[37], q[38];
cx q[55], q[56];
cx q[4], q[5];
cx q[28], q[29];
cx q[48], q[49];
cx q[64], q[65];
cx q[76], q[77];
cx q[18], q[19];
cx q[40], q[41];
cx q[58], q[59];
cx q[72], q[73];
cx q[82], q[83];
cx q[88], q[89];
cx q[7], q[8];
cx q[31], q[32];
cx q[51], q[52];
cx q[67], q[68];
cx q[79], q[80];
cx q[21], q[22];
cx q[43], q[44];
cx q[61], q[62];
cx q[10], q[11];
cx q[34], q[35];
barrier q;

cx q[14], q[2];
cx q[36], q[26];
cx q[27], q[16];
cx q[47], q[38];
cx q[63], q[56];
cx q[17], q[5];
cx q[39], q[29];
cx q[57], q[49];
cx q[71], q[65];
cx q[81], q[77];
cx q[30], q[19];
cx q[50], q[41];
cx q[66], q[59];
cx q[78], q[73];
cx q[86], q[83];
cx q[90], q[89];
cx q[20], q[8];
cx q[42], q[32];
cx q[60], q[52];
cx q[74], q[68];
cx q[84], q[80];
cx q[33], q[22];
cx q[53], q[44];
cx q[69], q[62];
cx q[23], q[11];
cx q[45], q[35];
barrier q;

cx q[0], q[13];
cx q[14], q[26];
cx q[36], q[46];
cx q[3], q[16];
cx q[27], q[38];
cx q[47], q[56];
cx q[63], q[70];
cx q[17], q[29];
cx q[39], q[49];
cx q[57], q[65];
cx q[71], q[77];
cx q[81], q[85];
cx q[6], q[19];
cx q[30], q[41];
cx q[50], q[59];
cx q[66], q[73];
cx q[78], q[83];
cx q[86], q[89];
cx q[20], q[32];
cx q[42], q[52];
cx q[60], q[68];
cx q[74], q[80];
cx q[9], q[22];
cx q[33], q[44];
cx q[53], q[62];
cx q[23], q[35];
barrier q;

measure q[2] -> rec[0]; reset q[2]; // decomposed MR
measure q[5] -> rec[1]; reset q[5]; // decomposed MR
measure q[8] -> rec[2]; reset q[8]; // decomposed MR
measure q[11] -> rec[3]; reset q[11]; // decomposed MR
measure q[13] -> rec[4]; reset q[13]; // decomposed MR
measure q[16] -> rec[5]; reset q[16]; // decomposed MR
measure q[19] -> rec[6]; reset q[19]; // decomposed MR
measure q[22] -> rec[7]; reset q[22]; // decomposed MR
measure q[26] -> rec[8]; reset q[26]; // decomposed MR
measure q[29] -> rec[9]; reset q[29]; // decomposed MR
measure q[32] -> rec[10]; reset q[32]; // decomposed MR
measure q[35] -> rec[11]; reset q[35]; // decomposed MR
measure q[38] -> rec[12]; reset q[38]; // decomposed MR
measure q[41] -> rec[13]; reset q[41]; // decomposed MR
measure q[44] -> rec[14]; reset q[44]; // decomposed MR
measure q[46] -> rec[15]; reset q[46]; // decomposed MR
measure q[49] -> rec[16]; reset q[49]; // decomposed MR
measure q[52] -> rec[17]; reset q[52]; // decomposed MR
measure q[56] -> rec[18]; reset q[56]; // decomposed MR
measure q[59] -> rec[19]; reset q[59]; // decomposed MR
measure q[62] -> rec[20]; reset q[62]; // decomposed MR
measure q[65] -> rec[21]; reset q[65]; // decomposed MR
measure q[68] -> rec[22]; reset q[68]; // decomposed MR
measure q[70] -> rec[23]; reset q[70]; // decomposed MR
measure q[73] -> rec[24]; reset q[73]; // decomposed MR
measure q[77] -> rec[25]; reset q[77]; // decomposed MR
measure q[80] -> rec[26]; reset q[80]; // decomposed MR
measure q[83] -> rec[27]; reset q[83]; // decomposed MR
measure q[85] -> rec[28]; reset q[85]; // decomposed MR
measure q[89] -> rec[29]; reset q[89]; // decomposed MR
barrier q;

cxyz q[0];
cxyz q[1];
cxyz q[3];
cxyz q[4];
cxyz q[6];
cxyz q[7];
cxyz q[9];
cxyz q[10];
cxyz q[12];
cxyz q[14];
cxyz q[15];
cxyz q[17];
cxyz q[18];
cxyz q[20];
cxyz q[21];
cxyz q[23];
cxyz q[24];
cxyz q[25];
cxyz q[27];
cxyz q[28];
cxyz q[30];
cxyz q[31];
cxyz q[33];
cxyz q[34];
cxyz q[36];
cxyz q[37];
cxyz q[39];
cxyz q[40];
cxyz q[42];
cxyz q[43];
cxyz q[45];
cxyz q[47];
cxyz q[48];
cxyz q[50];
cxyz q[51];
cxyz q[53];
cxyz q[54];
cxyz q[55];
cxyz q[57];
cxyz q[58];
cxyz q[60];
cxyz q[61];
cxyz q[63];
cxyz q[64];
cxyz q[66];
cxyz q[67];
cxyz q[69];
cxyz q[71];
cxyz q[72];
cxyz q[74];
cxyz q[75];
cxyz q[76];
cxyz q[78];
cxyz q[79];
cxyz q[81];
cxyz q[82];
cxyz q[84];
cxyz q[86];
cxyz q[87];
cxyz q[88];
cxyz q[90];
barrier q;

cx q[14], q[13];
cx q[3], q[2];
cx q[27], q[26];
cx q[47], q[46];
cx q[17], q[16];
cx q[39], q[38];
cx q[57], q[56];
cx q[71], q[70];
cx q[6], q[5];
cx q[30], q[29];
cx q[50], q[49];
cx q[66], q[65];
cx q[78], q[77];
cx q[86], q[85];
cx q[20], q[19];
cx q[42], q[41];
cx q[60], q[59];
cx q[74], q[73];
cx q[84], q[83];
cx q[9], q[8];
cx q[33], q[32];
cx q[53], q[52];
cx q[69], q[68];
cx q[23], q[22];
cx q[45], q[44];
cx q[12], q[11];
barrier q;

cx q[25], q[13];
cx q[15], q[2];
cx q[37], q[26];
cx q[55], q[46];
cx q[28], q[16];
cx q[48], q[38];
cx q[64], q[56];
cx q[76], q[70];
cx q[18], q[5];
cx q[40], q[29];
cx q[58], q[49];
cx q[72], q[65];
cx q[82], q[77];
cx q[88], q[85];
cx q[31], q[19];
cx q[51], q[41];
cx q[67], q[59];
cx q[79], q[73];
cx q[87], q[83];
cx q[21], q[8];
cx q[43], q[32];
cx q[61], q[52];
cx q[75], q[68];
cx q[34], q[22];
cx q[54], q[44];
cx q[24], q[11];
barrier q;

cx q[1], q[13];
cx q[15], q[26];
cx q[37], q[46];
cx q[4], q[16];
cx q[28], q[38];
cx q[48], q[56];
cx q[64], q[70];
cx q[18], q[29];
cx q[40], q[49];
cx q[58], q[65];
cx q[72], q[77];
cx q[82], q[85];
cx q[7], q[19];
cx q[31], q[41];
cx q[51], q[59];
cx q[67], q[73];
cx q[79], q[83];
cx q[87], q[89];
cx q[21], q[32];
cx q[43], q[52];
cx q[61], q[68];
cx q[75], q[80];
cx q[10], q[22];
cx q[34], q[44];
cx q[54], q[62];
cx q[24], q[35];
barrier q;

cx q[1], q[2];
cx q[25], q[26];
cx q[15], q[16];
cx q[37], q[38];
cx q[55], q[56];
cx q[4], q[5];
cx q[28], q[29];
cx q[48], q[49];
cx q[64], q[65];
cx q[76], q[77];
cx q[18], q[19];
cx q[40], q[41];
cx q[58], q[59];
cx q[72], q[73];
cx q[82], q[83];
cx q[88], q[89];
cx q[7], q[8];
cx q[31], q[32];
cx q[51], q[52];
cx q[67], q[68];
cx q[79], q[80];
cx q[21], q[22];
cx q[43], q[44];
cx q[61], q[62];
cx q[10], q[11];
cx q[34], q[35];
barrier q;

cx q[14], q[2];
cx q[36], q[26];
cx q[27], q[16];
cx q[47], q[38];
cx q[63], q[56];
cx q[17], q[5];
cx q[39], q[29];
cx q[57], q[49];
cx q[71], q[65];
cx q[81], q[77];
cx q[30], q[19];
cx q[50], q[41];
cx q[66], q[59];
cx q[78], q[73];
cx q[86], q[83];
cx q[90], q[89];
cx q[20], q[8];
cx q[42], q[32];
cx q[60], q[52];
cx q[74], q[68];
cx q[84], q[80];
cx q[33], q[22];
cx q[53], q[44];
cx q[69], q[62];
cx q[23], q[11];
cx q[45], q[35];
barrier q;

cx q[0], q[13];
cx q[14], q[26];
cx q[36], q[46];
cx q[3], q[16];
cx q[27], q[38];
cx q[47], q[56];
cx q[63], q[70];
cx q[17], q[29];
cx q[39], q[49];
cx q[57], q[65];
cx q[71], q[77];
cx q[81], q[85];
cx q[6], q[19];
cx q[30], q[41];
cx q[50], q[59];
cx q[66], q[73];
cx q[78], q[83];
cx q[86], q[89];
cx q[20], q[32];
cx q[42], q[52];
cx q[60], q[68];
cx q[74], q[80];
cx q[9], q[22];
cx q[33], q[44];
cx q[53], q[62];
cx q[23], q[35];
barrier q;

measure q[2] -> rec[30]; reset q[2]; // decomposed MR
measure q[5] -> rec[31]; reset q[5]; // decomposed MR
measure q[8] -> rec[32]; reset q[8]; // decomposed MR
measure q[11] -> rec[33]; reset q[11]; // decomposed MR
measure q[13] -> rec[34]; reset q[13]; // decomposed MR
measure q[16] -> rec[35]; reset q[16]; // decomposed MR
measure q[19] -> rec[36]; reset q[19]; // decomposed MR
measure q[22] -> rec[37]; reset q[22]; // decomposed MR
measure q[26] -> rec[38]; reset q[26]; // decomposed MR
measure q[29] -> rec[39]; reset q[29]; // decomposed MR
measure q[32] -> rec[40]; reset q[32]; // decomposed MR
measure q[35] -> rec[41]; reset q[35]; // decomposed MR
measure q[38] -> rec[42]; reset q[38]; // decomposed MR
measure q[41] -> rec[43]; reset q[41]; // decomposed MR
measure q[44] -> rec[44]; reset q[44]; // decomposed MR
measure q[46] -> rec[45]; reset q[46]; // decomposed MR
measure q[49] -> rec[46]; reset q[49]; // decomposed MR
measure q[52] -> rec[47]; reset q[52]; // decomposed MR
measure q[56] -> rec[48]; reset q[56]; // decomposed MR
measure q[59] -> rec[49]; reset q[59]; // decomposed MR
measure q[62] -> rec[50]; reset q[62]; // decomposed MR
measure q[65] -> rec[51]; reset q[65]; // decomposed MR
measure q[68] -> rec[52]; reset q[68]; // decomposed MR
measure q[70] -> rec[53]; reset q[70]; // decomposed MR
measure q[73] -> rec[54]; reset q[73]; // decomposed MR
measure q[77] -> rec[55]; reset q[77]; // decomposed MR
measure q[80] -> rec[56]; reset q[80]; // decomposed MR
measure q[83] -> rec[57]; reset q[83]; // decomposed MR
measure q[85] -> rec[58]; reset q[85]; // decomposed MR
measure q[89] -> rec[59]; reset q[89]; // decomposed MR
s q[0]; s q[0]; s q[0]; h q[0]; measure q[0] -> rec[60]; h q[0]; s q[0]; // decomposed MY
s q[1]; s q[1]; s q[1]; h q[1]; measure q[1] -> rec[61]; h q[1]; s q[1]; // decomposed MY
s q[3]; s q[3]; s q[3]; h q[3]; measure q[3] -> rec[62]; h q[3]; s q[3]; // decomposed MY
s q[4]; s q[4]; s q[4]; h q[4]; measure q[4] -> rec[63]; h q[4]; s q[4]; // decomposed MY
s q[6]; s q[6]; s q[6]; h q[6]; measure q[6] -> rec[64]; h q[6]; s q[6]; // decomposed MY
s q[7]; s q[7]; s q[7]; h q[7]; measure q[7] -> rec[65]; h q[7]; s q[7]; // decomposed MY
s q[9]; s q[9]; s q[9]; h q[9]; measure q[9] -> rec[66]; h q[9]; s q[9]; // decomposed MY
s q[10]; s q[10]; s q[10]; h q[10]; measure q[10] -> rec[67]; h q[10]; s q[10]; // decomposed MY
s q[12]; s q[12]; s q[12]; h q[12]; measure q[12] -> rec[68]; h q[12]; s q[12]; // decomposed MY
s q[14]; s q[14]; s q[14]; h q[14]; measure q[14] -> rec[69]; h q[14]; s q[14]; // decomposed MY
s q[15]; s q[15]; s q[15]; h q[15]; measure q[15] -> rec[70]; h q[15]; s q[15]; // decomposed MY
s q[17]; s q[17]; s q[17]; h q[17]; measure q[17] -> rec[71]; h q[17]; s q[17]; // decomposed MY
s q[18]; s q[18]; s q[18]; h q[18]; measure q[18] -> rec[72]; h q[18]; s q[18]; // decomposed MY
s q[20]; s q[20]; s q[20]; h q[20]; measure q[20] -> rec[73]; h q[20]; s q[20]; // decomposed MY
s q[21]; s q[21]; s q[21]; h q[21]; measure q[21] -> rec[74]; h q[21]; s q[21]; // decomposed MY
s q[23]; s q[23]; s q[23]; h q[23]; measure q[23] -> rec[75]; h q[23]; s q[23]; // decomposed MY
s q[24]; s q[24]; s q[24]; h q[24]; measure q[24] -> rec[76]; h q[24]; s q[24]; // decomposed MY
s q[25]; s q[25]; s q[25]; h q[25]; measure q[25] -> rec[77]; h q[25]; s q[25]; // decomposed MY
s q[27]; s q[27]; s q[27]; h q[27]; measure q[27] -> rec[78]; h q[27]; s q[27]; // decomposed MY
s q[28]; s q[28]; s q[28]; h q[28]; measure q[28] -> rec[79]; h q[28]; s q[28]; // decomposed MY
s q[30]; s q[30]; s q[30]; h q[30]; measure q[30] -> rec[80]; h q[30]; s q[30]; // decomposed MY
s q[31]; s q[31]; s q[31]; h q[31]; measure q[31] -> rec[81]; h q[31]; s q[31]; // decomposed MY
s q[33]; s q[33]; s q[33]; h q[33]; measure q[33] -> rec[82]; h q[33]; s q[33]; // decomposed MY
s q[34]; s q[34]; s q[34]; h q[34]; measure q[34] -> rec[83]; h q[34]; s q[34]; // decomposed MY
s q[36]; s q[36]; s q[36]; h q[36]; measure q[36] -> rec[84]; h q[36]; s q[36]; // decomposed MY
s q[37]; s q[37]; s q[37]; h q[37]; measure q[37] -> rec[85]; h q[37]; s q[37]; // decomposed MY
s q[39]; s q[39]; s q[39]; h q[39]; measure q[39] -> rec[86]; h q[39]; s q[39]; // decomposed MY
s q[40]; s q[40]; s q[40]; h q[40]; measure q[40] -> rec[87]; h q[40]; s q[40]; // decomposed MY
s q[42]; s q[42]; s q[42]; h q[42]; measure q[42] -> rec[88]; h q[42]; s q[42]; // decomposed MY
s q[43]; s q[43]; s q[43]; h q[43]; measure q[43] -> rec[89]; h q[43]; s q[43]; // decomposed MY
s q[45]; s q[45]; s q[45]; h q[45]; measure q[45] -> rec[90]; h q[45]; s q[45]; // decomposed MY
s q[47]; s q[47]; s q[47]; h q[47]; measure q[47] -> rec[91]; h q[47]; s q[47]; // decomposed MY
s q[48]; s q[48]; s q[48]; h q[48]; measure q[48] -> rec[92]; h q[48]; s q[48]; // decomposed MY
s q[50]; s q[50]; s q[50]; h q[50]; measure q[50] -> rec[93]; h q[50]; s q[50]; // decomposed MY
s q[51]; s q[51]; s q[51]; h q[51]; measure q[51] -> rec[94]; h q[51]; s q[51]; // decomposed MY
s q[53]; s q[53]; s q[53]; h q[53]; measure q[53] -> rec[95]; h q[53]; s q[53]; // decomposed MY
s q[54]; s q[54]; s q[54]; h q[54]; measure q[54] -> rec[96]; h q[54]; s q[54]; // decomposed MY
s q[55]; s q[55]; s q[55]; h q[55]; measure q[55] -> rec[97]; h q[55]; s q[55]; // decomposed MY
s q[57]; s q[57]; s q[57]; h q[57]; measure q[57] -> rec[98]; h q[57]; s q[57]; // decomposed MY
s q[58]; s q[58]; s q[58]; h q[58]; measure q[58] -> rec[99]; h q[58]; s q[58]; // decomposed MY
s q[60]; s q[60]; s q[60]; h q[60]; measure q[60] -> rec[100]; h q[60]; s q[60]; // decomposed MY
s q[61]; s q[61]; s q[61]; h q[61]; measure q[61] -> rec[101]; h q[61]; s q[61]; // decomposed MY
s q[63]; s q[63]; s q[63]; h q[63]; measure q[63] -> rec[102]; h q[63]; s q[63]; // decomposed MY
s q[64]; s q[64]; s q[64]; h q[64]; measure q[64] -> rec[103]; h q[64]; s q[64]; // decomposed MY
s q[66]; s q[66]; s q[66]; h q[66]; measure q[66] -> rec[104]; h q[66]; s q[66]; // decomposed MY
s q[67]; s q[67]; s q[67]; h q[67]; measure q[67] -> rec[105]; h q[67]; s q[67]; // decomposed MY
s q[69]; s q[69]; s q[69]; h q[69]; measure q[69] -> rec[106]; h q[69]; s q[69]; // decomposed MY
s q[71]; s q[71]; s q[71]; h q[71]; measure q[71] -> rec[107]; h q[71]; s q[71]; // decomposed MY
s q[72]; s q[72]; s q[72]; h q[72]; measure q[72] -> rec[108]; h q[72]; s q[72]; // decomposed MY
s q[74]; s q[74]; s q[74]; h q[74]; measure q[74] -> rec[109]; h q[74]; s q[74]; // decomposed MY
s q[75]; s q[75]; s q[75]; h q[75]; measure q[75] -> rec[110]; h q[75]; s q[75]; // decomposed MY
s q[76]; s q[76]; s q[76]; h q[76]; measure q[76] -> rec[111]; h q[76]; s q[76]; // decomposed MY
s q[78]; s q[78]; s q[78]; h q[78]; measure q[78] -> rec[112]; h q[78]; s q[78]; // decomposed MY
s q[79]; s q[79]; s q[79]; h q[79]; measure q[79] -> rec[113]; h q[79]; s q[79]; // decomposed MY
s q[81]; s q[81]; s q[81]; h q[81]; measure q[81] -> rec[114]; h q[81]; s q[81]; // decomposed MY
s q[82]; s q[82]; s q[82]; h q[82]; measure q[82] -> rec[115]; h q[82]; s q[82]; // decomposed MY
s q[84]; s q[84]; s q[84]; h q[84]; measure q[84] -> rec[116]; h q[84]; s q[84]; // decomposed MY
s q[86]; s q[86]; s q[86]; h q[86]; measure q[86] -> rec[117]; h q[86]; s q[86]; // decomposed MY
s q[87]; s q[87]; s q[87]; h q[87]; measure q[87] -> rec[118]; h q[87]; s q[87]; // decomposed MY
s q[88]; s q[88]; s q[88]; h q[88]; measure q[88] -> rec[119]; h q[88]; s q[88]; // decomposed MY
s q[90]; s q[90]; s q[90]; h q[90]; measure q[90] -> rec[120]; h q[90]; s q[90]; // decomposed MY

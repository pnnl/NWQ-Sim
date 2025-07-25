OPENQASM 2.0;
include "qelib1.inc";
gate cxyz q0 { U(pi/2, 0, pi/2) q0; }

qreg q[136];
creg rec[181];

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
reset q[91];
reset q[92];
reset q[93];
reset q[94];
reset q[95];
reset q[96];
reset q[97];
reset q[98];
reset q[99];
reset q[100];
reset q[101];
reset q[102];
reset q[103];
reset q[104];
reset q[105];
reset q[106];
reset q[107];
reset q[108];
reset q[109];
reset q[110];
reset q[111];
reset q[112];
reset q[113];
reset q[114];
reset q[115];
reset q[116];
reset q[117];
reset q[118];
reset q[119];
reset q[120];
reset q[121];
reset q[122];
reset q[123];
reset q[124];
reset q[125];
reset q[126];
reset q[127];
reset q[128];
reset q[129];
reset q[130];
reset q[131];
reset q[132];
reset q[133];
reset q[134];
reset q[135];
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
cxyz q[13];
cxyz q[15];
cxyz q[17];
cxyz q[18];
cxyz q[20];
cxyz q[21];
cxyz q[23];
cxyz q[24];
cxyz q[26];
cxyz q[27];
cxyz q[29];
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
cxyz q[46];
cxyz q[48];
cxyz q[49];
cxyz q[51];
cxyz q[52];
cxyz q[54];
cxyz q[55];
cxyz q[57];
cxyz q[59];
cxyz q[60];
cxyz q[62];
cxyz q[63];
cxyz q[65];
cxyz q[66];
cxyz q[68];
cxyz q[69];
cxyz q[70];
cxyz q[72];
cxyz q[73];
cxyz q[75];
cxyz q[76];
cxyz q[78];
cxyz q[79];
cxyz q[81];
cxyz q[82];
cxyz q[84];
cxyz q[85];
cxyz q[87];
cxyz q[88];
cxyz q[90];
cxyz q[92];
cxyz q[93];
cxyz q[95];
cxyz q[96];
cxyz q[98];
cxyz q[99];
cxyz q[100];
cxyz q[102];
cxyz q[103];
cxyz q[105];
cxyz q[106];
cxyz q[108];
cxyz q[109];
cxyz q[111];
cxyz q[112];
cxyz q[114];
cxyz q[116];
cxyz q[117];
cxyz q[119];
cxyz q[120];
cxyz q[121];
cxyz q[123];
cxyz q[124];
cxyz q[126];
cxyz q[127];
cxyz q[129];
cxyz q[131];
cxyz q[132];
cxyz q[133];
cxyz q[135];
barrier q;

cx q[17], q[16];
cx q[3], q[2];
cx q[33], q[32];
cx q[59], q[58];
cx q[20], q[19];
cx q[48], q[47];
cx q[72], q[71];
cx q[92], q[91];
cx q[6], q[5];
cx q[36], q[35];
cx q[62], q[61];
cx q[84], q[83];
cx q[102], q[101];
cx q[116], q[115];
cx q[23], q[22];
cx q[51], q[50];
cx q[75], q[74];
cx q[95], q[94];
cx q[111], q[110];
cx q[123], q[122];
cx q[131], q[130];
cx q[9], q[8];
cx q[39], q[38];
cx q[65], q[64];
cx q[87], q[86];
cx q[105], q[104];
cx q[119], q[118];
cx q[129], q[128];
cx q[26], q[25];
cx q[54], q[53];
cx q[78], q[77];
cx q[98], q[97];
cx q[114], q[113];
cx q[12], q[11];
cx q[42], q[41];
cx q[68], q[67];
cx q[90], q[89];
cx q[29], q[28];
cx q[57], q[56];
cx q[15], q[14];
barrier q;

cx q[31], q[16];
cx q[18], q[2];
cx q[46], q[32];
cx q[70], q[58];
cx q[34], q[19];
cx q[60], q[47];
cx q[82], q[71];
cx q[100], q[91];
cx q[21], q[5];
cx q[49], q[35];
cx q[73], q[61];
cx q[93], q[83];
cx q[109], q[101];
cx q[121], q[115];
cx q[37], q[22];
cx q[63], q[50];
cx q[85], q[74];
cx q[103], q[94];
cx q[117], q[110];
cx q[127], q[122];
cx q[133], q[130];
cx q[24], q[8];
cx q[52], q[38];
cx q[76], q[64];
cx q[96], q[86];
cx q[112], q[104];
cx q[124], q[118];
cx q[132], q[128];
cx q[40], q[25];
cx q[66], q[53];
cx q[88], q[77];
cx q[106], q[97];
cx q[120], q[113];
cx q[27], q[11];
cx q[55], q[41];
cx q[79], q[67];
cx q[99], q[89];
cx q[43], q[28];
cx q[69], q[56];
cx q[30], q[14];
barrier q;

cx q[1], q[16];
cx q[18], q[32];
cx q[46], q[58];
cx q[4], q[19];
cx q[34], q[47];
cx q[60], q[71];
cx q[82], q[91];
cx q[21], q[35];
cx q[49], q[61];
cx q[73], q[83];
cx q[93], q[101];
cx q[109], q[115];
cx q[7], q[22];
cx q[37], q[50];
cx q[63], q[74];
cx q[85], q[94];
cx q[103], q[110];
cx q[117], q[122];
cx q[127], q[130];
cx q[24], q[38];
cx q[52], q[64];
cx q[76], q[86];
cx q[96], q[104];
cx q[112], q[118];
cx q[124], q[128];
cx q[132], q[134];
cx q[10], q[25];
cx q[40], q[53];
cx q[66], q[77];
cx q[88], q[97];
cx q[106], q[113];
cx q[120], q[125];
cx q[27], q[41];
cx q[55], q[67];
cx q[79], q[89];
cx q[99], q[107];
cx q[13], q[28];
cx q[43], q[56];
cx q[69], q[80];
cx q[30], q[44];
barrier q;

cx q[1], q[2];
cx q[31], q[32];
cx q[18], q[19];
cx q[46], q[47];
cx q[70], q[71];
cx q[4], q[5];
cx q[34], q[35];
cx q[60], q[61];
cx q[82], q[83];
cx q[100], q[101];
cx q[21], q[22];
cx q[49], q[50];
cx q[73], q[74];
cx q[93], q[94];
cx q[109], q[110];
cx q[121], q[122];
cx q[7], q[8];
cx q[37], q[38];
cx q[63], q[64];
cx q[85], q[86];
cx q[103], q[104];
cx q[117], q[118];
cx q[127], q[128];
cx q[133], q[134];
cx q[24], q[25];
cx q[52], q[53];
cx q[76], q[77];
cx q[96], q[97];
cx q[112], q[113];
cx q[124], q[125];
cx q[10], q[11];
cx q[40], q[41];
cx q[66], q[67];
cx q[88], q[89];
cx q[106], q[107];
cx q[27], q[28];
cx q[55], q[56];
cx q[79], q[80];
cx q[13], q[14];
cx q[43], q[44];
barrier q;

cx q[17], q[2];
cx q[45], q[32];
cx q[33], q[19];
cx q[59], q[47];
cx q[81], q[71];
cx q[20], q[5];
cx q[48], q[35];
cx q[72], q[61];
cx q[92], q[83];
cx q[108], q[101];
cx q[36], q[22];
cx q[62], q[50];
cx q[84], q[74];
cx q[102], q[94];
cx q[116], q[110];
cx q[126], q[122];
cx q[23], q[8];
cx q[51], q[38];
cx q[75], q[64];
cx q[95], q[86];
cx q[111], q[104];
cx q[123], q[118];
cx q[131], q[128];
cx q[135], q[134];
cx q[39], q[25];
cx q[65], q[53];
cx q[87], q[77];
cx q[105], q[97];
cx q[119], q[113];
cx q[129], q[125];
cx q[26], q[11];
cx q[54], q[41];
cx q[78], q[67];
cx q[98], q[89];
cx q[114], q[107];
cx q[42], q[28];
cx q[68], q[56];
cx q[90], q[80];
cx q[29], q[14];
cx q[57], q[44];
barrier q;

cx q[0], q[16];
cx q[17], q[32];
cx q[45], q[58];
cx q[3], q[19];
cx q[33], q[47];
cx q[59], q[71];
cx q[81], q[91];
cx q[20], q[35];
cx q[48], q[61];
cx q[72], q[83];
cx q[92], q[101];
cx q[108], q[115];
cx q[6], q[22];
cx q[36], q[50];
cx q[62], q[74];
cx q[84], q[94];
cx q[102], q[110];
cx q[116], q[122];
cx q[126], q[130];
cx q[23], q[38];
cx q[51], q[64];
cx q[75], q[86];
cx q[95], q[104];
cx q[111], q[118];
cx q[123], q[128];
cx q[131], q[134];
cx q[9], q[25];
cx q[39], q[53];
cx q[65], q[77];
cx q[87], q[97];
cx q[105], q[113];
cx q[119], q[125];
cx q[26], q[41];
cx q[54], q[67];
cx q[78], q[89];
cx q[98], q[107];
cx q[12], q[28];
cx q[42], q[56];
cx q[68], q[80];
cx q[29], q[44];
barrier q;

measure q[2] -> rec[0]; reset q[2]; // decomposed MR
measure q[5] -> rec[1]; reset q[5]; // decomposed MR
measure q[8] -> rec[2]; reset q[8]; // decomposed MR
measure q[11] -> rec[3]; reset q[11]; // decomposed MR
measure q[14] -> rec[4]; reset q[14]; // decomposed MR
measure q[16] -> rec[5]; reset q[16]; // decomposed MR
measure q[19] -> rec[6]; reset q[19]; // decomposed MR
measure q[22] -> rec[7]; reset q[22]; // decomposed MR
measure q[25] -> rec[8]; reset q[25]; // decomposed MR
measure q[28] -> rec[9]; reset q[28]; // decomposed MR
measure q[32] -> rec[10]; reset q[32]; // decomposed MR
measure q[35] -> rec[11]; reset q[35]; // decomposed MR
measure q[38] -> rec[12]; reset q[38]; // decomposed MR
measure q[41] -> rec[13]; reset q[41]; // decomposed MR
measure q[44] -> rec[14]; reset q[44]; // decomposed MR
measure q[47] -> rec[15]; reset q[47]; // decomposed MR
measure q[50] -> rec[16]; reset q[50]; // decomposed MR
measure q[53] -> rec[17]; reset q[53]; // decomposed MR
measure q[56] -> rec[18]; reset q[56]; // decomposed MR
measure q[58] -> rec[19]; reset q[58]; // decomposed MR
measure q[61] -> rec[20]; reset q[61]; // decomposed MR
measure q[64] -> rec[21]; reset q[64]; // decomposed MR
measure q[67] -> rec[22]; reset q[67]; // decomposed MR
measure q[71] -> rec[23]; reset q[71]; // decomposed MR
measure q[74] -> rec[24]; reset q[74]; // decomposed MR
measure q[77] -> rec[25]; reset q[77]; // decomposed MR
measure q[80] -> rec[26]; reset q[80]; // decomposed MR
measure q[83] -> rec[27]; reset q[83]; // decomposed MR
measure q[86] -> rec[28]; reset q[86]; // decomposed MR
measure q[89] -> rec[29]; reset q[89]; // decomposed MR
measure q[91] -> rec[30]; reset q[91]; // decomposed MR
measure q[94] -> rec[31]; reset q[94]; // decomposed MR
measure q[97] -> rec[32]; reset q[97]; // decomposed MR
measure q[101] -> rec[33]; reset q[101]; // decomposed MR
measure q[104] -> rec[34]; reset q[104]; // decomposed MR
measure q[107] -> rec[35]; reset q[107]; // decomposed MR
measure q[110] -> rec[36]; reset q[110]; // decomposed MR
measure q[113] -> rec[37]; reset q[113]; // decomposed MR
measure q[115] -> rec[38]; reset q[115]; // decomposed MR
measure q[118] -> rec[39]; reset q[118]; // decomposed MR
measure q[122] -> rec[40]; reset q[122]; // decomposed MR
measure q[125] -> rec[41]; reset q[125]; // decomposed MR
measure q[128] -> rec[42]; reset q[128]; // decomposed MR
measure q[130] -> rec[43]; reset q[130]; // decomposed MR
measure q[134] -> rec[44]; reset q[134]; // decomposed MR
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
cxyz q[13];
cxyz q[15];
cxyz q[17];
cxyz q[18];
cxyz q[20];
cxyz q[21];
cxyz q[23];
cxyz q[24];
cxyz q[26];
cxyz q[27];
cxyz q[29];
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
cxyz q[46];
cxyz q[48];
cxyz q[49];
cxyz q[51];
cxyz q[52];
cxyz q[54];
cxyz q[55];
cxyz q[57];
cxyz q[59];
cxyz q[60];
cxyz q[62];
cxyz q[63];
cxyz q[65];
cxyz q[66];
cxyz q[68];
cxyz q[69];
cxyz q[70];
cxyz q[72];
cxyz q[73];
cxyz q[75];
cxyz q[76];
cxyz q[78];
cxyz q[79];
cxyz q[81];
cxyz q[82];
cxyz q[84];
cxyz q[85];
cxyz q[87];
cxyz q[88];
cxyz q[90];
cxyz q[92];
cxyz q[93];
cxyz q[95];
cxyz q[96];
cxyz q[98];
cxyz q[99];
cxyz q[100];
cxyz q[102];
cxyz q[103];
cxyz q[105];
cxyz q[106];
cxyz q[108];
cxyz q[109];
cxyz q[111];
cxyz q[112];
cxyz q[114];
cxyz q[116];
cxyz q[117];
cxyz q[119];
cxyz q[120];
cxyz q[121];
cxyz q[123];
cxyz q[124];
cxyz q[126];
cxyz q[127];
cxyz q[129];
cxyz q[131];
cxyz q[132];
cxyz q[133];
cxyz q[135];
barrier q;

cx q[17], q[16];
cx q[3], q[2];
cx q[33], q[32];
cx q[59], q[58];
cx q[20], q[19];
cx q[48], q[47];
cx q[72], q[71];
cx q[92], q[91];
cx q[6], q[5];
cx q[36], q[35];
cx q[62], q[61];
cx q[84], q[83];
cx q[102], q[101];
cx q[116], q[115];
cx q[23], q[22];
cx q[51], q[50];
cx q[75], q[74];
cx q[95], q[94];
cx q[111], q[110];
cx q[123], q[122];
cx q[131], q[130];
cx q[9], q[8];
cx q[39], q[38];
cx q[65], q[64];
cx q[87], q[86];
cx q[105], q[104];
cx q[119], q[118];
cx q[129], q[128];
cx q[26], q[25];
cx q[54], q[53];
cx q[78], q[77];
cx q[98], q[97];
cx q[114], q[113];
cx q[12], q[11];
cx q[42], q[41];
cx q[68], q[67];
cx q[90], q[89];
cx q[29], q[28];
cx q[57], q[56];
cx q[15], q[14];
barrier q;

cx q[31], q[16];
cx q[18], q[2];
cx q[46], q[32];
cx q[70], q[58];
cx q[34], q[19];
cx q[60], q[47];
cx q[82], q[71];
cx q[100], q[91];
cx q[21], q[5];
cx q[49], q[35];
cx q[73], q[61];
cx q[93], q[83];
cx q[109], q[101];
cx q[121], q[115];
cx q[37], q[22];
cx q[63], q[50];
cx q[85], q[74];
cx q[103], q[94];
cx q[117], q[110];
cx q[127], q[122];
cx q[133], q[130];
cx q[24], q[8];
cx q[52], q[38];
cx q[76], q[64];
cx q[96], q[86];
cx q[112], q[104];
cx q[124], q[118];
cx q[132], q[128];
cx q[40], q[25];
cx q[66], q[53];
cx q[88], q[77];
cx q[106], q[97];
cx q[120], q[113];
cx q[27], q[11];
cx q[55], q[41];
cx q[79], q[67];
cx q[99], q[89];
cx q[43], q[28];
cx q[69], q[56];
cx q[30], q[14];
barrier q;

cx q[1], q[16];
cx q[18], q[32];
cx q[46], q[58];
cx q[4], q[19];
cx q[34], q[47];
cx q[60], q[71];
cx q[82], q[91];
cx q[21], q[35];
cx q[49], q[61];
cx q[73], q[83];
cx q[93], q[101];
cx q[109], q[115];
cx q[7], q[22];
cx q[37], q[50];
cx q[63], q[74];
cx q[85], q[94];
cx q[103], q[110];
cx q[117], q[122];
cx q[127], q[130];
cx q[24], q[38];
cx q[52], q[64];
cx q[76], q[86];
cx q[96], q[104];
cx q[112], q[118];
cx q[124], q[128];
cx q[132], q[134];
cx q[10], q[25];
cx q[40], q[53];
cx q[66], q[77];
cx q[88], q[97];
cx q[106], q[113];
cx q[120], q[125];
cx q[27], q[41];
cx q[55], q[67];
cx q[79], q[89];
cx q[99], q[107];
cx q[13], q[28];
cx q[43], q[56];
cx q[69], q[80];
cx q[30], q[44];
barrier q;

cx q[1], q[2];
cx q[31], q[32];
cx q[18], q[19];
cx q[46], q[47];
cx q[70], q[71];
cx q[4], q[5];
cx q[34], q[35];
cx q[60], q[61];
cx q[82], q[83];
cx q[100], q[101];
cx q[21], q[22];
cx q[49], q[50];
cx q[73], q[74];
cx q[93], q[94];
cx q[109], q[110];
cx q[121], q[122];
cx q[7], q[8];
cx q[37], q[38];
cx q[63], q[64];
cx q[85], q[86];
cx q[103], q[104];
cx q[117], q[118];
cx q[127], q[128];
cx q[133], q[134];
cx q[24], q[25];
cx q[52], q[53];
cx q[76], q[77];
cx q[96], q[97];
cx q[112], q[113];
cx q[124], q[125];
cx q[10], q[11];
cx q[40], q[41];
cx q[66], q[67];
cx q[88], q[89];
cx q[106], q[107];
cx q[27], q[28];
cx q[55], q[56];
cx q[79], q[80];
cx q[13], q[14];
cx q[43], q[44];
barrier q;

cx q[17], q[2];
cx q[45], q[32];
cx q[33], q[19];
cx q[59], q[47];
cx q[81], q[71];
cx q[20], q[5];
cx q[48], q[35];
cx q[72], q[61];
cx q[92], q[83];
cx q[108], q[101];
cx q[36], q[22];
cx q[62], q[50];
cx q[84], q[74];
cx q[102], q[94];
cx q[116], q[110];
cx q[126], q[122];
cx q[23], q[8];
cx q[51], q[38];
cx q[75], q[64];
cx q[95], q[86];
cx q[111], q[104];
cx q[123], q[118];
cx q[131], q[128];
cx q[135], q[134];
cx q[39], q[25];
cx q[65], q[53];
cx q[87], q[77];
cx q[105], q[97];
cx q[119], q[113];
cx q[129], q[125];
cx q[26], q[11];
cx q[54], q[41];
cx q[78], q[67];
cx q[98], q[89];
cx q[114], q[107];
cx q[42], q[28];
cx q[68], q[56];
cx q[90], q[80];
cx q[29], q[14];
cx q[57], q[44];
barrier q;

cx q[0], q[16];
cx q[17], q[32];
cx q[45], q[58];
cx q[3], q[19];
cx q[33], q[47];
cx q[59], q[71];
cx q[81], q[91];
cx q[20], q[35];
cx q[48], q[61];
cx q[72], q[83];
cx q[92], q[101];
cx q[108], q[115];
cx q[6], q[22];
cx q[36], q[50];
cx q[62], q[74];
cx q[84], q[94];
cx q[102], q[110];
cx q[116], q[122];
cx q[126], q[130];
cx q[23], q[38];
cx q[51], q[64];
cx q[75], q[86];
cx q[95], q[104];
cx q[111], q[118];
cx q[123], q[128];
cx q[131], q[134];
cx q[9], q[25];
cx q[39], q[53];
cx q[65], q[77];
cx q[87], q[97];
cx q[105], q[113];
cx q[119], q[125];
cx q[26], q[41];
cx q[54], q[67];
cx q[78], q[89];
cx q[98], q[107];
cx q[12], q[28];
cx q[42], q[56];
cx q[68], q[80];
cx q[29], q[44];
barrier q;

measure q[2] -> rec[45]; reset q[2]; // decomposed MR
measure q[5] -> rec[46]; reset q[5]; // decomposed MR
measure q[8] -> rec[47]; reset q[8]; // decomposed MR
measure q[11] -> rec[48]; reset q[11]; // decomposed MR
measure q[14] -> rec[49]; reset q[14]; // decomposed MR
measure q[16] -> rec[50]; reset q[16]; // decomposed MR
measure q[19] -> rec[51]; reset q[19]; // decomposed MR
measure q[22] -> rec[52]; reset q[22]; // decomposed MR
measure q[25] -> rec[53]; reset q[25]; // decomposed MR
measure q[28] -> rec[54]; reset q[28]; // decomposed MR
measure q[32] -> rec[55]; reset q[32]; // decomposed MR
measure q[35] -> rec[56]; reset q[35]; // decomposed MR
measure q[38] -> rec[57]; reset q[38]; // decomposed MR
measure q[41] -> rec[58]; reset q[41]; // decomposed MR
measure q[44] -> rec[59]; reset q[44]; // decomposed MR
measure q[47] -> rec[60]; reset q[47]; // decomposed MR
measure q[50] -> rec[61]; reset q[50]; // decomposed MR
measure q[53] -> rec[62]; reset q[53]; // decomposed MR
measure q[56] -> rec[63]; reset q[56]; // decomposed MR
measure q[58] -> rec[64]; reset q[58]; // decomposed MR
measure q[61] -> rec[65]; reset q[61]; // decomposed MR
measure q[64] -> rec[66]; reset q[64]; // decomposed MR
measure q[67] -> rec[67]; reset q[67]; // decomposed MR
measure q[71] -> rec[68]; reset q[71]; // decomposed MR
measure q[74] -> rec[69]; reset q[74]; // decomposed MR
measure q[77] -> rec[70]; reset q[77]; // decomposed MR
measure q[80] -> rec[71]; reset q[80]; // decomposed MR
measure q[83] -> rec[72]; reset q[83]; // decomposed MR
measure q[86] -> rec[73]; reset q[86]; // decomposed MR
measure q[89] -> rec[74]; reset q[89]; // decomposed MR
measure q[91] -> rec[75]; reset q[91]; // decomposed MR
measure q[94] -> rec[76]; reset q[94]; // decomposed MR
measure q[97] -> rec[77]; reset q[97]; // decomposed MR
measure q[101] -> rec[78]; reset q[101]; // decomposed MR
measure q[104] -> rec[79]; reset q[104]; // decomposed MR
measure q[107] -> rec[80]; reset q[107]; // decomposed MR
measure q[110] -> rec[81]; reset q[110]; // decomposed MR
measure q[113] -> rec[82]; reset q[113]; // decomposed MR
measure q[115] -> rec[83]; reset q[115]; // decomposed MR
measure q[118] -> rec[84]; reset q[118]; // decomposed MR
measure q[122] -> rec[85]; reset q[122]; // decomposed MR
measure q[125] -> rec[86]; reset q[125]; // decomposed MR
measure q[128] -> rec[87]; reset q[128]; // decomposed MR
measure q[130] -> rec[88]; reset q[130]; // decomposed MR
measure q[134] -> rec[89]; reset q[134]; // decomposed MR
s q[0]; s q[0]; s q[0]; h q[0]; measure q[0] -> rec[90]; h q[0]; s q[0]; // decomposed MY
s q[1]; s q[1]; s q[1]; h q[1]; measure q[1] -> rec[91]; h q[1]; s q[1]; // decomposed MY
s q[3]; s q[3]; s q[3]; h q[3]; measure q[3] -> rec[92]; h q[3]; s q[3]; // decomposed MY
s q[4]; s q[4]; s q[4]; h q[4]; measure q[4] -> rec[93]; h q[4]; s q[4]; // decomposed MY
s q[6]; s q[6]; s q[6]; h q[6]; measure q[6] -> rec[94]; h q[6]; s q[6]; // decomposed MY
s q[7]; s q[7]; s q[7]; h q[7]; measure q[7] -> rec[95]; h q[7]; s q[7]; // decomposed MY
s q[9]; s q[9]; s q[9]; h q[9]; measure q[9] -> rec[96]; h q[9]; s q[9]; // decomposed MY
s q[10]; s q[10]; s q[10]; h q[10]; measure q[10] -> rec[97]; h q[10]; s q[10]; // decomposed MY
s q[12]; s q[12]; s q[12]; h q[12]; measure q[12] -> rec[98]; h q[12]; s q[12]; // decomposed MY
s q[13]; s q[13]; s q[13]; h q[13]; measure q[13] -> rec[99]; h q[13]; s q[13]; // decomposed MY
s q[15]; s q[15]; s q[15]; h q[15]; measure q[15] -> rec[100]; h q[15]; s q[15]; // decomposed MY
s q[17]; s q[17]; s q[17]; h q[17]; measure q[17] -> rec[101]; h q[17]; s q[17]; // decomposed MY
s q[18]; s q[18]; s q[18]; h q[18]; measure q[18] -> rec[102]; h q[18]; s q[18]; // decomposed MY
s q[20]; s q[20]; s q[20]; h q[20]; measure q[20] -> rec[103]; h q[20]; s q[20]; // decomposed MY
s q[21]; s q[21]; s q[21]; h q[21]; measure q[21] -> rec[104]; h q[21]; s q[21]; // decomposed MY
s q[23]; s q[23]; s q[23]; h q[23]; measure q[23] -> rec[105]; h q[23]; s q[23]; // decomposed MY
s q[24]; s q[24]; s q[24]; h q[24]; measure q[24] -> rec[106]; h q[24]; s q[24]; // decomposed MY
s q[26]; s q[26]; s q[26]; h q[26]; measure q[26] -> rec[107]; h q[26]; s q[26]; // decomposed MY
s q[27]; s q[27]; s q[27]; h q[27]; measure q[27] -> rec[108]; h q[27]; s q[27]; // decomposed MY
s q[29]; s q[29]; s q[29]; h q[29]; measure q[29] -> rec[109]; h q[29]; s q[29]; // decomposed MY
s q[30]; s q[30]; s q[30]; h q[30]; measure q[30] -> rec[110]; h q[30]; s q[30]; // decomposed MY
s q[31]; s q[31]; s q[31]; h q[31]; measure q[31] -> rec[111]; h q[31]; s q[31]; // decomposed MY
s q[33]; s q[33]; s q[33]; h q[33]; measure q[33] -> rec[112]; h q[33]; s q[33]; // decomposed MY
s q[34]; s q[34]; s q[34]; h q[34]; measure q[34] -> rec[113]; h q[34]; s q[34]; // decomposed MY
s q[36]; s q[36]; s q[36]; h q[36]; measure q[36] -> rec[114]; h q[36]; s q[36]; // decomposed MY
s q[37]; s q[37]; s q[37]; h q[37]; measure q[37] -> rec[115]; h q[37]; s q[37]; // decomposed MY
s q[39]; s q[39]; s q[39]; h q[39]; measure q[39] -> rec[116]; h q[39]; s q[39]; // decomposed MY
s q[40]; s q[40]; s q[40]; h q[40]; measure q[40] -> rec[117]; h q[40]; s q[40]; // decomposed MY
s q[42]; s q[42]; s q[42]; h q[42]; measure q[42] -> rec[118]; h q[42]; s q[42]; // decomposed MY
s q[43]; s q[43]; s q[43]; h q[43]; measure q[43] -> rec[119]; h q[43]; s q[43]; // decomposed MY
s q[45]; s q[45]; s q[45]; h q[45]; measure q[45] -> rec[120]; h q[45]; s q[45]; // decomposed MY
s q[46]; s q[46]; s q[46]; h q[46]; measure q[46] -> rec[121]; h q[46]; s q[46]; // decomposed MY
s q[48]; s q[48]; s q[48]; h q[48]; measure q[48] -> rec[122]; h q[48]; s q[48]; // decomposed MY
s q[49]; s q[49]; s q[49]; h q[49]; measure q[49] -> rec[123]; h q[49]; s q[49]; // decomposed MY
s q[51]; s q[51]; s q[51]; h q[51]; measure q[51] -> rec[124]; h q[51]; s q[51]; // decomposed MY
s q[52]; s q[52]; s q[52]; h q[52]; measure q[52] -> rec[125]; h q[52]; s q[52]; // decomposed MY
s q[54]; s q[54]; s q[54]; h q[54]; measure q[54] -> rec[126]; h q[54]; s q[54]; // decomposed MY
s q[55]; s q[55]; s q[55]; h q[55]; measure q[55] -> rec[127]; h q[55]; s q[55]; // decomposed MY
s q[57]; s q[57]; s q[57]; h q[57]; measure q[57] -> rec[128]; h q[57]; s q[57]; // decomposed MY
s q[59]; s q[59]; s q[59]; h q[59]; measure q[59] -> rec[129]; h q[59]; s q[59]; // decomposed MY
s q[60]; s q[60]; s q[60]; h q[60]; measure q[60] -> rec[130]; h q[60]; s q[60]; // decomposed MY
s q[62]; s q[62]; s q[62]; h q[62]; measure q[62] -> rec[131]; h q[62]; s q[62]; // decomposed MY
s q[63]; s q[63]; s q[63]; h q[63]; measure q[63] -> rec[132]; h q[63]; s q[63]; // decomposed MY
s q[65]; s q[65]; s q[65]; h q[65]; measure q[65] -> rec[133]; h q[65]; s q[65]; // decomposed MY
s q[66]; s q[66]; s q[66]; h q[66]; measure q[66] -> rec[134]; h q[66]; s q[66]; // decomposed MY
s q[68]; s q[68]; s q[68]; h q[68]; measure q[68] -> rec[135]; h q[68]; s q[68]; // decomposed MY
s q[69]; s q[69]; s q[69]; h q[69]; measure q[69] -> rec[136]; h q[69]; s q[69]; // decomposed MY
s q[70]; s q[70]; s q[70]; h q[70]; measure q[70] -> rec[137]; h q[70]; s q[70]; // decomposed MY
s q[72]; s q[72]; s q[72]; h q[72]; measure q[72] -> rec[138]; h q[72]; s q[72]; // decomposed MY
s q[73]; s q[73]; s q[73]; h q[73]; measure q[73] -> rec[139]; h q[73]; s q[73]; // decomposed MY
s q[75]; s q[75]; s q[75]; h q[75]; measure q[75] -> rec[140]; h q[75]; s q[75]; // decomposed MY
s q[76]; s q[76]; s q[76]; h q[76]; measure q[76] -> rec[141]; h q[76]; s q[76]; // decomposed MY
s q[78]; s q[78]; s q[78]; h q[78]; measure q[78] -> rec[142]; h q[78]; s q[78]; // decomposed MY
s q[79]; s q[79]; s q[79]; h q[79]; measure q[79] -> rec[143]; h q[79]; s q[79]; // decomposed MY
s q[81]; s q[81]; s q[81]; h q[81]; measure q[81] -> rec[144]; h q[81]; s q[81]; // decomposed MY
s q[82]; s q[82]; s q[82]; h q[82]; measure q[82] -> rec[145]; h q[82]; s q[82]; // decomposed MY
s q[84]; s q[84]; s q[84]; h q[84]; measure q[84] -> rec[146]; h q[84]; s q[84]; // decomposed MY
s q[85]; s q[85]; s q[85]; h q[85]; measure q[85] -> rec[147]; h q[85]; s q[85]; // decomposed MY
s q[87]; s q[87]; s q[87]; h q[87]; measure q[87] -> rec[148]; h q[87]; s q[87]; // decomposed MY
s q[88]; s q[88]; s q[88]; h q[88]; measure q[88] -> rec[149]; h q[88]; s q[88]; // decomposed MY
s q[90]; s q[90]; s q[90]; h q[90]; measure q[90] -> rec[150]; h q[90]; s q[90]; // decomposed MY
s q[92]; s q[92]; s q[92]; h q[92]; measure q[92] -> rec[151]; h q[92]; s q[92]; // decomposed MY
s q[93]; s q[93]; s q[93]; h q[93]; measure q[93] -> rec[152]; h q[93]; s q[93]; // decomposed MY
s q[95]; s q[95]; s q[95]; h q[95]; measure q[95] -> rec[153]; h q[95]; s q[95]; // decomposed MY
s q[96]; s q[96]; s q[96]; h q[96]; measure q[96] -> rec[154]; h q[96]; s q[96]; // decomposed MY
s q[98]; s q[98]; s q[98]; h q[98]; measure q[98] -> rec[155]; h q[98]; s q[98]; // decomposed MY
s q[99]; s q[99]; s q[99]; h q[99]; measure q[99] -> rec[156]; h q[99]; s q[99]; // decomposed MY
s q[100]; s q[100]; s q[100]; h q[100]; measure q[100] -> rec[157]; h q[100]; s q[100]; // decomposed MY
s q[102]; s q[102]; s q[102]; h q[102]; measure q[102] -> rec[158]; h q[102]; s q[102]; // decomposed MY
s q[103]; s q[103]; s q[103]; h q[103]; measure q[103] -> rec[159]; h q[103]; s q[103]; // decomposed MY
s q[105]; s q[105]; s q[105]; h q[105]; measure q[105] -> rec[160]; h q[105]; s q[105]; // decomposed MY
s q[106]; s q[106]; s q[106]; h q[106]; measure q[106] -> rec[161]; h q[106]; s q[106]; // decomposed MY
s q[108]; s q[108]; s q[108]; h q[108]; measure q[108] -> rec[162]; h q[108]; s q[108]; // decomposed MY
s q[109]; s q[109]; s q[109]; h q[109]; measure q[109] -> rec[163]; h q[109]; s q[109]; // decomposed MY
s q[111]; s q[111]; s q[111]; h q[111]; measure q[111] -> rec[164]; h q[111]; s q[111]; // decomposed MY
s q[112]; s q[112]; s q[112]; h q[112]; measure q[112] -> rec[165]; h q[112]; s q[112]; // decomposed MY
s q[114]; s q[114]; s q[114]; h q[114]; measure q[114] -> rec[166]; h q[114]; s q[114]; // decomposed MY
s q[116]; s q[116]; s q[116]; h q[116]; measure q[116] -> rec[167]; h q[116]; s q[116]; // decomposed MY
s q[117]; s q[117]; s q[117]; h q[117]; measure q[117] -> rec[168]; h q[117]; s q[117]; // decomposed MY
s q[119]; s q[119]; s q[119]; h q[119]; measure q[119] -> rec[169]; h q[119]; s q[119]; // decomposed MY
s q[120]; s q[120]; s q[120]; h q[120]; measure q[120] -> rec[170]; h q[120]; s q[120]; // decomposed MY
s q[121]; s q[121]; s q[121]; h q[121]; measure q[121] -> rec[171]; h q[121]; s q[121]; // decomposed MY
s q[123]; s q[123]; s q[123]; h q[123]; measure q[123] -> rec[172]; h q[123]; s q[123]; // decomposed MY
s q[124]; s q[124]; s q[124]; h q[124]; measure q[124] -> rec[173]; h q[124]; s q[124]; // decomposed MY
s q[126]; s q[126]; s q[126]; h q[126]; measure q[126] -> rec[174]; h q[126]; s q[126]; // decomposed MY
s q[127]; s q[127]; s q[127]; h q[127]; measure q[127] -> rec[175]; h q[127]; s q[127]; // decomposed MY
s q[129]; s q[129]; s q[129]; h q[129]; measure q[129] -> rec[176]; h q[129]; s q[129]; // decomposed MY
s q[131]; s q[131]; s q[131]; h q[131]; measure q[131] -> rec[177]; h q[131]; s q[131]; // decomposed MY
s q[132]; s q[132]; s q[132]; h q[132]; measure q[132] -> rec[178]; h q[132]; s q[132]; // decomposed MY
s q[133]; s q[133]; s q[133]; h q[133]; measure q[133] -> rec[179]; h q[133]; s q[133]; // decomposed MY
s q[135]; s q[135]; s q[135]; h q[135]; measure q[135] -> rec[180]; h q[135]; s q[135]; // decomposed MY

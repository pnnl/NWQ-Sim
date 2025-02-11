OPENQASM 2.0;
include "qelib1.inc";

qreg q[188];
creg rec[161];

reset q[1];
reset q[3];
reset q[5];
reset q[7];
reset q[9];
reset q[11];
reset q[13];
reset q[15];
reset q[17];
reset q[20];
reset q[22];
reset q[24];
reset q[26];
reset q[28];
reset q[30];
reset q[32];
reset q[34];
reset q[36];
reset q[39];
reset q[41];
reset q[43];
reset q[45];
reset q[47];
reset q[49];
reset q[51];
reset q[53];
reset q[55];
reset q[58];
reset q[60];
reset q[62];
reset q[64];
reset q[66];
reset q[68];
reset q[70];
reset q[72];
reset q[74];
reset q[77];
reset q[79];
reset q[81];
reset q[83];
reset q[85];
reset q[87];
reset q[89];
reset q[91];
reset q[93];
reset q[96];
reset q[98];
reset q[100];
reset q[102];
reset q[104];
reset q[106];
reset q[108];
reset q[110];
reset q[112];
reset q[115];
reset q[117];
reset q[119];
reset q[121];
reset q[123];
reset q[125];
reset q[127];
reset q[129];
reset q[131];
reset q[134];
reset q[136];
reset q[138];
reset q[140];
reset q[142];
reset q[144];
reset q[146];
reset q[148];
reset q[150];
reset q[153];
reset q[155];
reset q[157];
reset q[159];
reset q[161];
reset q[163];
reset q[165];
reset q[167];
reset q[169];
reset q[2];
reset q[6];
reset q[10];
reset q[14];
reset q[21];
reset q[23];
reset q[25];
reset q[27];
reset q[29];
reset q[31];
reset q[33];
reset q[35];
reset q[37];
reset q[38];
reset q[40];
reset q[42];
reset q[44];
reset q[46];
reset q[48];
reset q[50];
reset q[52];
reset q[54];
reset q[59];
reset q[61];
reset q[63];
reset q[65];
reset q[67];
reset q[69];
reset q[71];
reset q[73];
reset q[75];
reset q[76];
reset q[78];
reset q[80];
reset q[82];
reset q[84];
reset q[86];
reset q[88];
reset q[90];
reset q[92];
reset q[97];
reset q[99];
reset q[101];
reset q[103];
reset q[105];
reset q[107];
reset q[109];
reset q[111];
reset q[113];
reset q[114];
reset q[116];
reset q[118];
reset q[120];
reset q[122];
reset q[124];
reset q[126];
reset q[128];
reset q[130];
reset q[135];
reset q[137];
reset q[139];
reset q[141];
reset q[143];
reset q[145];
reset q[147];
reset q[149];
reset q[151];
reset q[152];
reset q[154];
reset q[156];
reset q[158];
reset q[160];
reset q[162];
reset q[164];
reset q[166];
reset q[168];
reset q[175];
reset q[179];
reset q[183];
reset q[187];
barrier q;

h q[2];
h q[6];
h q[10];
h q[14];
h q[23];
h q[27];
h q[31];
h q[35];
h q[40];
h q[44];
h q[48];
h q[52];
h q[61];
h q[65];
h q[69];
h q[73];
h q[78];
h q[82];
h q[86];
h q[90];
h q[99];
h q[103];
h q[107];
h q[111];
h q[116];
h q[120];
h q[124];
h q[128];
h q[137];
h q[141];
h q[145];
h q[149];
h q[154];
h q[158];
h q[162];
h q[166];
h q[175];
h q[179];
h q[183];
h q[187];
barrier q;

cx q[2], q[3];
cx q[40], q[41];
cx q[78], q[79];
cx q[116], q[117];
cx q[154], q[155];
cx q[23], q[24];
cx q[61], q[62];
cx q[99], q[100];
cx q[137], q[138];
cx q[6], q[7];
cx q[44], q[45];
cx q[82], q[83];
cx q[120], q[121];
cx q[158], q[159];
cx q[27], q[28];
cx q[65], q[66];
cx q[103], q[104];
cx q[141], q[142];
cx q[10], q[11];
cx q[48], q[49];
cx q[86], q[87];
cx q[124], q[125];
cx q[162], q[163];
cx q[31], q[32];
cx q[69], q[70];
cx q[107], q[108];
cx q[145], q[146];
cx q[14], q[15];
cx q[52], q[53];
cx q[90], q[91];
cx q[128], q[129];
cx q[166], q[167];
cx q[35], q[36];
cx q[73], q[74];
cx q[111], q[112];
cx q[149], q[150];
cx q[39], q[38];
cx q[77], q[76];
cx q[115], q[114];
cx q[153], q[152];
cx q[22], q[21];
cx q[60], q[59];
cx q[98], q[97];
cx q[136], q[135];
cx q[43], q[42];
cx q[81], q[80];
cx q[119], q[118];
cx q[157], q[156];
cx q[26], q[25];
cx q[64], q[63];
cx q[102], q[101];
cx q[140], q[139];
cx q[47], q[46];
cx q[85], q[84];
cx q[123], q[122];
cx q[161], q[160];
cx q[30], q[29];
cx q[68], q[67];
cx q[106], q[105];
cx q[144], q[143];
cx q[51], q[50];
cx q[89], q[88];
cx q[127], q[126];
cx q[165], q[164];
cx q[34], q[33];
cx q[72], q[71];
cx q[110], q[109];
cx q[148], q[147];
cx q[55], q[54];
cx q[93], q[92];
cx q[131], q[130];
cx q[169], q[168];
barrier q;

cx q[2], q[1];
cx q[40], q[39];
cx q[78], q[77];
cx q[116], q[115];
cx q[154], q[153];
cx q[23], q[22];
cx q[61], q[60];
cx q[99], q[98];
cx q[137], q[136];
cx q[6], q[5];
cx q[44], q[43];
cx q[82], q[81];
cx q[120], q[119];
cx q[158], q[157];
cx q[27], q[26];
cx q[65], q[64];
cx q[103], q[102];
cx q[141], q[140];
cx q[10], q[9];
cx q[48], q[47];
cx q[86], q[85];
cx q[124], q[123];
cx q[162], q[161];
cx q[31], q[30];
cx q[69], q[68];
cx q[107], q[106];
cx q[145], q[144];
cx q[14], q[13];
cx q[52], q[51];
cx q[90], q[89];
cx q[128], q[127];
cx q[166], q[165];
cx q[35], q[34];
cx q[73], q[72];
cx q[111], q[110];
cx q[149], q[148];
cx q[20], q[38];
cx q[58], q[76];
cx q[96], q[114];
cx q[134], q[152];
cx q[3], q[21];
cx q[41], q[59];
cx q[79], q[97];
cx q[117], q[135];
cx q[24], q[42];
cx q[62], q[80];
cx q[100], q[118];
cx q[138], q[156];
cx q[7], q[25];
cx q[45], q[63];
cx q[83], q[101];
cx q[121], q[139];
cx q[28], q[46];
cx q[66], q[84];
cx q[104], q[122];
cx q[142], q[160];
cx q[11], q[29];
cx q[49], q[67];
cx q[87], q[105];
cx q[125], q[143];
cx q[32], q[50];
cx q[70], q[88];
cx q[108], q[126];
cx q[146], q[164];
cx q[15], q[33];
cx q[53], q[71];
cx q[91], q[109];
cx q[129], q[147];
cx q[36], q[54];
cx q[74], q[92];
cx q[112], q[130];
cx q[150], q[168];
barrier q;

cx q[40], q[22];
cx q[78], q[60];
cx q[116], q[98];
cx q[154], q[136];
cx q[23], q[5];
cx q[61], q[43];
cx q[99], q[81];
cx q[137], q[119];
cx q[175], q[157];
cx q[44], q[26];
cx q[82], q[64];
cx q[120], q[102];
cx q[158], q[140];
cx q[27], q[9];
cx q[65], q[47];
cx q[103], q[85];
cx q[141], q[123];
cx q[179], q[161];
cx q[48], q[30];
cx q[86], q[68];
cx q[124], q[106];
cx q[162], q[144];
cx q[31], q[13];
cx q[69], q[51];
cx q[107], q[89];
cx q[145], q[127];
cx q[183], q[165];
cx q[52], q[34];
cx q[90], q[72];
cx q[128], q[110];
cx q[166], q[148];
cx q[35], q[17];
cx q[73], q[55];
cx q[111], q[93];
cx q[149], q[131];
cx q[187], q[169];
cx q[20], q[21];
cx q[58], q[59];
cx q[96], q[97];
cx q[134], q[135];
cx q[41], q[42];
cx q[79], q[80];
cx q[117], q[118];
cx q[155], q[156];
cx q[24], q[25];
cx q[62], q[63];
cx q[100], q[101];
cx q[138], q[139];
cx q[45], q[46];
cx q[83], q[84];
cx q[121], q[122];
cx q[159], q[160];
cx q[28], q[29];
cx q[66], q[67];
cx q[104], q[105];
cx q[142], q[143];
cx q[49], q[50];
cx q[87], q[88];
cx q[125], q[126];
cx q[163], q[164];
cx q[32], q[33];
cx q[70], q[71];
cx q[108], q[109];
cx q[146], q[147];
cx q[53], q[54];
cx q[91], q[92];
cx q[129], q[130];
cx q[167], q[168];
cx q[36], q[37];
cx q[74], q[75];
cx q[112], q[113];
cx q[150], q[151];
barrier q;

cx q[40], q[20];
cx q[78], q[58];
cx q[116], q[96];
cx q[154], q[134];
cx q[23], q[3];
cx q[61], q[41];
cx q[99], q[79];
cx q[137], q[117];
cx q[175], q[155];
cx q[44], q[24];
cx q[82], q[62];
cx q[120], q[100];
cx q[158], q[138];
cx q[27], q[7];
cx q[65], q[45];
cx q[103], q[83];
cx q[141], q[121];
cx q[179], q[159];
cx q[48], q[28];
cx q[86], q[66];
cx q[124], q[104];
cx q[162], q[142];
cx q[31], q[11];
cx q[69], q[49];
cx q[107], q[87];
cx q[145], q[125];
cx q[183], q[163];
cx q[52], q[32];
cx q[90], q[70];
cx q[128], q[108];
cx q[166], q[146];
cx q[35], q[15];
cx q[73], q[53];
cx q[111], q[91];
cx q[149], q[129];
cx q[187], q[167];
cx q[1], q[21];
cx q[39], q[59];
cx q[77], q[97];
cx q[115], q[135];
cx q[22], q[42];
cx q[60], q[80];
cx q[98], q[118];
cx q[136], q[156];
cx q[5], q[25];
cx q[43], q[63];
cx q[81], q[101];
cx q[119], q[139];
cx q[26], q[46];
cx q[64], q[84];
cx q[102], q[122];
cx q[140], q[160];
cx q[9], q[29];
cx q[47], q[67];
cx q[85], q[105];
cx q[123], q[143];
cx q[30], q[50];
cx q[68], q[88];
cx q[106], q[126];
cx q[144], q[164];
cx q[13], q[33];
cx q[51], q[71];
cx q[89], q[109];
cx q[127], q[147];
cx q[34], q[54];
cx q[72], q[92];
cx q[110], q[130];
cx q[148], q[168];
cx q[17], q[37];
cx q[55], q[75];
cx q[93], q[113];
cx q[131], q[151];
barrier q;

h q[2];
h q[6];
h q[10];
h q[14];
h q[23];
h q[27];
h q[31];
h q[35];
h q[40];
h q[44];
h q[48];
h q[52];
h q[61];
h q[65];
h q[69];
h q[73];
h q[78];
h q[82];
h q[86];
h q[90];
h q[99];
h q[103];
h q[107];
h q[111];
h q[116];
h q[120];
h q[124];
h q[128];
h q[137];
h q[141];
h q[145];
h q[149];
h q[154];
h q[158];
h q[162];
h q[166];
h q[175];
h q[179];
h q[183];
h q[187];
barrier q;

measure q[2] -> rec[0]; reset q[2]; // decomposed MR
measure q[6] -> rec[1]; reset q[6]; // decomposed MR
measure q[10] -> rec[2]; reset q[10]; // decomposed MR
measure q[14] -> rec[3]; reset q[14]; // decomposed MR
measure q[21] -> rec[4]; reset q[21]; // decomposed MR
measure q[23] -> rec[5]; reset q[23]; // decomposed MR
measure q[25] -> rec[6]; reset q[25]; // decomposed MR
measure q[27] -> rec[7]; reset q[27]; // decomposed MR
measure q[29] -> rec[8]; reset q[29]; // decomposed MR
measure q[31] -> rec[9]; reset q[31]; // decomposed MR
measure q[33] -> rec[10]; reset q[33]; // decomposed MR
measure q[35] -> rec[11]; reset q[35]; // decomposed MR
measure q[37] -> rec[12]; reset q[37]; // decomposed MR
measure q[38] -> rec[13]; reset q[38]; // decomposed MR
measure q[40] -> rec[14]; reset q[40]; // decomposed MR
measure q[42] -> rec[15]; reset q[42]; // decomposed MR
measure q[44] -> rec[16]; reset q[44]; // decomposed MR
measure q[46] -> rec[17]; reset q[46]; // decomposed MR
measure q[48] -> rec[18]; reset q[48]; // decomposed MR
measure q[50] -> rec[19]; reset q[50]; // decomposed MR
measure q[52] -> rec[20]; reset q[52]; // decomposed MR
measure q[54] -> rec[21]; reset q[54]; // decomposed MR
measure q[59] -> rec[22]; reset q[59]; // decomposed MR
measure q[61] -> rec[23]; reset q[61]; // decomposed MR
measure q[63] -> rec[24]; reset q[63]; // decomposed MR
measure q[65] -> rec[25]; reset q[65]; // decomposed MR
measure q[67] -> rec[26]; reset q[67]; // decomposed MR
measure q[69] -> rec[27]; reset q[69]; // decomposed MR
measure q[71] -> rec[28]; reset q[71]; // decomposed MR
measure q[73] -> rec[29]; reset q[73]; // decomposed MR
measure q[75] -> rec[30]; reset q[75]; // decomposed MR
measure q[76] -> rec[31]; reset q[76]; // decomposed MR
measure q[78] -> rec[32]; reset q[78]; // decomposed MR
measure q[80] -> rec[33]; reset q[80]; // decomposed MR
measure q[82] -> rec[34]; reset q[82]; // decomposed MR
measure q[84] -> rec[35]; reset q[84]; // decomposed MR
measure q[86] -> rec[36]; reset q[86]; // decomposed MR
measure q[88] -> rec[37]; reset q[88]; // decomposed MR
measure q[90] -> rec[38]; reset q[90]; // decomposed MR
measure q[92] -> rec[39]; reset q[92]; // decomposed MR
measure q[97] -> rec[40]; reset q[97]; // decomposed MR
measure q[99] -> rec[41]; reset q[99]; // decomposed MR
measure q[101] -> rec[42]; reset q[101]; // decomposed MR
measure q[103] -> rec[43]; reset q[103]; // decomposed MR
measure q[105] -> rec[44]; reset q[105]; // decomposed MR
measure q[107] -> rec[45]; reset q[107]; // decomposed MR
measure q[109] -> rec[46]; reset q[109]; // decomposed MR
measure q[111] -> rec[47]; reset q[111]; // decomposed MR
measure q[113] -> rec[48]; reset q[113]; // decomposed MR
measure q[114] -> rec[49]; reset q[114]; // decomposed MR
measure q[116] -> rec[50]; reset q[116]; // decomposed MR
measure q[118] -> rec[51]; reset q[118]; // decomposed MR
measure q[120] -> rec[52]; reset q[120]; // decomposed MR
measure q[122] -> rec[53]; reset q[122]; // decomposed MR
measure q[124] -> rec[54]; reset q[124]; // decomposed MR
measure q[126] -> rec[55]; reset q[126]; // decomposed MR
measure q[128] -> rec[56]; reset q[128]; // decomposed MR
measure q[130] -> rec[57]; reset q[130]; // decomposed MR
measure q[135] -> rec[58]; reset q[135]; // decomposed MR
measure q[137] -> rec[59]; reset q[137]; // decomposed MR
measure q[139] -> rec[60]; reset q[139]; // decomposed MR
measure q[141] -> rec[61]; reset q[141]; // decomposed MR
measure q[143] -> rec[62]; reset q[143]; // decomposed MR
measure q[145] -> rec[63]; reset q[145]; // decomposed MR
measure q[147] -> rec[64]; reset q[147]; // decomposed MR
measure q[149] -> rec[65]; reset q[149]; // decomposed MR
measure q[151] -> rec[66]; reset q[151]; // decomposed MR
measure q[152] -> rec[67]; reset q[152]; // decomposed MR
measure q[154] -> rec[68]; reset q[154]; // decomposed MR
measure q[156] -> rec[69]; reset q[156]; // decomposed MR
measure q[158] -> rec[70]; reset q[158]; // decomposed MR
measure q[160] -> rec[71]; reset q[160]; // decomposed MR
measure q[162] -> rec[72]; reset q[162]; // decomposed MR
measure q[164] -> rec[73]; reset q[164]; // decomposed MR
measure q[166] -> rec[74]; reset q[166]; // decomposed MR
measure q[168] -> rec[75]; reset q[168]; // decomposed MR
measure q[175] -> rec[76]; reset q[175]; // decomposed MR
measure q[179] -> rec[77]; reset q[179]; // decomposed MR
measure q[183] -> rec[78]; reset q[183]; // decomposed MR
measure q[187] -> rec[79]; reset q[187]; // decomposed MR
measure q[1] -> rec[80];
measure q[3] -> rec[81];
measure q[5] -> rec[82];
measure q[7] -> rec[83];
measure q[9] -> rec[84];
measure q[11] -> rec[85];
measure q[13] -> rec[86];
measure q[15] -> rec[87];
measure q[17] -> rec[88];
measure q[20] -> rec[89];
measure q[22] -> rec[90];
measure q[24] -> rec[91];
measure q[26] -> rec[92];
measure q[28] -> rec[93];
measure q[30] -> rec[94];
measure q[32] -> rec[95];
measure q[34] -> rec[96];
measure q[36] -> rec[97];
measure q[39] -> rec[98];
measure q[41] -> rec[99];
measure q[43] -> rec[100];
measure q[45] -> rec[101];
measure q[47] -> rec[102];
measure q[49] -> rec[103];
measure q[51] -> rec[104];
measure q[53] -> rec[105];
measure q[55] -> rec[106];
measure q[58] -> rec[107];
measure q[60] -> rec[108];
measure q[62] -> rec[109];
measure q[64] -> rec[110];
measure q[66] -> rec[111];
measure q[68] -> rec[112];
measure q[70] -> rec[113];
measure q[72] -> rec[114];
measure q[74] -> rec[115];
measure q[77] -> rec[116];
measure q[79] -> rec[117];
measure q[81] -> rec[118];
measure q[83] -> rec[119];
measure q[85] -> rec[120];
measure q[87] -> rec[121];
measure q[89] -> rec[122];
measure q[91] -> rec[123];
measure q[93] -> rec[124];
measure q[96] -> rec[125];
measure q[98] -> rec[126];
measure q[100] -> rec[127];
measure q[102] -> rec[128];
measure q[104] -> rec[129];
measure q[106] -> rec[130];
measure q[108] -> rec[131];
measure q[110] -> rec[132];
measure q[112] -> rec[133];
measure q[115] -> rec[134];
measure q[117] -> rec[135];
measure q[119] -> rec[136];
measure q[121] -> rec[137];
measure q[123] -> rec[138];
measure q[125] -> rec[139];
measure q[127] -> rec[140];
measure q[129] -> rec[141];
measure q[131] -> rec[142];
measure q[134] -> rec[143];
measure q[136] -> rec[144];
measure q[138] -> rec[145];
measure q[140] -> rec[146];
measure q[142] -> rec[147];
measure q[144] -> rec[148];
measure q[146] -> rec[149];
measure q[148] -> rec[150];
measure q[150] -> rec[151];
measure q[153] -> rec[152];
measure q[155] -> rec[153];
measure q[157] -> rec[154];
measure q[159] -> rec[155];
measure q[161] -> rec[156];
measure q[163] -> rec[157];
measure q[165] -> rec[158];
measure q[167] -> rec[159];
measure q[169] -> rec[160];

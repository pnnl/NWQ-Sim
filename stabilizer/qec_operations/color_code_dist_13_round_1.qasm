OPENQASM 2.0;
include "qelib1.inc";
gate cxyz q0 { U(pi/2, 0, pi/2) q0; }

qreg q[190];
creg rec[253];

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
reset q[136];
reset q[137];
reset q[138];
reset q[139];
reset q[140];
reset q[141];
reset q[142];
reset q[143];
reset q[144];
reset q[145];
reset q[146];
reset q[147];
reset q[148];
reset q[149];
reset q[150];
reset q[151];
reset q[152];
reset q[153];
reset q[154];
reset q[155];
reset q[156];
reset q[157];
reset q[158];
reset q[159];
reset q[160];
reset q[161];
reset q[162];
reset q[163];
reset q[164];
reset q[165];
reset q[166];
reset q[167];
reset q[168];
reset q[169];
reset q[170];
reset q[171];
reset q[172];
reset q[173];
reset q[174];
reset q[175];
reset q[176];
reset q[177];
reset q[178];
reset q[179];
reset q[180];
reset q[181];
reset q[182];
reset q[183];
reset q[184];
reset q[185];
reset q[186];
reset q[187];
reset q[188];
reset q[189];
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
cxyz q[16];
cxyz q[18];
cxyz q[20];
cxyz q[21];
cxyz q[23];
cxyz q[24];
cxyz q[26];
cxyz q[27];
cxyz q[29];
cxyz q[30];
cxyz q[32];
cxyz q[33];
cxyz q[35];
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
cxyz q[77];
cxyz q[78];
cxyz q[80];
cxyz q[81];
cxyz q[83];
cxyz q[84];
cxyz q[85];
cxyz q[87];
cxyz q[88];
cxyz q[90];
cxyz q[91];
cxyz q[93];
cxyz q[94];
cxyz q[96];
cxyz q[97];
cxyz q[99];
cxyz q[100];
cxyz q[102];
cxyz q[103];
cxyz q[105];
cxyz q[106];
cxyz q[108];
cxyz q[109];
cxyz q[111];
cxyz q[113];
cxyz q[114];
cxyz q[116];
cxyz q[117];
cxyz q[119];
cxyz q[120];
cxyz q[122];
cxyz q[123];
cxyz q[124];
cxyz q[126];
cxyz q[127];
cxyz q[129];
cxyz q[130];
cxyz q[132];
cxyz q[133];
cxyz q[135];
cxyz q[136];
cxyz q[138];
cxyz q[139];
cxyz q[141];
cxyz q[142];
cxyz q[144];
cxyz q[146];
cxyz q[147];
cxyz q[149];
cxyz q[150];
cxyz q[152];
cxyz q[153];
cxyz q[154];
cxyz q[156];
cxyz q[157];
cxyz q[159];
cxyz q[160];
cxyz q[162];
cxyz q[163];
cxyz q[165];
cxyz q[166];
cxyz q[168];
cxyz q[170];
cxyz q[171];
cxyz q[173];
cxyz q[174];
cxyz q[175];
cxyz q[177];
cxyz q[178];
cxyz q[180];
cxyz q[181];
cxyz q[183];
cxyz q[185];
cxyz q[186];
cxyz q[187];
cxyz q[189];
barrier q;

cx q[20], q[19];
cx q[3], q[2];
cx q[39], q[38];
cx q[71], q[70];
cx q[23], q[22];
cx q[57], q[56];
cx q[87], q[86];
cx q[113], q[112];
cx q[6], q[5];
cx q[42], q[41];
cx q[74], q[73];
cx q[102], q[101];
cx q[126], q[125];
cx q[146], q[145];
cx q[26], q[25];
cx q[60], q[59];
cx q[90], q[89];
cx q[116], q[115];
cx q[138], q[137];
cx q[156], q[155];
cx q[170], q[169];
cx q[9], q[8];
cx q[45], q[44];
cx q[77], q[76];
cx q[105], q[104];
cx q[129], q[128];
cx q[149], q[148];
cx q[165], q[164];
cx q[177], q[176];
cx q[185], q[184];
cx q[29], q[28];
cx q[63], q[62];
cx q[93], q[92];
cx q[119], q[118];
cx q[141], q[140];
cx q[159], q[158];
cx q[173], q[172];
cx q[183], q[182];
cx q[12], q[11];
cx q[48], q[47];
cx q[80], q[79];
cx q[108], q[107];
cx q[132], q[131];
cx q[152], q[151];
cx q[168], q[167];
cx q[32], q[31];
cx q[66], q[65];
cx q[96], q[95];
cx q[122], q[121];
cx q[144], q[143];
cx q[15], q[14];
cx q[51], q[50];
cx q[83], q[82];
cx q[111], q[110];
cx q[35], q[34];
cx q[69], q[68];
cx q[18], q[17];
barrier q;

cx q[37], q[19];
cx q[21], q[2];
cx q[55], q[38];
cx q[85], q[70];
cx q[40], q[22];
cx q[72], q[56];
cx q[100], q[86];
cx q[124], q[112];
cx q[24], q[5];
cx q[58], q[41];
cx q[88], q[73];
cx q[114], q[101];
cx q[136], q[125];
cx q[154], q[145];
cx q[43], q[25];
cx q[75], q[59];
cx q[103], q[89];
cx q[127], q[115];
cx q[147], q[137];
cx q[163], q[155];
cx q[175], q[169];
cx q[27], q[8];
cx q[61], q[44];
cx q[91], q[76];
cx q[117], q[104];
cx q[139], q[128];
cx q[157], q[148];
cx q[171], q[164];
cx q[181], q[176];
cx q[187], q[184];
cx q[46], q[28];
cx q[78], q[62];
cx q[106], q[92];
cx q[130], q[118];
cx q[150], q[140];
cx q[166], q[158];
cx q[178], q[172];
cx q[186], q[182];
cx q[30], q[11];
cx q[64], q[47];
cx q[94], q[79];
cx q[120], q[107];
cx q[142], q[131];
cx q[160], q[151];
cx q[174], q[167];
cx q[49], q[31];
cx q[81], q[65];
cx q[109], q[95];
cx q[133], q[121];
cx q[153], q[143];
cx q[33], q[14];
cx q[67], q[50];
cx q[97], q[82];
cx q[123], q[110];
cx q[52], q[34];
cx q[84], q[68];
cx q[36], q[17];
barrier q;

cx q[1], q[19];
cx q[21], q[38];
cx q[55], q[70];
cx q[4], q[22];
cx q[40], q[56];
cx q[72], q[86];
cx q[100], q[112];
cx q[24], q[41];
cx q[58], q[73];
cx q[88], q[101];
cx q[114], q[125];
cx q[136], q[145];
cx q[7], q[25];
cx q[43], q[59];
cx q[75], q[89];
cx q[103], q[115];
cx q[127], q[137];
cx q[147], q[155];
cx q[163], q[169];
cx q[27], q[44];
cx q[61], q[76];
cx q[91], q[104];
cx q[117], q[128];
cx q[139], q[148];
cx q[157], q[164];
cx q[171], q[176];
cx q[181], q[184];
cx q[10], q[28];
cx q[46], q[62];
cx q[78], q[92];
cx q[106], q[118];
cx q[130], q[140];
cx q[150], q[158];
cx q[166], q[172];
cx q[178], q[182];
cx q[186], q[188];
cx q[30], q[47];
cx q[64], q[79];
cx q[94], q[107];
cx q[120], q[131];
cx q[142], q[151];
cx q[160], q[167];
cx q[174], q[179];
cx q[13], q[31];
cx q[49], q[65];
cx q[81], q[95];
cx q[109], q[121];
cx q[133], q[143];
cx q[153], q[161];
cx q[33], q[50];
cx q[67], q[82];
cx q[97], q[110];
cx q[123], q[134];
cx q[16], q[34];
cx q[52], q[68];
cx q[84], q[98];
cx q[36], q[53];
barrier q;

cx q[1], q[2];
cx q[37], q[38];
cx q[21], q[22];
cx q[55], q[56];
cx q[85], q[86];
cx q[4], q[5];
cx q[40], q[41];
cx q[72], q[73];
cx q[100], q[101];
cx q[124], q[125];
cx q[24], q[25];
cx q[58], q[59];
cx q[88], q[89];
cx q[114], q[115];
cx q[136], q[137];
cx q[154], q[155];
cx q[7], q[8];
cx q[43], q[44];
cx q[75], q[76];
cx q[103], q[104];
cx q[127], q[128];
cx q[147], q[148];
cx q[163], q[164];
cx q[175], q[176];
cx q[27], q[28];
cx q[61], q[62];
cx q[91], q[92];
cx q[117], q[118];
cx q[139], q[140];
cx q[157], q[158];
cx q[171], q[172];
cx q[181], q[182];
cx q[187], q[188];
cx q[10], q[11];
cx q[46], q[47];
cx q[78], q[79];
cx q[106], q[107];
cx q[130], q[131];
cx q[150], q[151];
cx q[166], q[167];
cx q[178], q[179];
cx q[30], q[31];
cx q[64], q[65];
cx q[94], q[95];
cx q[120], q[121];
cx q[142], q[143];
cx q[160], q[161];
cx q[13], q[14];
cx q[49], q[50];
cx q[81], q[82];
cx q[109], q[110];
cx q[133], q[134];
cx q[33], q[34];
cx q[67], q[68];
cx q[97], q[98];
cx q[16], q[17];
cx q[52], q[53];
barrier q;

cx q[20], q[2];
cx q[54], q[38];
cx q[39], q[22];
cx q[71], q[56];
cx q[99], q[86];
cx q[23], q[5];
cx q[57], q[41];
cx q[87], q[73];
cx q[113], q[101];
cx q[135], q[125];
cx q[42], q[25];
cx q[74], q[59];
cx q[102], q[89];
cx q[126], q[115];
cx q[146], q[137];
cx q[162], q[155];
cx q[26], q[8];
cx q[60], q[44];
cx q[90], q[76];
cx q[116], q[104];
cx q[138], q[128];
cx q[156], q[148];
cx q[170], q[164];
cx q[180], q[176];
cx q[45], q[28];
cx q[77], q[62];
cx q[105], q[92];
cx q[129], q[118];
cx q[149], q[140];
cx q[165], q[158];
cx q[177], q[172];
cx q[185], q[182];
cx q[189], q[188];
cx q[29], q[11];
cx q[63], q[47];
cx q[93], q[79];
cx q[119], q[107];
cx q[141], q[131];
cx q[159], q[151];
cx q[173], q[167];
cx q[183], q[179];
cx q[48], q[31];
cx q[80], q[65];
cx q[108], q[95];
cx q[132], q[121];
cx q[152], q[143];
cx q[168], q[161];
cx q[32], q[14];
cx q[66], q[50];
cx q[96], q[82];
cx q[122], q[110];
cx q[144], q[134];
cx q[51], q[34];
cx q[83], q[68];
cx q[111], q[98];
cx q[35], q[17];
cx q[69], q[53];
barrier q;

cx q[0], q[19];
cx q[20], q[38];
cx q[54], q[70];
cx q[3], q[22];
cx q[39], q[56];
cx q[71], q[86];
cx q[99], q[112];
cx q[23], q[41];
cx q[57], q[73];
cx q[87], q[101];
cx q[113], q[125];
cx q[135], q[145];
cx q[6], q[25];
cx q[42], q[59];
cx q[74], q[89];
cx q[102], q[115];
cx q[126], q[137];
cx q[146], q[155];
cx q[162], q[169];
cx q[26], q[44];
cx q[60], q[76];
cx q[90], q[104];
cx q[116], q[128];
cx q[138], q[148];
cx q[156], q[164];
cx q[170], q[176];
cx q[180], q[184];
cx q[9], q[28];
cx q[45], q[62];
cx q[77], q[92];
cx q[105], q[118];
cx q[129], q[140];
cx q[149], q[158];
cx q[165], q[172];
cx q[177], q[182];
cx q[185], q[188];
cx q[29], q[47];
cx q[63], q[79];
cx q[93], q[107];
cx q[119], q[131];
cx q[141], q[151];
cx q[159], q[167];
cx q[173], q[179];
cx q[12], q[31];
cx q[48], q[65];
cx q[80], q[95];
cx q[108], q[121];
cx q[132], q[143];
cx q[152], q[161];
cx q[32], q[50];
cx q[66], q[82];
cx q[96], q[110];
cx q[122], q[134];
cx q[15], q[34];
cx q[51], q[68];
cx q[83], q[98];
cx q[35], q[53];
barrier q;

measure q[2] -> rec[0]; reset q[2]; // decomposed MR
measure q[5] -> rec[1]; reset q[5]; // decomposed MR
measure q[8] -> rec[2]; reset q[8]; // decomposed MR
measure q[11] -> rec[3]; reset q[11]; // decomposed MR
measure q[14] -> rec[4]; reset q[14]; // decomposed MR
measure q[17] -> rec[5]; reset q[17]; // decomposed MR
measure q[19] -> rec[6]; reset q[19]; // decomposed MR
measure q[22] -> rec[7]; reset q[22]; // decomposed MR
measure q[25] -> rec[8]; reset q[25]; // decomposed MR
measure q[28] -> rec[9]; reset q[28]; // decomposed MR
measure q[31] -> rec[10]; reset q[31]; // decomposed MR
measure q[34] -> rec[11]; reset q[34]; // decomposed MR
measure q[38] -> rec[12]; reset q[38]; // decomposed MR
measure q[41] -> rec[13]; reset q[41]; // decomposed MR
measure q[44] -> rec[14]; reset q[44]; // decomposed MR
measure q[47] -> rec[15]; reset q[47]; // decomposed MR
measure q[50] -> rec[16]; reset q[50]; // decomposed MR
measure q[53] -> rec[17]; reset q[53]; // decomposed MR
measure q[56] -> rec[18]; reset q[56]; // decomposed MR
measure q[59] -> rec[19]; reset q[59]; // decomposed MR
measure q[62] -> rec[20]; reset q[62]; // decomposed MR
measure q[65] -> rec[21]; reset q[65]; // decomposed MR
measure q[68] -> rec[22]; reset q[68]; // decomposed MR
measure q[70] -> rec[23]; reset q[70]; // decomposed MR
measure q[73] -> rec[24]; reset q[73]; // decomposed MR
measure q[76] -> rec[25]; reset q[76]; // decomposed MR
measure q[79] -> rec[26]; reset q[79]; // decomposed MR
measure q[82] -> rec[27]; reset q[82]; // decomposed MR
measure q[86] -> rec[28]; reset q[86]; // decomposed MR
measure q[89] -> rec[29]; reset q[89]; // decomposed MR
measure q[92] -> rec[30]; reset q[92]; // decomposed MR
measure q[95] -> rec[31]; reset q[95]; // decomposed MR
measure q[98] -> rec[32]; reset q[98]; // decomposed MR
measure q[101] -> rec[33]; reset q[101]; // decomposed MR
measure q[104] -> rec[34]; reset q[104]; // decomposed MR
measure q[107] -> rec[35]; reset q[107]; // decomposed MR
measure q[110] -> rec[36]; reset q[110]; // decomposed MR
measure q[112] -> rec[37]; reset q[112]; // decomposed MR
measure q[115] -> rec[38]; reset q[115]; // decomposed MR
measure q[118] -> rec[39]; reset q[118]; // decomposed MR
measure q[121] -> rec[40]; reset q[121]; // decomposed MR
measure q[125] -> rec[41]; reset q[125]; // decomposed MR
measure q[128] -> rec[42]; reset q[128]; // decomposed MR
measure q[131] -> rec[43]; reset q[131]; // decomposed MR
measure q[134] -> rec[44]; reset q[134]; // decomposed MR
measure q[137] -> rec[45]; reset q[137]; // decomposed MR
measure q[140] -> rec[46]; reset q[140]; // decomposed MR
measure q[143] -> rec[47]; reset q[143]; // decomposed MR
measure q[145] -> rec[48]; reset q[145]; // decomposed MR
measure q[148] -> rec[49]; reset q[148]; // decomposed MR
measure q[151] -> rec[50]; reset q[151]; // decomposed MR
measure q[155] -> rec[51]; reset q[155]; // decomposed MR
measure q[158] -> rec[52]; reset q[158]; // decomposed MR
measure q[161] -> rec[53]; reset q[161]; // decomposed MR
measure q[164] -> rec[54]; reset q[164]; // decomposed MR
measure q[167] -> rec[55]; reset q[167]; // decomposed MR
measure q[169] -> rec[56]; reset q[169]; // decomposed MR
measure q[172] -> rec[57]; reset q[172]; // decomposed MR
measure q[176] -> rec[58]; reset q[176]; // decomposed MR
measure q[179] -> rec[59]; reset q[179]; // decomposed MR
measure q[182] -> rec[60]; reset q[182]; // decomposed MR
measure q[184] -> rec[61]; reset q[184]; // decomposed MR
measure q[188] -> rec[62]; reset q[188]; // decomposed MR
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
cxyz q[16];
cxyz q[18];
cxyz q[20];
cxyz q[21];
cxyz q[23];
cxyz q[24];
cxyz q[26];
cxyz q[27];
cxyz q[29];
cxyz q[30];
cxyz q[32];
cxyz q[33];
cxyz q[35];
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
cxyz q[77];
cxyz q[78];
cxyz q[80];
cxyz q[81];
cxyz q[83];
cxyz q[84];
cxyz q[85];
cxyz q[87];
cxyz q[88];
cxyz q[90];
cxyz q[91];
cxyz q[93];
cxyz q[94];
cxyz q[96];
cxyz q[97];
cxyz q[99];
cxyz q[100];
cxyz q[102];
cxyz q[103];
cxyz q[105];
cxyz q[106];
cxyz q[108];
cxyz q[109];
cxyz q[111];
cxyz q[113];
cxyz q[114];
cxyz q[116];
cxyz q[117];
cxyz q[119];
cxyz q[120];
cxyz q[122];
cxyz q[123];
cxyz q[124];
cxyz q[126];
cxyz q[127];
cxyz q[129];
cxyz q[130];
cxyz q[132];
cxyz q[133];
cxyz q[135];
cxyz q[136];
cxyz q[138];
cxyz q[139];
cxyz q[141];
cxyz q[142];
cxyz q[144];
cxyz q[146];
cxyz q[147];
cxyz q[149];
cxyz q[150];
cxyz q[152];
cxyz q[153];
cxyz q[154];
cxyz q[156];
cxyz q[157];
cxyz q[159];
cxyz q[160];
cxyz q[162];
cxyz q[163];
cxyz q[165];
cxyz q[166];
cxyz q[168];
cxyz q[170];
cxyz q[171];
cxyz q[173];
cxyz q[174];
cxyz q[175];
cxyz q[177];
cxyz q[178];
cxyz q[180];
cxyz q[181];
cxyz q[183];
cxyz q[185];
cxyz q[186];
cxyz q[187];
cxyz q[189];
barrier q;

cx q[20], q[19];
cx q[3], q[2];
cx q[39], q[38];
cx q[71], q[70];
cx q[23], q[22];
cx q[57], q[56];
cx q[87], q[86];
cx q[113], q[112];
cx q[6], q[5];
cx q[42], q[41];
cx q[74], q[73];
cx q[102], q[101];
cx q[126], q[125];
cx q[146], q[145];
cx q[26], q[25];
cx q[60], q[59];
cx q[90], q[89];
cx q[116], q[115];
cx q[138], q[137];
cx q[156], q[155];
cx q[170], q[169];
cx q[9], q[8];
cx q[45], q[44];
cx q[77], q[76];
cx q[105], q[104];
cx q[129], q[128];
cx q[149], q[148];
cx q[165], q[164];
cx q[177], q[176];
cx q[185], q[184];
cx q[29], q[28];
cx q[63], q[62];
cx q[93], q[92];
cx q[119], q[118];
cx q[141], q[140];
cx q[159], q[158];
cx q[173], q[172];
cx q[183], q[182];
cx q[12], q[11];
cx q[48], q[47];
cx q[80], q[79];
cx q[108], q[107];
cx q[132], q[131];
cx q[152], q[151];
cx q[168], q[167];
cx q[32], q[31];
cx q[66], q[65];
cx q[96], q[95];
cx q[122], q[121];
cx q[144], q[143];
cx q[15], q[14];
cx q[51], q[50];
cx q[83], q[82];
cx q[111], q[110];
cx q[35], q[34];
cx q[69], q[68];
cx q[18], q[17];
barrier q;

cx q[37], q[19];
cx q[21], q[2];
cx q[55], q[38];
cx q[85], q[70];
cx q[40], q[22];
cx q[72], q[56];
cx q[100], q[86];
cx q[124], q[112];
cx q[24], q[5];
cx q[58], q[41];
cx q[88], q[73];
cx q[114], q[101];
cx q[136], q[125];
cx q[154], q[145];
cx q[43], q[25];
cx q[75], q[59];
cx q[103], q[89];
cx q[127], q[115];
cx q[147], q[137];
cx q[163], q[155];
cx q[175], q[169];
cx q[27], q[8];
cx q[61], q[44];
cx q[91], q[76];
cx q[117], q[104];
cx q[139], q[128];
cx q[157], q[148];
cx q[171], q[164];
cx q[181], q[176];
cx q[187], q[184];
cx q[46], q[28];
cx q[78], q[62];
cx q[106], q[92];
cx q[130], q[118];
cx q[150], q[140];
cx q[166], q[158];
cx q[178], q[172];
cx q[186], q[182];
cx q[30], q[11];
cx q[64], q[47];
cx q[94], q[79];
cx q[120], q[107];
cx q[142], q[131];
cx q[160], q[151];
cx q[174], q[167];
cx q[49], q[31];
cx q[81], q[65];
cx q[109], q[95];
cx q[133], q[121];
cx q[153], q[143];
cx q[33], q[14];
cx q[67], q[50];
cx q[97], q[82];
cx q[123], q[110];
cx q[52], q[34];
cx q[84], q[68];
cx q[36], q[17];
barrier q;

cx q[1], q[19];
cx q[21], q[38];
cx q[55], q[70];
cx q[4], q[22];
cx q[40], q[56];
cx q[72], q[86];
cx q[100], q[112];
cx q[24], q[41];
cx q[58], q[73];
cx q[88], q[101];
cx q[114], q[125];
cx q[136], q[145];
cx q[7], q[25];
cx q[43], q[59];
cx q[75], q[89];
cx q[103], q[115];
cx q[127], q[137];
cx q[147], q[155];
cx q[163], q[169];
cx q[27], q[44];
cx q[61], q[76];
cx q[91], q[104];
cx q[117], q[128];
cx q[139], q[148];
cx q[157], q[164];
cx q[171], q[176];
cx q[181], q[184];
cx q[10], q[28];
cx q[46], q[62];
cx q[78], q[92];
cx q[106], q[118];
cx q[130], q[140];
cx q[150], q[158];
cx q[166], q[172];
cx q[178], q[182];
cx q[186], q[188];
cx q[30], q[47];
cx q[64], q[79];
cx q[94], q[107];
cx q[120], q[131];
cx q[142], q[151];
cx q[160], q[167];
cx q[174], q[179];
cx q[13], q[31];
cx q[49], q[65];
cx q[81], q[95];
cx q[109], q[121];
cx q[133], q[143];
cx q[153], q[161];
cx q[33], q[50];
cx q[67], q[82];
cx q[97], q[110];
cx q[123], q[134];
cx q[16], q[34];
cx q[52], q[68];
cx q[84], q[98];
cx q[36], q[53];
barrier q;

cx q[1], q[2];
cx q[37], q[38];
cx q[21], q[22];
cx q[55], q[56];
cx q[85], q[86];
cx q[4], q[5];
cx q[40], q[41];
cx q[72], q[73];
cx q[100], q[101];
cx q[124], q[125];
cx q[24], q[25];
cx q[58], q[59];
cx q[88], q[89];
cx q[114], q[115];
cx q[136], q[137];
cx q[154], q[155];
cx q[7], q[8];
cx q[43], q[44];
cx q[75], q[76];
cx q[103], q[104];
cx q[127], q[128];
cx q[147], q[148];
cx q[163], q[164];
cx q[175], q[176];
cx q[27], q[28];
cx q[61], q[62];
cx q[91], q[92];
cx q[117], q[118];
cx q[139], q[140];
cx q[157], q[158];
cx q[171], q[172];
cx q[181], q[182];
cx q[187], q[188];
cx q[10], q[11];
cx q[46], q[47];
cx q[78], q[79];
cx q[106], q[107];
cx q[130], q[131];
cx q[150], q[151];
cx q[166], q[167];
cx q[178], q[179];
cx q[30], q[31];
cx q[64], q[65];
cx q[94], q[95];
cx q[120], q[121];
cx q[142], q[143];
cx q[160], q[161];
cx q[13], q[14];
cx q[49], q[50];
cx q[81], q[82];
cx q[109], q[110];
cx q[133], q[134];
cx q[33], q[34];
cx q[67], q[68];
cx q[97], q[98];
cx q[16], q[17];
cx q[52], q[53];
barrier q;

cx q[20], q[2];
cx q[54], q[38];
cx q[39], q[22];
cx q[71], q[56];
cx q[99], q[86];
cx q[23], q[5];
cx q[57], q[41];
cx q[87], q[73];
cx q[113], q[101];
cx q[135], q[125];
cx q[42], q[25];
cx q[74], q[59];
cx q[102], q[89];
cx q[126], q[115];
cx q[146], q[137];
cx q[162], q[155];
cx q[26], q[8];
cx q[60], q[44];
cx q[90], q[76];
cx q[116], q[104];
cx q[138], q[128];
cx q[156], q[148];
cx q[170], q[164];
cx q[180], q[176];
cx q[45], q[28];
cx q[77], q[62];
cx q[105], q[92];
cx q[129], q[118];
cx q[149], q[140];
cx q[165], q[158];
cx q[177], q[172];
cx q[185], q[182];
cx q[189], q[188];
cx q[29], q[11];
cx q[63], q[47];
cx q[93], q[79];
cx q[119], q[107];
cx q[141], q[131];
cx q[159], q[151];
cx q[173], q[167];
cx q[183], q[179];
cx q[48], q[31];
cx q[80], q[65];
cx q[108], q[95];
cx q[132], q[121];
cx q[152], q[143];
cx q[168], q[161];
cx q[32], q[14];
cx q[66], q[50];
cx q[96], q[82];
cx q[122], q[110];
cx q[144], q[134];
cx q[51], q[34];
cx q[83], q[68];
cx q[111], q[98];
cx q[35], q[17];
cx q[69], q[53];
barrier q;

cx q[0], q[19];
cx q[20], q[38];
cx q[54], q[70];
cx q[3], q[22];
cx q[39], q[56];
cx q[71], q[86];
cx q[99], q[112];
cx q[23], q[41];
cx q[57], q[73];
cx q[87], q[101];
cx q[113], q[125];
cx q[135], q[145];
cx q[6], q[25];
cx q[42], q[59];
cx q[74], q[89];
cx q[102], q[115];
cx q[126], q[137];
cx q[146], q[155];
cx q[162], q[169];
cx q[26], q[44];
cx q[60], q[76];
cx q[90], q[104];
cx q[116], q[128];
cx q[138], q[148];
cx q[156], q[164];
cx q[170], q[176];
cx q[180], q[184];
cx q[9], q[28];
cx q[45], q[62];
cx q[77], q[92];
cx q[105], q[118];
cx q[129], q[140];
cx q[149], q[158];
cx q[165], q[172];
cx q[177], q[182];
cx q[185], q[188];
cx q[29], q[47];
cx q[63], q[79];
cx q[93], q[107];
cx q[119], q[131];
cx q[141], q[151];
cx q[159], q[167];
cx q[173], q[179];
cx q[12], q[31];
cx q[48], q[65];
cx q[80], q[95];
cx q[108], q[121];
cx q[132], q[143];
cx q[152], q[161];
cx q[32], q[50];
cx q[66], q[82];
cx q[96], q[110];
cx q[122], q[134];
cx q[15], q[34];
cx q[51], q[68];
cx q[83], q[98];
cx q[35], q[53];
barrier q;

measure q[2] -> rec[63]; reset q[2]; // decomposed MR
measure q[5] -> rec[64]; reset q[5]; // decomposed MR
measure q[8] -> rec[65]; reset q[8]; // decomposed MR
measure q[11] -> rec[66]; reset q[11]; // decomposed MR
measure q[14] -> rec[67]; reset q[14]; // decomposed MR
measure q[17] -> rec[68]; reset q[17]; // decomposed MR
measure q[19] -> rec[69]; reset q[19]; // decomposed MR
measure q[22] -> rec[70]; reset q[22]; // decomposed MR
measure q[25] -> rec[71]; reset q[25]; // decomposed MR
measure q[28] -> rec[72]; reset q[28]; // decomposed MR
measure q[31] -> rec[73]; reset q[31]; // decomposed MR
measure q[34] -> rec[74]; reset q[34]; // decomposed MR
measure q[38] -> rec[75]; reset q[38]; // decomposed MR
measure q[41] -> rec[76]; reset q[41]; // decomposed MR
measure q[44] -> rec[77]; reset q[44]; // decomposed MR
measure q[47] -> rec[78]; reset q[47]; // decomposed MR
measure q[50] -> rec[79]; reset q[50]; // decomposed MR
measure q[53] -> rec[80]; reset q[53]; // decomposed MR
measure q[56] -> rec[81]; reset q[56]; // decomposed MR
measure q[59] -> rec[82]; reset q[59]; // decomposed MR
measure q[62] -> rec[83]; reset q[62]; // decomposed MR
measure q[65] -> rec[84]; reset q[65]; // decomposed MR
measure q[68] -> rec[85]; reset q[68]; // decomposed MR
measure q[70] -> rec[86]; reset q[70]; // decomposed MR
measure q[73] -> rec[87]; reset q[73]; // decomposed MR
measure q[76] -> rec[88]; reset q[76]; // decomposed MR
measure q[79] -> rec[89]; reset q[79]; // decomposed MR
measure q[82] -> rec[90]; reset q[82]; // decomposed MR
measure q[86] -> rec[91]; reset q[86]; // decomposed MR
measure q[89] -> rec[92]; reset q[89]; // decomposed MR
measure q[92] -> rec[93]; reset q[92]; // decomposed MR
measure q[95] -> rec[94]; reset q[95]; // decomposed MR
measure q[98] -> rec[95]; reset q[98]; // decomposed MR
measure q[101] -> rec[96]; reset q[101]; // decomposed MR
measure q[104] -> rec[97]; reset q[104]; // decomposed MR
measure q[107] -> rec[98]; reset q[107]; // decomposed MR
measure q[110] -> rec[99]; reset q[110]; // decomposed MR
measure q[112] -> rec[100]; reset q[112]; // decomposed MR
measure q[115] -> rec[101]; reset q[115]; // decomposed MR
measure q[118] -> rec[102]; reset q[118]; // decomposed MR
measure q[121] -> rec[103]; reset q[121]; // decomposed MR
measure q[125] -> rec[104]; reset q[125]; // decomposed MR
measure q[128] -> rec[105]; reset q[128]; // decomposed MR
measure q[131] -> rec[106]; reset q[131]; // decomposed MR
measure q[134] -> rec[107]; reset q[134]; // decomposed MR
measure q[137] -> rec[108]; reset q[137]; // decomposed MR
measure q[140] -> rec[109]; reset q[140]; // decomposed MR
measure q[143] -> rec[110]; reset q[143]; // decomposed MR
measure q[145] -> rec[111]; reset q[145]; // decomposed MR
measure q[148] -> rec[112]; reset q[148]; // decomposed MR
measure q[151] -> rec[113]; reset q[151]; // decomposed MR
measure q[155] -> rec[114]; reset q[155]; // decomposed MR
measure q[158] -> rec[115]; reset q[158]; // decomposed MR
measure q[161] -> rec[116]; reset q[161]; // decomposed MR
measure q[164] -> rec[117]; reset q[164]; // decomposed MR
measure q[167] -> rec[118]; reset q[167]; // decomposed MR
measure q[169] -> rec[119]; reset q[169]; // decomposed MR
measure q[172] -> rec[120]; reset q[172]; // decomposed MR
measure q[176] -> rec[121]; reset q[176]; // decomposed MR
measure q[179] -> rec[122]; reset q[179]; // decomposed MR
measure q[182] -> rec[123]; reset q[182]; // decomposed MR
measure q[184] -> rec[124]; reset q[184]; // decomposed MR
measure q[188] -> rec[125]; reset q[188]; // decomposed MR
s q[0]; s q[0]; s q[0]; h q[0]; measure q[0] -> rec[126]; h q[0]; s q[0]; // decomposed MY
s q[1]; s q[1]; s q[1]; h q[1]; measure q[1] -> rec[127]; h q[1]; s q[1]; // decomposed MY
s q[3]; s q[3]; s q[3]; h q[3]; measure q[3] -> rec[128]; h q[3]; s q[3]; // decomposed MY
s q[4]; s q[4]; s q[4]; h q[4]; measure q[4] -> rec[129]; h q[4]; s q[4]; // decomposed MY
s q[6]; s q[6]; s q[6]; h q[6]; measure q[6] -> rec[130]; h q[6]; s q[6]; // decomposed MY
s q[7]; s q[7]; s q[7]; h q[7]; measure q[7] -> rec[131]; h q[7]; s q[7]; // decomposed MY
s q[9]; s q[9]; s q[9]; h q[9]; measure q[9] -> rec[132]; h q[9]; s q[9]; // decomposed MY
s q[10]; s q[10]; s q[10]; h q[10]; measure q[10] -> rec[133]; h q[10]; s q[10]; // decomposed MY
s q[12]; s q[12]; s q[12]; h q[12]; measure q[12] -> rec[134]; h q[12]; s q[12]; // decomposed MY
s q[13]; s q[13]; s q[13]; h q[13]; measure q[13] -> rec[135]; h q[13]; s q[13]; // decomposed MY
s q[15]; s q[15]; s q[15]; h q[15]; measure q[15] -> rec[136]; h q[15]; s q[15]; // decomposed MY
s q[16]; s q[16]; s q[16]; h q[16]; measure q[16] -> rec[137]; h q[16]; s q[16]; // decomposed MY
s q[18]; s q[18]; s q[18]; h q[18]; measure q[18] -> rec[138]; h q[18]; s q[18]; // decomposed MY
s q[20]; s q[20]; s q[20]; h q[20]; measure q[20] -> rec[139]; h q[20]; s q[20]; // decomposed MY
s q[21]; s q[21]; s q[21]; h q[21]; measure q[21] -> rec[140]; h q[21]; s q[21]; // decomposed MY
s q[23]; s q[23]; s q[23]; h q[23]; measure q[23] -> rec[141]; h q[23]; s q[23]; // decomposed MY
s q[24]; s q[24]; s q[24]; h q[24]; measure q[24] -> rec[142]; h q[24]; s q[24]; // decomposed MY
s q[26]; s q[26]; s q[26]; h q[26]; measure q[26] -> rec[143]; h q[26]; s q[26]; // decomposed MY
s q[27]; s q[27]; s q[27]; h q[27]; measure q[27] -> rec[144]; h q[27]; s q[27]; // decomposed MY
s q[29]; s q[29]; s q[29]; h q[29]; measure q[29] -> rec[145]; h q[29]; s q[29]; // decomposed MY
s q[30]; s q[30]; s q[30]; h q[30]; measure q[30] -> rec[146]; h q[30]; s q[30]; // decomposed MY
s q[32]; s q[32]; s q[32]; h q[32]; measure q[32] -> rec[147]; h q[32]; s q[32]; // decomposed MY
s q[33]; s q[33]; s q[33]; h q[33]; measure q[33] -> rec[148]; h q[33]; s q[33]; // decomposed MY
s q[35]; s q[35]; s q[35]; h q[35]; measure q[35] -> rec[149]; h q[35]; s q[35]; // decomposed MY
s q[36]; s q[36]; s q[36]; h q[36]; measure q[36] -> rec[150]; h q[36]; s q[36]; // decomposed MY
s q[37]; s q[37]; s q[37]; h q[37]; measure q[37] -> rec[151]; h q[37]; s q[37]; // decomposed MY
s q[39]; s q[39]; s q[39]; h q[39]; measure q[39] -> rec[152]; h q[39]; s q[39]; // decomposed MY
s q[40]; s q[40]; s q[40]; h q[40]; measure q[40] -> rec[153]; h q[40]; s q[40]; // decomposed MY
s q[42]; s q[42]; s q[42]; h q[42]; measure q[42] -> rec[154]; h q[42]; s q[42]; // decomposed MY
s q[43]; s q[43]; s q[43]; h q[43]; measure q[43] -> rec[155]; h q[43]; s q[43]; // decomposed MY
s q[45]; s q[45]; s q[45]; h q[45]; measure q[45] -> rec[156]; h q[45]; s q[45]; // decomposed MY
s q[46]; s q[46]; s q[46]; h q[46]; measure q[46] -> rec[157]; h q[46]; s q[46]; // decomposed MY
s q[48]; s q[48]; s q[48]; h q[48]; measure q[48] -> rec[158]; h q[48]; s q[48]; // decomposed MY
s q[49]; s q[49]; s q[49]; h q[49]; measure q[49] -> rec[159]; h q[49]; s q[49]; // decomposed MY
s q[51]; s q[51]; s q[51]; h q[51]; measure q[51] -> rec[160]; h q[51]; s q[51]; // decomposed MY
s q[52]; s q[52]; s q[52]; h q[52]; measure q[52] -> rec[161]; h q[52]; s q[52]; // decomposed MY
s q[54]; s q[54]; s q[54]; h q[54]; measure q[54] -> rec[162]; h q[54]; s q[54]; // decomposed MY
s q[55]; s q[55]; s q[55]; h q[55]; measure q[55] -> rec[163]; h q[55]; s q[55]; // decomposed MY
s q[57]; s q[57]; s q[57]; h q[57]; measure q[57] -> rec[164]; h q[57]; s q[57]; // decomposed MY
s q[58]; s q[58]; s q[58]; h q[58]; measure q[58] -> rec[165]; h q[58]; s q[58]; // decomposed MY
s q[60]; s q[60]; s q[60]; h q[60]; measure q[60] -> rec[166]; h q[60]; s q[60]; // decomposed MY
s q[61]; s q[61]; s q[61]; h q[61]; measure q[61] -> rec[167]; h q[61]; s q[61]; // decomposed MY
s q[63]; s q[63]; s q[63]; h q[63]; measure q[63] -> rec[168]; h q[63]; s q[63]; // decomposed MY
s q[64]; s q[64]; s q[64]; h q[64]; measure q[64] -> rec[169]; h q[64]; s q[64]; // decomposed MY
s q[66]; s q[66]; s q[66]; h q[66]; measure q[66] -> rec[170]; h q[66]; s q[66]; // decomposed MY
s q[67]; s q[67]; s q[67]; h q[67]; measure q[67] -> rec[171]; h q[67]; s q[67]; // decomposed MY
s q[69]; s q[69]; s q[69]; h q[69]; measure q[69] -> rec[172]; h q[69]; s q[69]; // decomposed MY
s q[71]; s q[71]; s q[71]; h q[71]; measure q[71] -> rec[173]; h q[71]; s q[71]; // decomposed MY
s q[72]; s q[72]; s q[72]; h q[72]; measure q[72] -> rec[174]; h q[72]; s q[72]; // decomposed MY
s q[74]; s q[74]; s q[74]; h q[74]; measure q[74] -> rec[175]; h q[74]; s q[74]; // decomposed MY
s q[75]; s q[75]; s q[75]; h q[75]; measure q[75] -> rec[176]; h q[75]; s q[75]; // decomposed MY
s q[77]; s q[77]; s q[77]; h q[77]; measure q[77] -> rec[177]; h q[77]; s q[77]; // decomposed MY
s q[78]; s q[78]; s q[78]; h q[78]; measure q[78] -> rec[178]; h q[78]; s q[78]; // decomposed MY
s q[80]; s q[80]; s q[80]; h q[80]; measure q[80] -> rec[179]; h q[80]; s q[80]; // decomposed MY
s q[81]; s q[81]; s q[81]; h q[81]; measure q[81] -> rec[180]; h q[81]; s q[81]; // decomposed MY
s q[83]; s q[83]; s q[83]; h q[83]; measure q[83] -> rec[181]; h q[83]; s q[83]; // decomposed MY
s q[84]; s q[84]; s q[84]; h q[84]; measure q[84] -> rec[182]; h q[84]; s q[84]; // decomposed MY
s q[85]; s q[85]; s q[85]; h q[85]; measure q[85] -> rec[183]; h q[85]; s q[85]; // decomposed MY
s q[87]; s q[87]; s q[87]; h q[87]; measure q[87] -> rec[184]; h q[87]; s q[87]; // decomposed MY
s q[88]; s q[88]; s q[88]; h q[88]; measure q[88] -> rec[185]; h q[88]; s q[88]; // decomposed MY
s q[90]; s q[90]; s q[90]; h q[90]; measure q[90] -> rec[186]; h q[90]; s q[90]; // decomposed MY
s q[91]; s q[91]; s q[91]; h q[91]; measure q[91] -> rec[187]; h q[91]; s q[91]; // decomposed MY
s q[93]; s q[93]; s q[93]; h q[93]; measure q[93] -> rec[188]; h q[93]; s q[93]; // decomposed MY
s q[94]; s q[94]; s q[94]; h q[94]; measure q[94] -> rec[189]; h q[94]; s q[94]; // decomposed MY
s q[96]; s q[96]; s q[96]; h q[96]; measure q[96] -> rec[190]; h q[96]; s q[96]; // decomposed MY
s q[97]; s q[97]; s q[97]; h q[97]; measure q[97] -> rec[191]; h q[97]; s q[97]; // decomposed MY
s q[99]; s q[99]; s q[99]; h q[99]; measure q[99] -> rec[192]; h q[99]; s q[99]; // decomposed MY
s q[100]; s q[100]; s q[100]; h q[100]; measure q[100] -> rec[193]; h q[100]; s q[100]; // decomposed MY
s q[102]; s q[102]; s q[102]; h q[102]; measure q[102] -> rec[194]; h q[102]; s q[102]; // decomposed MY
s q[103]; s q[103]; s q[103]; h q[103]; measure q[103] -> rec[195]; h q[103]; s q[103]; // decomposed MY
s q[105]; s q[105]; s q[105]; h q[105]; measure q[105] -> rec[196]; h q[105]; s q[105]; // decomposed MY
s q[106]; s q[106]; s q[106]; h q[106]; measure q[106] -> rec[197]; h q[106]; s q[106]; // decomposed MY
s q[108]; s q[108]; s q[108]; h q[108]; measure q[108] -> rec[198]; h q[108]; s q[108]; // decomposed MY
s q[109]; s q[109]; s q[109]; h q[109]; measure q[109] -> rec[199]; h q[109]; s q[109]; // decomposed MY
s q[111]; s q[111]; s q[111]; h q[111]; measure q[111] -> rec[200]; h q[111]; s q[111]; // decomposed MY
s q[113]; s q[113]; s q[113]; h q[113]; measure q[113] -> rec[201]; h q[113]; s q[113]; // decomposed MY
s q[114]; s q[114]; s q[114]; h q[114]; measure q[114] -> rec[202]; h q[114]; s q[114]; // decomposed MY
s q[116]; s q[116]; s q[116]; h q[116]; measure q[116] -> rec[203]; h q[116]; s q[116]; // decomposed MY
s q[117]; s q[117]; s q[117]; h q[117]; measure q[117] -> rec[204]; h q[117]; s q[117]; // decomposed MY
s q[119]; s q[119]; s q[119]; h q[119]; measure q[119] -> rec[205]; h q[119]; s q[119]; // decomposed MY
s q[120]; s q[120]; s q[120]; h q[120]; measure q[120] -> rec[206]; h q[120]; s q[120]; // decomposed MY
s q[122]; s q[122]; s q[122]; h q[122]; measure q[122] -> rec[207]; h q[122]; s q[122]; // decomposed MY
s q[123]; s q[123]; s q[123]; h q[123]; measure q[123] -> rec[208]; h q[123]; s q[123]; // decomposed MY
s q[124]; s q[124]; s q[124]; h q[124]; measure q[124] -> rec[209]; h q[124]; s q[124]; // decomposed MY
s q[126]; s q[126]; s q[126]; h q[126]; measure q[126] -> rec[210]; h q[126]; s q[126]; // decomposed MY
s q[127]; s q[127]; s q[127]; h q[127]; measure q[127] -> rec[211]; h q[127]; s q[127]; // decomposed MY
s q[129]; s q[129]; s q[129]; h q[129]; measure q[129] -> rec[212]; h q[129]; s q[129]; // decomposed MY
s q[130]; s q[130]; s q[130]; h q[130]; measure q[130] -> rec[213]; h q[130]; s q[130]; // decomposed MY
s q[132]; s q[132]; s q[132]; h q[132]; measure q[132] -> rec[214]; h q[132]; s q[132]; // decomposed MY
s q[133]; s q[133]; s q[133]; h q[133]; measure q[133] -> rec[215]; h q[133]; s q[133]; // decomposed MY
s q[135]; s q[135]; s q[135]; h q[135]; measure q[135] -> rec[216]; h q[135]; s q[135]; // decomposed MY
s q[136]; s q[136]; s q[136]; h q[136]; measure q[136] -> rec[217]; h q[136]; s q[136]; // decomposed MY
s q[138]; s q[138]; s q[138]; h q[138]; measure q[138] -> rec[218]; h q[138]; s q[138]; // decomposed MY
s q[139]; s q[139]; s q[139]; h q[139]; measure q[139] -> rec[219]; h q[139]; s q[139]; // decomposed MY
s q[141]; s q[141]; s q[141]; h q[141]; measure q[141] -> rec[220]; h q[141]; s q[141]; // decomposed MY
s q[142]; s q[142]; s q[142]; h q[142]; measure q[142] -> rec[221]; h q[142]; s q[142]; // decomposed MY
s q[144]; s q[144]; s q[144]; h q[144]; measure q[144] -> rec[222]; h q[144]; s q[144]; // decomposed MY
s q[146]; s q[146]; s q[146]; h q[146]; measure q[146] -> rec[223]; h q[146]; s q[146]; // decomposed MY
s q[147]; s q[147]; s q[147]; h q[147]; measure q[147] -> rec[224]; h q[147]; s q[147]; // decomposed MY
s q[149]; s q[149]; s q[149]; h q[149]; measure q[149] -> rec[225]; h q[149]; s q[149]; // decomposed MY
s q[150]; s q[150]; s q[150]; h q[150]; measure q[150] -> rec[226]; h q[150]; s q[150]; // decomposed MY
s q[152]; s q[152]; s q[152]; h q[152]; measure q[152] -> rec[227]; h q[152]; s q[152]; // decomposed MY
s q[153]; s q[153]; s q[153]; h q[153]; measure q[153] -> rec[228]; h q[153]; s q[153]; // decomposed MY
s q[154]; s q[154]; s q[154]; h q[154]; measure q[154] -> rec[229]; h q[154]; s q[154]; // decomposed MY
s q[156]; s q[156]; s q[156]; h q[156]; measure q[156] -> rec[230]; h q[156]; s q[156]; // decomposed MY
s q[157]; s q[157]; s q[157]; h q[157]; measure q[157] -> rec[231]; h q[157]; s q[157]; // decomposed MY
s q[159]; s q[159]; s q[159]; h q[159]; measure q[159] -> rec[232]; h q[159]; s q[159]; // decomposed MY
s q[160]; s q[160]; s q[160]; h q[160]; measure q[160] -> rec[233]; h q[160]; s q[160]; // decomposed MY
s q[162]; s q[162]; s q[162]; h q[162]; measure q[162] -> rec[234]; h q[162]; s q[162]; // decomposed MY
s q[163]; s q[163]; s q[163]; h q[163]; measure q[163] -> rec[235]; h q[163]; s q[163]; // decomposed MY
s q[165]; s q[165]; s q[165]; h q[165]; measure q[165] -> rec[236]; h q[165]; s q[165]; // decomposed MY
s q[166]; s q[166]; s q[166]; h q[166]; measure q[166] -> rec[237]; h q[166]; s q[166]; // decomposed MY
s q[168]; s q[168]; s q[168]; h q[168]; measure q[168] -> rec[238]; h q[168]; s q[168]; // decomposed MY
s q[170]; s q[170]; s q[170]; h q[170]; measure q[170] -> rec[239]; h q[170]; s q[170]; // decomposed MY
s q[171]; s q[171]; s q[171]; h q[171]; measure q[171] -> rec[240]; h q[171]; s q[171]; // decomposed MY
s q[173]; s q[173]; s q[173]; h q[173]; measure q[173] -> rec[241]; h q[173]; s q[173]; // decomposed MY
s q[174]; s q[174]; s q[174]; h q[174]; measure q[174] -> rec[242]; h q[174]; s q[174]; // decomposed MY
s q[175]; s q[175]; s q[175]; h q[175]; measure q[175] -> rec[243]; h q[175]; s q[175]; // decomposed MY
s q[177]; s q[177]; s q[177]; h q[177]; measure q[177] -> rec[244]; h q[177]; s q[177]; // decomposed MY
s q[178]; s q[178]; s q[178]; h q[178]; measure q[178] -> rec[245]; h q[178]; s q[178]; // decomposed MY
s q[180]; s q[180]; s q[180]; h q[180]; measure q[180] -> rec[246]; h q[180]; s q[180]; // decomposed MY
s q[181]; s q[181]; s q[181]; h q[181]; measure q[181] -> rec[247]; h q[181]; s q[181]; // decomposed MY
s q[183]; s q[183]; s q[183]; h q[183]; measure q[183] -> rec[248]; h q[183]; s q[183]; // decomposed MY
s q[185]; s q[185]; s q[185]; h q[185]; measure q[185] -> rec[249]; h q[185]; s q[185]; // decomposed MY
s q[186]; s q[186]; s q[186]; h q[186]; measure q[186] -> rec[250]; h q[186]; s q[186]; // decomposed MY
s q[187]; s q[187]; s q[187]; h q[187]; measure q[187] -> rec[251]; h q[187]; s q[187]; // decomposed MY
s q[189]; s q[189]; s q[189]; h q[189]; measure q[189] -> rec[252]; h q[189]; s q[189]; // decomposed MY

OPENQASM 2.0;
include "qelib1.inc";
gate cxyz q0 { U(pi/2, 0, pi/2) q0; }

qreg q[253];
creg rec[337];

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
reset q[190];
reset q[191];
reset q[192];
reset q[193];
reset q[194];
reset q[195];
reset q[196];
reset q[197];
reset q[198];
reset q[199];
reset q[200];
reset q[201];
reset q[202];
reset q[203];
reset q[204];
reset q[205];
reset q[206];
reset q[207];
reset q[208];
reset q[209];
reset q[210];
reset q[211];
reset q[212];
reset q[213];
reset q[214];
reset q[215];
reset q[216];
reset q[217];
reset q[218];
reset q[219];
reset q[220];
reset q[221];
reset q[222];
reset q[223];
reset q[224];
reset q[225];
reset q[226];
reset q[227];
reset q[228];
reset q[229];
reset q[230];
reset q[231];
reset q[232];
reset q[233];
reset q[234];
reset q[235];
reset q[236];
reset q[237];
reset q[238];
reset q[239];
reset q[240];
reset q[241];
reset q[242];
reset q[243];
reset q[244];
reset q[245];
reset q[246];
reset q[247];
reset q[248];
reset q[249];
reset q[250];
reset q[251];
reset q[252];
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
cxyz q[19];
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
cxyz q[38];
cxyz q[39];
cxyz q[41];
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
cxyz q[70];
cxyz q[72];
cxyz q[73];
cxyz q[75];
cxyz q[76];
cxyz q[78];
cxyz q[79];
cxyz q[81];
cxyz q[83];
cxyz q[84];
cxyz q[86];
cxyz q[87];
cxyz q[89];
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
cxyz q[115];
cxyz q[117];
cxyz q[118];
cxyz q[120];
cxyz q[121];
cxyz q[123];
cxyz q[124];
cxyz q[126];
cxyz q[127];
cxyz q[129];
cxyz q[130];
cxyz q[132];
cxyz q[134];
cxyz q[135];
cxyz q[137];
cxyz q[138];
cxyz q[140];
cxyz q[141];
cxyz q[143];
cxyz q[144];
cxyz q[146];
cxyz q[147];
cxyz q[148];
cxyz q[150];
cxyz q[151];
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
cxyz q[169];
cxyz q[171];
cxyz q[172];
cxyz q[174];
cxyz q[176];
cxyz q[177];
cxyz q[179];
cxyz q[180];
cxyz q[182];
cxyz q[183];
cxyz q[185];
cxyz q[186];
cxyz q[187];
cxyz q[189];
cxyz q[190];
cxyz q[192];
cxyz q[193];
cxyz q[195];
cxyz q[196];
cxyz q[198];
cxyz q[199];
cxyz q[201];
cxyz q[202];
cxyz q[204];
cxyz q[205];
cxyz q[207];
cxyz q[209];
cxyz q[210];
cxyz q[212];
cxyz q[213];
cxyz q[215];
cxyz q[216];
cxyz q[217];
cxyz q[219];
cxyz q[220];
cxyz q[222];
cxyz q[223];
cxyz q[225];
cxyz q[226];
cxyz q[228];
cxyz q[229];
cxyz q[231];
cxyz q[233];
cxyz q[234];
cxyz q[236];
cxyz q[237];
cxyz q[238];
cxyz q[240];
cxyz q[241];
cxyz q[243];
cxyz q[244];
cxyz q[246];
cxyz q[248];
cxyz q[249];
cxyz q[250];
cxyz q[252];
barrier q;

cx q[23], q[22];
cx q[3], q[2];
cx q[45], q[44];
cx q[83], q[82];
cx q[26], q[25];
cx q[66], q[65];
cx q[102], q[101];
cx q[134], q[133];
cx q[6], q[5];
cx q[48], q[47];
cx q[86], q[85];
cx q[120], q[119];
cx q[150], q[149];
cx q[176], q[175];
cx q[29], q[28];
cx q[69], q[68];
cx q[105], q[104];
cx q[137], q[136];
cx q[165], q[164];
cx q[189], q[188];
cx q[209], q[208];
cx q[9], q[8];
cx q[51], q[50];
cx q[89], q[88];
cx q[123], q[122];
cx q[153], q[152];
cx q[179], q[178];
cx q[201], q[200];
cx q[219], q[218];
cx q[233], q[232];
cx q[32], q[31];
cx q[72], q[71];
cx q[108], q[107];
cx q[140], q[139];
cx q[168], q[167];
cx q[192], q[191];
cx q[212], q[211];
cx q[228], q[227];
cx q[240], q[239];
cx q[248], q[247];
cx q[12], q[11];
cx q[54], q[53];
cx q[92], q[91];
cx q[126], q[125];
cx q[156], q[155];
cx q[182], q[181];
cx q[204], q[203];
cx q[222], q[221];
cx q[236], q[235];
cx q[246], q[245];
cx q[35], q[34];
cx q[75], q[74];
cx q[111], q[110];
cx q[143], q[142];
cx q[171], q[170];
cx q[195], q[194];
cx q[215], q[214];
cx q[231], q[230];
cx q[15], q[14];
cx q[57], q[56];
cx q[95], q[94];
cx q[129], q[128];
cx q[159], q[158];
cx q[185], q[184];
cx q[207], q[206];
cx q[38], q[37];
cx q[78], q[77];
cx q[114], q[113];
cx q[146], q[145];
cx q[174], q[173];
cx q[18], q[17];
cx q[60], q[59];
cx q[98], q[97];
cx q[132], q[131];
cx q[41], q[40];
cx q[81], q[80];
cx q[21], q[20];
barrier q;

cx q[43], q[22];
cx q[24], q[2];
cx q[64], q[44];
cx q[100], q[82];
cx q[46], q[25];
cx q[84], q[65];
cx q[118], q[101];
cx q[148], q[133];
cx q[27], q[5];
cx q[67], q[47];
cx q[103], q[85];
cx q[135], q[119];
cx q[163], q[149];
cx q[187], q[175];
cx q[49], q[28];
cx q[87], q[68];
cx q[121], q[104];
cx q[151], q[136];
cx q[177], q[164];
cx q[199], q[188];
cx q[217], q[208];
cx q[30], q[8];
cx q[70], q[50];
cx q[106], q[88];
cx q[138], q[122];
cx q[166], q[152];
cx q[190], q[178];
cx q[210], q[200];
cx q[226], q[218];
cx q[238], q[232];
cx q[52], q[31];
cx q[90], q[71];
cx q[124], q[107];
cx q[154], q[139];
cx q[180], q[167];
cx q[202], q[191];
cx q[220], q[211];
cx q[234], q[227];
cx q[244], q[239];
cx q[250], q[247];
cx q[33], q[11];
cx q[73], q[53];
cx q[109], q[91];
cx q[141], q[125];
cx q[169], q[155];
cx q[193], q[181];
cx q[213], q[203];
cx q[229], q[221];
cx q[241], q[235];
cx q[249], q[245];
cx q[55], q[34];
cx q[93], q[74];
cx q[127], q[110];
cx q[157], q[142];
cx q[183], q[170];
cx q[205], q[194];
cx q[223], q[214];
cx q[237], q[230];
cx q[36], q[14];
cx q[76], q[56];
cx q[112], q[94];
cx q[144], q[128];
cx q[172], q[158];
cx q[196], q[184];
cx q[216], q[206];
cx q[58], q[37];
cx q[96], q[77];
cx q[130], q[113];
cx q[160], q[145];
cx q[186], q[173];
cx q[39], q[17];
cx q[79], q[59];
cx q[115], q[97];
cx q[147], q[131];
cx q[61], q[40];
cx q[99], q[80];
cx q[42], q[20];
barrier q;

cx q[1], q[22];
cx q[24], q[44];
cx q[64], q[82];
cx q[4], q[25];
cx q[46], q[65];
cx q[84], q[101];
cx q[118], q[133];
cx q[27], q[47];
cx q[67], q[85];
cx q[103], q[119];
cx q[135], q[149];
cx q[163], q[175];
cx q[7], q[28];
cx q[49], q[68];
cx q[87], q[104];
cx q[121], q[136];
cx q[151], q[164];
cx q[177], q[188];
cx q[199], q[208];
cx q[30], q[50];
cx q[70], q[88];
cx q[106], q[122];
cx q[138], q[152];
cx q[166], q[178];
cx q[190], q[200];
cx q[210], q[218];
cx q[226], q[232];
cx q[10], q[31];
cx q[52], q[71];
cx q[90], q[107];
cx q[124], q[139];
cx q[154], q[167];
cx q[180], q[191];
cx q[202], q[211];
cx q[220], q[227];
cx q[234], q[239];
cx q[244], q[247];
cx q[33], q[53];
cx q[73], q[91];
cx q[109], q[125];
cx q[141], q[155];
cx q[169], q[181];
cx q[193], q[203];
cx q[213], q[221];
cx q[229], q[235];
cx q[241], q[245];
cx q[249], q[251];
cx q[13], q[34];
cx q[55], q[74];
cx q[93], q[110];
cx q[127], q[142];
cx q[157], q[170];
cx q[183], q[194];
cx q[205], q[214];
cx q[223], q[230];
cx q[237], q[242];
cx q[36], q[56];
cx q[76], q[94];
cx q[112], q[128];
cx q[144], q[158];
cx q[172], q[184];
cx q[196], q[206];
cx q[216], q[224];
cx q[16], q[37];
cx q[58], q[77];
cx q[96], q[113];
cx q[130], q[145];
cx q[160], q[173];
cx q[186], q[197];
cx q[39], q[59];
cx q[79], q[97];
cx q[115], q[131];
cx q[147], q[161];
cx q[19], q[40];
cx q[61], q[80];
cx q[99], q[116];
cx q[42], q[62];
barrier q;

cx q[1], q[2];
cx q[43], q[44];
cx q[24], q[25];
cx q[64], q[65];
cx q[100], q[101];
cx q[4], q[5];
cx q[46], q[47];
cx q[84], q[85];
cx q[118], q[119];
cx q[148], q[149];
cx q[27], q[28];
cx q[67], q[68];
cx q[103], q[104];
cx q[135], q[136];
cx q[163], q[164];
cx q[187], q[188];
cx q[7], q[8];
cx q[49], q[50];
cx q[87], q[88];
cx q[121], q[122];
cx q[151], q[152];
cx q[177], q[178];
cx q[199], q[200];
cx q[217], q[218];
cx q[30], q[31];
cx q[70], q[71];
cx q[106], q[107];
cx q[138], q[139];
cx q[166], q[167];
cx q[190], q[191];
cx q[210], q[211];
cx q[226], q[227];
cx q[238], q[239];
cx q[10], q[11];
cx q[52], q[53];
cx q[90], q[91];
cx q[124], q[125];
cx q[154], q[155];
cx q[180], q[181];
cx q[202], q[203];
cx q[220], q[221];
cx q[234], q[235];
cx q[244], q[245];
cx q[250], q[251];
cx q[33], q[34];
cx q[73], q[74];
cx q[109], q[110];
cx q[141], q[142];
cx q[169], q[170];
cx q[193], q[194];
cx q[213], q[214];
cx q[229], q[230];
cx q[241], q[242];
cx q[13], q[14];
cx q[55], q[56];
cx q[93], q[94];
cx q[127], q[128];
cx q[157], q[158];
cx q[183], q[184];
cx q[205], q[206];
cx q[223], q[224];
cx q[36], q[37];
cx q[76], q[77];
cx q[112], q[113];
cx q[144], q[145];
cx q[172], q[173];
cx q[196], q[197];
cx q[16], q[17];
cx q[58], q[59];
cx q[96], q[97];
cx q[130], q[131];
cx q[160], q[161];
cx q[39], q[40];
cx q[79], q[80];
cx q[115], q[116];
cx q[19], q[20];
cx q[61], q[62];
barrier q;

cx q[23], q[2];
cx q[63], q[44];
cx q[45], q[25];
cx q[83], q[65];
cx q[117], q[101];
cx q[26], q[5];
cx q[66], q[47];
cx q[102], q[85];
cx q[134], q[119];
cx q[162], q[149];
cx q[48], q[28];
cx q[86], q[68];
cx q[120], q[104];
cx q[150], q[136];
cx q[176], q[164];
cx q[198], q[188];
cx q[29], q[8];
cx q[69], q[50];
cx q[105], q[88];
cx q[137], q[122];
cx q[165], q[152];
cx q[189], q[178];
cx q[209], q[200];
cx q[225], q[218];
cx q[51], q[31];
cx q[89], q[71];
cx q[123], q[107];
cx q[153], q[139];
cx q[179], q[167];
cx q[201], q[191];
cx q[219], q[211];
cx q[233], q[227];
cx q[243], q[239];
cx q[32], q[11];
cx q[72], q[53];
cx q[108], q[91];
cx q[140], q[125];
cx q[168], q[155];
cx q[192], q[181];
cx q[212], q[203];
cx q[228], q[221];
cx q[240], q[235];
cx q[248], q[245];
cx q[252], q[251];
cx q[54], q[34];
cx q[92], q[74];
cx q[126], q[110];
cx q[156], q[142];
cx q[182], q[170];
cx q[204], q[194];
cx q[222], q[214];
cx q[236], q[230];
cx q[246], q[242];
cx q[35], q[14];
cx q[75], q[56];
cx q[111], q[94];
cx q[143], q[128];
cx q[171], q[158];
cx q[195], q[184];
cx q[215], q[206];
cx q[231], q[224];
cx q[57], q[37];
cx q[95], q[77];
cx q[129], q[113];
cx q[159], q[145];
cx q[185], q[173];
cx q[207], q[197];
cx q[38], q[17];
cx q[78], q[59];
cx q[114], q[97];
cx q[146], q[131];
cx q[174], q[161];
cx q[60], q[40];
cx q[98], q[80];
cx q[132], q[116];
cx q[41], q[20];
cx q[81], q[62];
barrier q;

cx q[0], q[22];
cx q[23], q[44];
cx q[63], q[82];
cx q[3], q[25];
cx q[45], q[65];
cx q[83], q[101];
cx q[117], q[133];
cx q[26], q[47];
cx q[66], q[85];
cx q[102], q[119];
cx q[134], q[149];
cx q[162], q[175];
cx q[6], q[28];
cx q[48], q[68];
cx q[86], q[104];
cx q[120], q[136];
cx q[150], q[164];
cx q[176], q[188];
cx q[198], q[208];
cx q[29], q[50];
cx q[69], q[88];
cx q[105], q[122];
cx q[137], q[152];
cx q[165], q[178];
cx q[189], q[200];
cx q[209], q[218];
cx q[225], q[232];
cx q[9], q[31];
cx q[51], q[71];
cx q[89], q[107];
cx q[123], q[139];
cx q[153], q[167];
cx q[179], q[191];
cx q[201], q[211];
cx q[219], q[227];
cx q[233], q[239];
cx q[243], q[247];
cx q[32], q[53];
cx q[72], q[91];
cx q[108], q[125];
cx q[140], q[155];
cx q[168], q[181];
cx q[192], q[203];
cx q[212], q[221];
cx q[228], q[235];
cx q[240], q[245];
cx q[248], q[251];
cx q[12], q[34];
cx q[54], q[74];
cx q[92], q[110];
cx q[126], q[142];
cx q[156], q[170];
cx q[182], q[194];
cx q[204], q[214];
cx q[222], q[230];
cx q[236], q[242];
cx q[35], q[56];
cx q[75], q[94];
cx q[111], q[128];
cx q[143], q[158];
cx q[171], q[184];
cx q[195], q[206];
cx q[215], q[224];
cx q[15], q[37];
cx q[57], q[77];
cx q[95], q[113];
cx q[129], q[145];
cx q[159], q[173];
cx q[185], q[197];
cx q[38], q[59];
cx q[78], q[97];
cx q[114], q[131];
cx q[146], q[161];
cx q[18], q[40];
cx q[60], q[80];
cx q[98], q[116];
cx q[41], q[62];
barrier q;

measure q[2] -> rec[0]; reset q[2]; // decomposed MR
measure q[5] -> rec[1]; reset q[5]; // decomposed MR
measure q[8] -> rec[2]; reset q[8]; // decomposed MR
measure q[11] -> rec[3]; reset q[11]; // decomposed MR
measure q[14] -> rec[4]; reset q[14]; // decomposed MR
measure q[17] -> rec[5]; reset q[17]; // decomposed MR
measure q[20] -> rec[6]; reset q[20]; // decomposed MR
measure q[22] -> rec[7]; reset q[22]; // decomposed MR
measure q[25] -> rec[8]; reset q[25]; // decomposed MR
measure q[28] -> rec[9]; reset q[28]; // decomposed MR
measure q[31] -> rec[10]; reset q[31]; // decomposed MR
measure q[34] -> rec[11]; reset q[34]; // decomposed MR
measure q[37] -> rec[12]; reset q[37]; // decomposed MR
measure q[40] -> rec[13]; reset q[40]; // decomposed MR
measure q[44] -> rec[14]; reset q[44]; // decomposed MR
measure q[47] -> rec[15]; reset q[47]; // decomposed MR
measure q[50] -> rec[16]; reset q[50]; // decomposed MR
measure q[53] -> rec[17]; reset q[53]; // decomposed MR
measure q[56] -> rec[18]; reset q[56]; // decomposed MR
measure q[59] -> rec[19]; reset q[59]; // decomposed MR
measure q[62] -> rec[20]; reset q[62]; // decomposed MR
measure q[65] -> rec[21]; reset q[65]; // decomposed MR
measure q[68] -> rec[22]; reset q[68]; // decomposed MR
measure q[71] -> rec[23]; reset q[71]; // decomposed MR
measure q[74] -> rec[24]; reset q[74]; // decomposed MR
measure q[77] -> rec[25]; reset q[77]; // decomposed MR
measure q[80] -> rec[26]; reset q[80]; // decomposed MR
measure q[82] -> rec[27]; reset q[82]; // decomposed MR
measure q[85] -> rec[28]; reset q[85]; // decomposed MR
measure q[88] -> rec[29]; reset q[88]; // decomposed MR
measure q[91] -> rec[30]; reset q[91]; // decomposed MR
measure q[94] -> rec[31]; reset q[94]; // decomposed MR
measure q[97] -> rec[32]; reset q[97]; // decomposed MR
measure q[101] -> rec[33]; reset q[101]; // decomposed MR
measure q[104] -> rec[34]; reset q[104]; // decomposed MR
measure q[107] -> rec[35]; reset q[107]; // decomposed MR
measure q[110] -> rec[36]; reset q[110]; // decomposed MR
measure q[113] -> rec[37]; reset q[113]; // decomposed MR
measure q[116] -> rec[38]; reset q[116]; // decomposed MR
measure q[119] -> rec[39]; reset q[119]; // decomposed MR
measure q[122] -> rec[40]; reset q[122]; // decomposed MR
measure q[125] -> rec[41]; reset q[125]; // decomposed MR
measure q[128] -> rec[42]; reset q[128]; // decomposed MR
measure q[131] -> rec[43]; reset q[131]; // decomposed MR
measure q[133] -> rec[44]; reset q[133]; // decomposed MR
measure q[136] -> rec[45]; reset q[136]; // decomposed MR
measure q[139] -> rec[46]; reset q[139]; // decomposed MR
measure q[142] -> rec[47]; reset q[142]; // decomposed MR
measure q[145] -> rec[48]; reset q[145]; // decomposed MR
measure q[149] -> rec[49]; reset q[149]; // decomposed MR
measure q[152] -> rec[50]; reset q[152]; // decomposed MR
measure q[155] -> rec[51]; reset q[155]; // decomposed MR
measure q[158] -> rec[52]; reset q[158]; // decomposed MR
measure q[161] -> rec[53]; reset q[161]; // decomposed MR
measure q[164] -> rec[54]; reset q[164]; // decomposed MR
measure q[167] -> rec[55]; reset q[167]; // decomposed MR
measure q[170] -> rec[56]; reset q[170]; // decomposed MR
measure q[173] -> rec[57]; reset q[173]; // decomposed MR
measure q[175] -> rec[58]; reset q[175]; // decomposed MR
measure q[178] -> rec[59]; reset q[178]; // decomposed MR
measure q[181] -> rec[60]; reset q[181]; // decomposed MR
measure q[184] -> rec[61]; reset q[184]; // decomposed MR
measure q[188] -> rec[62]; reset q[188]; // decomposed MR
measure q[191] -> rec[63]; reset q[191]; // decomposed MR
measure q[194] -> rec[64]; reset q[194]; // decomposed MR
measure q[197] -> rec[65]; reset q[197]; // decomposed MR
measure q[200] -> rec[66]; reset q[200]; // decomposed MR
measure q[203] -> rec[67]; reset q[203]; // decomposed MR
measure q[206] -> rec[68]; reset q[206]; // decomposed MR
measure q[208] -> rec[69]; reset q[208]; // decomposed MR
measure q[211] -> rec[70]; reset q[211]; // decomposed MR
measure q[214] -> rec[71]; reset q[214]; // decomposed MR
measure q[218] -> rec[72]; reset q[218]; // decomposed MR
measure q[221] -> rec[73]; reset q[221]; // decomposed MR
measure q[224] -> rec[74]; reset q[224]; // decomposed MR
measure q[227] -> rec[75]; reset q[227]; // decomposed MR
measure q[230] -> rec[76]; reset q[230]; // decomposed MR
measure q[232] -> rec[77]; reset q[232]; // decomposed MR
measure q[235] -> rec[78]; reset q[235]; // decomposed MR
measure q[239] -> rec[79]; reset q[239]; // decomposed MR
measure q[242] -> rec[80]; reset q[242]; // decomposed MR
measure q[245] -> rec[81]; reset q[245]; // decomposed MR
measure q[247] -> rec[82]; reset q[247]; // decomposed MR
measure q[251] -> rec[83]; reset q[251]; // decomposed MR
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
cxyz q[19];
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
cxyz q[38];
cxyz q[39];
cxyz q[41];
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
cxyz q[70];
cxyz q[72];
cxyz q[73];
cxyz q[75];
cxyz q[76];
cxyz q[78];
cxyz q[79];
cxyz q[81];
cxyz q[83];
cxyz q[84];
cxyz q[86];
cxyz q[87];
cxyz q[89];
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
cxyz q[115];
cxyz q[117];
cxyz q[118];
cxyz q[120];
cxyz q[121];
cxyz q[123];
cxyz q[124];
cxyz q[126];
cxyz q[127];
cxyz q[129];
cxyz q[130];
cxyz q[132];
cxyz q[134];
cxyz q[135];
cxyz q[137];
cxyz q[138];
cxyz q[140];
cxyz q[141];
cxyz q[143];
cxyz q[144];
cxyz q[146];
cxyz q[147];
cxyz q[148];
cxyz q[150];
cxyz q[151];
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
cxyz q[169];
cxyz q[171];
cxyz q[172];
cxyz q[174];
cxyz q[176];
cxyz q[177];
cxyz q[179];
cxyz q[180];
cxyz q[182];
cxyz q[183];
cxyz q[185];
cxyz q[186];
cxyz q[187];
cxyz q[189];
cxyz q[190];
cxyz q[192];
cxyz q[193];
cxyz q[195];
cxyz q[196];
cxyz q[198];
cxyz q[199];
cxyz q[201];
cxyz q[202];
cxyz q[204];
cxyz q[205];
cxyz q[207];
cxyz q[209];
cxyz q[210];
cxyz q[212];
cxyz q[213];
cxyz q[215];
cxyz q[216];
cxyz q[217];
cxyz q[219];
cxyz q[220];
cxyz q[222];
cxyz q[223];
cxyz q[225];
cxyz q[226];
cxyz q[228];
cxyz q[229];
cxyz q[231];
cxyz q[233];
cxyz q[234];
cxyz q[236];
cxyz q[237];
cxyz q[238];
cxyz q[240];
cxyz q[241];
cxyz q[243];
cxyz q[244];
cxyz q[246];
cxyz q[248];
cxyz q[249];
cxyz q[250];
cxyz q[252];
barrier q;

cx q[23], q[22];
cx q[3], q[2];
cx q[45], q[44];
cx q[83], q[82];
cx q[26], q[25];
cx q[66], q[65];
cx q[102], q[101];
cx q[134], q[133];
cx q[6], q[5];
cx q[48], q[47];
cx q[86], q[85];
cx q[120], q[119];
cx q[150], q[149];
cx q[176], q[175];
cx q[29], q[28];
cx q[69], q[68];
cx q[105], q[104];
cx q[137], q[136];
cx q[165], q[164];
cx q[189], q[188];
cx q[209], q[208];
cx q[9], q[8];
cx q[51], q[50];
cx q[89], q[88];
cx q[123], q[122];
cx q[153], q[152];
cx q[179], q[178];
cx q[201], q[200];
cx q[219], q[218];
cx q[233], q[232];
cx q[32], q[31];
cx q[72], q[71];
cx q[108], q[107];
cx q[140], q[139];
cx q[168], q[167];
cx q[192], q[191];
cx q[212], q[211];
cx q[228], q[227];
cx q[240], q[239];
cx q[248], q[247];
cx q[12], q[11];
cx q[54], q[53];
cx q[92], q[91];
cx q[126], q[125];
cx q[156], q[155];
cx q[182], q[181];
cx q[204], q[203];
cx q[222], q[221];
cx q[236], q[235];
cx q[246], q[245];
cx q[35], q[34];
cx q[75], q[74];
cx q[111], q[110];
cx q[143], q[142];
cx q[171], q[170];
cx q[195], q[194];
cx q[215], q[214];
cx q[231], q[230];
cx q[15], q[14];
cx q[57], q[56];
cx q[95], q[94];
cx q[129], q[128];
cx q[159], q[158];
cx q[185], q[184];
cx q[207], q[206];
cx q[38], q[37];
cx q[78], q[77];
cx q[114], q[113];
cx q[146], q[145];
cx q[174], q[173];
cx q[18], q[17];
cx q[60], q[59];
cx q[98], q[97];
cx q[132], q[131];
cx q[41], q[40];
cx q[81], q[80];
cx q[21], q[20];
barrier q;

cx q[43], q[22];
cx q[24], q[2];
cx q[64], q[44];
cx q[100], q[82];
cx q[46], q[25];
cx q[84], q[65];
cx q[118], q[101];
cx q[148], q[133];
cx q[27], q[5];
cx q[67], q[47];
cx q[103], q[85];
cx q[135], q[119];
cx q[163], q[149];
cx q[187], q[175];
cx q[49], q[28];
cx q[87], q[68];
cx q[121], q[104];
cx q[151], q[136];
cx q[177], q[164];
cx q[199], q[188];
cx q[217], q[208];
cx q[30], q[8];
cx q[70], q[50];
cx q[106], q[88];
cx q[138], q[122];
cx q[166], q[152];
cx q[190], q[178];
cx q[210], q[200];
cx q[226], q[218];
cx q[238], q[232];
cx q[52], q[31];
cx q[90], q[71];
cx q[124], q[107];
cx q[154], q[139];
cx q[180], q[167];
cx q[202], q[191];
cx q[220], q[211];
cx q[234], q[227];
cx q[244], q[239];
cx q[250], q[247];
cx q[33], q[11];
cx q[73], q[53];
cx q[109], q[91];
cx q[141], q[125];
cx q[169], q[155];
cx q[193], q[181];
cx q[213], q[203];
cx q[229], q[221];
cx q[241], q[235];
cx q[249], q[245];
cx q[55], q[34];
cx q[93], q[74];
cx q[127], q[110];
cx q[157], q[142];
cx q[183], q[170];
cx q[205], q[194];
cx q[223], q[214];
cx q[237], q[230];
cx q[36], q[14];
cx q[76], q[56];
cx q[112], q[94];
cx q[144], q[128];
cx q[172], q[158];
cx q[196], q[184];
cx q[216], q[206];
cx q[58], q[37];
cx q[96], q[77];
cx q[130], q[113];
cx q[160], q[145];
cx q[186], q[173];
cx q[39], q[17];
cx q[79], q[59];
cx q[115], q[97];
cx q[147], q[131];
cx q[61], q[40];
cx q[99], q[80];
cx q[42], q[20];
barrier q;

cx q[1], q[22];
cx q[24], q[44];
cx q[64], q[82];
cx q[4], q[25];
cx q[46], q[65];
cx q[84], q[101];
cx q[118], q[133];
cx q[27], q[47];
cx q[67], q[85];
cx q[103], q[119];
cx q[135], q[149];
cx q[163], q[175];
cx q[7], q[28];
cx q[49], q[68];
cx q[87], q[104];
cx q[121], q[136];
cx q[151], q[164];
cx q[177], q[188];
cx q[199], q[208];
cx q[30], q[50];
cx q[70], q[88];
cx q[106], q[122];
cx q[138], q[152];
cx q[166], q[178];
cx q[190], q[200];
cx q[210], q[218];
cx q[226], q[232];
cx q[10], q[31];
cx q[52], q[71];
cx q[90], q[107];
cx q[124], q[139];
cx q[154], q[167];
cx q[180], q[191];
cx q[202], q[211];
cx q[220], q[227];
cx q[234], q[239];
cx q[244], q[247];
cx q[33], q[53];
cx q[73], q[91];
cx q[109], q[125];
cx q[141], q[155];
cx q[169], q[181];
cx q[193], q[203];
cx q[213], q[221];
cx q[229], q[235];
cx q[241], q[245];
cx q[249], q[251];
cx q[13], q[34];
cx q[55], q[74];
cx q[93], q[110];
cx q[127], q[142];
cx q[157], q[170];
cx q[183], q[194];
cx q[205], q[214];
cx q[223], q[230];
cx q[237], q[242];
cx q[36], q[56];
cx q[76], q[94];
cx q[112], q[128];
cx q[144], q[158];
cx q[172], q[184];
cx q[196], q[206];
cx q[216], q[224];
cx q[16], q[37];
cx q[58], q[77];
cx q[96], q[113];
cx q[130], q[145];
cx q[160], q[173];
cx q[186], q[197];
cx q[39], q[59];
cx q[79], q[97];
cx q[115], q[131];
cx q[147], q[161];
cx q[19], q[40];
cx q[61], q[80];
cx q[99], q[116];
cx q[42], q[62];
barrier q;

cx q[1], q[2];
cx q[43], q[44];
cx q[24], q[25];
cx q[64], q[65];
cx q[100], q[101];
cx q[4], q[5];
cx q[46], q[47];
cx q[84], q[85];
cx q[118], q[119];
cx q[148], q[149];
cx q[27], q[28];
cx q[67], q[68];
cx q[103], q[104];
cx q[135], q[136];
cx q[163], q[164];
cx q[187], q[188];
cx q[7], q[8];
cx q[49], q[50];
cx q[87], q[88];
cx q[121], q[122];
cx q[151], q[152];
cx q[177], q[178];
cx q[199], q[200];
cx q[217], q[218];
cx q[30], q[31];
cx q[70], q[71];
cx q[106], q[107];
cx q[138], q[139];
cx q[166], q[167];
cx q[190], q[191];
cx q[210], q[211];
cx q[226], q[227];
cx q[238], q[239];
cx q[10], q[11];
cx q[52], q[53];
cx q[90], q[91];
cx q[124], q[125];
cx q[154], q[155];
cx q[180], q[181];
cx q[202], q[203];
cx q[220], q[221];
cx q[234], q[235];
cx q[244], q[245];
cx q[250], q[251];
cx q[33], q[34];
cx q[73], q[74];
cx q[109], q[110];
cx q[141], q[142];
cx q[169], q[170];
cx q[193], q[194];
cx q[213], q[214];
cx q[229], q[230];
cx q[241], q[242];
cx q[13], q[14];
cx q[55], q[56];
cx q[93], q[94];
cx q[127], q[128];
cx q[157], q[158];
cx q[183], q[184];
cx q[205], q[206];
cx q[223], q[224];
cx q[36], q[37];
cx q[76], q[77];
cx q[112], q[113];
cx q[144], q[145];
cx q[172], q[173];
cx q[196], q[197];
cx q[16], q[17];
cx q[58], q[59];
cx q[96], q[97];
cx q[130], q[131];
cx q[160], q[161];
cx q[39], q[40];
cx q[79], q[80];
cx q[115], q[116];
cx q[19], q[20];
cx q[61], q[62];
barrier q;

cx q[23], q[2];
cx q[63], q[44];
cx q[45], q[25];
cx q[83], q[65];
cx q[117], q[101];
cx q[26], q[5];
cx q[66], q[47];
cx q[102], q[85];
cx q[134], q[119];
cx q[162], q[149];
cx q[48], q[28];
cx q[86], q[68];
cx q[120], q[104];
cx q[150], q[136];
cx q[176], q[164];
cx q[198], q[188];
cx q[29], q[8];
cx q[69], q[50];
cx q[105], q[88];
cx q[137], q[122];
cx q[165], q[152];
cx q[189], q[178];
cx q[209], q[200];
cx q[225], q[218];
cx q[51], q[31];
cx q[89], q[71];
cx q[123], q[107];
cx q[153], q[139];
cx q[179], q[167];
cx q[201], q[191];
cx q[219], q[211];
cx q[233], q[227];
cx q[243], q[239];
cx q[32], q[11];
cx q[72], q[53];
cx q[108], q[91];
cx q[140], q[125];
cx q[168], q[155];
cx q[192], q[181];
cx q[212], q[203];
cx q[228], q[221];
cx q[240], q[235];
cx q[248], q[245];
cx q[252], q[251];
cx q[54], q[34];
cx q[92], q[74];
cx q[126], q[110];
cx q[156], q[142];
cx q[182], q[170];
cx q[204], q[194];
cx q[222], q[214];
cx q[236], q[230];
cx q[246], q[242];
cx q[35], q[14];
cx q[75], q[56];
cx q[111], q[94];
cx q[143], q[128];
cx q[171], q[158];
cx q[195], q[184];
cx q[215], q[206];
cx q[231], q[224];
cx q[57], q[37];
cx q[95], q[77];
cx q[129], q[113];
cx q[159], q[145];
cx q[185], q[173];
cx q[207], q[197];
cx q[38], q[17];
cx q[78], q[59];
cx q[114], q[97];
cx q[146], q[131];
cx q[174], q[161];
cx q[60], q[40];
cx q[98], q[80];
cx q[132], q[116];
cx q[41], q[20];
cx q[81], q[62];
barrier q;

cx q[0], q[22];
cx q[23], q[44];
cx q[63], q[82];
cx q[3], q[25];
cx q[45], q[65];
cx q[83], q[101];
cx q[117], q[133];
cx q[26], q[47];
cx q[66], q[85];
cx q[102], q[119];
cx q[134], q[149];
cx q[162], q[175];
cx q[6], q[28];
cx q[48], q[68];
cx q[86], q[104];
cx q[120], q[136];
cx q[150], q[164];
cx q[176], q[188];
cx q[198], q[208];
cx q[29], q[50];
cx q[69], q[88];
cx q[105], q[122];
cx q[137], q[152];
cx q[165], q[178];
cx q[189], q[200];
cx q[209], q[218];
cx q[225], q[232];
cx q[9], q[31];
cx q[51], q[71];
cx q[89], q[107];
cx q[123], q[139];
cx q[153], q[167];
cx q[179], q[191];
cx q[201], q[211];
cx q[219], q[227];
cx q[233], q[239];
cx q[243], q[247];
cx q[32], q[53];
cx q[72], q[91];
cx q[108], q[125];
cx q[140], q[155];
cx q[168], q[181];
cx q[192], q[203];
cx q[212], q[221];
cx q[228], q[235];
cx q[240], q[245];
cx q[248], q[251];
cx q[12], q[34];
cx q[54], q[74];
cx q[92], q[110];
cx q[126], q[142];
cx q[156], q[170];
cx q[182], q[194];
cx q[204], q[214];
cx q[222], q[230];
cx q[236], q[242];
cx q[35], q[56];
cx q[75], q[94];
cx q[111], q[128];
cx q[143], q[158];
cx q[171], q[184];
cx q[195], q[206];
cx q[215], q[224];
cx q[15], q[37];
cx q[57], q[77];
cx q[95], q[113];
cx q[129], q[145];
cx q[159], q[173];
cx q[185], q[197];
cx q[38], q[59];
cx q[78], q[97];
cx q[114], q[131];
cx q[146], q[161];
cx q[18], q[40];
cx q[60], q[80];
cx q[98], q[116];
cx q[41], q[62];
barrier q;

measure q[2] -> rec[84]; reset q[2]; // decomposed MR
measure q[5] -> rec[85]; reset q[5]; // decomposed MR
measure q[8] -> rec[86]; reset q[8]; // decomposed MR
measure q[11] -> rec[87]; reset q[11]; // decomposed MR
measure q[14] -> rec[88]; reset q[14]; // decomposed MR
measure q[17] -> rec[89]; reset q[17]; // decomposed MR
measure q[20] -> rec[90]; reset q[20]; // decomposed MR
measure q[22] -> rec[91]; reset q[22]; // decomposed MR
measure q[25] -> rec[92]; reset q[25]; // decomposed MR
measure q[28] -> rec[93]; reset q[28]; // decomposed MR
measure q[31] -> rec[94]; reset q[31]; // decomposed MR
measure q[34] -> rec[95]; reset q[34]; // decomposed MR
measure q[37] -> rec[96]; reset q[37]; // decomposed MR
measure q[40] -> rec[97]; reset q[40]; // decomposed MR
measure q[44] -> rec[98]; reset q[44]; // decomposed MR
measure q[47] -> rec[99]; reset q[47]; // decomposed MR
measure q[50] -> rec[100]; reset q[50]; // decomposed MR
measure q[53] -> rec[101]; reset q[53]; // decomposed MR
measure q[56] -> rec[102]; reset q[56]; // decomposed MR
measure q[59] -> rec[103]; reset q[59]; // decomposed MR
measure q[62] -> rec[104]; reset q[62]; // decomposed MR
measure q[65] -> rec[105]; reset q[65]; // decomposed MR
measure q[68] -> rec[106]; reset q[68]; // decomposed MR
measure q[71] -> rec[107]; reset q[71]; // decomposed MR
measure q[74] -> rec[108]; reset q[74]; // decomposed MR
measure q[77] -> rec[109]; reset q[77]; // decomposed MR
measure q[80] -> rec[110]; reset q[80]; // decomposed MR
measure q[82] -> rec[111]; reset q[82]; // decomposed MR
measure q[85] -> rec[112]; reset q[85]; // decomposed MR
measure q[88] -> rec[113]; reset q[88]; // decomposed MR
measure q[91] -> rec[114]; reset q[91]; // decomposed MR
measure q[94] -> rec[115]; reset q[94]; // decomposed MR
measure q[97] -> rec[116]; reset q[97]; // decomposed MR
measure q[101] -> rec[117]; reset q[101]; // decomposed MR
measure q[104] -> rec[118]; reset q[104]; // decomposed MR
measure q[107] -> rec[119]; reset q[107]; // decomposed MR
measure q[110] -> rec[120]; reset q[110]; // decomposed MR
measure q[113] -> rec[121]; reset q[113]; // decomposed MR
measure q[116] -> rec[122]; reset q[116]; // decomposed MR
measure q[119] -> rec[123]; reset q[119]; // decomposed MR
measure q[122] -> rec[124]; reset q[122]; // decomposed MR
measure q[125] -> rec[125]; reset q[125]; // decomposed MR
measure q[128] -> rec[126]; reset q[128]; // decomposed MR
measure q[131] -> rec[127]; reset q[131]; // decomposed MR
measure q[133] -> rec[128]; reset q[133]; // decomposed MR
measure q[136] -> rec[129]; reset q[136]; // decomposed MR
measure q[139] -> rec[130]; reset q[139]; // decomposed MR
measure q[142] -> rec[131]; reset q[142]; // decomposed MR
measure q[145] -> rec[132]; reset q[145]; // decomposed MR
measure q[149] -> rec[133]; reset q[149]; // decomposed MR
measure q[152] -> rec[134]; reset q[152]; // decomposed MR
measure q[155] -> rec[135]; reset q[155]; // decomposed MR
measure q[158] -> rec[136]; reset q[158]; // decomposed MR
measure q[161] -> rec[137]; reset q[161]; // decomposed MR
measure q[164] -> rec[138]; reset q[164]; // decomposed MR
measure q[167] -> rec[139]; reset q[167]; // decomposed MR
measure q[170] -> rec[140]; reset q[170]; // decomposed MR
measure q[173] -> rec[141]; reset q[173]; // decomposed MR
measure q[175] -> rec[142]; reset q[175]; // decomposed MR
measure q[178] -> rec[143]; reset q[178]; // decomposed MR
measure q[181] -> rec[144]; reset q[181]; // decomposed MR
measure q[184] -> rec[145]; reset q[184]; // decomposed MR
measure q[188] -> rec[146]; reset q[188]; // decomposed MR
measure q[191] -> rec[147]; reset q[191]; // decomposed MR
measure q[194] -> rec[148]; reset q[194]; // decomposed MR
measure q[197] -> rec[149]; reset q[197]; // decomposed MR
measure q[200] -> rec[150]; reset q[200]; // decomposed MR
measure q[203] -> rec[151]; reset q[203]; // decomposed MR
measure q[206] -> rec[152]; reset q[206]; // decomposed MR
measure q[208] -> rec[153]; reset q[208]; // decomposed MR
measure q[211] -> rec[154]; reset q[211]; // decomposed MR
measure q[214] -> rec[155]; reset q[214]; // decomposed MR
measure q[218] -> rec[156]; reset q[218]; // decomposed MR
measure q[221] -> rec[157]; reset q[221]; // decomposed MR
measure q[224] -> rec[158]; reset q[224]; // decomposed MR
measure q[227] -> rec[159]; reset q[227]; // decomposed MR
measure q[230] -> rec[160]; reset q[230]; // decomposed MR
measure q[232] -> rec[161]; reset q[232]; // decomposed MR
measure q[235] -> rec[162]; reset q[235]; // decomposed MR
measure q[239] -> rec[163]; reset q[239]; // decomposed MR
measure q[242] -> rec[164]; reset q[242]; // decomposed MR
measure q[245] -> rec[165]; reset q[245]; // decomposed MR
measure q[247] -> rec[166]; reset q[247]; // decomposed MR
measure q[251] -> rec[167]; reset q[251]; // decomposed MR
s q[0]; s q[0]; s q[0]; h q[0]; measure q[0] -> rec[168]; h q[0]; s q[0]; // decomposed MY
s q[1]; s q[1]; s q[1]; h q[1]; measure q[1] -> rec[169]; h q[1]; s q[1]; // decomposed MY
s q[3]; s q[3]; s q[3]; h q[3]; measure q[3] -> rec[170]; h q[3]; s q[3]; // decomposed MY
s q[4]; s q[4]; s q[4]; h q[4]; measure q[4] -> rec[171]; h q[4]; s q[4]; // decomposed MY
s q[6]; s q[6]; s q[6]; h q[6]; measure q[6] -> rec[172]; h q[6]; s q[6]; // decomposed MY
s q[7]; s q[7]; s q[7]; h q[7]; measure q[7] -> rec[173]; h q[7]; s q[7]; // decomposed MY
s q[9]; s q[9]; s q[9]; h q[9]; measure q[9] -> rec[174]; h q[9]; s q[9]; // decomposed MY
s q[10]; s q[10]; s q[10]; h q[10]; measure q[10] -> rec[175]; h q[10]; s q[10]; // decomposed MY
s q[12]; s q[12]; s q[12]; h q[12]; measure q[12] -> rec[176]; h q[12]; s q[12]; // decomposed MY
s q[13]; s q[13]; s q[13]; h q[13]; measure q[13] -> rec[177]; h q[13]; s q[13]; // decomposed MY
s q[15]; s q[15]; s q[15]; h q[15]; measure q[15] -> rec[178]; h q[15]; s q[15]; // decomposed MY
s q[16]; s q[16]; s q[16]; h q[16]; measure q[16] -> rec[179]; h q[16]; s q[16]; // decomposed MY
s q[18]; s q[18]; s q[18]; h q[18]; measure q[18] -> rec[180]; h q[18]; s q[18]; // decomposed MY
s q[19]; s q[19]; s q[19]; h q[19]; measure q[19] -> rec[181]; h q[19]; s q[19]; // decomposed MY
s q[21]; s q[21]; s q[21]; h q[21]; measure q[21] -> rec[182]; h q[21]; s q[21]; // decomposed MY
s q[23]; s q[23]; s q[23]; h q[23]; measure q[23] -> rec[183]; h q[23]; s q[23]; // decomposed MY
s q[24]; s q[24]; s q[24]; h q[24]; measure q[24] -> rec[184]; h q[24]; s q[24]; // decomposed MY
s q[26]; s q[26]; s q[26]; h q[26]; measure q[26] -> rec[185]; h q[26]; s q[26]; // decomposed MY
s q[27]; s q[27]; s q[27]; h q[27]; measure q[27] -> rec[186]; h q[27]; s q[27]; // decomposed MY
s q[29]; s q[29]; s q[29]; h q[29]; measure q[29] -> rec[187]; h q[29]; s q[29]; // decomposed MY
s q[30]; s q[30]; s q[30]; h q[30]; measure q[30] -> rec[188]; h q[30]; s q[30]; // decomposed MY
s q[32]; s q[32]; s q[32]; h q[32]; measure q[32] -> rec[189]; h q[32]; s q[32]; // decomposed MY
s q[33]; s q[33]; s q[33]; h q[33]; measure q[33] -> rec[190]; h q[33]; s q[33]; // decomposed MY
s q[35]; s q[35]; s q[35]; h q[35]; measure q[35] -> rec[191]; h q[35]; s q[35]; // decomposed MY
s q[36]; s q[36]; s q[36]; h q[36]; measure q[36] -> rec[192]; h q[36]; s q[36]; // decomposed MY
s q[38]; s q[38]; s q[38]; h q[38]; measure q[38] -> rec[193]; h q[38]; s q[38]; // decomposed MY
s q[39]; s q[39]; s q[39]; h q[39]; measure q[39] -> rec[194]; h q[39]; s q[39]; // decomposed MY
s q[41]; s q[41]; s q[41]; h q[41]; measure q[41] -> rec[195]; h q[41]; s q[41]; // decomposed MY
s q[42]; s q[42]; s q[42]; h q[42]; measure q[42] -> rec[196]; h q[42]; s q[42]; // decomposed MY
s q[43]; s q[43]; s q[43]; h q[43]; measure q[43] -> rec[197]; h q[43]; s q[43]; // decomposed MY
s q[45]; s q[45]; s q[45]; h q[45]; measure q[45] -> rec[198]; h q[45]; s q[45]; // decomposed MY
s q[46]; s q[46]; s q[46]; h q[46]; measure q[46] -> rec[199]; h q[46]; s q[46]; // decomposed MY
s q[48]; s q[48]; s q[48]; h q[48]; measure q[48] -> rec[200]; h q[48]; s q[48]; // decomposed MY
s q[49]; s q[49]; s q[49]; h q[49]; measure q[49] -> rec[201]; h q[49]; s q[49]; // decomposed MY
s q[51]; s q[51]; s q[51]; h q[51]; measure q[51] -> rec[202]; h q[51]; s q[51]; // decomposed MY
s q[52]; s q[52]; s q[52]; h q[52]; measure q[52] -> rec[203]; h q[52]; s q[52]; // decomposed MY
s q[54]; s q[54]; s q[54]; h q[54]; measure q[54] -> rec[204]; h q[54]; s q[54]; // decomposed MY
s q[55]; s q[55]; s q[55]; h q[55]; measure q[55] -> rec[205]; h q[55]; s q[55]; // decomposed MY
s q[57]; s q[57]; s q[57]; h q[57]; measure q[57] -> rec[206]; h q[57]; s q[57]; // decomposed MY
s q[58]; s q[58]; s q[58]; h q[58]; measure q[58] -> rec[207]; h q[58]; s q[58]; // decomposed MY
s q[60]; s q[60]; s q[60]; h q[60]; measure q[60] -> rec[208]; h q[60]; s q[60]; // decomposed MY
s q[61]; s q[61]; s q[61]; h q[61]; measure q[61] -> rec[209]; h q[61]; s q[61]; // decomposed MY
s q[63]; s q[63]; s q[63]; h q[63]; measure q[63] -> rec[210]; h q[63]; s q[63]; // decomposed MY
s q[64]; s q[64]; s q[64]; h q[64]; measure q[64] -> rec[211]; h q[64]; s q[64]; // decomposed MY
s q[66]; s q[66]; s q[66]; h q[66]; measure q[66] -> rec[212]; h q[66]; s q[66]; // decomposed MY
s q[67]; s q[67]; s q[67]; h q[67]; measure q[67] -> rec[213]; h q[67]; s q[67]; // decomposed MY
s q[69]; s q[69]; s q[69]; h q[69]; measure q[69] -> rec[214]; h q[69]; s q[69]; // decomposed MY
s q[70]; s q[70]; s q[70]; h q[70]; measure q[70] -> rec[215]; h q[70]; s q[70]; // decomposed MY
s q[72]; s q[72]; s q[72]; h q[72]; measure q[72] -> rec[216]; h q[72]; s q[72]; // decomposed MY
s q[73]; s q[73]; s q[73]; h q[73]; measure q[73] -> rec[217]; h q[73]; s q[73]; // decomposed MY
s q[75]; s q[75]; s q[75]; h q[75]; measure q[75] -> rec[218]; h q[75]; s q[75]; // decomposed MY
s q[76]; s q[76]; s q[76]; h q[76]; measure q[76] -> rec[219]; h q[76]; s q[76]; // decomposed MY
s q[78]; s q[78]; s q[78]; h q[78]; measure q[78] -> rec[220]; h q[78]; s q[78]; // decomposed MY
s q[79]; s q[79]; s q[79]; h q[79]; measure q[79] -> rec[221]; h q[79]; s q[79]; // decomposed MY
s q[81]; s q[81]; s q[81]; h q[81]; measure q[81] -> rec[222]; h q[81]; s q[81]; // decomposed MY
s q[83]; s q[83]; s q[83]; h q[83]; measure q[83] -> rec[223]; h q[83]; s q[83]; // decomposed MY
s q[84]; s q[84]; s q[84]; h q[84]; measure q[84] -> rec[224]; h q[84]; s q[84]; // decomposed MY
s q[86]; s q[86]; s q[86]; h q[86]; measure q[86] -> rec[225]; h q[86]; s q[86]; // decomposed MY
s q[87]; s q[87]; s q[87]; h q[87]; measure q[87] -> rec[226]; h q[87]; s q[87]; // decomposed MY
s q[89]; s q[89]; s q[89]; h q[89]; measure q[89] -> rec[227]; h q[89]; s q[89]; // decomposed MY
s q[90]; s q[90]; s q[90]; h q[90]; measure q[90] -> rec[228]; h q[90]; s q[90]; // decomposed MY
s q[92]; s q[92]; s q[92]; h q[92]; measure q[92] -> rec[229]; h q[92]; s q[92]; // decomposed MY
s q[93]; s q[93]; s q[93]; h q[93]; measure q[93] -> rec[230]; h q[93]; s q[93]; // decomposed MY
s q[95]; s q[95]; s q[95]; h q[95]; measure q[95] -> rec[231]; h q[95]; s q[95]; // decomposed MY
s q[96]; s q[96]; s q[96]; h q[96]; measure q[96] -> rec[232]; h q[96]; s q[96]; // decomposed MY
s q[98]; s q[98]; s q[98]; h q[98]; measure q[98] -> rec[233]; h q[98]; s q[98]; // decomposed MY
s q[99]; s q[99]; s q[99]; h q[99]; measure q[99] -> rec[234]; h q[99]; s q[99]; // decomposed MY
s q[100]; s q[100]; s q[100]; h q[100]; measure q[100] -> rec[235]; h q[100]; s q[100]; // decomposed MY
s q[102]; s q[102]; s q[102]; h q[102]; measure q[102] -> rec[236]; h q[102]; s q[102]; // decomposed MY
s q[103]; s q[103]; s q[103]; h q[103]; measure q[103] -> rec[237]; h q[103]; s q[103]; // decomposed MY
s q[105]; s q[105]; s q[105]; h q[105]; measure q[105] -> rec[238]; h q[105]; s q[105]; // decomposed MY
s q[106]; s q[106]; s q[106]; h q[106]; measure q[106] -> rec[239]; h q[106]; s q[106]; // decomposed MY
s q[108]; s q[108]; s q[108]; h q[108]; measure q[108] -> rec[240]; h q[108]; s q[108]; // decomposed MY
s q[109]; s q[109]; s q[109]; h q[109]; measure q[109] -> rec[241]; h q[109]; s q[109]; // decomposed MY
s q[111]; s q[111]; s q[111]; h q[111]; measure q[111] -> rec[242]; h q[111]; s q[111]; // decomposed MY
s q[112]; s q[112]; s q[112]; h q[112]; measure q[112] -> rec[243]; h q[112]; s q[112]; // decomposed MY
s q[114]; s q[114]; s q[114]; h q[114]; measure q[114] -> rec[244]; h q[114]; s q[114]; // decomposed MY
s q[115]; s q[115]; s q[115]; h q[115]; measure q[115] -> rec[245]; h q[115]; s q[115]; // decomposed MY
s q[117]; s q[117]; s q[117]; h q[117]; measure q[117] -> rec[246]; h q[117]; s q[117]; // decomposed MY
s q[118]; s q[118]; s q[118]; h q[118]; measure q[118] -> rec[247]; h q[118]; s q[118]; // decomposed MY
s q[120]; s q[120]; s q[120]; h q[120]; measure q[120] -> rec[248]; h q[120]; s q[120]; // decomposed MY
s q[121]; s q[121]; s q[121]; h q[121]; measure q[121] -> rec[249]; h q[121]; s q[121]; // decomposed MY
s q[123]; s q[123]; s q[123]; h q[123]; measure q[123] -> rec[250]; h q[123]; s q[123]; // decomposed MY
s q[124]; s q[124]; s q[124]; h q[124]; measure q[124] -> rec[251]; h q[124]; s q[124]; // decomposed MY
s q[126]; s q[126]; s q[126]; h q[126]; measure q[126] -> rec[252]; h q[126]; s q[126]; // decomposed MY
s q[127]; s q[127]; s q[127]; h q[127]; measure q[127] -> rec[253]; h q[127]; s q[127]; // decomposed MY
s q[129]; s q[129]; s q[129]; h q[129]; measure q[129] -> rec[254]; h q[129]; s q[129]; // decomposed MY
s q[130]; s q[130]; s q[130]; h q[130]; measure q[130] -> rec[255]; h q[130]; s q[130]; // decomposed MY
s q[132]; s q[132]; s q[132]; h q[132]; measure q[132] -> rec[256]; h q[132]; s q[132]; // decomposed MY
s q[134]; s q[134]; s q[134]; h q[134]; measure q[134] -> rec[257]; h q[134]; s q[134]; // decomposed MY
s q[135]; s q[135]; s q[135]; h q[135]; measure q[135] -> rec[258]; h q[135]; s q[135]; // decomposed MY
s q[137]; s q[137]; s q[137]; h q[137]; measure q[137] -> rec[259]; h q[137]; s q[137]; // decomposed MY
s q[138]; s q[138]; s q[138]; h q[138]; measure q[138] -> rec[260]; h q[138]; s q[138]; // decomposed MY
s q[140]; s q[140]; s q[140]; h q[140]; measure q[140] -> rec[261]; h q[140]; s q[140]; // decomposed MY
s q[141]; s q[141]; s q[141]; h q[141]; measure q[141] -> rec[262]; h q[141]; s q[141]; // decomposed MY
s q[143]; s q[143]; s q[143]; h q[143]; measure q[143] -> rec[263]; h q[143]; s q[143]; // decomposed MY
s q[144]; s q[144]; s q[144]; h q[144]; measure q[144] -> rec[264]; h q[144]; s q[144]; // decomposed MY
s q[146]; s q[146]; s q[146]; h q[146]; measure q[146] -> rec[265]; h q[146]; s q[146]; // decomposed MY
s q[147]; s q[147]; s q[147]; h q[147]; measure q[147] -> rec[266]; h q[147]; s q[147]; // decomposed MY
s q[148]; s q[148]; s q[148]; h q[148]; measure q[148] -> rec[267]; h q[148]; s q[148]; // decomposed MY
s q[150]; s q[150]; s q[150]; h q[150]; measure q[150] -> rec[268]; h q[150]; s q[150]; // decomposed MY
s q[151]; s q[151]; s q[151]; h q[151]; measure q[151] -> rec[269]; h q[151]; s q[151]; // decomposed MY
s q[153]; s q[153]; s q[153]; h q[153]; measure q[153] -> rec[270]; h q[153]; s q[153]; // decomposed MY
s q[154]; s q[154]; s q[154]; h q[154]; measure q[154] -> rec[271]; h q[154]; s q[154]; // decomposed MY
s q[156]; s q[156]; s q[156]; h q[156]; measure q[156] -> rec[272]; h q[156]; s q[156]; // decomposed MY
s q[157]; s q[157]; s q[157]; h q[157]; measure q[157] -> rec[273]; h q[157]; s q[157]; // decomposed MY
s q[159]; s q[159]; s q[159]; h q[159]; measure q[159] -> rec[274]; h q[159]; s q[159]; // decomposed MY
s q[160]; s q[160]; s q[160]; h q[160]; measure q[160] -> rec[275]; h q[160]; s q[160]; // decomposed MY
s q[162]; s q[162]; s q[162]; h q[162]; measure q[162] -> rec[276]; h q[162]; s q[162]; // decomposed MY
s q[163]; s q[163]; s q[163]; h q[163]; measure q[163] -> rec[277]; h q[163]; s q[163]; // decomposed MY
s q[165]; s q[165]; s q[165]; h q[165]; measure q[165] -> rec[278]; h q[165]; s q[165]; // decomposed MY
s q[166]; s q[166]; s q[166]; h q[166]; measure q[166] -> rec[279]; h q[166]; s q[166]; // decomposed MY
s q[168]; s q[168]; s q[168]; h q[168]; measure q[168] -> rec[280]; h q[168]; s q[168]; // decomposed MY
s q[169]; s q[169]; s q[169]; h q[169]; measure q[169] -> rec[281]; h q[169]; s q[169]; // decomposed MY
s q[171]; s q[171]; s q[171]; h q[171]; measure q[171] -> rec[282]; h q[171]; s q[171]; // decomposed MY
s q[172]; s q[172]; s q[172]; h q[172]; measure q[172] -> rec[283]; h q[172]; s q[172]; // decomposed MY
s q[174]; s q[174]; s q[174]; h q[174]; measure q[174] -> rec[284]; h q[174]; s q[174]; // decomposed MY
s q[176]; s q[176]; s q[176]; h q[176]; measure q[176] -> rec[285]; h q[176]; s q[176]; // decomposed MY
s q[177]; s q[177]; s q[177]; h q[177]; measure q[177] -> rec[286]; h q[177]; s q[177]; // decomposed MY
s q[179]; s q[179]; s q[179]; h q[179]; measure q[179] -> rec[287]; h q[179]; s q[179]; // decomposed MY
s q[180]; s q[180]; s q[180]; h q[180]; measure q[180] -> rec[288]; h q[180]; s q[180]; // decomposed MY
s q[182]; s q[182]; s q[182]; h q[182]; measure q[182] -> rec[289]; h q[182]; s q[182]; // decomposed MY
s q[183]; s q[183]; s q[183]; h q[183]; measure q[183] -> rec[290]; h q[183]; s q[183]; // decomposed MY
s q[185]; s q[185]; s q[185]; h q[185]; measure q[185] -> rec[291]; h q[185]; s q[185]; // decomposed MY
s q[186]; s q[186]; s q[186]; h q[186]; measure q[186] -> rec[292]; h q[186]; s q[186]; // decomposed MY
s q[187]; s q[187]; s q[187]; h q[187]; measure q[187] -> rec[293]; h q[187]; s q[187]; // decomposed MY
s q[189]; s q[189]; s q[189]; h q[189]; measure q[189] -> rec[294]; h q[189]; s q[189]; // decomposed MY
s q[190]; s q[190]; s q[190]; h q[190]; measure q[190] -> rec[295]; h q[190]; s q[190]; // decomposed MY
s q[192]; s q[192]; s q[192]; h q[192]; measure q[192] -> rec[296]; h q[192]; s q[192]; // decomposed MY
s q[193]; s q[193]; s q[193]; h q[193]; measure q[193] -> rec[297]; h q[193]; s q[193]; // decomposed MY
s q[195]; s q[195]; s q[195]; h q[195]; measure q[195] -> rec[298]; h q[195]; s q[195]; // decomposed MY
s q[196]; s q[196]; s q[196]; h q[196]; measure q[196] -> rec[299]; h q[196]; s q[196]; // decomposed MY
s q[198]; s q[198]; s q[198]; h q[198]; measure q[198] -> rec[300]; h q[198]; s q[198]; // decomposed MY
s q[199]; s q[199]; s q[199]; h q[199]; measure q[199] -> rec[301]; h q[199]; s q[199]; // decomposed MY
s q[201]; s q[201]; s q[201]; h q[201]; measure q[201] -> rec[302]; h q[201]; s q[201]; // decomposed MY
s q[202]; s q[202]; s q[202]; h q[202]; measure q[202] -> rec[303]; h q[202]; s q[202]; // decomposed MY
s q[204]; s q[204]; s q[204]; h q[204]; measure q[204] -> rec[304]; h q[204]; s q[204]; // decomposed MY
s q[205]; s q[205]; s q[205]; h q[205]; measure q[205] -> rec[305]; h q[205]; s q[205]; // decomposed MY
s q[207]; s q[207]; s q[207]; h q[207]; measure q[207] -> rec[306]; h q[207]; s q[207]; // decomposed MY
s q[209]; s q[209]; s q[209]; h q[209]; measure q[209] -> rec[307]; h q[209]; s q[209]; // decomposed MY
s q[210]; s q[210]; s q[210]; h q[210]; measure q[210] -> rec[308]; h q[210]; s q[210]; // decomposed MY
s q[212]; s q[212]; s q[212]; h q[212]; measure q[212] -> rec[309]; h q[212]; s q[212]; // decomposed MY
s q[213]; s q[213]; s q[213]; h q[213]; measure q[213] -> rec[310]; h q[213]; s q[213]; // decomposed MY
s q[215]; s q[215]; s q[215]; h q[215]; measure q[215] -> rec[311]; h q[215]; s q[215]; // decomposed MY
s q[216]; s q[216]; s q[216]; h q[216]; measure q[216] -> rec[312]; h q[216]; s q[216]; // decomposed MY
s q[217]; s q[217]; s q[217]; h q[217]; measure q[217] -> rec[313]; h q[217]; s q[217]; // decomposed MY
s q[219]; s q[219]; s q[219]; h q[219]; measure q[219] -> rec[314]; h q[219]; s q[219]; // decomposed MY
s q[220]; s q[220]; s q[220]; h q[220]; measure q[220] -> rec[315]; h q[220]; s q[220]; // decomposed MY
s q[222]; s q[222]; s q[222]; h q[222]; measure q[222] -> rec[316]; h q[222]; s q[222]; // decomposed MY
s q[223]; s q[223]; s q[223]; h q[223]; measure q[223] -> rec[317]; h q[223]; s q[223]; // decomposed MY
s q[225]; s q[225]; s q[225]; h q[225]; measure q[225] -> rec[318]; h q[225]; s q[225]; // decomposed MY
s q[226]; s q[226]; s q[226]; h q[226]; measure q[226] -> rec[319]; h q[226]; s q[226]; // decomposed MY
s q[228]; s q[228]; s q[228]; h q[228]; measure q[228] -> rec[320]; h q[228]; s q[228]; // decomposed MY
s q[229]; s q[229]; s q[229]; h q[229]; measure q[229] -> rec[321]; h q[229]; s q[229]; // decomposed MY
s q[231]; s q[231]; s q[231]; h q[231]; measure q[231] -> rec[322]; h q[231]; s q[231]; // decomposed MY
s q[233]; s q[233]; s q[233]; h q[233]; measure q[233] -> rec[323]; h q[233]; s q[233]; // decomposed MY
s q[234]; s q[234]; s q[234]; h q[234]; measure q[234] -> rec[324]; h q[234]; s q[234]; // decomposed MY
s q[236]; s q[236]; s q[236]; h q[236]; measure q[236] -> rec[325]; h q[236]; s q[236]; // decomposed MY
s q[237]; s q[237]; s q[237]; h q[237]; measure q[237] -> rec[326]; h q[237]; s q[237]; // decomposed MY
s q[238]; s q[238]; s q[238]; h q[238]; measure q[238] -> rec[327]; h q[238]; s q[238]; // decomposed MY
s q[240]; s q[240]; s q[240]; h q[240]; measure q[240] -> rec[328]; h q[240]; s q[240]; // decomposed MY
s q[241]; s q[241]; s q[241]; h q[241]; measure q[241] -> rec[329]; h q[241]; s q[241]; // decomposed MY
s q[243]; s q[243]; s q[243]; h q[243]; measure q[243] -> rec[330]; h q[243]; s q[243]; // decomposed MY
s q[244]; s q[244]; s q[244]; h q[244]; measure q[244] -> rec[331]; h q[244]; s q[244]; // decomposed MY
s q[246]; s q[246]; s q[246]; h q[246]; measure q[246] -> rec[332]; h q[246]; s q[246]; // decomposed MY
s q[248]; s q[248]; s q[248]; h q[248]; measure q[248] -> rec[333]; h q[248]; s q[248]; // decomposed MY
s q[249]; s q[249]; s q[249]; h q[249]; measure q[249] -> rec[334]; h q[249]; s q[249]; // decomposed MY
s q[250]; s q[250]; s q[250]; h q[250]; measure q[250] -> rec[335]; h q[250]; s q[250]; // decomposed MY
s q[252]; s q[252]; s q[252]; h q[252]; measure q[252] -> rec[336]; h q[252]; s q[252]; // decomposed MY

OPENQASM 2.0;
include "qelib1.inc";

qreg q[274];
creg rec[241];

reset q[1];
reset q[3];
reset q[5];
reset q[7];
reset q[9];
reset q[11];
reset q[13];
reset q[15];
reset q[17];
reset q[19];
reset q[21];
reset q[24];
reset q[26];
reset q[28];
reset q[30];
reset q[32];
reset q[34];
reset q[36];
reset q[38];
reset q[40];
reset q[42];
reset q[44];
reset q[47];
reset q[49];
reset q[51];
reset q[53];
reset q[55];
reset q[57];
reset q[59];
reset q[61];
reset q[63];
reset q[65];
reset q[67];
reset q[70];
reset q[72];
reset q[74];
reset q[76];
reset q[78];
reset q[80];
reset q[82];
reset q[84];
reset q[86];
reset q[88];
reset q[90];
reset q[93];
reset q[95];
reset q[97];
reset q[99];
reset q[101];
reset q[103];
reset q[105];
reset q[107];
reset q[109];
reset q[111];
reset q[113];
reset q[116];
reset q[118];
reset q[120];
reset q[122];
reset q[124];
reset q[126];
reset q[128];
reset q[130];
reset q[132];
reset q[134];
reset q[136];
reset q[139];
reset q[141];
reset q[143];
reset q[145];
reset q[147];
reset q[149];
reset q[151];
reset q[153];
reset q[155];
reset q[157];
reset q[159];
reset q[162];
reset q[164];
reset q[166];
reset q[168];
reset q[170];
reset q[172];
reset q[174];
reset q[176];
reset q[178];
reset q[180];
reset q[182];
reset q[185];
reset q[187];
reset q[189];
reset q[191];
reset q[193];
reset q[195];
reset q[197];
reset q[199];
reset q[201];
reset q[203];
reset q[205];
reset q[208];
reset q[210];
reset q[212];
reset q[214];
reset q[216];
reset q[218];
reset q[220];
reset q[222];
reset q[224];
reset q[226];
reset q[228];
reset q[231];
reset q[233];
reset q[235];
reset q[237];
reset q[239];
reset q[241];
reset q[243];
reset q[245];
reset q[247];
reset q[249];
reset q[251];
reset q[2];
reset q[6];
reset q[10];
reset q[14];
reset q[18];
reset q[25];
reset q[27];
reset q[29];
reset q[31];
reset q[33];
reset q[35];
reset q[37];
reset q[39];
reset q[41];
reset q[43];
reset q[45];
reset q[46];
reset q[48];
reset q[50];
reset q[52];
reset q[54];
reset q[56];
reset q[58];
reset q[60];
reset q[62];
reset q[64];
reset q[66];
reset q[71];
reset q[73];
reset q[75];
reset q[77];
reset q[79];
reset q[81];
reset q[83];
reset q[85];
reset q[87];
reset q[89];
reset q[91];
reset q[92];
reset q[94];
reset q[96];
reset q[98];
reset q[100];
reset q[102];
reset q[104];
reset q[106];
reset q[108];
reset q[110];
reset q[112];
reset q[117];
reset q[119];
reset q[121];
reset q[123];
reset q[125];
reset q[127];
reset q[129];
reset q[131];
reset q[133];
reset q[135];
reset q[137];
reset q[138];
reset q[140];
reset q[142];
reset q[144];
reset q[146];
reset q[148];
reset q[150];
reset q[152];
reset q[154];
reset q[156];
reset q[158];
reset q[163];
reset q[165];
reset q[167];
reset q[169];
reset q[171];
reset q[173];
reset q[175];
reset q[177];
reset q[179];
reset q[181];
reset q[183];
reset q[184];
reset q[186];
reset q[188];
reset q[190];
reset q[192];
reset q[194];
reset q[196];
reset q[198];
reset q[200];
reset q[202];
reset q[204];
reset q[209];
reset q[211];
reset q[213];
reset q[215];
reset q[217];
reset q[219];
reset q[221];
reset q[223];
reset q[225];
reset q[227];
reset q[229];
reset q[230];
reset q[232];
reset q[234];
reset q[236];
reset q[238];
reset q[240];
reset q[242];
reset q[244];
reset q[246];
reset q[248];
reset q[250];
reset q[257];
reset q[261];
reset q[265];
reset q[269];
reset q[273];
barrier q;

h q[2];
h q[6];
h q[10];
h q[14];
h q[18];
h q[27];
h q[31];
h q[35];
h q[39];
h q[43];
h q[48];
h q[52];
h q[56];
h q[60];
h q[64];
h q[73];
h q[77];
h q[81];
h q[85];
h q[89];
h q[94];
h q[98];
h q[102];
h q[106];
h q[110];
h q[119];
h q[123];
h q[127];
h q[131];
h q[135];
h q[140];
h q[144];
h q[148];
h q[152];
h q[156];
h q[165];
h q[169];
h q[173];
h q[177];
h q[181];
h q[186];
h q[190];
h q[194];
h q[198];
h q[202];
h q[211];
h q[215];
h q[219];
h q[223];
h q[227];
h q[232];
h q[236];
h q[240];
h q[244];
h q[248];
h q[257];
h q[261];
h q[265];
h q[269];
h q[273];
barrier q;

cx q[2], q[3];
cx q[48], q[49];
cx q[94], q[95];
cx q[140], q[141];
cx q[186], q[187];
cx q[232], q[233];
cx q[27], q[28];
cx q[73], q[74];
cx q[119], q[120];
cx q[165], q[166];
cx q[211], q[212];
cx q[6], q[7];
cx q[52], q[53];
cx q[98], q[99];
cx q[144], q[145];
cx q[190], q[191];
cx q[236], q[237];
cx q[31], q[32];
cx q[77], q[78];
cx q[123], q[124];
cx q[169], q[170];
cx q[215], q[216];
cx q[10], q[11];
cx q[56], q[57];
cx q[102], q[103];
cx q[148], q[149];
cx q[194], q[195];
cx q[240], q[241];
cx q[35], q[36];
cx q[81], q[82];
cx q[127], q[128];
cx q[173], q[174];
cx q[219], q[220];
cx q[14], q[15];
cx q[60], q[61];
cx q[106], q[107];
cx q[152], q[153];
cx q[198], q[199];
cx q[244], q[245];
cx q[39], q[40];
cx q[85], q[86];
cx q[131], q[132];
cx q[177], q[178];
cx q[223], q[224];
cx q[18], q[19];
cx q[64], q[65];
cx q[110], q[111];
cx q[156], q[157];
cx q[202], q[203];
cx q[248], q[249];
cx q[43], q[44];
cx q[89], q[90];
cx q[135], q[136];
cx q[181], q[182];
cx q[227], q[228];
cx q[47], q[46];
cx q[93], q[92];
cx q[139], q[138];
cx q[185], q[184];
cx q[231], q[230];
cx q[26], q[25];
cx q[72], q[71];
cx q[118], q[117];
cx q[164], q[163];
cx q[210], q[209];
cx q[51], q[50];
cx q[97], q[96];
cx q[143], q[142];
cx q[189], q[188];
cx q[235], q[234];
cx q[30], q[29];
cx q[76], q[75];
cx q[122], q[121];
cx q[168], q[167];
cx q[214], q[213];
cx q[55], q[54];
cx q[101], q[100];
cx q[147], q[146];
cx q[193], q[192];
cx q[239], q[238];
cx q[34], q[33];
cx q[80], q[79];
cx q[126], q[125];
cx q[172], q[171];
cx q[218], q[217];
cx q[59], q[58];
cx q[105], q[104];
cx q[151], q[150];
cx q[197], q[196];
cx q[243], q[242];
cx q[38], q[37];
cx q[84], q[83];
cx q[130], q[129];
cx q[176], q[175];
cx q[222], q[221];
cx q[63], q[62];
cx q[109], q[108];
cx q[155], q[154];
cx q[201], q[200];
cx q[247], q[246];
cx q[42], q[41];
cx q[88], q[87];
cx q[134], q[133];
cx q[180], q[179];
cx q[226], q[225];
cx q[67], q[66];
cx q[113], q[112];
cx q[159], q[158];
cx q[205], q[204];
cx q[251], q[250];
barrier q;

cx q[2], q[1];
cx q[48], q[47];
cx q[94], q[93];
cx q[140], q[139];
cx q[186], q[185];
cx q[232], q[231];
cx q[27], q[26];
cx q[73], q[72];
cx q[119], q[118];
cx q[165], q[164];
cx q[211], q[210];
cx q[6], q[5];
cx q[52], q[51];
cx q[98], q[97];
cx q[144], q[143];
cx q[190], q[189];
cx q[236], q[235];
cx q[31], q[30];
cx q[77], q[76];
cx q[123], q[122];
cx q[169], q[168];
cx q[215], q[214];
cx q[10], q[9];
cx q[56], q[55];
cx q[102], q[101];
cx q[148], q[147];
cx q[194], q[193];
cx q[240], q[239];
cx q[35], q[34];
cx q[81], q[80];
cx q[127], q[126];
cx q[173], q[172];
cx q[219], q[218];
cx q[14], q[13];
cx q[60], q[59];
cx q[106], q[105];
cx q[152], q[151];
cx q[198], q[197];
cx q[244], q[243];
cx q[39], q[38];
cx q[85], q[84];
cx q[131], q[130];
cx q[177], q[176];
cx q[223], q[222];
cx q[18], q[17];
cx q[64], q[63];
cx q[110], q[109];
cx q[156], q[155];
cx q[202], q[201];
cx q[248], q[247];
cx q[43], q[42];
cx q[89], q[88];
cx q[135], q[134];
cx q[181], q[180];
cx q[227], q[226];
cx q[24], q[46];
cx q[70], q[92];
cx q[116], q[138];
cx q[162], q[184];
cx q[208], q[230];
cx q[3], q[25];
cx q[49], q[71];
cx q[95], q[117];
cx q[141], q[163];
cx q[187], q[209];
cx q[28], q[50];
cx q[74], q[96];
cx q[120], q[142];
cx q[166], q[188];
cx q[212], q[234];
cx q[7], q[29];
cx q[53], q[75];
cx q[99], q[121];
cx q[145], q[167];
cx q[191], q[213];
cx q[32], q[54];
cx q[78], q[100];
cx q[124], q[146];
cx q[170], q[192];
cx q[216], q[238];
cx q[11], q[33];
cx q[57], q[79];
cx q[103], q[125];
cx q[149], q[171];
cx q[195], q[217];
cx q[36], q[58];
cx q[82], q[104];
cx q[128], q[150];
cx q[174], q[196];
cx q[220], q[242];
cx q[15], q[37];
cx q[61], q[83];
cx q[107], q[129];
cx q[153], q[175];
cx q[199], q[221];
cx q[40], q[62];
cx q[86], q[108];
cx q[132], q[154];
cx q[178], q[200];
cx q[224], q[246];
cx q[19], q[41];
cx q[65], q[87];
cx q[111], q[133];
cx q[157], q[179];
cx q[203], q[225];
cx q[44], q[66];
cx q[90], q[112];
cx q[136], q[158];
cx q[182], q[204];
cx q[228], q[250];
barrier q;

cx q[48], q[26];
cx q[94], q[72];
cx q[140], q[118];
cx q[186], q[164];
cx q[232], q[210];
cx q[27], q[5];
cx q[73], q[51];
cx q[119], q[97];
cx q[165], q[143];
cx q[211], q[189];
cx q[257], q[235];
cx q[52], q[30];
cx q[98], q[76];
cx q[144], q[122];
cx q[190], q[168];
cx q[236], q[214];
cx q[31], q[9];
cx q[77], q[55];
cx q[123], q[101];
cx q[169], q[147];
cx q[215], q[193];
cx q[261], q[239];
cx q[56], q[34];
cx q[102], q[80];
cx q[148], q[126];
cx q[194], q[172];
cx q[240], q[218];
cx q[35], q[13];
cx q[81], q[59];
cx q[127], q[105];
cx q[173], q[151];
cx q[219], q[197];
cx q[265], q[243];
cx q[60], q[38];
cx q[106], q[84];
cx q[152], q[130];
cx q[198], q[176];
cx q[244], q[222];
cx q[39], q[17];
cx q[85], q[63];
cx q[131], q[109];
cx q[177], q[155];
cx q[223], q[201];
cx q[269], q[247];
cx q[64], q[42];
cx q[110], q[88];
cx q[156], q[134];
cx q[202], q[180];
cx q[248], q[226];
cx q[43], q[21];
cx q[89], q[67];
cx q[135], q[113];
cx q[181], q[159];
cx q[227], q[205];
cx q[273], q[251];
cx q[24], q[25];
cx q[70], q[71];
cx q[116], q[117];
cx q[162], q[163];
cx q[208], q[209];
cx q[49], q[50];
cx q[95], q[96];
cx q[141], q[142];
cx q[187], q[188];
cx q[233], q[234];
cx q[28], q[29];
cx q[74], q[75];
cx q[120], q[121];
cx q[166], q[167];
cx q[212], q[213];
cx q[53], q[54];
cx q[99], q[100];
cx q[145], q[146];
cx q[191], q[192];
cx q[237], q[238];
cx q[32], q[33];
cx q[78], q[79];
cx q[124], q[125];
cx q[170], q[171];
cx q[216], q[217];
cx q[57], q[58];
cx q[103], q[104];
cx q[149], q[150];
cx q[195], q[196];
cx q[241], q[242];
cx q[36], q[37];
cx q[82], q[83];
cx q[128], q[129];
cx q[174], q[175];
cx q[220], q[221];
cx q[61], q[62];
cx q[107], q[108];
cx q[153], q[154];
cx q[199], q[200];
cx q[245], q[246];
cx q[40], q[41];
cx q[86], q[87];
cx q[132], q[133];
cx q[178], q[179];
cx q[224], q[225];
cx q[65], q[66];
cx q[111], q[112];
cx q[157], q[158];
cx q[203], q[204];
cx q[249], q[250];
cx q[44], q[45];
cx q[90], q[91];
cx q[136], q[137];
cx q[182], q[183];
cx q[228], q[229];
barrier q;

cx q[48], q[24];
cx q[94], q[70];
cx q[140], q[116];
cx q[186], q[162];
cx q[232], q[208];
cx q[27], q[3];
cx q[73], q[49];
cx q[119], q[95];
cx q[165], q[141];
cx q[211], q[187];
cx q[257], q[233];
cx q[52], q[28];
cx q[98], q[74];
cx q[144], q[120];
cx q[190], q[166];
cx q[236], q[212];
cx q[31], q[7];
cx q[77], q[53];
cx q[123], q[99];
cx q[169], q[145];
cx q[215], q[191];
cx q[261], q[237];
cx q[56], q[32];
cx q[102], q[78];
cx q[148], q[124];
cx q[194], q[170];
cx q[240], q[216];
cx q[35], q[11];
cx q[81], q[57];
cx q[127], q[103];
cx q[173], q[149];
cx q[219], q[195];
cx q[265], q[241];
cx q[60], q[36];
cx q[106], q[82];
cx q[152], q[128];
cx q[198], q[174];
cx q[244], q[220];
cx q[39], q[15];
cx q[85], q[61];
cx q[131], q[107];
cx q[177], q[153];
cx q[223], q[199];
cx q[269], q[245];
cx q[64], q[40];
cx q[110], q[86];
cx q[156], q[132];
cx q[202], q[178];
cx q[248], q[224];
cx q[43], q[19];
cx q[89], q[65];
cx q[135], q[111];
cx q[181], q[157];
cx q[227], q[203];
cx q[273], q[249];
cx q[1], q[25];
cx q[47], q[71];
cx q[93], q[117];
cx q[139], q[163];
cx q[185], q[209];
cx q[26], q[50];
cx q[72], q[96];
cx q[118], q[142];
cx q[164], q[188];
cx q[210], q[234];
cx q[5], q[29];
cx q[51], q[75];
cx q[97], q[121];
cx q[143], q[167];
cx q[189], q[213];
cx q[30], q[54];
cx q[76], q[100];
cx q[122], q[146];
cx q[168], q[192];
cx q[214], q[238];
cx q[9], q[33];
cx q[55], q[79];
cx q[101], q[125];
cx q[147], q[171];
cx q[193], q[217];
cx q[34], q[58];
cx q[80], q[104];
cx q[126], q[150];
cx q[172], q[196];
cx q[218], q[242];
cx q[13], q[37];
cx q[59], q[83];
cx q[105], q[129];
cx q[151], q[175];
cx q[197], q[221];
cx q[38], q[62];
cx q[84], q[108];
cx q[130], q[154];
cx q[176], q[200];
cx q[222], q[246];
cx q[17], q[41];
cx q[63], q[87];
cx q[109], q[133];
cx q[155], q[179];
cx q[201], q[225];
cx q[42], q[66];
cx q[88], q[112];
cx q[134], q[158];
cx q[180], q[204];
cx q[226], q[250];
cx q[21], q[45];
cx q[67], q[91];
cx q[113], q[137];
cx q[159], q[183];
cx q[205], q[229];
barrier q;

h q[2];
h q[6];
h q[10];
h q[14];
h q[18];
h q[27];
h q[31];
h q[35];
h q[39];
h q[43];
h q[48];
h q[52];
h q[56];
h q[60];
h q[64];
h q[73];
h q[77];
h q[81];
h q[85];
h q[89];
h q[94];
h q[98];
h q[102];
h q[106];
h q[110];
h q[119];
h q[123];
h q[127];
h q[131];
h q[135];
h q[140];
h q[144];
h q[148];
h q[152];
h q[156];
h q[165];
h q[169];
h q[173];
h q[177];
h q[181];
h q[186];
h q[190];
h q[194];
h q[198];
h q[202];
h q[211];
h q[215];
h q[219];
h q[223];
h q[227];
h q[232];
h q[236];
h q[240];
h q[244];
h q[248];
h q[257];
h q[261];
h q[265];
h q[269];
h q[273];
barrier q;

measure q[2] -> rec[0]; reset q[2]; // decomposed MR
measure q[6] -> rec[1]; reset q[6]; // decomposed MR
measure q[10] -> rec[2]; reset q[10]; // decomposed MR
measure q[14] -> rec[3]; reset q[14]; // decomposed MR
measure q[18] -> rec[4]; reset q[18]; // decomposed MR
measure q[25] -> rec[5]; reset q[25]; // decomposed MR
measure q[27] -> rec[6]; reset q[27]; // decomposed MR
measure q[29] -> rec[7]; reset q[29]; // decomposed MR
measure q[31] -> rec[8]; reset q[31]; // decomposed MR
measure q[33] -> rec[9]; reset q[33]; // decomposed MR
measure q[35] -> rec[10]; reset q[35]; // decomposed MR
measure q[37] -> rec[11]; reset q[37]; // decomposed MR
measure q[39] -> rec[12]; reset q[39]; // decomposed MR
measure q[41] -> rec[13]; reset q[41]; // decomposed MR
measure q[43] -> rec[14]; reset q[43]; // decomposed MR
measure q[45] -> rec[15]; reset q[45]; // decomposed MR
measure q[46] -> rec[16]; reset q[46]; // decomposed MR
measure q[48] -> rec[17]; reset q[48]; // decomposed MR
measure q[50] -> rec[18]; reset q[50]; // decomposed MR
measure q[52] -> rec[19]; reset q[52]; // decomposed MR
measure q[54] -> rec[20]; reset q[54]; // decomposed MR
measure q[56] -> rec[21]; reset q[56]; // decomposed MR
measure q[58] -> rec[22]; reset q[58]; // decomposed MR
measure q[60] -> rec[23]; reset q[60]; // decomposed MR
measure q[62] -> rec[24]; reset q[62]; // decomposed MR
measure q[64] -> rec[25]; reset q[64]; // decomposed MR
measure q[66] -> rec[26]; reset q[66]; // decomposed MR
measure q[71] -> rec[27]; reset q[71]; // decomposed MR
measure q[73] -> rec[28]; reset q[73]; // decomposed MR
measure q[75] -> rec[29]; reset q[75]; // decomposed MR
measure q[77] -> rec[30]; reset q[77]; // decomposed MR
measure q[79] -> rec[31]; reset q[79]; // decomposed MR
measure q[81] -> rec[32]; reset q[81]; // decomposed MR
measure q[83] -> rec[33]; reset q[83]; // decomposed MR
measure q[85] -> rec[34]; reset q[85]; // decomposed MR
measure q[87] -> rec[35]; reset q[87]; // decomposed MR
measure q[89] -> rec[36]; reset q[89]; // decomposed MR
measure q[91] -> rec[37]; reset q[91]; // decomposed MR
measure q[92] -> rec[38]; reset q[92]; // decomposed MR
measure q[94] -> rec[39]; reset q[94]; // decomposed MR
measure q[96] -> rec[40]; reset q[96]; // decomposed MR
measure q[98] -> rec[41]; reset q[98]; // decomposed MR
measure q[100] -> rec[42]; reset q[100]; // decomposed MR
measure q[102] -> rec[43]; reset q[102]; // decomposed MR
measure q[104] -> rec[44]; reset q[104]; // decomposed MR
measure q[106] -> rec[45]; reset q[106]; // decomposed MR
measure q[108] -> rec[46]; reset q[108]; // decomposed MR
measure q[110] -> rec[47]; reset q[110]; // decomposed MR
measure q[112] -> rec[48]; reset q[112]; // decomposed MR
measure q[117] -> rec[49]; reset q[117]; // decomposed MR
measure q[119] -> rec[50]; reset q[119]; // decomposed MR
measure q[121] -> rec[51]; reset q[121]; // decomposed MR
measure q[123] -> rec[52]; reset q[123]; // decomposed MR
measure q[125] -> rec[53]; reset q[125]; // decomposed MR
measure q[127] -> rec[54]; reset q[127]; // decomposed MR
measure q[129] -> rec[55]; reset q[129]; // decomposed MR
measure q[131] -> rec[56]; reset q[131]; // decomposed MR
measure q[133] -> rec[57]; reset q[133]; // decomposed MR
measure q[135] -> rec[58]; reset q[135]; // decomposed MR
measure q[137] -> rec[59]; reset q[137]; // decomposed MR
measure q[138] -> rec[60]; reset q[138]; // decomposed MR
measure q[140] -> rec[61]; reset q[140]; // decomposed MR
measure q[142] -> rec[62]; reset q[142]; // decomposed MR
measure q[144] -> rec[63]; reset q[144]; // decomposed MR
measure q[146] -> rec[64]; reset q[146]; // decomposed MR
measure q[148] -> rec[65]; reset q[148]; // decomposed MR
measure q[150] -> rec[66]; reset q[150]; // decomposed MR
measure q[152] -> rec[67]; reset q[152]; // decomposed MR
measure q[154] -> rec[68]; reset q[154]; // decomposed MR
measure q[156] -> rec[69]; reset q[156]; // decomposed MR
measure q[158] -> rec[70]; reset q[158]; // decomposed MR
measure q[163] -> rec[71]; reset q[163]; // decomposed MR
measure q[165] -> rec[72]; reset q[165]; // decomposed MR
measure q[167] -> rec[73]; reset q[167]; // decomposed MR
measure q[169] -> rec[74]; reset q[169]; // decomposed MR
measure q[171] -> rec[75]; reset q[171]; // decomposed MR
measure q[173] -> rec[76]; reset q[173]; // decomposed MR
measure q[175] -> rec[77]; reset q[175]; // decomposed MR
measure q[177] -> rec[78]; reset q[177]; // decomposed MR
measure q[179] -> rec[79]; reset q[179]; // decomposed MR
measure q[181] -> rec[80]; reset q[181]; // decomposed MR
measure q[183] -> rec[81]; reset q[183]; // decomposed MR
measure q[184] -> rec[82]; reset q[184]; // decomposed MR
measure q[186] -> rec[83]; reset q[186]; // decomposed MR
measure q[188] -> rec[84]; reset q[188]; // decomposed MR
measure q[190] -> rec[85]; reset q[190]; // decomposed MR
measure q[192] -> rec[86]; reset q[192]; // decomposed MR
measure q[194] -> rec[87]; reset q[194]; // decomposed MR
measure q[196] -> rec[88]; reset q[196]; // decomposed MR
measure q[198] -> rec[89]; reset q[198]; // decomposed MR
measure q[200] -> rec[90]; reset q[200]; // decomposed MR
measure q[202] -> rec[91]; reset q[202]; // decomposed MR
measure q[204] -> rec[92]; reset q[204]; // decomposed MR
measure q[209] -> rec[93]; reset q[209]; // decomposed MR
measure q[211] -> rec[94]; reset q[211]; // decomposed MR
measure q[213] -> rec[95]; reset q[213]; // decomposed MR
measure q[215] -> rec[96]; reset q[215]; // decomposed MR
measure q[217] -> rec[97]; reset q[217]; // decomposed MR
measure q[219] -> rec[98]; reset q[219]; // decomposed MR
measure q[221] -> rec[99]; reset q[221]; // decomposed MR
measure q[223] -> rec[100]; reset q[223]; // decomposed MR
measure q[225] -> rec[101]; reset q[225]; // decomposed MR
measure q[227] -> rec[102]; reset q[227]; // decomposed MR
measure q[229] -> rec[103]; reset q[229]; // decomposed MR
measure q[230] -> rec[104]; reset q[230]; // decomposed MR
measure q[232] -> rec[105]; reset q[232]; // decomposed MR
measure q[234] -> rec[106]; reset q[234]; // decomposed MR
measure q[236] -> rec[107]; reset q[236]; // decomposed MR
measure q[238] -> rec[108]; reset q[238]; // decomposed MR
measure q[240] -> rec[109]; reset q[240]; // decomposed MR
measure q[242] -> rec[110]; reset q[242]; // decomposed MR
measure q[244] -> rec[111]; reset q[244]; // decomposed MR
measure q[246] -> rec[112]; reset q[246]; // decomposed MR
measure q[248] -> rec[113]; reset q[248]; // decomposed MR
measure q[250] -> rec[114]; reset q[250]; // decomposed MR
measure q[257] -> rec[115]; reset q[257]; // decomposed MR
measure q[261] -> rec[116]; reset q[261]; // decomposed MR
measure q[265] -> rec[117]; reset q[265]; // decomposed MR
measure q[269] -> rec[118]; reset q[269]; // decomposed MR
measure q[273] -> rec[119]; reset q[273]; // decomposed MR
measure q[1] -> rec[120];
measure q[3] -> rec[121];
measure q[5] -> rec[122];
measure q[7] -> rec[123];
measure q[9] -> rec[124];
measure q[11] -> rec[125];
measure q[13] -> rec[126];
measure q[15] -> rec[127];
measure q[17] -> rec[128];
measure q[19] -> rec[129];
measure q[21] -> rec[130];
measure q[24] -> rec[131];
measure q[26] -> rec[132];
measure q[28] -> rec[133];
measure q[30] -> rec[134];
measure q[32] -> rec[135];
measure q[34] -> rec[136];
measure q[36] -> rec[137];
measure q[38] -> rec[138];
measure q[40] -> rec[139];
measure q[42] -> rec[140];
measure q[44] -> rec[141];
measure q[47] -> rec[142];
measure q[49] -> rec[143];
measure q[51] -> rec[144];
measure q[53] -> rec[145];
measure q[55] -> rec[146];
measure q[57] -> rec[147];
measure q[59] -> rec[148];
measure q[61] -> rec[149];
measure q[63] -> rec[150];
measure q[65] -> rec[151];
measure q[67] -> rec[152];
measure q[70] -> rec[153];
measure q[72] -> rec[154];
measure q[74] -> rec[155];
measure q[76] -> rec[156];
measure q[78] -> rec[157];
measure q[80] -> rec[158];
measure q[82] -> rec[159];
measure q[84] -> rec[160];
measure q[86] -> rec[161];
measure q[88] -> rec[162];
measure q[90] -> rec[163];
measure q[93] -> rec[164];
measure q[95] -> rec[165];
measure q[97] -> rec[166];
measure q[99] -> rec[167];
measure q[101] -> rec[168];
measure q[103] -> rec[169];
measure q[105] -> rec[170];
measure q[107] -> rec[171];
measure q[109] -> rec[172];
measure q[111] -> rec[173];
measure q[113] -> rec[174];
measure q[116] -> rec[175];
measure q[118] -> rec[176];
measure q[120] -> rec[177];
measure q[122] -> rec[178];
measure q[124] -> rec[179];
measure q[126] -> rec[180];
measure q[128] -> rec[181];
measure q[130] -> rec[182];
measure q[132] -> rec[183];
measure q[134] -> rec[184];
measure q[136] -> rec[185];
measure q[139] -> rec[186];
measure q[141] -> rec[187];
measure q[143] -> rec[188];
measure q[145] -> rec[189];
measure q[147] -> rec[190];
measure q[149] -> rec[191];
measure q[151] -> rec[192];
measure q[153] -> rec[193];
measure q[155] -> rec[194];
measure q[157] -> rec[195];
measure q[159] -> rec[196];
measure q[162] -> rec[197];
measure q[164] -> rec[198];
measure q[166] -> rec[199];
measure q[168] -> rec[200];
measure q[170] -> rec[201];
measure q[172] -> rec[202];
measure q[174] -> rec[203];
measure q[176] -> rec[204];
measure q[178] -> rec[205];
measure q[180] -> rec[206];
measure q[182] -> rec[207];
measure q[185] -> rec[208];
measure q[187] -> rec[209];
measure q[189] -> rec[210];
measure q[191] -> rec[211];
measure q[193] -> rec[212];
measure q[195] -> rec[213];
measure q[197] -> rec[214];
measure q[199] -> rec[215];
measure q[201] -> rec[216];
measure q[203] -> rec[217];
measure q[205] -> rec[218];
measure q[208] -> rec[219];
measure q[210] -> rec[220];
measure q[212] -> rec[221];
measure q[214] -> rec[222];
measure q[216] -> rec[223];
measure q[218] -> rec[224];
measure q[220] -> rec[225];
measure q[222] -> rec[226];
measure q[224] -> rec[227];
measure q[226] -> rec[228];
measure q[228] -> rec[229];
measure q[231] -> rec[230];
measure q[233] -> rec[231];
measure q[235] -> rec[232];
measure q[237] -> rec[233];
measure q[239] -> rec[234];
measure q[241] -> rec[235];
measure q[243] -> rec[236];
measure q[245] -> rec[237];
measure q[247] -> rec[238];
measure q[249] -> rec[239];
measure q[251] -> rec[240];

OPENQASM 2.0;
include "qelib1.inc";

qreg q[376];
creg rec[337];

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
reset q[23];
reset q[25];
reset q[28];
reset q[30];
reset q[32];
reset q[34];
reset q[36];
reset q[38];
reset q[40];
reset q[42];
reset q[44];
reset q[46];
reset q[48];
reset q[50];
reset q[52];
reset q[55];
reset q[57];
reset q[59];
reset q[61];
reset q[63];
reset q[65];
reset q[67];
reset q[69];
reset q[71];
reset q[73];
reset q[75];
reset q[77];
reset q[79];
reset q[82];
reset q[84];
reset q[86];
reset q[88];
reset q[90];
reset q[92];
reset q[94];
reset q[96];
reset q[98];
reset q[100];
reset q[102];
reset q[104];
reset q[106];
reset q[109];
reset q[111];
reset q[113];
reset q[115];
reset q[117];
reset q[119];
reset q[121];
reset q[123];
reset q[125];
reset q[127];
reset q[129];
reset q[131];
reset q[133];
reset q[136];
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
reset q[160];
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
reset q[185];
reset q[187];
reset q[190];
reset q[192];
reset q[194];
reset q[196];
reset q[198];
reset q[200];
reset q[202];
reset q[204];
reset q[206];
reset q[208];
reset q[210];
reset q[212];
reset q[214];
reset q[217];
reset q[219];
reset q[221];
reset q[223];
reset q[225];
reset q[227];
reset q[229];
reset q[231];
reset q[233];
reset q[235];
reset q[237];
reset q[239];
reset q[241];
reset q[244];
reset q[246];
reset q[248];
reset q[250];
reset q[252];
reset q[254];
reset q[256];
reset q[258];
reset q[260];
reset q[262];
reset q[264];
reset q[266];
reset q[268];
reset q[271];
reset q[273];
reset q[275];
reset q[277];
reset q[279];
reset q[281];
reset q[283];
reset q[285];
reset q[287];
reset q[289];
reset q[291];
reset q[293];
reset q[295];
reset q[298];
reset q[300];
reset q[302];
reset q[304];
reset q[306];
reset q[308];
reset q[310];
reset q[312];
reset q[314];
reset q[316];
reset q[318];
reset q[320];
reset q[322];
reset q[325];
reset q[327];
reset q[329];
reset q[331];
reset q[333];
reset q[335];
reset q[337];
reset q[339];
reset q[341];
reset q[343];
reset q[345];
reset q[347];
reset q[349];
reset q[2];
reset q[6];
reset q[10];
reset q[14];
reset q[18];
reset q[22];
reset q[29];
reset q[31];
reset q[33];
reset q[35];
reset q[37];
reset q[39];
reset q[41];
reset q[43];
reset q[45];
reset q[47];
reset q[49];
reset q[51];
reset q[53];
reset q[54];
reset q[56];
reset q[58];
reset q[60];
reset q[62];
reset q[64];
reset q[66];
reset q[68];
reset q[70];
reset q[72];
reset q[74];
reset q[76];
reset q[78];
reset q[83];
reset q[85];
reset q[87];
reset q[89];
reset q[91];
reset q[93];
reset q[95];
reset q[97];
reset q[99];
reset q[101];
reset q[103];
reset q[105];
reset q[107];
reset q[108];
reset q[110];
reset q[112];
reset q[114];
reset q[116];
reset q[118];
reset q[120];
reset q[122];
reset q[124];
reset q[126];
reset q[128];
reset q[130];
reset q[132];
reset q[137];
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
reset q[161];
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
reset q[184];
reset q[186];
reset q[191];
reset q[193];
reset q[195];
reset q[197];
reset q[199];
reset q[201];
reset q[203];
reset q[205];
reset q[207];
reset q[209];
reset q[211];
reset q[213];
reset q[215];
reset q[216];
reset q[218];
reset q[220];
reset q[222];
reset q[224];
reset q[226];
reset q[228];
reset q[230];
reset q[232];
reset q[234];
reset q[236];
reset q[238];
reset q[240];
reset q[245];
reset q[247];
reset q[249];
reset q[251];
reset q[253];
reset q[255];
reset q[257];
reset q[259];
reset q[261];
reset q[263];
reset q[265];
reset q[267];
reset q[269];
reset q[270];
reset q[272];
reset q[274];
reset q[276];
reset q[278];
reset q[280];
reset q[282];
reset q[284];
reset q[286];
reset q[288];
reset q[290];
reset q[292];
reset q[294];
reset q[299];
reset q[301];
reset q[303];
reset q[305];
reset q[307];
reset q[309];
reset q[311];
reset q[313];
reset q[315];
reset q[317];
reset q[319];
reset q[321];
reset q[323];
reset q[324];
reset q[326];
reset q[328];
reset q[330];
reset q[332];
reset q[334];
reset q[336];
reset q[338];
reset q[340];
reset q[342];
reset q[344];
reset q[346];
reset q[348];
reset q[355];
reset q[359];
reset q[363];
reset q[367];
reset q[371];
reset q[375];
barrier q;

h q[2];
h q[6];
h q[10];
h q[14];
h q[18];
h q[22];
h q[31];
h q[35];
h q[39];
h q[43];
h q[47];
h q[51];
h q[56];
h q[60];
h q[64];
h q[68];
h q[72];
h q[76];
h q[85];
h q[89];
h q[93];
h q[97];
h q[101];
h q[105];
h q[110];
h q[114];
h q[118];
h q[122];
h q[126];
h q[130];
h q[139];
h q[143];
h q[147];
h q[151];
h q[155];
h q[159];
h q[164];
h q[168];
h q[172];
h q[176];
h q[180];
h q[184];
h q[193];
h q[197];
h q[201];
h q[205];
h q[209];
h q[213];
h q[218];
h q[222];
h q[226];
h q[230];
h q[234];
h q[238];
h q[247];
h q[251];
h q[255];
h q[259];
h q[263];
h q[267];
h q[272];
h q[276];
h q[280];
h q[284];
h q[288];
h q[292];
h q[301];
h q[305];
h q[309];
h q[313];
h q[317];
h q[321];
h q[326];
h q[330];
h q[334];
h q[338];
h q[342];
h q[346];
h q[355];
h q[359];
h q[363];
h q[367];
h q[371];
h q[375];
barrier q;

cx q[2], q[3];
cx q[56], q[57];
cx q[110], q[111];
cx q[164], q[165];
cx q[218], q[219];
cx q[272], q[273];
cx q[326], q[327];
cx q[31], q[32];
cx q[85], q[86];
cx q[139], q[140];
cx q[193], q[194];
cx q[247], q[248];
cx q[301], q[302];
cx q[6], q[7];
cx q[60], q[61];
cx q[114], q[115];
cx q[168], q[169];
cx q[222], q[223];
cx q[276], q[277];
cx q[330], q[331];
cx q[35], q[36];
cx q[89], q[90];
cx q[143], q[144];
cx q[197], q[198];
cx q[251], q[252];
cx q[305], q[306];
cx q[10], q[11];
cx q[64], q[65];
cx q[118], q[119];
cx q[172], q[173];
cx q[226], q[227];
cx q[280], q[281];
cx q[334], q[335];
cx q[39], q[40];
cx q[93], q[94];
cx q[147], q[148];
cx q[201], q[202];
cx q[255], q[256];
cx q[309], q[310];
cx q[14], q[15];
cx q[68], q[69];
cx q[122], q[123];
cx q[176], q[177];
cx q[230], q[231];
cx q[284], q[285];
cx q[338], q[339];
cx q[43], q[44];
cx q[97], q[98];
cx q[151], q[152];
cx q[205], q[206];
cx q[259], q[260];
cx q[313], q[314];
cx q[18], q[19];
cx q[72], q[73];
cx q[126], q[127];
cx q[180], q[181];
cx q[234], q[235];
cx q[288], q[289];
cx q[342], q[343];
cx q[47], q[48];
cx q[101], q[102];
cx q[155], q[156];
cx q[209], q[210];
cx q[263], q[264];
cx q[317], q[318];
cx q[22], q[23];
cx q[76], q[77];
cx q[130], q[131];
cx q[184], q[185];
cx q[238], q[239];
cx q[292], q[293];
cx q[346], q[347];
cx q[51], q[52];
cx q[105], q[106];
cx q[159], q[160];
cx q[213], q[214];
cx q[267], q[268];
cx q[321], q[322];
cx q[55], q[54];
cx q[109], q[108];
cx q[163], q[162];
cx q[217], q[216];
cx q[271], q[270];
cx q[325], q[324];
cx q[30], q[29];
cx q[84], q[83];
cx q[138], q[137];
cx q[192], q[191];
cx q[246], q[245];
cx q[300], q[299];
cx q[59], q[58];
cx q[113], q[112];
cx q[167], q[166];
cx q[221], q[220];
cx q[275], q[274];
cx q[329], q[328];
cx q[34], q[33];
cx q[88], q[87];
cx q[142], q[141];
cx q[196], q[195];
cx q[250], q[249];
cx q[304], q[303];
cx q[63], q[62];
cx q[117], q[116];
cx q[171], q[170];
cx q[225], q[224];
cx q[279], q[278];
cx q[333], q[332];
cx q[38], q[37];
cx q[92], q[91];
cx q[146], q[145];
cx q[200], q[199];
cx q[254], q[253];
cx q[308], q[307];
cx q[67], q[66];
cx q[121], q[120];
cx q[175], q[174];
cx q[229], q[228];
cx q[283], q[282];
cx q[337], q[336];
cx q[42], q[41];
cx q[96], q[95];
cx q[150], q[149];
cx q[204], q[203];
cx q[258], q[257];
cx q[312], q[311];
cx q[71], q[70];
cx q[125], q[124];
cx q[179], q[178];
cx q[233], q[232];
cx q[287], q[286];
cx q[341], q[340];
cx q[46], q[45];
cx q[100], q[99];
cx q[154], q[153];
cx q[208], q[207];
cx q[262], q[261];
cx q[316], q[315];
cx q[75], q[74];
cx q[129], q[128];
cx q[183], q[182];
cx q[237], q[236];
cx q[291], q[290];
cx q[345], q[344];
cx q[50], q[49];
cx q[104], q[103];
cx q[158], q[157];
cx q[212], q[211];
cx q[266], q[265];
cx q[320], q[319];
cx q[79], q[78];
cx q[133], q[132];
cx q[187], q[186];
cx q[241], q[240];
cx q[295], q[294];
cx q[349], q[348];
barrier q;

cx q[2], q[1];
cx q[56], q[55];
cx q[110], q[109];
cx q[164], q[163];
cx q[218], q[217];
cx q[272], q[271];
cx q[326], q[325];
cx q[31], q[30];
cx q[85], q[84];
cx q[139], q[138];
cx q[193], q[192];
cx q[247], q[246];
cx q[301], q[300];
cx q[6], q[5];
cx q[60], q[59];
cx q[114], q[113];
cx q[168], q[167];
cx q[222], q[221];
cx q[276], q[275];
cx q[330], q[329];
cx q[35], q[34];
cx q[89], q[88];
cx q[143], q[142];
cx q[197], q[196];
cx q[251], q[250];
cx q[305], q[304];
cx q[10], q[9];
cx q[64], q[63];
cx q[118], q[117];
cx q[172], q[171];
cx q[226], q[225];
cx q[280], q[279];
cx q[334], q[333];
cx q[39], q[38];
cx q[93], q[92];
cx q[147], q[146];
cx q[201], q[200];
cx q[255], q[254];
cx q[309], q[308];
cx q[14], q[13];
cx q[68], q[67];
cx q[122], q[121];
cx q[176], q[175];
cx q[230], q[229];
cx q[284], q[283];
cx q[338], q[337];
cx q[43], q[42];
cx q[97], q[96];
cx q[151], q[150];
cx q[205], q[204];
cx q[259], q[258];
cx q[313], q[312];
cx q[18], q[17];
cx q[72], q[71];
cx q[126], q[125];
cx q[180], q[179];
cx q[234], q[233];
cx q[288], q[287];
cx q[342], q[341];
cx q[47], q[46];
cx q[101], q[100];
cx q[155], q[154];
cx q[209], q[208];
cx q[263], q[262];
cx q[317], q[316];
cx q[22], q[21];
cx q[76], q[75];
cx q[130], q[129];
cx q[184], q[183];
cx q[238], q[237];
cx q[292], q[291];
cx q[346], q[345];
cx q[51], q[50];
cx q[105], q[104];
cx q[159], q[158];
cx q[213], q[212];
cx q[267], q[266];
cx q[321], q[320];
cx q[28], q[54];
cx q[82], q[108];
cx q[136], q[162];
cx q[190], q[216];
cx q[244], q[270];
cx q[298], q[324];
cx q[3], q[29];
cx q[57], q[83];
cx q[111], q[137];
cx q[165], q[191];
cx q[219], q[245];
cx q[273], q[299];
cx q[32], q[58];
cx q[86], q[112];
cx q[140], q[166];
cx q[194], q[220];
cx q[248], q[274];
cx q[302], q[328];
cx q[7], q[33];
cx q[61], q[87];
cx q[115], q[141];
cx q[169], q[195];
cx q[223], q[249];
cx q[277], q[303];
cx q[36], q[62];
cx q[90], q[116];
cx q[144], q[170];
cx q[198], q[224];
cx q[252], q[278];
cx q[306], q[332];
cx q[11], q[37];
cx q[65], q[91];
cx q[119], q[145];
cx q[173], q[199];
cx q[227], q[253];
cx q[281], q[307];
cx q[40], q[66];
cx q[94], q[120];
cx q[148], q[174];
cx q[202], q[228];
cx q[256], q[282];
cx q[310], q[336];
cx q[15], q[41];
cx q[69], q[95];
cx q[123], q[149];
cx q[177], q[203];
cx q[231], q[257];
cx q[285], q[311];
cx q[44], q[70];
cx q[98], q[124];
cx q[152], q[178];
cx q[206], q[232];
cx q[260], q[286];
cx q[314], q[340];
cx q[19], q[45];
cx q[73], q[99];
cx q[127], q[153];
cx q[181], q[207];
cx q[235], q[261];
cx q[289], q[315];
cx q[48], q[74];
cx q[102], q[128];
cx q[156], q[182];
cx q[210], q[236];
cx q[264], q[290];
cx q[318], q[344];
cx q[23], q[49];
cx q[77], q[103];
cx q[131], q[157];
cx q[185], q[211];
cx q[239], q[265];
cx q[293], q[319];
cx q[52], q[78];
cx q[106], q[132];
cx q[160], q[186];
cx q[214], q[240];
cx q[268], q[294];
cx q[322], q[348];
barrier q;

cx q[56], q[30];
cx q[110], q[84];
cx q[164], q[138];
cx q[218], q[192];
cx q[272], q[246];
cx q[326], q[300];
cx q[31], q[5];
cx q[85], q[59];
cx q[139], q[113];
cx q[193], q[167];
cx q[247], q[221];
cx q[301], q[275];
cx q[355], q[329];
cx q[60], q[34];
cx q[114], q[88];
cx q[168], q[142];
cx q[222], q[196];
cx q[276], q[250];
cx q[330], q[304];
cx q[35], q[9];
cx q[89], q[63];
cx q[143], q[117];
cx q[197], q[171];
cx q[251], q[225];
cx q[305], q[279];
cx q[359], q[333];
cx q[64], q[38];
cx q[118], q[92];
cx q[172], q[146];
cx q[226], q[200];
cx q[280], q[254];
cx q[334], q[308];
cx q[39], q[13];
cx q[93], q[67];
cx q[147], q[121];
cx q[201], q[175];
cx q[255], q[229];
cx q[309], q[283];
cx q[363], q[337];
cx q[68], q[42];
cx q[122], q[96];
cx q[176], q[150];
cx q[230], q[204];
cx q[284], q[258];
cx q[338], q[312];
cx q[43], q[17];
cx q[97], q[71];
cx q[151], q[125];
cx q[205], q[179];
cx q[259], q[233];
cx q[313], q[287];
cx q[367], q[341];
cx q[72], q[46];
cx q[126], q[100];
cx q[180], q[154];
cx q[234], q[208];
cx q[288], q[262];
cx q[342], q[316];
cx q[47], q[21];
cx q[101], q[75];
cx q[155], q[129];
cx q[209], q[183];
cx q[263], q[237];
cx q[317], q[291];
cx q[371], q[345];
cx q[76], q[50];
cx q[130], q[104];
cx q[184], q[158];
cx q[238], q[212];
cx q[292], q[266];
cx q[346], q[320];
cx q[51], q[25];
cx q[105], q[79];
cx q[159], q[133];
cx q[213], q[187];
cx q[267], q[241];
cx q[321], q[295];
cx q[375], q[349];
cx q[28], q[29];
cx q[82], q[83];
cx q[136], q[137];
cx q[190], q[191];
cx q[244], q[245];
cx q[298], q[299];
cx q[57], q[58];
cx q[111], q[112];
cx q[165], q[166];
cx q[219], q[220];
cx q[273], q[274];
cx q[327], q[328];
cx q[32], q[33];
cx q[86], q[87];
cx q[140], q[141];
cx q[194], q[195];
cx q[248], q[249];
cx q[302], q[303];
cx q[61], q[62];
cx q[115], q[116];
cx q[169], q[170];
cx q[223], q[224];
cx q[277], q[278];
cx q[331], q[332];
cx q[36], q[37];
cx q[90], q[91];
cx q[144], q[145];
cx q[198], q[199];
cx q[252], q[253];
cx q[306], q[307];
cx q[65], q[66];
cx q[119], q[120];
cx q[173], q[174];
cx q[227], q[228];
cx q[281], q[282];
cx q[335], q[336];
cx q[40], q[41];
cx q[94], q[95];
cx q[148], q[149];
cx q[202], q[203];
cx q[256], q[257];
cx q[310], q[311];
cx q[69], q[70];
cx q[123], q[124];
cx q[177], q[178];
cx q[231], q[232];
cx q[285], q[286];
cx q[339], q[340];
cx q[44], q[45];
cx q[98], q[99];
cx q[152], q[153];
cx q[206], q[207];
cx q[260], q[261];
cx q[314], q[315];
cx q[73], q[74];
cx q[127], q[128];
cx q[181], q[182];
cx q[235], q[236];
cx q[289], q[290];
cx q[343], q[344];
cx q[48], q[49];
cx q[102], q[103];
cx q[156], q[157];
cx q[210], q[211];
cx q[264], q[265];
cx q[318], q[319];
cx q[77], q[78];
cx q[131], q[132];
cx q[185], q[186];
cx q[239], q[240];
cx q[293], q[294];
cx q[347], q[348];
cx q[52], q[53];
cx q[106], q[107];
cx q[160], q[161];
cx q[214], q[215];
cx q[268], q[269];
cx q[322], q[323];
barrier q;

cx q[56], q[28];
cx q[110], q[82];
cx q[164], q[136];
cx q[218], q[190];
cx q[272], q[244];
cx q[326], q[298];
cx q[31], q[3];
cx q[85], q[57];
cx q[139], q[111];
cx q[193], q[165];
cx q[247], q[219];
cx q[301], q[273];
cx q[355], q[327];
cx q[60], q[32];
cx q[114], q[86];
cx q[168], q[140];
cx q[222], q[194];
cx q[276], q[248];
cx q[330], q[302];
cx q[35], q[7];
cx q[89], q[61];
cx q[143], q[115];
cx q[197], q[169];
cx q[251], q[223];
cx q[305], q[277];
cx q[359], q[331];
cx q[64], q[36];
cx q[118], q[90];
cx q[172], q[144];
cx q[226], q[198];
cx q[280], q[252];
cx q[334], q[306];
cx q[39], q[11];
cx q[93], q[65];
cx q[147], q[119];
cx q[201], q[173];
cx q[255], q[227];
cx q[309], q[281];
cx q[363], q[335];
cx q[68], q[40];
cx q[122], q[94];
cx q[176], q[148];
cx q[230], q[202];
cx q[284], q[256];
cx q[338], q[310];
cx q[43], q[15];
cx q[97], q[69];
cx q[151], q[123];
cx q[205], q[177];
cx q[259], q[231];
cx q[313], q[285];
cx q[367], q[339];
cx q[72], q[44];
cx q[126], q[98];
cx q[180], q[152];
cx q[234], q[206];
cx q[288], q[260];
cx q[342], q[314];
cx q[47], q[19];
cx q[101], q[73];
cx q[155], q[127];
cx q[209], q[181];
cx q[263], q[235];
cx q[317], q[289];
cx q[371], q[343];
cx q[76], q[48];
cx q[130], q[102];
cx q[184], q[156];
cx q[238], q[210];
cx q[292], q[264];
cx q[346], q[318];
cx q[51], q[23];
cx q[105], q[77];
cx q[159], q[131];
cx q[213], q[185];
cx q[267], q[239];
cx q[321], q[293];
cx q[375], q[347];
cx q[1], q[29];
cx q[55], q[83];
cx q[109], q[137];
cx q[163], q[191];
cx q[217], q[245];
cx q[271], q[299];
cx q[30], q[58];
cx q[84], q[112];
cx q[138], q[166];
cx q[192], q[220];
cx q[246], q[274];
cx q[300], q[328];
cx q[5], q[33];
cx q[59], q[87];
cx q[113], q[141];
cx q[167], q[195];
cx q[221], q[249];
cx q[275], q[303];
cx q[34], q[62];
cx q[88], q[116];
cx q[142], q[170];
cx q[196], q[224];
cx q[250], q[278];
cx q[304], q[332];
cx q[9], q[37];
cx q[63], q[91];
cx q[117], q[145];
cx q[171], q[199];
cx q[225], q[253];
cx q[279], q[307];
cx q[38], q[66];
cx q[92], q[120];
cx q[146], q[174];
cx q[200], q[228];
cx q[254], q[282];
cx q[308], q[336];
cx q[13], q[41];
cx q[67], q[95];
cx q[121], q[149];
cx q[175], q[203];
cx q[229], q[257];
cx q[283], q[311];
cx q[42], q[70];
cx q[96], q[124];
cx q[150], q[178];
cx q[204], q[232];
cx q[258], q[286];
cx q[312], q[340];
cx q[17], q[45];
cx q[71], q[99];
cx q[125], q[153];
cx q[179], q[207];
cx q[233], q[261];
cx q[287], q[315];
cx q[46], q[74];
cx q[100], q[128];
cx q[154], q[182];
cx q[208], q[236];
cx q[262], q[290];
cx q[316], q[344];
cx q[21], q[49];
cx q[75], q[103];
cx q[129], q[157];
cx q[183], q[211];
cx q[237], q[265];
cx q[291], q[319];
cx q[50], q[78];
cx q[104], q[132];
cx q[158], q[186];
cx q[212], q[240];
cx q[266], q[294];
cx q[320], q[348];
cx q[25], q[53];
cx q[79], q[107];
cx q[133], q[161];
cx q[187], q[215];
cx q[241], q[269];
cx q[295], q[323];
barrier q;

h q[2];
h q[6];
h q[10];
h q[14];
h q[18];
h q[22];
h q[31];
h q[35];
h q[39];
h q[43];
h q[47];
h q[51];
h q[56];
h q[60];
h q[64];
h q[68];
h q[72];
h q[76];
h q[85];
h q[89];
h q[93];
h q[97];
h q[101];
h q[105];
h q[110];
h q[114];
h q[118];
h q[122];
h q[126];
h q[130];
h q[139];
h q[143];
h q[147];
h q[151];
h q[155];
h q[159];
h q[164];
h q[168];
h q[172];
h q[176];
h q[180];
h q[184];
h q[193];
h q[197];
h q[201];
h q[205];
h q[209];
h q[213];
h q[218];
h q[222];
h q[226];
h q[230];
h q[234];
h q[238];
h q[247];
h q[251];
h q[255];
h q[259];
h q[263];
h q[267];
h q[272];
h q[276];
h q[280];
h q[284];
h q[288];
h q[292];
h q[301];
h q[305];
h q[309];
h q[313];
h q[317];
h q[321];
h q[326];
h q[330];
h q[334];
h q[338];
h q[342];
h q[346];
h q[355];
h q[359];
h q[363];
h q[367];
h q[371];
h q[375];
barrier q;

measure q[2] -> rec[0]; reset q[2]; // decomposed MR
measure q[6] -> rec[1]; reset q[6]; // decomposed MR
measure q[10] -> rec[2]; reset q[10]; // decomposed MR
measure q[14] -> rec[3]; reset q[14]; // decomposed MR
measure q[18] -> rec[4]; reset q[18]; // decomposed MR
measure q[22] -> rec[5]; reset q[22]; // decomposed MR
measure q[29] -> rec[6]; reset q[29]; // decomposed MR
measure q[31] -> rec[7]; reset q[31]; // decomposed MR
measure q[33] -> rec[8]; reset q[33]; // decomposed MR
measure q[35] -> rec[9]; reset q[35]; // decomposed MR
measure q[37] -> rec[10]; reset q[37]; // decomposed MR
measure q[39] -> rec[11]; reset q[39]; // decomposed MR
measure q[41] -> rec[12]; reset q[41]; // decomposed MR
measure q[43] -> rec[13]; reset q[43]; // decomposed MR
measure q[45] -> rec[14]; reset q[45]; // decomposed MR
measure q[47] -> rec[15]; reset q[47]; // decomposed MR
measure q[49] -> rec[16]; reset q[49]; // decomposed MR
measure q[51] -> rec[17]; reset q[51]; // decomposed MR
measure q[53] -> rec[18]; reset q[53]; // decomposed MR
measure q[54] -> rec[19]; reset q[54]; // decomposed MR
measure q[56] -> rec[20]; reset q[56]; // decomposed MR
measure q[58] -> rec[21]; reset q[58]; // decomposed MR
measure q[60] -> rec[22]; reset q[60]; // decomposed MR
measure q[62] -> rec[23]; reset q[62]; // decomposed MR
measure q[64] -> rec[24]; reset q[64]; // decomposed MR
measure q[66] -> rec[25]; reset q[66]; // decomposed MR
measure q[68] -> rec[26]; reset q[68]; // decomposed MR
measure q[70] -> rec[27]; reset q[70]; // decomposed MR
measure q[72] -> rec[28]; reset q[72]; // decomposed MR
measure q[74] -> rec[29]; reset q[74]; // decomposed MR
measure q[76] -> rec[30]; reset q[76]; // decomposed MR
measure q[78] -> rec[31]; reset q[78]; // decomposed MR
measure q[83] -> rec[32]; reset q[83]; // decomposed MR
measure q[85] -> rec[33]; reset q[85]; // decomposed MR
measure q[87] -> rec[34]; reset q[87]; // decomposed MR
measure q[89] -> rec[35]; reset q[89]; // decomposed MR
measure q[91] -> rec[36]; reset q[91]; // decomposed MR
measure q[93] -> rec[37]; reset q[93]; // decomposed MR
measure q[95] -> rec[38]; reset q[95]; // decomposed MR
measure q[97] -> rec[39]; reset q[97]; // decomposed MR
measure q[99] -> rec[40]; reset q[99]; // decomposed MR
measure q[101] -> rec[41]; reset q[101]; // decomposed MR
measure q[103] -> rec[42]; reset q[103]; // decomposed MR
measure q[105] -> rec[43]; reset q[105]; // decomposed MR
measure q[107] -> rec[44]; reset q[107]; // decomposed MR
measure q[108] -> rec[45]; reset q[108]; // decomposed MR
measure q[110] -> rec[46]; reset q[110]; // decomposed MR
measure q[112] -> rec[47]; reset q[112]; // decomposed MR
measure q[114] -> rec[48]; reset q[114]; // decomposed MR
measure q[116] -> rec[49]; reset q[116]; // decomposed MR
measure q[118] -> rec[50]; reset q[118]; // decomposed MR
measure q[120] -> rec[51]; reset q[120]; // decomposed MR
measure q[122] -> rec[52]; reset q[122]; // decomposed MR
measure q[124] -> rec[53]; reset q[124]; // decomposed MR
measure q[126] -> rec[54]; reset q[126]; // decomposed MR
measure q[128] -> rec[55]; reset q[128]; // decomposed MR
measure q[130] -> rec[56]; reset q[130]; // decomposed MR
measure q[132] -> rec[57]; reset q[132]; // decomposed MR
measure q[137] -> rec[58]; reset q[137]; // decomposed MR
measure q[139] -> rec[59]; reset q[139]; // decomposed MR
measure q[141] -> rec[60]; reset q[141]; // decomposed MR
measure q[143] -> rec[61]; reset q[143]; // decomposed MR
measure q[145] -> rec[62]; reset q[145]; // decomposed MR
measure q[147] -> rec[63]; reset q[147]; // decomposed MR
measure q[149] -> rec[64]; reset q[149]; // decomposed MR
measure q[151] -> rec[65]; reset q[151]; // decomposed MR
measure q[153] -> rec[66]; reset q[153]; // decomposed MR
measure q[155] -> rec[67]; reset q[155]; // decomposed MR
measure q[157] -> rec[68]; reset q[157]; // decomposed MR
measure q[159] -> rec[69]; reset q[159]; // decomposed MR
measure q[161] -> rec[70]; reset q[161]; // decomposed MR
measure q[162] -> rec[71]; reset q[162]; // decomposed MR
measure q[164] -> rec[72]; reset q[164]; // decomposed MR
measure q[166] -> rec[73]; reset q[166]; // decomposed MR
measure q[168] -> rec[74]; reset q[168]; // decomposed MR
measure q[170] -> rec[75]; reset q[170]; // decomposed MR
measure q[172] -> rec[76]; reset q[172]; // decomposed MR
measure q[174] -> rec[77]; reset q[174]; // decomposed MR
measure q[176] -> rec[78]; reset q[176]; // decomposed MR
measure q[178] -> rec[79]; reset q[178]; // decomposed MR
measure q[180] -> rec[80]; reset q[180]; // decomposed MR
measure q[182] -> rec[81]; reset q[182]; // decomposed MR
measure q[184] -> rec[82]; reset q[184]; // decomposed MR
measure q[186] -> rec[83]; reset q[186]; // decomposed MR
measure q[191] -> rec[84]; reset q[191]; // decomposed MR
measure q[193] -> rec[85]; reset q[193]; // decomposed MR
measure q[195] -> rec[86]; reset q[195]; // decomposed MR
measure q[197] -> rec[87]; reset q[197]; // decomposed MR
measure q[199] -> rec[88]; reset q[199]; // decomposed MR
measure q[201] -> rec[89]; reset q[201]; // decomposed MR
measure q[203] -> rec[90]; reset q[203]; // decomposed MR
measure q[205] -> rec[91]; reset q[205]; // decomposed MR
measure q[207] -> rec[92]; reset q[207]; // decomposed MR
measure q[209] -> rec[93]; reset q[209]; // decomposed MR
measure q[211] -> rec[94]; reset q[211]; // decomposed MR
measure q[213] -> rec[95]; reset q[213]; // decomposed MR
measure q[215] -> rec[96]; reset q[215]; // decomposed MR
measure q[216] -> rec[97]; reset q[216]; // decomposed MR
measure q[218] -> rec[98]; reset q[218]; // decomposed MR
measure q[220] -> rec[99]; reset q[220]; // decomposed MR
measure q[222] -> rec[100]; reset q[222]; // decomposed MR
measure q[224] -> rec[101]; reset q[224]; // decomposed MR
measure q[226] -> rec[102]; reset q[226]; // decomposed MR
measure q[228] -> rec[103]; reset q[228]; // decomposed MR
measure q[230] -> rec[104]; reset q[230]; // decomposed MR
measure q[232] -> rec[105]; reset q[232]; // decomposed MR
measure q[234] -> rec[106]; reset q[234]; // decomposed MR
measure q[236] -> rec[107]; reset q[236]; // decomposed MR
measure q[238] -> rec[108]; reset q[238]; // decomposed MR
measure q[240] -> rec[109]; reset q[240]; // decomposed MR
measure q[245] -> rec[110]; reset q[245]; // decomposed MR
measure q[247] -> rec[111]; reset q[247]; // decomposed MR
measure q[249] -> rec[112]; reset q[249]; // decomposed MR
measure q[251] -> rec[113]; reset q[251]; // decomposed MR
measure q[253] -> rec[114]; reset q[253]; // decomposed MR
measure q[255] -> rec[115]; reset q[255]; // decomposed MR
measure q[257] -> rec[116]; reset q[257]; // decomposed MR
measure q[259] -> rec[117]; reset q[259]; // decomposed MR
measure q[261] -> rec[118]; reset q[261]; // decomposed MR
measure q[263] -> rec[119]; reset q[263]; // decomposed MR
measure q[265] -> rec[120]; reset q[265]; // decomposed MR
measure q[267] -> rec[121]; reset q[267]; // decomposed MR
measure q[269] -> rec[122]; reset q[269]; // decomposed MR
measure q[270] -> rec[123]; reset q[270]; // decomposed MR
measure q[272] -> rec[124]; reset q[272]; // decomposed MR
measure q[274] -> rec[125]; reset q[274]; // decomposed MR
measure q[276] -> rec[126]; reset q[276]; // decomposed MR
measure q[278] -> rec[127]; reset q[278]; // decomposed MR
measure q[280] -> rec[128]; reset q[280]; // decomposed MR
measure q[282] -> rec[129]; reset q[282]; // decomposed MR
measure q[284] -> rec[130]; reset q[284]; // decomposed MR
measure q[286] -> rec[131]; reset q[286]; // decomposed MR
measure q[288] -> rec[132]; reset q[288]; // decomposed MR
measure q[290] -> rec[133]; reset q[290]; // decomposed MR
measure q[292] -> rec[134]; reset q[292]; // decomposed MR
measure q[294] -> rec[135]; reset q[294]; // decomposed MR
measure q[299] -> rec[136]; reset q[299]; // decomposed MR
measure q[301] -> rec[137]; reset q[301]; // decomposed MR
measure q[303] -> rec[138]; reset q[303]; // decomposed MR
measure q[305] -> rec[139]; reset q[305]; // decomposed MR
measure q[307] -> rec[140]; reset q[307]; // decomposed MR
measure q[309] -> rec[141]; reset q[309]; // decomposed MR
measure q[311] -> rec[142]; reset q[311]; // decomposed MR
measure q[313] -> rec[143]; reset q[313]; // decomposed MR
measure q[315] -> rec[144]; reset q[315]; // decomposed MR
measure q[317] -> rec[145]; reset q[317]; // decomposed MR
measure q[319] -> rec[146]; reset q[319]; // decomposed MR
measure q[321] -> rec[147]; reset q[321]; // decomposed MR
measure q[323] -> rec[148]; reset q[323]; // decomposed MR
measure q[324] -> rec[149]; reset q[324]; // decomposed MR
measure q[326] -> rec[150]; reset q[326]; // decomposed MR
measure q[328] -> rec[151]; reset q[328]; // decomposed MR
measure q[330] -> rec[152]; reset q[330]; // decomposed MR
measure q[332] -> rec[153]; reset q[332]; // decomposed MR
measure q[334] -> rec[154]; reset q[334]; // decomposed MR
measure q[336] -> rec[155]; reset q[336]; // decomposed MR
measure q[338] -> rec[156]; reset q[338]; // decomposed MR
measure q[340] -> rec[157]; reset q[340]; // decomposed MR
measure q[342] -> rec[158]; reset q[342]; // decomposed MR
measure q[344] -> rec[159]; reset q[344]; // decomposed MR
measure q[346] -> rec[160]; reset q[346]; // decomposed MR
measure q[348] -> rec[161]; reset q[348]; // decomposed MR
measure q[355] -> rec[162]; reset q[355]; // decomposed MR
measure q[359] -> rec[163]; reset q[359]; // decomposed MR
measure q[363] -> rec[164]; reset q[363]; // decomposed MR
measure q[367] -> rec[165]; reset q[367]; // decomposed MR
measure q[371] -> rec[166]; reset q[371]; // decomposed MR
measure q[375] -> rec[167]; reset q[375]; // decomposed MR
measure q[1] -> rec[168];
measure q[3] -> rec[169];
measure q[5] -> rec[170];
measure q[7] -> rec[171];
measure q[9] -> rec[172];
measure q[11] -> rec[173];
measure q[13] -> rec[174];
measure q[15] -> rec[175];
measure q[17] -> rec[176];
measure q[19] -> rec[177];
measure q[21] -> rec[178];
measure q[23] -> rec[179];
measure q[25] -> rec[180];
measure q[28] -> rec[181];
measure q[30] -> rec[182];
measure q[32] -> rec[183];
measure q[34] -> rec[184];
measure q[36] -> rec[185];
measure q[38] -> rec[186];
measure q[40] -> rec[187];
measure q[42] -> rec[188];
measure q[44] -> rec[189];
measure q[46] -> rec[190];
measure q[48] -> rec[191];
measure q[50] -> rec[192];
measure q[52] -> rec[193];
measure q[55] -> rec[194];
measure q[57] -> rec[195];
measure q[59] -> rec[196];
measure q[61] -> rec[197];
measure q[63] -> rec[198];
measure q[65] -> rec[199];
measure q[67] -> rec[200];
measure q[69] -> rec[201];
measure q[71] -> rec[202];
measure q[73] -> rec[203];
measure q[75] -> rec[204];
measure q[77] -> rec[205];
measure q[79] -> rec[206];
measure q[82] -> rec[207];
measure q[84] -> rec[208];
measure q[86] -> rec[209];
measure q[88] -> rec[210];
measure q[90] -> rec[211];
measure q[92] -> rec[212];
measure q[94] -> rec[213];
measure q[96] -> rec[214];
measure q[98] -> rec[215];
measure q[100] -> rec[216];
measure q[102] -> rec[217];
measure q[104] -> rec[218];
measure q[106] -> rec[219];
measure q[109] -> rec[220];
measure q[111] -> rec[221];
measure q[113] -> rec[222];
measure q[115] -> rec[223];
measure q[117] -> rec[224];
measure q[119] -> rec[225];
measure q[121] -> rec[226];
measure q[123] -> rec[227];
measure q[125] -> rec[228];
measure q[127] -> rec[229];
measure q[129] -> rec[230];
measure q[131] -> rec[231];
measure q[133] -> rec[232];
measure q[136] -> rec[233];
measure q[138] -> rec[234];
measure q[140] -> rec[235];
measure q[142] -> rec[236];
measure q[144] -> rec[237];
measure q[146] -> rec[238];
measure q[148] -> rec[239];
measure q[150] -> rec[240];
measure q[152] -> rec[241];
measure q[154] -> rec[242];
measure q[156] -> rec[243];
measure q[158] -> rec[244];
measure q[160] -> rec[245];
measure q[163] -> rec[246];
measure q[165] -> rec[247];
measure q[167] -> rec[248];
measure q[169] -> rec[249];
measure q[171] -> rec[250];
measure q[173] -> rec[251];
measure q[175] -> rec[252];
measure q[177] -> rec[253];
measure q[179] -> rec[254];
measure q[181] -> rec[255];
measure q[183] -> rec[256];
measure q[185] -> rec[257];
measure q[187] -> rec[258];
measure q[190] -> rec[259];
measure q[192] -> rec[260];
measure q[194] -> rec[261];
measure q[196] -> rec[262];
measure q[198] -> rec[263];
measure q[200] -> rec[264];
measure q[202] -> rec[265];
measure q[204] -> rec[266];
measure q[206] -> rec[267];
measure q[208] -> rec[268];
measure q[210] -> rec[269];
measure q[212] -> rec[270];
measure q[214] -> rec[271];
measure q[217] -> rec[272];
measure q[219] -> rec[273];
measure q[221] -> rec[274];
measure q[223] -> rec[275];
measure q[225] -> rec[276];
measure q[227] -> rec[277];
measure q[229] -> rec[278];
measure q[231] -> rec[279];
measure q[233] -> rec[280];
measure q[235] -> rec[281];
measure q[237] -> rec[282];
measure q[239] -> rec[283];
measure q[241] -> rec[284];
measure q[244] -> rec[285];
measure q[246] -> rec[286];
measure q[248] -> rec[287];
measure q[250] -> rec[288];
measure q[252] -> rec[289];
measure q[254] -> rec[290];
measure q[256] -> rec[291];
measure q[258] -> rec[292];
measure q[260] -> rec[293];
measure q[262] -> rec[294];
measure q[264] -> rec[295];
measure q[266] -> rec[296];
measure q[268] -> rec[297];
measure q[271] -> rec[298];
measure q[273] -> rec[299];
measure q[275] -> rec[300];
measure q[277] -> rec[301];
measure q[279] -> rec[302];
measure q[281] -> rec[303];
measure q[283] -> rec[304];
measure q[285] -> rec[305];
measure q[287] -> rec[306];
measure q[289] -> rec[307];
measure q[291] -> rec[308];
measure q[293] -> rec[309];
measure q[295] -> rec[310];
measure q[298] -> rec[311];
measure q[300] -> rec[312];
measure q[302] -> rec[313];
measure q[304] -> rec[314];
measure q[306] -> rec[315];
measure q[308] -> rec[316];
measure q[310] -> rec[317];
measure q[312] -> rec[318];
measure q[314] -> rec[319];
measure q[316] -> rec[320];
measure q[318] -> rec[321];
measure q[320] -> rec[322];
measure q[322] -> rec[323];
measure q[325] -> rec[324];
measure q[327] -> rec[325];
measure q[329] -> rec[326];
measure q[331] -> rec[327];
measure q[333] -> rec[328];
measure q[335] -> rec[329];
measure q[337] -> rec[330];
measure q[339] -> rec[331];
measure q[341] -> rec[332];
measure q[343] -> rec[333];
measure q[345] -> rec[334];
measure q[347] -> rec[335];
measure q[349] -> rec[336];

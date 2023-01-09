OPENQASM 2.0;
include "qelib1.inc";
gate multiplex1_reverse_reverse_dg q0 { ry(-pi/2) q0; }
gate multiplex1_reverse_dg q0 { ry(pi/2) q0; }
gate multiplex2_reverse_dg q0,q1 { multiplex1_reverse_dg q0; cx q1,q0; multiplex1_reverse_reverse_dg q0; cx q1,q0; }
gate multiplex1_reverse_dg_v2 q0 { ry(pi) q0; }
gate disentangler_dg q0,q1,q2,q3,q4,q5 { multiplex1_reverse_dg_v2 q5; multiplex2_reverse_dg q4,q5; }
gate state_preparation(param0,param1,param2,param3,param4,param5,param6,param7,param8,param9,param10,param11,param12,param13,param14,param15,param16,param17,param18,param19,param20,param21,param22,param23,param24,param25,param26,param27,param28,param29,param30,param31,param32,param33,param34,param35,param36,param37,param38,param39,param40,param41,param42,param43,param44,param45,param46,param47,param48,param49,param50,param51,param52,param53,param54,param55,param56,param57,param58,param59,param60,param61,param62,param63) q0,q1,q2,q3,q4,q5 { disentangler_dg q0,q1,q2,q3,q4,q5; }
gate initialize(param0,param1,param2,param3,param4,param5,param6,param7,param8,param9,param10,param11,param12,param13,param14,param15,param16,param17,param18,param19,param20,param21,param22,param23,param24,param25,param26,param27,param28,param29,param30,param31,param32,param33,param34,param35,param36,param37,param38,param39,param40,param41,param42,param43,param44,param45,param46,param47,param48,param49,param50,param51,param52,param53,param54,param55,param56,param57,param58,param59,param60,param61,param62,param63) q0,q1,q2,q3,q4,q5 { reset q0; reset q1; reset q2; reset q3; reset q4; reset q5; state_preparation(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) q0,q1,q2,q3,q4,q5; }
qreg q517[11];
creg c2[5];
ry(1.4744939) q517[0];
ry(1.5288237) q517[1];
ry(1.5348263) q517[2];
ry(1.5349545) q517[3];
ry(1.5280419) q517[4];
cx q517[4],q517[3];
ry(0.26362996) q517[3];
cx q517[4],q517[3];
cx q517[3],q517[2];
ry(0.09089708) q517[2];
cx q517[4],q517[2];
ry(-0.028078987) q517[2];
cx q517[3],q517[2];
ry(0.154098) q517[2];
cx q517[4],q517[2];
cx q517[2],q517[1];
ry(0.053909553) q517[1];
cx q517[3],q517[1];
ry(-0.037992821) q517[1];
cx q517[2],q517[1];
ry(0.069829301) q517[1];
cx q517[4],q517[1];
ry(-0.037989791) q517[1];
cx q517[2],q517[1];
ry(0.037977681) q517[1];
cx q517[3],q517[1];
ry(-0.03799181) q517[1];
cx q517[2],q517[1];
ry(0.10169709) q517[1];
cx q517[4],q517[1];
cx q517[1],q517[0];
ry(0.098299744) q517[0];
cx q517[2],q517[0];
ry(-0.094305097) q517[0];
cx q517[1],q517[0];
ry(0.10229433) q517[0];
cx q517[3],q517[0];
ry(-0.094305289) q517[0];
cx q517[1],q517[0];
ry(0.094305544) q517[0];
cx q517[2],q517[0];
ry(-0.094305161) q517[0];
cx q517[1],q517[0];
ry(0.1102826) q517[0];
cx q517[4],q517[0];
ry(-0.094306053) q517[0];
cx q517[1],q517[0];
ry(0.094307073) q517[0];
cx q517[2],q517[0];
ry(-0.094305033) q517[0];
cx q517[1],q517[0];
ry(0.094309113) q517[0];
cx q517[3],q517[0];
ry(-0.094305543) q517[0];
cx q517[1],q517[0];
ry(0.094306054) q517[0];
cx q517[2],q517[0];
ry(-0.094305288) q517[0];
cx q517[1],q517[0];
ry(0.12625201) q517[0];
cx q517[4],q517[0];
initialize(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) q517[5],q517[6],q517[7],q517[8],q517[9],q517[10];
OPENQASM 2.0;
include "qelib1.inc";
gate multiplex1_reverse_reverse_dg q0 
{ 
    ry(-pi/2) q0; 
}
gate multiplex1_reverse_dg q0 
{ 
    ry(pi/2) q0; 
}
gate multiplex2_reverse_dg q0,q1 
{ 
    multiplex1_reverse_dg q0; 
    cx q1,q0; 
    multiplex1_reverse_reverse_dg q0; 
    cx q1,q0; 
}
gate multiplex1_reverse_dg_v2 q0 
{ 
    ry(pi) q0; 
}
gate disentangler_dg q0,q1,q2,q3,q4,q5 
{ 
    multiplex1_reverse_dg_v2 q5; 
    multiplex2_reverse_dg q4,q5; 
}
gate state_preparation(param0,param1,param2,param3,param4,param5,param6,param7,param8,param9,param10,param11,param12,param13,param14,param15,param16,param17,param18,param19,param20,param21,param22,param23,param24,param25,param26,param27,param28,param29,param30,param31,param32,param33,param34,param35,param36,param37,param38,param39,param40,param41,param42,param43,param44,param45,param46,param47,param48,param49,param50,param51,param52,param53,param54,param55,param56,param57,param58,param59,param60,param61,param62,param63) q0,q1,q2,q3,q4,q5 
{ 
    disentangler_dg q0,q1,q2,q3,q4,q5; 
}

gate initialize(param0,param1,param2,param3,param4,param5,param6,param7,param8,param9,param10,param11,param12,param13,param14,param15,param16,param17,param18,param19,param20,param21,param22,param23,param24,param25,param26,param27,param28,param29,param30,param31,param32,param33,param34,param35,param36,param37,param38,param39,param40,param41,param42,param43,param44,param45,param46,param47,param48,param49,param50,param51,param52,param53,param54,param55,param56,param57,param58,param59,param60,param61,param62,param63) q0,q1,q2,q3,q4,q5 
{ 
    reset q0; 
    reset q1; 
    reset q2; 
    reset q3; 
    reset q4; 
    reset q5; 
    state_preparation(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) q0,q1,q2,q3,q4,q5; 
}
qreg q322[11];
creg c1[5];
ry(1.4744939) q322[0];
ry(1.5288237) q322[1];
ry(1.5348263) q322[2];
ry(1.5349545) q322[3];
ry(1.5280419) q322[4];
cx q322[4],q322[3];
ry(0.26362996) q322[3];
cx q322[4],q322[3];
cx q322[3],q322[2];

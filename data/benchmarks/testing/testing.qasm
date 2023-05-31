OPENQASM 2.0;
//include "qelib1.inc;
gate majority a,b,c 
{ 
  cx c,b; 
  cx c,a; 
  ccx a,b,c; 
}
gate unmaj a,b,c 
{ 
  ccx a,b,c; 
  cx c,a; 
  cx a,b; 
}

// add a to b, storing result in b
gate add4 a0,a1,a2,a3,b0,b1,b2,b3,cin,cout 
{
  majority cin,b0,a0;
  cx a3,cout;
  unmaj a2,b3,a3;
}

qreg carry[2];
qreg a[8];
qreg b[8];
creg ans[8];
creg carryout[1];
// set input states
x a[0]; // a = 00000001
x b;
x b[6]; // b = 10111111
// output should be 11000000 0
majority a[0],b[1],carry[0];
add4 a[0],a[1],a[2],a[3],b[0],b[1],b[2],b[3],carry[0],carry[1];

rx(pi*-0.5) a[1];
rz(pi*0.5) b[1];
cx a[0],b[1];
u3(pi*1.0,0,pi*0.5) carry[0];
u3(pi*1.0,0,pi*1.5) carry[1];
measure a -> ans;
if(ans==0) add4 a[0],a[1],a[2],a[3],b[0],b[1],b[2],b[3],carry[0],carry[1];
// Gate: XX**1.1
u3(pi*0.5,pi*0.5,pi/4*3) a[1];
measure a[0] -> ans[0];
measure carry[0] -> carryout[0];
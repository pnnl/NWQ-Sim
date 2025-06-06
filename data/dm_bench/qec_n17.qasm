OPENQASM 2.0;
include "qelib1.inc";
qreg data[9];
qreg s_x[4];
qreg s_z[4];
creg meas_s_x[4];
creg meas_s_z[4];
rz(pi/2) s_x[0];
sx s_x[0];
rz(pi/2) s_x[0];
cx s_x[0],data[0];
cx s_x[0],data[1];
cx s_x[0],data[3];
cx s_x[0],data[4];
rz(pi/2) s_x[0];
sx s_x[0];
rz(pi/2) s_x[0];
rz(pi/2) s_x[1];
sx s_x[1];
rz(pi/2) s_x[1];
cx s_x[1],data[1];
cx s_x[1],data[2];
cx s_x[1],data[4];
cx s_x[1],data[5];
rz(pi/2) s_x[1];
sx s_x[1];
rz(pi/2) s_x[1];
rz(pi/2) s_x[2];
sx s_x[2];
rz(pi/2) s_x[2];
cx s_x[2],data[3];
cx s_x[2],data[4];
cx s_x[2],data[6];
cx s_x[2],data[7];
rz(pi/2) s_x[2];
sx s_x[2];
rz(pi/2) s_x[2];
rz(pi/2) s_x[3];
sx s_x[3];
rz(pi/2) s_x[3];
cx s_x[3],data[4];
cx s_x[3],data[5];
cx s_x[3],data[7];
cx s_x[3],data[8];
rz(pi/2) s_x[3];
sx s_x[3];
rz(pi/2) s_x[3];
cx data[0],s_z[0];
cx data[1],s_z[0];
cx data[3],s_z[0];
cx data[4],s_z[0];
cx data[1],s_z[1];
cx data[2],s_z[1];
cx data[4],s_z[1];
cx data[5],s_z[1];
cx data[3],s_z[2];
cx data[4],s_z[2];
cx data[6],s_z[2];
cx data[7],s_z[2];
cx data[4],s_z[3];
cx data[5],s_z[3];
cx data[7],s_z[3];
cx data[8],s_z[3];
measure s_x[0] -> meas_s_x[0];
measure s_x[1] -> meas_s_x[1];
measure s_x[2] -> meas_s_x[2];
measure s_x[3] -> meas_s_x[3];
measure s_z[0] -> meas_s_z[0];
measure s_z[1] -> meas_s_z[1];
measure s_z[2] -> meas_s_z[2];
measure s_z[3] -> meas_s_z[3];

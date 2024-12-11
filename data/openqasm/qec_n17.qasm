// Surface Code (9,1,3) - One Round of Syndrome Measurement
OPENQASM 2.0;
include "qelib1.inc";

// Define quantum registers
qreg data[9];    // Data qubits D0-D8
qreg s_x[4];      // X-type syndrome qubits Sx0-Sx3
qreg s_z[4];      // Z-type syndrome qubits Sz0-Sz3

// Define classical registers for measurements
creg meas_s_x[4]; // Measurement outcomes for X-type syndrome qubits
creg meas_s_z[4]; // Measurement outcomes for Z-type syndrome qubits

// Initialize X-type syndrome qubits to |+⟩ state
h s_x[0];
h s_x[1];
h s_x[2];
h s_x[3];

// **X-type Stabilizer Measurements**

// Sx0 stabilizer involving data qubits D0, D1, D3, D4
// Apply CNOT gates: syndrome qubit is control, data qubits are targets
cx s_x[0], data[0];
cx s_x[0], data[1];
cx s_x[0], data[3];
cx s_x[0], data[4];
h s_x[0];                // Rotate back to computational basis

// Sx1 stabilizer involving data qubits D1, D2, D4, D5
cx s_x[1], data[1];
cx s_x[1], data[2];
cx s_x[1], data[4];
cx s_x[1], data[5];
h s_x[1];

// Sx2 stabilizer involving data qubits D3, D4, D6, D7
cx s_x[2], data[3];
cx s_x[2], data[4];
cx s_x[2], data[6];
cx s_x[2], data[7];
h s_x[2];

// Sx3 stabilizer involving data qubits D4, D5, D7, D8
cx s_x[3], data[4];
cx s_x[3], data[5];
cx s_x[3], data[7];
cx s_x[3], data[8];
h s_x[3];

// **Z-type Stabilizer Measurements**

// Sz0 stabilizer involving data qubits D0, D1, D3, D4
// Apply CNOT gates: data qubits are controls, syndrome qubit is target
cx data[0], s_z[0];
cx data[1], s_z[0];
cx data[3], s_z[0];
cx data[4], s_z[0];

// Sz1 stabilizer involving data qubits D1, D2, D4, D5
cx data[1], s_z[1];
cx data[2], s_z[1];
cx data[4], s_z[1];
cx data[5], s_z[1];

// Sz2 stabilizer involving data qubits D3, D4, D6, D7
cx data[3], s_z[2];
cx data[4], s_z[2];
cx data[6], s_z[2];
cx data[7], s_z[2];

// Sz3 stabilizer involving data qubits D4, D5, D7, D8
cx data[4], s_z[3];
cx data[5], s_z[3];
cx data[7], s_z[3];
cx data[8], s_z[3];

measure s_x[0] -> meas_s_x[0];
measure s_x[1] -> meas_s_x[1];
measure s_x[2] -> meas_s_x[2];
measure s_x[3] -> meas_s_x[3];
measure s_z[0] -> meas_s_z[0];
measure s_z[1] -> meas_s_z[1];
measure s_z[2] -> meas_s_z[2];
measure s_z[3] -> meas_s_z[3];
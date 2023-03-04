// ---------------------------------------------------------------------------
// NWQsim: Northwest Quantum Circuit Simulation Environment
// ---------------------------------------------------------------------------
// Ang Li, Senior Computer Scientist
// Pacific Northwest National Laboratory(PNNL), U.S.
// Homepage: http://www.angliphd.com
// GitHub repo: http://www.github.com/pnnl/DM-Sim
// PNNL-IPID: 31919-E, ECCN: EAR99, IR: PNNL-SA-143160
// GitHub repo: http://www.github.com/pnnl/DM-Sim
// PNNL-IPID: 32166, ECCN: EAR99, IR: PNNL-SA-161181
// BSD Lincese.
// ---------------------------------------------------------------------------
// File: metric.h
// Define metric profiling functions for SV-Sim
// ---------------------------------------------------------------------------
#ifndef METRIC_H_
#define METRIC_H_

#include <cmath>

namespace NWQSim
{

    /**
    Calculates various metrics of a given quantum circuit and prints them to the console.
    @param circuit The quantum circuit represented as a vector of Gates.
    @param n_qubits The number of qubits in the circuit.
    @return None.
    */
    void circuit_metrics(std::vector<Gate> &circuit, int n_qubits)
    {
        // Initialize a vector to track the current depth of each qubit
        std::vector<int> qubit_depth(n_qubits, 0);
        // Initialize a vector to track the number of two-qubit gates applied to each qubit
        std::vector<int> qubit_g2_gates(n_qubits, 0);

        // Track the total number of single-qubit gates and two-qubit gates in the circuit
        IdxType g1_gates = 0;
        IdxType g2_gates = 0;

        // Loop over each gate in the circuit and update the depth of the corresponding qubits
        int max_depth = 0;
        for (int i = 0; i < circuit.size(); i++)
        {
            if (circuit[i].op_name == OP::MA)
                continue;
            int ctrl = circuit[i].ctrl;
            int target = circuit[i].qubit;

            // Calculate the depth of this gate based on the depths of the control and target qubits
            int depth;
            if (ctrl == -1)
            {
                // Single-qubit gate
                depth = qubit_depth[target] + 1;
                g1_gates++;
            }
            else
            {
                // Two-qubit gate
                depth = max(qubit_depth[ctrl], qubit_depth[target]) + 1;
                g2_gates++;

                // Increment the number of two-qubit gates applied to each qubit
                qubit_g2_gates[ctrl]++;
                qubit_g2_gates[target]++;

                if (ctrl == target)
                    printf("Wrong\n");
            }

            // Update the depth of the control and target qubits
            if (ctrl != -1)
                qubit_depth[ctrl] = depth;
            qubit_depth[target] = depth;

            // Update the maximum depth if the current depth is greater than the previous maximum
            if (depth > max_depth)
            {
                max_depth = depth;
            }
        }

        // Calculate the gate density, retention lifespan, and entanglement variance of the circuit
        ValType gate_density = (g1_gates + 2 * g2_gates) / (ValType)(max_depth * n_qubits);
        ValType retention_lifespan = log(max_depth);

        IdxType sum_g2_gates = 0;
        for (auto val : qubit_g2_gates)
        {
            sum_g2_gates += val;
        }
        ValType average_g2_gates = sum_g2_gates / (ValType)n_qubits;

        ValType entanglement_var = 0;
        for (auto val : qubit_g2_gates)
        {
            entanglement_var += log(pow(val - average_g2_gates, 2) + 1);
        }
        entanglement_var /= n_qubits;

        // Print the results to the console
        printf("Circuit Depth: %lld; Gate Density: %.4f; Retention Lifespan: %.4f; Entanglement Variance: %.4f\n", max_depth, gate_density, retention_lifespan, entanglement_var);
    }

} // end of NWQSim
#endif

#include "../include/NWQSim.hpp"
#include "../include/util.hpp"
#include <vector>
#define SVSIM

namespace NWQSim
{
    class svsim_cpu_circuit : public Circuit
    {
    public:
     
        // fused gate sequence
        std::vector<Gate> fused_circuit;
    

        void clear()
        {
            circuit.clear();
            n_gates = 0;
            // SAFE_FREE_HOST(circuit_cpu);
        }

    public:
        svsim_cpu_circuit(IdxType _n_qubits = 0) : n_qubits(_n_qubits), n_gates(0), circuit_cpu(NULL){};
        ~svsim_cpu_circuit(){clear();};


    };

    class svsim_cpu : public Simulation
    {
    };

} // namespace NWQSim
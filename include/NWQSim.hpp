// ---------------------------------------------------------------------------
// NWQsim: Northwest Quantum Circuit Simulation Environment
// ---------------------------------------------------------------------------
// Ang Li, Senior Computer Scientist
// Pacific Northwest National Laboratory(PNNL), U.S.
// Homepage: http://www.angliphd.com
// GitHub repo: http://www.github.com/pnnl/NWQ-Sim
// PNNL-IPID: 32166, ECCN: EAR99, IR: PNNL-SA-161181
// BSD Lincese.
// ---------------------------------------------------------------------------
// File:

// ---------------------------------------------------------------------------

#ifndef NWQSIM_HPP
#define NWQSIM_HPP

#include "util.hpp"
#include "circuit.hpp"
#include "gate.hpp"

#include <cmath>
#include <string>

namespace NWQSim
{
    class Simulation
    {
    public:
        Simulation(IdxType _n_qubits);
        ~Simulation();

        void AllocateQubit();
        void ReleaseQubit();

        void reset_sim();
        IdxType get_n_qubits();

        void set_seed(IdxType seed);

        void update(const IdxType _n_qubits, const IdxType _n_gates);

        void sim(Circuit circuit);

    private:
        IdxType measure(IdxType qubit);
        IdxType *measure_all(IdxType repetition = DEFAULT_REPETITIONS);
    };

} // namespace NWQSim

#endif // _NWQSIM_HPP_
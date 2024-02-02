// ---------------------------------------------------------------------------
// NWQSim: Northwest Quantum Simulation Environment
// ---------------------------------------------------------------------------
// Ang Li, Senior Computer Scientist
// Pacific Northwest National Laboratory(PNNL), U.S.
// Homepage: http://www.angliphd.com
// GitHub repo: http://www.github.com/pnnl/NWQ-Sim
// PNNL-IPID: 32166, ECCN: EAR99, IR: PNNL-SA-161181
// BSD Lincese.
// ---------------------------------------------------------------------------
// File: qshrap_wrapper.cpp
// NWQSim interface to Q# and QIR.
// ---------------------------------------------------------------------------

#include <memory>
#include <exception>
#include <iostream>
#include <stdexcept>

#include "../src/config.h"
#include "QirRuntimeApi_I.hpp"
#include "QSharpSimApi_I.hpp"
#include "QirContext.hpp"
#include "SimFactory.hpp"

#include "include/backendManager.hpp"
#include "include/state.hpp"
#include "include/circuit.hpp"
#include "include/nwq_util.hpp"

using namespace NWQSim;

// class NWQSimSimulator final : public Microsoft::Quantum::IRuntimeDriver,
// public Microsoft::Quantum::IQuantumGateSet

class NWQSimSimulator final : public Microsoft::Quantum::IRuntimeDriver
{
public:
    Result zero = reinterpret_cast<Result>(0xface0000);
    Result one = reinterpret_cast<Result>(0xface1000);
    char *qbase = reinterpret_cast<char *>(0xfce20000);

    //================= Basic ==================//
    IdxType to_qubit(QubitIdType Q)
    {
        CHECK_NULL_POINTER(Q);
        IdxType q = static_cast<IdxType>(reinterpret_cast<char *>(Q) - qbase);
        return q;
    }
    QubitIdType from_qubit(IdxType qubit)
    {
        return reinterpret_cast<QubitIdType>(qbase + qubit);
    }
    NWQSimSimulator() : circit(NULL), state(NULL), n_qubits(0)
    {
        std::string backend="CPU";
        std::string sim_method="sv";
        circuit = std::make_shared<Circuit>(_n_qubit);
        state = BackendManager::create_state(backend, _n_qubits, sim_method);
    }
    ~NWQSimSimulator()
    {
        n_qubits = 0;
    }
    std::string QubitToString(QubitIdType qubit) override
    {
        throw std::logic_error("QubitToString not_implemented");
    }

    //================= RuntimeDriver ==================//
    QubitIdType AllocateQubit() override
    {
        ++n_qubits;
        return from_qubit(n_qubits - 1);
    }
    void ReleaseQubit(QubitIdType Q) override
    {
        --n_qubits;
        if (n_qubits == 0)
            state->reset_state();
    }
    void ReleaseResult(Result result) override {}

    bool AreEqualResults(Result r1, Result r2) override
    {
        return (r1 == r2);
    }
    ResultValue GetResultValue(Result result) override
    {
        return (result == one) ? Result_One : Result_Zero;
    }
    Result UseZero() override { return zero; }
    Result UseOne() override { return one; }
    //================= Basis Gates ==================//
    void X(QubitIdType Qtarget)
    {
        Idxtype target = to_qubit(Qtarget);
#ifdef DUMP_QASM
        push_qasmstr(DUMP_QASM, "x", target);
#endif
        circuit->X(target);
    }
    void ID(QubitIdType Qtarget)
    {
        Idxtype target = to_qubit(Qtarget);
#ifdef DUMP_QASM
        push_qasmstr(DUMP_QASM, "id", target);
#endif
        circuit->ID(target);
    }
    void RZ(double theta, QubitIdType Qtarget)
    {
        Idxtype target = to_qubit(Qtarget);
#ifdef DUMP_QASM
        push_qasmstr(DUMP_QASM, "rz", target, -1, 1, theta);
#endif
        circuit->RZ(theta, target);
    }
    void SX(QubitIdType Qtarget)
    {
        Idxtype target = to_qubit(Qtarget);
#ifdef DUMP_QASM
        push_qasmstr(DUMP_QASM, "sx", target);
#endif
        circuit->SX(target);
    }
    void CX(QubitIdType Qcontrol, QubitIdType Qtarget)
    {
        Idxtype ctrl = to_qubit(Qcontrol);
        Idxtype target = to_qubit(Qtarget);
#ifdef DUMP_QASM
        push_qasmstr(DUMP_QASM, "cx", target, ctrl);
#endif
        circuit->CX(ctrl, target);
    }
    Result measure(QubitIdType Qtarget)
    {
        Idxtype target = to_qubit(Qtarget);
        circuit->M(target);
        circuit->set_num_qubits(n_qubits);
#ifdef DUMP_QASM
        push_qasmstr(DUMP_QASM, "measure", target, 0, 0, 0, n_qubits);
#endif
        state->sim(circuit);
        p_res = state->get_results();
        Idxtype res = p_res[0];
        return (res == (IdxType)1) ? UseOne() : UseZero();
    }

private:
    Circuit* circuit;
    QuantumState* state;
    IdxType n_qubits;
};

extern "C" Microsoft::Quantum::IRuntimeDriver *GetNWQSim()
{
    static Microsoft::Quantum::IRuntimeDriver *g_iqa = nullptr;
    if (!g_iqa)
    {
        g_iqa = new NWQSimSimulator{};
    }
    return g_iqa;
}

static NWQSimSimulator *getSim()
{
    return dynamic_cast<NWQSimSimulator *>(Microsoft::Quantum::GlobalContext()->GetDriver());
}

extern "C"
{
    void __quantum__qis__x__body(QubitIdType qubit) // NOLINT
    {
        // printf("Called X!\n");
        getSim()->X(qubit);
    }
    void __quantum__qis__rz__body(double theta, QubitIdType qubit)
    {
        // printf("Called RZ!\n");
        getSim()->RZ(theta, qubit);
    }
    void __quantum__qis__sqrtx__body(QubitIdType qubit)
    {
        // printf("Called SX!\n");
        getSim()->SX(qubit);
    }
    void __quantum__qis__cnot__body(QubitIdType ctrl, QubitIdType qubit)
    {
        // printf("Called CX!\n");
        getSim()->CX(ctrl, qubit);
    }
    Result __quantum__qis__m__body(QubitIdType qubit)
    {
        // printf("Called M!\n");
        return getSim()->measure(qubit);
    }
}

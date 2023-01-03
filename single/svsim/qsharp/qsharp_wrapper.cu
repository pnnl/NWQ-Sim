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
// File: svsim_qshrap_wrapper.cpp
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

#include "../src/util.h"
#include "../src/svsim_nvgpu_sin.cuh"

using namespace NWQSim;

//class SVSimSimulator final : public Microsoft::Quantum::IRuntimeDriver, 
//public Microsoft::Quantum::IQuantumGateSet

class SVSimSimulator final : public Microsoft::Quantum::IRuntimeDriver
{
public:
    Result zero = reinterpret_cast<Result>(0xface0000);
    Result one = reinterpret_cast<Result>(0xface1000);
    char* qbase = reinterpret_cast<char*>(0xfce20000);

    //================= Basic ==================//
    IdxType to_qubit(QubitIdType Q) 
    {
        CHECK_NULL_POINTER(Q);
        IdxType q = static_cast<IdxType>(reinterpret_cast<char*>(Q) - qbase);
        return q;
    }
    QubitIdType from_qubit(IdxType qubit) 
    {
        return reinterpret_cast<QubitIdType>(qbase + qubit);
    }
    SVSimSimulator():sim(NULL)
    {
        sim = new Simulation();
        n_qubits = 0;
        std::srand(time(0));
    }
    ~SVSimSimulator()
    {
        n_qubits = 0;
        delete sim;
    }
    std::string QubitToString(QubitIdType qubit) override
    {
        throw std::logic_error("QubitToString not_implemented");
    }
    
    //================= RuntimeDriver ==================//
    QubitIdType AllocateQubit() override 
    {
        ++n_qubits;
        if (n_qubits > N_QUBIT_SLOT)
            throw std::logic_error("Request qubits more than slots!");
        sim->AllocateQubit();
        return from_qubit(n_qubits-1);
    }
    void ReleaseQubit(QubitIdType Q) override
    {
        sim->ReleaseQubit();
        --n_qubits;
        if (n_qubits == 0) sim->reset_sim();
    }
    void ReleaseResult(Result result) override {} 

    bool AreEqualResults(Result r1, Result r2) override 
    {
        return (r1 == r2);
    }
    ResultValue GetResultValue(Result result) override {
        return (result == one) ? Result_One : Result_Zero;
    }
    Result UseZero() override { return zero; }
    Result UseOne() override { return one; }

    IdxType ControlledToMask(long numControls, QubitIdType controls[])
    {
        IdxType mask = 0;
        for (long i=0; i<numControls; i++)
        {
            mask = mask | ((IdxType)1<<to_qubit(controls[i]));
        }
        return mask;
    }

    //================= Basis Gates ==================//
    void X(QubitIdType Qtarget) 
    {
        sim->X(to_qubit(Qtarget));
    }
    void ID(QubitIdType Qtarget) 
    {
        sim->ID(to_qubit(Qtarget));
    }
    void RZ(double theta, QubitIdType Qtarget) 
    {
        sim->RZ(theta, to_qubit(Qtarget));
    }
    void SX(QubitIdType Qtarget) 
    {
        sim->SX(to_qubit(Qtarget));
    }
    void CX(QubitIdType Qcontrol, QubitIdType Qtarget) 
    {
        sim->CX(to_qubit(Qcontrol), to_qubit(Qtarget));
    }
    Result measure(QubitIdType Qtarget)
    { 
        IdxType res = sim->measure(to_qubit(Qtarget));
        return (res==(IdxType)1) ? UseOne() : UseZero();
    }
private:
    Simulation* sim;
    IdxType n_qubits;
};


extern "C" Microsoft::Quantum::IRuntimeDriver* GetSVSim() 
{
    static Microsoft::Quantum::IRuntimeDriver* g_iqa = nullptr;
    if(!g_iqa) 
    {
        g_iqa = new SVSimSimulator{};
    }
    return g_iqa;
}

static SVSimSimulator* getSim()
{
    return dynamic_cast<SVSimSimulator*>(Microsoft::Quantum::GlobalContext()->GetDriver());
}

extern "C"
{
    void __quantum__qis__x__body(QubitIdType qubit) // NOLINT
    {
        //printf("Called X!\n");
        getSim()->X(qubit);
    }
    void __quantum__qis__rz__body(double theta, QubitIdType qubit)
    {
        //printf("Called RZ!\n");
        getSim()->RZ(theta, qubit);
    }
    void __quantum__qis__sqrtx__body(QubitIdType qubit)
    {
        //printf("Called SX!\n");
        getSim()->SX(qubit);
    }
    void __quantum__qis__cnot__body(QubitIdType ctrl, QubitIdType qubit)
    {
        //printf("Called CX!\n");
        getSim()->CX(ctrl, qubit);
    }
    Result __quantum__qis__m__body(QubitIdType qubit)
    {
        //printf("Called M!\n");
        return getSim()->measure(qubit);
    }

}



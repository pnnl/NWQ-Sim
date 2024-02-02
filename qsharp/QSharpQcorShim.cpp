// Copyright (c) Microsoft Corporation.

// Shim for using Type1.Core gate set with a QSharp.Core simulator.

#include <stdexcept>
#include "QirTypes.hpp"
#include "QirContext.hpp"
#include "QirRuntime.hpp"

// #include "SimFactory.hpp"
extern "C" Microsoft::Quantum::IRuntimeDriver* GetSVSim(); 

extern "C"
{
    QIR_SHARED_API void __quantum__qis__x__ctl(QirArray* controls, QubitIdType qubit); // NOLINT
    void __quantum__qis__cnot__body(QubitIdType control, QubitIdType target) // NOLINT
    {
        QirArray* controls = __quantum__rt__array_create_1d(sizeof(QubitIdType), 1);
        *(QubitIdType*)__quantum__rt__array_get_element_ptr_1d(controls, 0) = control;
        __quantum__qis__x__ctl(controls, target);
        __quantum__rt__array_update_reference_count(controls, -1);
    }
    void __quantum__qis__cnot(QubitIdType control, QubitIdType target) // NOLINT
    {
        __quantum__qis__cnot__body(control, target);
    }
    void __quantum__qis__swap(QubitIdType q1, QubitIdType q2) // NOLINT
    {
        __quantum__qis__cnot__body(q1, q2);
        __quantum__qis__cnot__body(q2, q1);
        __quantum__qis__cnot__body(q1, q2);
    }

    QIR_SHARED_API void __quantum__qis__z__ctl(void* controls, QubitIdType qubit); // NOLINT
    void __quantum__qis__cz__body(QubitIdType control, QubitIdType target) // NOLINT
    {
        QirArray* controls = __quantum__rt__array_create_1d(sizeof(QubitIdType), 1);
        *(QubitIdType*)__quantum__rt__array_get_element_ptr_1d(controls, 0) = control;
        __quantum__qis__z__ctl(controls, target);
        __quantum__rt__array_update_reference_count(controls, -1);
    }
    void __quantum__qis__cz(QubitIdType control, QubitIdType target) // NOLINT
    {
        __quantum__qis__cz__body(control, target);
    }

    QIR_SHARED_API void __quantum__qis__r__body(PauliId pauli, double theta, QubitIdType qubit); // NOLINT
    void __quantum__qis__rx__body(double theta, QubitIdType qubit) // NOLINT
    {
        __quantum__qis__r__body(PauliId_X, theta, qubit);
    }
    void __quantum__qis__ry__body(double theta, QubitIdType qubit) // NOLINT
    {
        __quantum__qis__r__body(PauliId_Y, theta, qubit);
    }
    void __quantum__qis__rz__body(double theta, QubitIdType qubit) // NOLINT
    {
        __quantum__qis__r__body(PauliId_Z, theta, qubit);
    }
    void __quantum__qis__rz(double theta, QubitIdType qubit) // NOLINT
    {
        __quantum__qis__rz__body(theta, qubit);
    }

    void __quantum__qis__cphase(double theta, QubitIdType control, QubitIdType target) // NOLINT
    {
        double updatedTheta = theta / 2.0;
        __quantum__qis__rz__body(updatedTheta, target);
        __quantum__qis__rz__body(updatedTheta, control);
        __quantum__qis__cnot__body(control, target);
        __quantum__qis__rz__body(-updatedTheta, target);
        __quantum__qis__cnot__body(control, target);
    }

    QIR_SHARED_API Result __quantum__qis__measure__body(QirArray* paulis, QirArray* qubits); // NOLINT
    Result __quantum__qis__m__body(QubitIdType qubit) // NOLINT
    {
        QirArray* target = __quantum__rt__array_create_1d(sizeof(QubitIdType), 1);
        QirArray* paulis = __quantum__rt__array_create_1d(sizeof(PauliId), 1);
        *(QubitIdType*)__quantum__rt__array_get_element_ptr_1d(target, 0) = qubit;
        *__quantum__rt__array_get_element_ptr_1d(paulis, 0) = PauliId_Z;
        auto res = __quantum__qis__measure__body(paulis, target);
        __quantum__rt__array_update_reference_count(target, -1);
        __quantum__rt__array_update_reference_count(paulis, -1);
        return res;
    }
    Result __quantum__qis__mz(QubitIdType qubit) // NOLINT
    {
        auto res = __quantum__qis__m__body(qubit);
        bool* bit = new bool();
        *bit = __quantum__rt__result_equal(res, __quantum__rt__result_get_one());
        return reinterpret_cast<Result>(bit);
    }

    QIR_SHARED_API void __quantum__qis__x__body(QubitIdType); // NOLINT
    QIR_SHARED_API bool __quantum__rt__result_equal(Result, Result); // NOLINT
    QIR_SHARED_API Result __quantum__rt__result_get_one(); // NOLINT
    void __quantum__qis__reset__body(QubitIdType qubit) // NOLINT
    {
        auto result = __quantum__qis__m__body(qubit);
        if (__quantum__rt__result_equal(result, __quantum__rt__result_get_one()))
        {
            __quantum__qis__x__body(qubit);
        }
    }
    void __quantum__qis__reset(QubitIdType qubit) // NOLINT
    {
        __quantum__qis__reset__body(qubit);
    }

    QIR_SHARED_API void __quantum__qis__h__body(QubitIdType); // NOLINT
    void __quantum__qis__h(QubitIdType qubit) // NOLINT
    {
        __quantum__qis__h__body(qubit);
    }

    QIR_SHARED_API void __quantum__qis__y__body(QubitIdType); // NOLINT
    void __quantum__qis__y(QubitIdType qubit) // NOLINT
    {
        __quantum__qis__y__body(qubit);
    }

    QIR_SHARED_API void __quantum__qis__z__body(QubitIdType); // NOLINT
    void __quantum__qis__z(QubitIdType qubit) // NOLINT
    {
        __quantum__qis__z__body(qubit);
    }

    QIR_SHARED_API void __quantum__qis__t__body(QubitIdType); // NOLINT
    QIR_SHARED_API void __quantum__qis__t__adj(QubitIdType); // NOLINT
    void __quantum__qis__t(QubitIdType qubit) // NOLINT
    {
        __quantum__qis__t__body(qubit);
    }

    QIR_SHARED_API void __quantum__qis__x__body(QubitIdType); // NOLINT
    void __quantum__qis__x(QubitIdType qubit) // NOLINT
    {
        __quantum__qis__x__body(qubit);
    }

    void __quantum__qis__ccx(QubitIdType control1, QubitIdType control2, QubitIdType target) // NOLINT
    {
        __quantum__qis__h__body(target);
        __quantum__qis__t__adj(control1);
        __quantum__qis__t__adj(control2);
        __quantum__qis__cnot__body(target, control1);
        __quantum__qis__t__body(control1);
        __quantum__qis__cnot__body(control2, target);
        __quantum__qis__cnot__body(control2, control1);
        __quantum__qis__t__body(target);
        __quantum__qis__t__adj(control1);
        __quantum__qis__cnot__body(control2, target);
        __quantum__qis__cnot__body(target, control1);
        __quantum__qis__t__adj(target);
        __quantum__qis__t__body(control1);
        __quantum__qis__cnot__body(control2, control1);
        __quantum__qis__h__body(target);
    }

    void __quantum__rt__set_external_qreg(void*) {} // NOLINT
    void __quantum__rt__finalize() 
    {
    } // NOLINT
    void __quantum__rt__set_config_parameter(char*, char*) {} // NOLINT
    int32_t __quantum__rt__initialize(int32_t argc, char** argv) // NOLINT
    {
        Microsoft::Quantum::IRuntimeDriver* sim = GetSVSim();
        Microsoft::Quantum::InitializeQirContext(sim, false);
        return 0;
    }
}

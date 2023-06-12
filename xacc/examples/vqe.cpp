/*******************************************************************************
 * Copyright (c) 2019 UT-Battelle, LLC.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the Eclipse Public License v1.0
 * and Eclipse Distribution License v1.0 which accompanies this
 * distribution. The Eclipse Public License is available at
 * http://www.eclipse.org/legal/epl-v10.html and the Eclipse Distribution
 *License is available at https://eclipse.org/org/documents/edl-v10.php
 *
 * Contributors:
 *   Alexander J. McCaskey - initial API and implementation
 *******************************************************************************/
#include "xacc.hpp"
#include "Optimizer.hpp"
#include "xacc_observable.hpp"
#include "xacc_service.hpp"

#include "../nwq_accelerator.hpp"
#include <typeinfo>

int main(int argc, char **argv)
{
    xacc::Initialize(argc, argv);

    // Get reference to the Accelerator
    // specified by --accelerator argument
    // auto accelerator = xacc::getAccelerator("aer");

    std::shared_ptr<xacc::Accelerator> accelerator = std::make_shared<xacc::quantum::NWQAccelerator>();
    accelerator->updateConfiguration({std::make_pair("shots", 4096)});

    // Create the N=2 deuteron Hamiltonian
    auto H_N_2 = xacc::quantum::getObservable(
        "pauli", std::string("5.907 - 2.1433 X0X1 "
                             "- 2.1433 Y0Y1"
                             "+ .21829 Z0 - 6.125 Z1"));

    auto optimizer = xacc::getOptimizer("mlpack");

    // JIT map Quil QASM Ansatz to IR
    xacc::qasm(R"(
.compiler xasm
.circuit deuteron_ansatz
.parameters theta
.qbit q
X(q[0]);
Ry(q[1], theta);
CNOT(q[1],q[0]);
)");
    auto ansatz = xacc::getCompiled("deuteron_ansatz");

    // Get the VQE Algorithm and initialize it
    auto vqe = xacc::getAlgorithm("vqe");
    vqe->initialize({std::make_pair("ansatz", ansatz),
                     std::make_pair("observable", H_N_2),
                     std::make_pair("accelerator", accelerator),
                     std::make_pair("optimizer", optimizer)});

    // Allocate some qubits and execute
    auto buffer = xacc::qalloc(2);
    vqe->execute(buffer);

    std::cout << "XACC VQE (H_N_2 deuteron Hamiltonian) with " << accelerator->name() << " results:\n\n";

    std::cout << "opt_val: " << buffer->getInformation("opt-val").as<double>() << std::endl;

    auto keys = buffer->listExtraInfoKeys();

    for (auto k : keys)
    {
        std::cout << k << ": " << (*buffer)[k].toString() << "\n\n";
    }

    // buffer->print();
}
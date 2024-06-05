#include <iostream>
#include <mpi.h>

#include "xacc.hpp"
#include "nwq_accelerator.hpp"

std::shared_ptr<xacc::CompositeInstruction> generateBellCircuit(int n_qubits) {

    // Get a reference to the quantum IR service
    auto ir = xacc::getIRProvider("quantum");

    // Create a new quantum program
    auto program = ir->createComposite("bell_circuit");

    // Apply Hadamard gate to the first qubit
    auto h_gate = ir->createInstruction("H", {0});
    program->addInstruction(h_gate);

    // Apply CNOT gates across the chain of qubits
    for (unsigned long i = 1; i < n_qubits; ++i) {
        auto cnot = ir->createInstruction("CNOT", {0, i});
        program->addInstruction(cnot);
    }

    // Optionally, apply measurements to all qubits
    for (unsigned long i = 0; i < n_qubits; ++i) {
        auto measure = ir->createInstruction("Measure", {i});
        program->addInstruction(measure);
    }

    // Return the quantum circuit
    return program;
}

int main(int argc, char **argv)
{
    xacc::Initialize(argc, argv);
    // Accelerator:
    
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    auto nwq_acc = std::make_shared<xacc::quantum::NWQAccelerator>();

    int n_qubits = std::stoi(argv[1]);

    if (argc == 3)
        nwq_acc->updateConfiguration({std::make_pair("backend", std::string(argv[2])), std::make_pair("shots", 0)});

    auto bell = generateBellCircuit(n_qubits);

    // Allocate some qubits and execute
    auto buffer_nwq = xacc::qalloc(n_qubits);
    
    nwq_acc->execute(buffer_nwq, bell);
    
    if (rank == 0) {
        std::cout << "NWQ (EXP-Z): "
            << "\n";
        buffer_nwq->print();
    }
    
    xacc::Finalize();
    // MPI_Finalize();
    return 0;
}
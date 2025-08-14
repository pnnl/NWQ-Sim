#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include <memory>

#include "../include/backendManager.hpp"
#include "../include/state.hpp"
#include "../include/circuit.hpp"
// Adjust include path if needed:
#include "../include/stabsim/src/qasm_extraction.hpp"

static void usage(const char* prog) {
    std::cerr << "Usage: " << prog
              << " --backend cpu|nvgpu --qasm PATH --seed SEED --meas-count N\n";
}

int main(int argc, char** argv) {
    std::string backend;
    std::string qasm_path;
    long long seed = 0;
    int meas_count = -1;
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "--backend") && i + 1 < argc) backend = argv[++i];
        else if (!strcmp(argv[i], "--qasm") && i + 1 < argc) qasm_path = argv[++i];
        else if (!strcmp(argv[i], "--seed") && i + 1 < argc) seed = std::stoll(argv[++i]);
        else if (!strcmp(argv[i], "--meas-count") && i + 1 < argc) meas_count = std::stoi(argv[++i]);
        else { usage(argv[0]); return 2; }
    }
    if (backend.empty() || qasm_path.empty() || meas_count < 0) {
        usage(argv[0]);
        return 2;
    }

    // Load circuit from QASM
    std::shared_ptr<NWQSim::Circuit> circuit;
    try {
        // TODO: Replace with the actual helper from qasm_extraction.hpp.
        // Example placeholders (choose the one that matches your project):
        // circuit = NWQSim::load_qasm_circuit(qasm_path);
        // circuit = NWQSim::Circuit::from_qasm_file(qasm_path);
        // circuit = NWQSim::parse_qasm_to_circuit(qasm_path);
        circuit = appendQASMToCircuit(circuit, qasm_path, n_qubits);
    } catch (const std::exception& e) {
        std::cerr << "Failed to load QASM: " << e.what() << std::endl;
        return 3;
    }
    if (!circuit) {
        std::cerr << "Failed to construct circuit from QASM.\n";
        return 3;
    }

    // Create backend
    int n_qubits = circuit->get_qubits();
    auto state = BackendManager::create_state(backend.c_str(), n_qubits, "stab");
    if (!state) {
        std::cerr << "Failed to create backend '" << backend << "'.\n";
        return 4;
    }

    // Allocate measurement buffers (GPU needs this)
    if (meas_count > 0) {
        state->allocate_measurement_buffers(meas_count);
    }

    // Set seed for deterministic PRNG in your simulators
    state->set_seed(static_cast<NWQSim::IdxType>(seed));

    double sim_time_ms = 0;
    state->sim(circuit, sim_time_ms);

    auto results = state->get_measurement_results();

    // Print as a single line for the Python parser
    std::cout << "MEASURE:";
    for (int v : results) std::cout << " " << v;
    std::cout << std::endl;

    // Free GPU buffers if allocated
    state->free_measurement_buffers();
    return 0;
}
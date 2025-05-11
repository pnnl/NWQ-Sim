#include <memory>
#include <string>
#include <algorithm>
#include <cctype>
#include <iomanip>

/********** LOAD THE SIMULATION BACKEND ************/
#include "config.hpp"
#include "backendManager.hpp"
#include "state.hpp"
#include "nwq_util.hpp"
/**************************************************************************/

/********** UPDATE THE INCLUDE PATH FOR LOCAL HEADER FILE HERE ************/
#include "src/qasm_parser.hpp"
#include "src/parser_util.hpp"
/**************************************************************************/
/********** LOAD THE COMMAND LINE ARGUMENT PARSER ************/
#include "src/cl_parser.hpp"
/**************************************************************************/

#include <tamm/tamm.hpp>

using namespace NWQSim;
ValType pass_threshold = 0.98;
ValType run_brnchmark(std::string backend, IdxType index, IdxType total_shots, std::string simulation_method, bool is_basis);
void print_result(bool is_qobj, map<string, IdxType> *counts, std::shared_ptr<NWQSim::QuantumState> state, IdxType total_shots, std::string prefix = "");

int main(int argc, char **argv)
{
    // Create a Config object to parse command-line arguments
    ConfigParser config_parser;
    config_parser.parse_arguments(argc, argv);

    if (config_parser.is_flag_set("help") || argc == 1)
    {
        config_parser.print_help();
        return 0;
    }

    // Fetch and set configuration values
    IdxType total_shots = std::stoll(config_parser.get_value("shots"));
    std::string backend = config_parser.get_value("backend");
    std::transform(backend.begin(), backend.end(), backend.begin(),
                       [](unsigned char c)
                       { return std::toupper(c); });
        
    std::string simulation_method = config_parser.get_value("sim");
    bool run_with_basis = config_parser.is_flag_set("basis");

    std::string init_file = config_parser.get_value("init_file");
    std::string dump_file = config_parser.get_value("dump_file");
    bool report_fidelity = config_parser.is_flag_set("fidelity");

    bool print_metrics = config_parser.is_flag_set("metrics");

    // Handle init_format, fallback to simulation_method if empty
    std::string init_format = config_parser.get_value("init_format").empty() ? simulation_method : config_parser.get_value("init_format");

    // Enable specific features based on command-line options
    Config::ENABLE_TENSOR_CORE = config_parser.is_flag_set("hw_tensorcore");
    Config::ENABLE_MATRIX_CORE = config_parser.is_flag_set("hw_matrixcore");
    Config::ENABLE_AVX512 = config_parser.is_flag_set("hw_avx512");
    Config::OMP_NUM_THREADS = std::stoi(config_parser.get_value("hw_threads"));
    Config::ENABLE_FUSION = config_parser.is_flag_set("disable_fusion") ? false : true;
    Config::PRINT_SIM_TRACE = config_parser.is_flag_set("verbose");

    // Handle noise model related configurations
    Config::device_noise_file = config_parser.get_value("device");
    Config::ENABLE_NOISE = !Config::device_noise_file.empty();
    simulation_method = Config::ENABLE_NOISE ? "DM" : simulation_method;
    Config::device_layout_file = config_parser.get_value("layout");
    Config::device_layout_str = config_parser.get_value("layout_str");

    if (config_parser.is_flag_set("backend_list"))

#include <tamm/tamm.hpp>
    {
        BackendManager::print_available_backends();
        return 0;
    }

    if (config_parser.get_value("random_seed") != "")
    {
        Config::RANDOM_SEED = std::stoi(config_parser.get_value("random_seed"));
    }

// If MPI or NVSHMEM backend, initialize MPI
#ifdef MPI_ENABLED
    if (backend == "MPI" || backend == "NVGPU_MPI")
    {
        MPI_Init(&argc, &argv);
    }
    if (backend == "NVGPU_TAMM" || backend == "CPU_TAMM")
    {
        tamm::initialize(argc, argv);
    }
#endif

    if (Config::PRINT_SIM_TRACE)
    {
        config_parser.print_configurations();
    }

    qasm_parser parser;
    bool is_qobj = false;

    // Load the circuit from either a QASM file, QAQM string, JSON file, or JSON string
    if (config_parser.get_value("qasm_file") != "")
    {
        std::string qasmFileStr = config_parser.get_value("qasm_file");
        parser.load_qasm_file(qasmFileStr);
    }
    else if (config_parser.get_value("qasm_string") != "")
    {
        std::string qasmString = config_parser.get_value("qasm_string");
        parser.load_qasm_string(qasmString + ";\n");
    }
    else if (config_parser.get_value("json_file") != "")
    {
        std::string qobjFileStr = config_parser.get_value("json_file");
        parser.load_qobj_file(qobjFileStr);
        is_qobj = true;
    }
    else if (config_parser.get_value("json_string") != "")
    {
        std::string qobjString = config_parser.get_value("json_string");
        parser.load_qobj_string(qobjString);
        is_qobj = true;
    }
    else
    { // no input file, check if benchmark is set
        if (config_parser.get_value("test") != "")
        {
            total_shots = 16384; // for verification
            int benchmark_index = stoi(config_parser.get_value("test"));
            ValType fidelity = run_brnchmark(backend, benchmark_index, total_shots, simulation_method, run_with_basis);
            safe_print("Fidelity between NWQSim and Qiskit Execution: %.4f\n", fidelity);
        }
        else if (config_parser.is_flag_set("all_tests"))
        {
            bool passed = true;
            for (int benchmark_index = 12; benchmark_index < 36; benchmark_index++)
            {
                ValType fidelity = run_brnchmark(backend, benchmark_index, total_shots, simulation_method, run_with_basis);
                if (fidelity < pass_threshold)
                {
                    safe_print("Benchmark %d fidelity: %.4f Failed!\n", benchmark_index, fidelity);
                    passed = false;
                }
                else
                {
                    safe_print("Benchmark %d fidelity: %.4f Passed!\n", benchmark_index, fidelity);
                }
            }
            if (passed)
                safe_print("All benchmarks passed!\n");
            else
                safe_print("TESTING FAILED!\n");
        }

        return 0;
    }

    // Print metrics of the circuit, or execute the circuit by default

    if (print_metrics) // Print metrics of the circuit
    {
        parser.print_metrics();
    }
    else // Execute the circuit by default
    {
        if (!report_fidelity)
        {
            // Create the backend
            std::shared_ptr<NWQSim::QuantumState> state = BackendManager::create_state(backend, parser.num_qubits(), simulation_method);
            if (!state)
            {
                std::cerr << "Failed to create backend\n";
                return 1;
            }
            state->print_config(simulation_method);
            map<string, IdxType> *counts = parser.execute(state, init_file, init_format, total_shots);
            if (state->i_proc == 0)
            {
                print_result(is_qobj, counts, state, total_shots);
            }
            if (dump_file != "")
            {
                state->dump_res_state(dump_file);
            }
            delete counts;
        }
        else
        {
            safe_print("Starting Statevector Simulation...");
            std::string sv_backend = (backend == "NVGPU_MPI" ? "NVGPU" : backend);
            std::shared_ptr<NWQSim::QuantumState> sv_state = BackendManager::create_state(sv_backend, parser.num_qubits(), "sv");
            if (!sv_state)
            {
                std::cerr << "Failed to create backend\n";
                return 1;
            }
            if (sv_state->i_proc == 0)
            {
                sv_state->print_config(simulation_method);
                map<string, IdxType> *counts_sv = parser.execute(sv_state, init_file, init_format, total_shots);

                print_result(is_qobj, counts_sv, sv_state, total_shots, "sv_");
                delete counts_sv;
            }

            safe_print("Starting Density Matrix Simulation...");
            std::shared_ptr<NWQSim::QuantumState> dm_state = BackendManager::create_state(backend, parser.num_qubits(), "dm");
            if (!dm_state)
            {
                std::cerr << "Failed to create backend\n";
                return 1;
            }
            dm_state->print_config(simulation_method);
            map<string, IdxType> *counts_dm = parser.execute(dm_state, init_file, init_format, total_shots);
            if (dm_state->i_proc == 0)
            {
                print_result(is_qobj, counts_dm, dm_state, total_shots, "dm_");
            }
            if (dump_file != "")
            {
                dm_state->dump_res_state(dump_file);
            }
            ValType fidelity = dm_state->fidelity(sv_state);
            safe_print("State Fidelity: %e\n", fidelity);
            delete counts_dm;
        }
    }

// Finalize MPI if necessary
#ifdef MPI_ENABLED
    if (backend == "MPI" || backend == "NVGPU_MPI")
    {
        MPI_Finalize();
    }

    if (backend == "NVGPU_TAMM" || backend == "CPU_TAMM")
    {
        tamm::finalize();
    }
#endif
    return 0;
}

void print_result(bool is_qobj, map<string, IdxType> *counts, std::shared_ptr<NWQSim::QuantumState> state, IdxType total_shots, std::string prefix)
{
    if (is_qobj)
    {
        json result_count_json;
        for (const auto &r : (*counts))
        {
            result_count_json[(r.first)] = r.second;
        }
        cout << prefix << "nwq_sim_counts=" << result_count_json.dump() << endl;

        cout << "----------" << endl;
        state->print_res_state();
        cout << "----------" << endl;
    }
    else
    {
        print_counts(counts, total_shots);
    }
}

ValType run_brnchmark(std::string backend, IdxType index, IdxType total_shots, std::string simulation_method, bool is_basis)
{
    stringstream ss_file, ss_result;
    if (is_basis)
    {
        ss_file << "../../data/benchmarks_basis/circuits/" << index << ".qasm";
        ss_result << "../../data/benchmarks_basis/results/" << index << "_result.txt";
    }
    else
    {
        ss_file << "../../data/benchmarks/circuits/" << index << ".qasm";
        ss_result << "../../data/benchmarks/results/" << index << "_result.txt";
    }
    ifstream resultFile(ss_result.str().c_str());
    if (!resultFile)
    {
        printf("%s\n", "Could not open result file\n");
        return -1;
    }
    qasm_parser parser;
    parser.load_qasm_file(ss_file.str().c_str());
    // Create the backend
    std::shared_ptr<NWQSim::QuantumState> state = BackendManager::create_state(backend, parser.num_qubits(), simulation_method);
    if (!state)
    {
        std::cerr << "Failed to create backend\n";
        return 1;
    }
    state->print_config(simulation_method);
    map<string, IdxType> *svsim_counts = parser.execute(state, "", "", total_shots);
    map<string, ValType> ref_probs;
    map<string, ValType> svsim_probs;
    string line;
    ValType ref_norm = 0, svsim_norm = 0;

    while (getline(resultFile, line))
    {
        vector<string> ref_result = split(line, ' ');
        string outcome = ref_result[0];
        ValType ref_prob = stod(ref_result[1]);
        ref_probs.insert({outcome, ref_prob});
        ref_norm += ref_prob * ref_prob;
    }
    for (auto it = svsim_counts->begin(); it != svsim_counts->end(); ++it)
    {
        const auto &key = it->first;
        const auto &val = it->second;

        ValType out_prob = val / static_cast<ValType>(total_shots);
        svsim_norm += out_prob * out_prob;
        svsim_probs.insert(std::make_pair(key, out_prob));
    }
    ref_norm = sqrt(ref_norm);
    svsim_norm = sqrt(svsim_norm);
    ValType fidelity = 0;
    for (auto it = svsim_probs.begin(); it != svsim_probs.end(); ++it)
    {
        const auto &key = it->first;
        const auto &val = it->second;

        if (ref_probs.count(key))
            fidelity += (val / svsim_norm) * (ref_probs.at(key) / ref_norm);
    }
    return fidelity;
}

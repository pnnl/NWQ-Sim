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

using namespace NWQSim;
ValType pass_threshold = 0.98;
ValType run_brnchmark(std::string backend, IdxType index, IdxType total_shots, std::string simulation_method, bool is_basis);

int main(int argc, char **argv)
{
    // Create a Config object to parse command-line arguments
    ConfigParser config_parser;
    config_parser.parse_command_line_arguments(argc, argv);

    // Helper function to retrieve boolean values
    auto get_bool_value = [&](const std::string &key)
    {
        return config_parser.get_value(key) == "true";
    };

    // Fetch and set configuration values
    IdxType total_shots = std::stoll(config_parser.get_value("shots"));
    std::string backend = config_parser.get_value("backend");
    std::string simulation_method = config_parser.get_value("sim_method");
    bool run_with_basis = get_bool_value("basis");
    bool print_metrics = get_bool_value("metrics");
    std::string init_file = config_parser.get_value("init_file");
    std::string dumpfile = config_parser.get_value("dump_file");
    bool report_fidelity = get_bool_value("fidelity");

    // Handle init_format, fallback to simulation_method if empty
    std::string init_format = config_parser.get_value("init_format").empty() ? simulation_method : config_parser.get_value("init_format");

    // Enable specific features based on command-line options
    Config::ENABLE_TENSOR_CORE = get_bool_value("tensorcore");
    Config::ENABLE_MATRIX_CORE = get_bool_value("matrixcore");
    Config::ENABLE_AVX512 = get_bool_value("AVX512");
    Config::OMP_NUM_THREADS = std::stoi(config_parser.get_value("threads"));
    Config::ENABLE_FUSION = get_bool_value("disable_fusion") ? false : true;
    Config::PRINT_SIM_TRACE = get_bool_value("verbose");

    // Handle noise model related configurations
    Config::ENABLE_NOISE = !config_parser.get_value("noise_model").empty();
    if (Config::ENABLE_NOISE)
    {
        Config::readDeviceConfig(config_parser.get_value("noise_model"));
        simulation_method = "dm";
    }
    std::string layoutfile = config_parser.get_value("layout");
    if (!layoutfile.empty())
    {
        Config::loadLayoutFile(layoutfile);
    }
    std::string layoutstring = config_parser.get_value("layout_str");
    if (!layoutstring.empty())
    {
        Config::loadLayoutString(layoutstring);
    }

    if (config_parser.get_value("backend_list") == "true")
    {
        BackendManager::print_available_backends();
        return 0;
    }

    if (config_parser.get_value("h") == "true")
    {
        config_parser.print_help();
        return 0;
    }

    if (Config::PRINT_SIM_TRACE)
    {
        config_parser.print_configs();
    }

// If MPI or NVSHMEM backend, initialize MPI
#ifdef MPI_ENABLED
    if (backend == "MPI" || backend == "NVGPU_MPI")
    {
        MPI_Init(&argc, &argv);
    }
#endif
    if (config_parser.get_value("q") != "")
    {
        std::string qasmFileStr = config_parser.get_value("q");

        qasm_parser parser;

        parser.load_qasm_file(qasmFileStr.c_str());

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
            map<string, IdxType> *counts = parser.execute(state, init_file, init_format, total_shots, print_metrics);

            if (dumpfile != "")
            {
                state->dump_res_state(dumpfile);
            }
            if (state->i_proc == 0)
            {
                print_counts(counts, total_shots);
            }
            delete counts;
        }
        else
        {
            BackendManager::safe_print("Starting Statevector Simulation...");
            std::string sv_backend = (backend == "NVGPU_MPI" ? "NVGPU" : backend);
            std::shared_ptr<NWQSim::QuantumState> sv_state = BackendManager::create_state(sv_backend, parser.num_qubits(), "sv");
            std::shared_ptr<NWQSim::QuantumState> dm_state = BackendManager::create_state(backend, parser.num_qubits(), "dm");
            if (!sv_state)
            {
                std::cerr << "Failed to create backend\n";
                return 1;
            }
            if (dm_state->i_proc == 0)
            {
                sv_state->print_config(simulation_method);
                map<string, IdxType> *counts_sv = parser.execute(sv_state, init_file, init_format, total_shots, print_metrics);
                print_counts(counts_sv, total_shots);
                delete counts_sv;
            }

            BackendManager::safe_print("Starting Density Matrix Simulation...");
            if (!dm_state)
            {
                std::cerr << "Failed to create backend\n";
                return 1;
            }
            dm_state->print_config(simulation_method);
            map<string, IdxType> *counts_dm = parser.execute(dm_state, init_file, init_format, total_shots, print_metrics);
            if (dm_state->i_proc == 0)
            {
                print_counts(counts_dm, total_shots);
            }
            if (dumpfile != "")
            {
                dm_state->dump_res_state(dumpfile);
            }

            // if (dumpfile != "") {
            //     dm_state->dump_res_state(dumpfile);
            // }
            ValType fidelity = dm_state->fidelity(sv_state);
            BackendManager::safe_print("State Fidelity: %e\n", fidelity);
            delete counts_dm;
        }
    }

    if (config_parser.get_value("qs") != "")
    {
        std::string qasmString = config_parser.get_value("qs");

        qasm_parser parser;
        parser.load_qasm_string(qasmString + ";\n");
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
            map<string, IdxType> *counts = parser.execute(state, init_file, init_format, total_shots, print_metrics);
            if (state->i_proc == 0)
            {
                print_counts(counts, total_shots);
            }
            if (dumpfile != "")
            {
                state->dump_res_state(dumpfile);
            }
            delete counts;
        }
        else
        {
            BackendManager::safe_print("Starting Statevector Simulation...");
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
                map<string, IdxType> *counts_sv = parser.execute(sv_state, init_file, init_format, total_shots, print_metrics);
                print_counts(counts_sv, total_shots);
                delete counts_sv;
            }

            BackendManager::safe_print("Starting Density Matrix Simulation...");
            std::shared_ptr<NWQSim::QuantumState> dm_state = BackendManager::create_state(backend, parser.num_qubits(), "dm");
            if (!dm_state)
            {
                std::cerr << "Failed to create backend\n";
                return 1;
            }
            dm_state->print_config(simulation_method);
            map<string, IdxType> *counts_dm = parser.execute(dm_state, init_file, init_format, total_shots, print_metrics);
            if (dm_state->i_proc == 0)
            {
                print_counts(counts_dm, total_shots);
            }
            if (dumpfile != "")
            {
                dm_state->dump_res_state(dumpfile);
            }
            ValType fidelity = dm_state->fidelity(sv_state);
            BackendManager::safe_print("State Fidelity: %e\n", fidelity);
            delete counts_dm;
        }
    }
    if (config_parser.get_value("j") != "")
    {
        std::string qobjFileStr = config_parser.get_value("j");
        qasm_parser parser;
        parser.load_qobj_file(qobjFileStr.c_str());
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
            map<string, IdxType> *counts = parser.execute(state, init_file, init_format, total_shots, print_metrics);
            if (state->i_proc == 0)
            {
                json result_count_json;
                for (const auto &r : (*counts))
                {
                    result_count_json[(r.first)] = r.second;
                }
                cout << "nwq_sim_counts=" << result_count_json.dump() << endl;

                cout << "----------" << endl;
                state->print_res_state();
                cout << "----------" << endl;
            }
            if (dumpfile != "")
            {
                state->dump_res_state(dumpfile);
            }
            delete counts;
        }
        else
        {
            BackendManager::safe_print("Starting Statevector Simulation...");
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
                map<string, IdxType> *counts_sv = parser.execute(sv_state, init_file, init_format, total_shots, print_metrics);
                json result_count_json;
                for (const auto &r : (*counts_sv))
                {
                    result_count_json[(r.first)] = r.second;
                }
                cout << "sv_nwq_sim_counts=" << result_count_json.dump() << endl;

                cout << "----------" << endl;
                sv_state->print_res_state();
                cout << "----------" << endl;
                delete counts_sv;
            }

            BackendManager::safe_print("Starting Density Matrix Simulation...");
            std::shared_ptr<NWQSim::QuantumState> dm_state = BackendManager::create_state(backend, parser.num_qubits(), "dm");
            if (!dm_state)
            {
                std::cerr << "Failed to create backend\n";
                return 1;
            }
            dm_state->print_config(simulation_method);
            map<string, IdxType> *counts_dm = parser.execute(dm_state, init_file, init_format, total_shots, print_metrics);
            if (dm_state->i_proc == 0)
            {
                json result_count_json;
                for (const auto &r : (*counts_dm))
                {
                    result_count_json[(r.first)] = r.second;
                }
                cout << "dm_nwq_sim_counts=" << result_count_json.dump() << endl;

                cout << "----------" << endl;
                dm_state->print_res_state();
                cout << "----------" << endl;
            }
            ValType fidelity = dm_state->fidelity(sv_state);
            BackendManager::safe_print("State Fidelity: %e\n", fidelity);
            delete counts_dm;
        }
    }
    if (config_parser.get_value("js") != "")
    {
        std::string qobjString = config_parser.get_value("js");
        qasm_parser parser;
        parser.load_qobj_string(qobjString);
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
            map<string, IdxType> *counts = parser.execute(state, init_file, init_format, total_shots, print_metrics);
            if (state->i_proc == 0)
            {
                json result_count_json;
                for (const auto &r : (*counts))
                {
                    result_count_json[(r.first)] = r.second;
                }
                cout << "nwq_sim_counts=" << result_count_json.dump() << endl;

                cout << "----------" << endl;
                state->print_res_state();
                cout << "----------" << endl;
            }
            if (dumpfile != "")
            {
                state->dump_res_state(dumpfile);
            }
            delete counts;
        }
        else
        {
            BackendManager::safe_print("Starting Statevector Simulation...");
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
                map<string, IdxType> *counts_sv = parser.execute(sv_state, init_file, init_format, total_shots, print_metrics);
                json result_count_json;
                for (const auto &r : (*counts_sv))
                {
                    result_count_json[(r.first)] = r.second;
                }
                cout << "sv_nwq_sim_counts=" << result_count_json.dump() << endl;

                cout << "----------" << endl;
                sv_state->print_res_state();
                cout << "----------" << endl;
                delete counts_sv;
            }

            BackendManager::safe_print("Starting Density Matrix Simulation...");
            std::shared_ptr<NWQSim::QuantumState> dm_state = BackendManager::create_state(backend, parser.num_qubits(), "dm");
            if (!dm_state)
            {
                std::cerr << "Failed to create backend\n";
                return 1;
            }
            dm_state->print_config(simulation_method);
            map<string, IdxType> *counts_dm = parser.execute(dm_state, init_file, init_format, total_shots, print_metrics);
            if (dm_state->i_proc == 0)
            {
                json result_count_json;
                for (const auto &r : (*counts_dm))
                {
                    result_count_json[(r.first)] = r.second;
                }
                cout << "dm_nwq_sim_counts=" << result_count_json.dump() << endl;

                cout << "----------" << endl;
                dm_state->print_res_state();
                cout << "----------" << endl;
            }
            ValType fidelity = dm_state->fidelity(sv_state);
            BackendManager::safe_print("State Fidelity: %e\n", fidelity);
            delete counts_dm;
        }
    }

    if (config_parser.get_value("t") != "")
    {
        total_shots = 16384; // for verification
        int benchmark_index = stoi(config_parser.get_value("t"));
        ValType fidelity = run_brnchmark(backend, benchmark_index, total_shots, simulation_method, run_with_basis);
        BackendManager::safe_print("Fidelity between NWQSim and Qiskit Execution: %.4f\n", fidelity);
    }

    if (config_parser.get_value("a") == "true")
    {
        bool passed = true;
        for (int benchmark_index = 12; benchmark_index < 36; benchmark_index++)
        {
            ValType fidelity = run_brnchmark(backend, benchmark_index, total_shots, simulation_method, run_with_basis);
            if (fidelity < pass_threshold)
            {
                BackendManager::safe_print("Benchmark %d fidelity: %.4f Failed!\n", benchmark_index, fidelity);
                passed = false;
            }
        }
        if (passed)
            BackendManager::safe_print("All benchmarks passed!\n");
        else
            BackendManager::safe_print("TESTING FAILED!\n");
    }
// Finalize MPI if necessary
#ifdef MPI_ENABLED
    if (backend == "MPI" || backend == "NVGPU_MPI")
    {
        MPI_Finalize();
    }
#endif
    return 0;
}

ValType run_brnchmark(std::string backend, IdxType index, IdxType total_shots, std::string simulation_method, bool is_basis)
{
    stringstream ss_file, ss_result;
    if (is_basis)
    {
        ss_file << "../data/benchmarks_basis/circuits/" << index << ".qasm";
        ss_result << "../data/benchmarks_basis/results/" << index << "_result.txt";
    }
    else
    {
        ss_file << "../data/benchmarks/circuits/" << index << ".qasm";
        ss_result << "../data/benchmarks/results/" << index << "_result.txt";
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

#include <memory>
#include <string>
#include <algorithm>
#include <cctype>
#include <iomanip>

/********** UPDATE THE INCLUDE PATH FOR LOCAL HEADER FILE HERE ************/
#include "private/config.hpp"
#include "src/qasm_parser.hpp"
#include "src/parser_util.hpp"
/**************************************************************************/

/********** LOAD THE SIMULATION BACKEND ************/
#include "backendManager.hpp"
#include "state.hpp"
#include "nwq_util.hpp"
/**************************************************************************/

using namespace NWQSim;
ValType pass_threshold = 0.98;
ValType run_brnchmark(std::string backend, const std::string& config, IdxType index, IdxType total_shots, std::string simulation_method, bool is_basis);

int main(int argc, char **argv)
{
    IdxType total_shots = 1024;
    bool run_with_basis = false;
    bool print_metrics = false;
    std::string config = "../default_config.json";
    std::string initfile = "";
    std::string dumpfile = "";
    std::string layoutfile = "";
    bool report_fidelity = false;
    std::string backend = "CPU";
    std::string simulation_method = "sv";

    if (cmdOptionExists(argv, argv + argc, "-h") || argc == 1)
    {
        std::cout << "\033[1;33m"
                  << "Usage: ./nwq_qasm [options]"
                  << "\033[0m" << std::endl
                  << std::endl;
        std::cout << std::left << std::setw(20) << "\033[1;32mOption\033[0m"
                  << "           Description" << std::endl;
        std::cout << std::setw(20) << "-q"
                  << "Executes a simulation with the given QASM file." << std::endl;
        std::cout << std::setw(20) << "-qs"
                  << "Executes a simulation with the given QASM string." << std::endl;
        std::cout << std::setw(20) << "-j"
                  << "Executes a simulation with the given json file with Qiskit Experiment Qobj." << std::endl;
        std::cout << std::setw(20) << "-js"
                  << "Executes a simulation with the given json string." << std::endl;
        std::cout << std::setw(20) << "-t <index>"
                  << "Runs the testing benchmarks for the specific index provided." << std::endl;
        std::cout << std::setw(20) << "-a"
                  << "Runs all testing benchmarks." << std::endl;
        std::cout << std::setw(20) << "-device"
                  << "Path to device model JSON." << std::endl;
        std::cout << std::setw(20) << "-backend_list"
                  << "Lists all the available backends." << std::endl;
        std::cout << std::setw(20) << "-metrics"
                  << "Print the metrics of the circuit." << std::endl;
        std::cout << std::setw(20) << "-backend <name>"
                  << "Sets the backend for your program to the specified one (default: " << backend << "). The backend name string is case-insensitive." << std::endl;
        std::cout << std::setw(20) << "-shots <value>"
                  << "Configures the total number of shots (default: " << total_shots << ")." << std::endl;
        std::cout << std::setw(20) << "-sim <method>"
                  << "Select the simulation method: sv (state vector, default), dm (density matrix). (default: " << simulation_method << ")." << std::endl;
        std::cout << std::setw(20) << "-basis"
                  << "Run the transpiled benchmark circuits which only contain basis gates." << std::endl;
        std::cout << std::setw(20) << "-config"
                  << "Path to configuration JSON file. (default: \"../default_config.json\")" << std::endl;
        std::cout << std::setw(20) << "-initial"
                  << "Initial statevector/density matrix to enable simulation reuse/checkpointing. Must match qubit dimension of input circuit." 
                  << "Expected format is binary sv/dm_real.concat(sv/dm_imag)" << std::endl;
        std::cout << std::setw(20) << "-dump"
                  << "Path to dump binary statevector/density matrix result. (default: \"\", no output dump). Can be reused as a later -initial file" << std::endl;
        std::cout << std::setw(20) << "-layout"
                  << "Path to json mapping logical qubits to device qubits, used when constructing the DM-Sim noise gates. (default: \"\", no mapping)." << std::endl;
        std::cout << std::setw(20) << "-layout-string"
                  << "String denoting qubit layout. Format maps logical qubits (lq, 0...n_qubits) to physical qubits."
                  << " Format is lq0=pq0;lq1=pq1..." << std::endl;
        std::cout << std::setw(20) << "-fidelity"
                  << "Run both DM-Sim and SV-Sim and report the state fidelity." << std::endl;
    }
    if (cmdOptionExists(argv, argv + argc, "-shots"))
    {
        const char *shots_str = getCmdOption(argv, argv + argc, "-shots");
        total_shots = stoi(shots_str);
    }
    if (cmdOptionExists(argv, argv + argc, "-basis"))
    {
        run_with_basis = true;
    }
    if (cmdOptionExists(argv, argv + argc, "-metrics"))
    {
        print_metrics = true;
    }
    if (cmdOptionExists(argv, argv + argc, "-backend_list"))
    {
        BackendManager::print_available_backends();
    }
    if (cmdOptionExists(argv, argv + argc, "-backend"))
    {
        backend = std::string(getCmdOption(argv, argv + argc, "-backend"));
    }
    if (cmdOptionExists(argv, argv + argc, "-device"))
    {
        std::string device = std::string(getCmdOption(argv, argv + argc, "-device"));
        Config::readDeviceConfig(device);
    }
    if (cmdOptionExists(argv, argv + argc, "-sim"))
    {
        simulation_method = std::string(getCmdOption(argv, argv + argc, "-sim"));
    }
    if (cmdOptionExists(argv, argv + argc, "-initial"))
    {
        initfile = std::string(getCmdOption(argv, argv + argc, "-initial"));
    }
    if (cmdOptionExists(argv, argv + argc, "-dump"))
    {
        dumpfile = std::string(getCmdOption(argv, argv + argc, "-dump"));
    }
    if (cmdOptionExists(argv, argv + argc, "-config"))
    {
        config = std::string(getCmdOption(argv, argv + argc, "-config"));
    }
    if (cmdOptionExists(argv, argv + argc, "-layout"))
    {
        layoutfile = std::string(getCmdOption(argv, argv + argc, "-layout"));
        Config::LoadLayoutFromFile(layoutfile);
    }
    if (cmdOptionExists(argv, argv + argc, "-layout-string"))
    {
        std::string layoutstring = std::string(getCmdOption(argv, argv + argc, "-layout-string"));
        Config::readLayoutString(layoutstring);
    }
    if (cmdOptionExists(argv, argv + argc, "-fidelity"))
    {
        report_fidelity = true;
    }
// If MPI or NVSHMEM backend, initialize MPI
#ifdef MPI_ENABLED
    if (backend == "MPI" || backend == "NVGPU_MPI")
    {
        MPI_Init(&argc, &argv);
    }
#endif
    if (cmdOptionExists(argv, argv + argc, "-q"))
    {
        const char *filename = getCmdOption(argv, argv + argc, "-q");
        qasm_parser parser;
        parser.load_qasm_file(filename);
        
        if (!report_fidelity) {
            // Create the backend
            std::shared_ptr<NWQSim::QuantumState> state = BackendManager::create_state(backend, config, parser.num_qubits(), simulation_method);
            if (!state)
            {
                std::cerr << "Failed to create backend\n";
                return 1;
            }
            state->print_config(simulation_method);
            map<string, IdxType> *counts = parser.execute(state, initfile, total_shots, print_metrics);

            // std::vector<size_t> in_bits;

            // for (int i = 0; i < parser.num_qubits(); ++i)
            // {
            //     in_bits.push_back(i);
            // }

            // ValType exp_val_z = state->get_exp_z(in_bits);
            // ValType exp_val_z_all = state->get_exp_z();
            // for (size_t i = 0; i < parser.num_qubits(); ++i)
            // {
            //     in_bits[i] = i;
            // }
            if (dumpfile != "") {
                state->dump_res_state(dumpfile);
            }
            if (state->i_proc == 0)
            {
                print_counts(counts, total_shots);
                // fflush(stdout);
                // printf("exp-val-z: %f\n", exp_val_z);
                // printf("exp-val-z all: %f\n", exp_val_z_all);
            }
            delete counts;
           
        } else {
            
            BackendManager::safe_print("Starting Statevector Simulation...");
            std::string sv_backend = (backend == "NVGPU_MPI" ? "NVGPU" : backend);
            std::shared_ptr<NWQSim::QuantumState> sv_state = BackendManager::create_state(sv_backend, config, parser.num_qubits(), "sv");
            std::shared_ptr<NWQSim::QuantumState> dm_state = BackendManager::create_state(backend, config, parser.num_qubits(), "dm");
            if (!sv_state)
            {
                std::cerr << "Failed to create backend\n";
                return 1;
            }
            if (dm_state->i_proc == 0)
            {
                sv_state->print_config(simulation_method);
                map<string, IdxType> *counts_sv = parser.execute(sv_state, "", total_shots, print_metrics);
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
            map<string, IdxType> *counts_dm = parser.execute(dm_state, "", total_shots, print_metrics);
            if (dm_state->i_proc == 0)
            {
                print_counts(counts_dm, total_shots);
            }

            // if (dumpfile != "") {
            //     dm_state->dump_res_state(dumpfile);
            // }
            ValType fidelity = dm_state->fidelity(sv_state);
            BackendManager::safe_print("State Fidelity: %e\n", fidelity);
            delete counts_dm;
        }
        
    }
    if (cmdOptionExists(argv, argv + argc, "-qs"))
    {
        const char *qasmString = getCmdOption(argv, argv + argc, "-qs");
        qasm_parser parser;
        parser.load_qasm_string(std::string(qasmString)+";\n");
        if (!report_fidelity) {
            // Create the backend
            std::shared_ptr<NWQSim::QuantumState> state = BackendManager::create_state(backend, config, parser.num_qubits(), simulation_method);
            if (!state)
            {
                std::cerr << "Failed to create backend\n";
                return 1;
            }
            state->print_config(simulation_method);
            map<string, IdxType> *counts = parser.execute(state, initfile, total_shots, print_metrics);
            if (state->i_proc == 0)
            {
                print_counts(counts, total_shots);
            }
            if (dumpfile != "") {
                state->dump_res_state(dumpfile);
            }
            delete counts;
        } else {
            BackendManager::safe_print("Starting Statevector Simulation...");
            std::string sv_backend = (backend == "NVGPU_MPI" ? "NVGPU" : backend);
            std::shared_ptr<NWQSim::QuantumState> sv_state = BackendManager::create_state(sv_backend, config, parser.num_qubits(), "sv");
            if (!sv_state)
            {
                std::cerr << "Failed to create backend\n";
                return 1;
            }
            if (sv_state->i_proc == 0)
            {
                sv_state->print_config(simulation_method);
                map<string, IdxType> *counts_sv = parser.execute(sv_state, "", total_shots, print_metrics);
                print_counts(counts_sv, total_shots);
                delete counts_sv;
            }

            BackendManager::safe_print("Starting Density Matrix Simulation...");
            std::shared_ptr<NWQSim::QuantumState> dm_state = BackendManager::create_state(backend, config, parser.num_qubits(), "dm");
            if (!dm_state)
            {
                std::cerr << "Failed to create backend\n";
                return 1;
            }
            dm_state->print_config(simulation_method);
            map<string, IdxType> *counts_dm = parser.execute(dm_state, "", total_shots, print_metrics);
            if (dm_state->i_proc == 0)
            {
                print_counts(counts_dm, total_shots);
            }
            ValType fidelity = dm_state->fidelity(sv_state);
            BackendManager::safe_print("State Fidelity: %e\n", fidelity);
            delete counts_dm;
        }
        
    }
    if (cmdOptionExists(argv, argv + argc, "-j"))
    {
        const char *qobjFile = getCmdOption(argv, argv + argc, "-j");
        qasm_parser parser;
        parser.load_qobj_file(qobjFile);
        if (!report_fidelity) {
            // Create the backend
            std::shared_ptr<NWQSim::QuantumState> state = BackendManager::create_state(backend, config, parser.num_qubits(), simulation_method);
            if (!state)
            {
                std::cerr << "Failed to create backend\n";
                return 1;
            }
            state->print_config(simulation_method);
            map<string, IdxType> *counts = parser.execute(state, initfile, total_shots, print_metrics);
            if (state->i_proc == 0)
            {
                json result_count_json;
                for (const auto& r : (*counts))
                {
                    result_count_json[(r.first)] = r.second;
                }
                cout << "nwq_sim_counts=" << result_count_json.dump() << endl;

                cout << "----------" << endl;
                state->print_res_state();
                cout << "----------" << endl;
            }
            if (dumpfile != "") {
                state->dump_res_state(dumpfile);
            }
            delete counts;
            
        } else {
            BackendManager::safe_print("Starting Statevector Simulation...");
            std::string sv_backend = (backend == "NVGPU_MPI" ? "NVGPU" : backend);
            std::shared_ptr<NWQSim::QuantumState> sv_state = BackendManager::create_state(sv_backend, config, parser.num_qubits(), "sv");
            if (!sv_state)
            {
                std::cerr << "Failed to create backend\n";
                return 1;
            }
            if (sv_state->i_proc == 0)
            {
                sv_state->print_config(simulation_method);
                map<string, IdxType> *counts_sv = parser.execute(sv_state,  "",total_shots, print_metrics);
                json result_count_json;
                for (const auto& r : (*counts_sv))
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
            std::shared_ptr<NWQSim::QuantumState> dm_state = BackendManager::create_state(backend, config, parser.num_qubits(), "dm");
            if (!dm_state)
            {
                std::cerr << "Failed to create backend\n";
                return 1;
            }
            dm_state->print_config(simulation_method);
            map<string, IdxType> *counts_dm = parser.execute(dm_state, "", total_shots, print_metrics);
            if (dm_state->i_proc == 0)
            {
                json result_count_json;
                for (const auto& r : (*counts_dm))
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
    if (cmdOptionExists(argv, argv + argc, "-js"))
    {
        const char *qobjString = getCmdOption(argv, argv + argc, "-js");
        qasm_parser parser;
        parser.load_qobj_string(std::string(qobjString));
        if (!report_fidelity) {
            // Create the backend
            std::shared_ptr<NWQSim::QuantumState> state = BackendManager::create_state(backend, config, parser.num_qubits(), simulation_method);
            if (!state)
            {
                std::cerr << "Failed to create backend\n";
                return 1;
            }
            state->print_config(simulation_method);
            map<string, IdxType> *counts = parser.execute(state, initfile, total_shots, print_metrics);
            if (state->i_proc == 0)
            {
                json result_count_json;
                for (const auto& r : (*counts))
                {
                    result_count_json[(r.first)] = r.second;
                }
                cout << "nwq_sim_counts=" << result_count_json.dump() << endl;

                cout << "----------" << endl;
                state->print_res_state();
                cout << "----------" << endl;
            }
            if (dumpfile != "") {
                state->dump_res_state(dumpfile);
            }
            delete counts;
            
        } else {
            BackendManager::safe_print("Starting Statevector Simulation...");
            std::string sv_backend = (backend == "NVGPU_MPI" ? "NVGPU" : backend);
            std::shared_ptr<NWQSim::QuantumState> sv_state = BackendManager::create_state(sv_backend, config, parser.num_qubits(), "sv");
            if (!sv_state)
            {
                std::cerr << "Failed to create backend\n";
                return 1;
            }
            
            if (sv_state->i_proc == 0)
            {
                sv_state->print_config(simulation_method);
                map<string, IdxType> *counts_sv = parser.execute(sv_state, "", total_shots, print_metrics);
                json result_count_json;
                for (const auto& r : (*counts_sv))
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
            std::shared_ptr<NWQSim::QuantumState> dm_state = BackendManager::create_state(backend, config, parser.num_qubits(), "dm");
            if (!dm_state)
            {
                std::cerr << "Failed to create backend\n";
                return 1;
            }
            dm_state->print_config(simulation_method);
            map<string, IdxType> *counts_dm = parser.execute(dm_state, "", total_shots, print_metrics);
            if (dm_state->i_proc == 0)
            {
                json result_count_json;
                for (const auto& r : (*counts_dm))
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

    if (cmdOptionExists(argv, argv + argc, "-t"))
    {
        total_shots = 16384; //for verification
        int benchmark_index = stoi(getCmdOption(argv, argv + argc, "-t"));
        ValType fidelity = run_brnchmark(backend, config, benchmark_index, total_shots, simulation_method, run_with_basis);
        BackendManager::safe_print("Fidelity between NWQSim and Qiskit Execution: %.4f\n", fidelity);
    }
    if (cmdOptionExists(argv, argv + argc, "-a"))
    {
        bool passed = true;
        for (int benchmark_index = 12; benchmark_index < 36; benchmark_index++)
        {
            ValType fidelity = run_brnchmark(backend, config, benchmark_index, total_shots, simulation_method, run_with_basis);
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

ValType run_brnchmark(std::string backend, const std::string& config, IdxType index, IdxType total_shots, std::string simulation_method, bool is_basis)
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
    std::shared_ptr<NWQSim::QuantumState> state = BackendManager::create_state(backend, config, parser.num_qubits(), simulation_method);
    if (!state)
    {
        std::cerr << "Failed to create backend\n";
        return 1;
    }
    state->print_config(simulation_method);
    map<string, IdxType> *svsim_counts = parser.execute(state, "", total_shots);
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

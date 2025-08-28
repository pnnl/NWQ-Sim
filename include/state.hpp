#pragma once

#include "nwq_util.hpp"
#include "circuit.hpp"
#include "config.hpp"

#ifdef EIGEN
#include "../stabilizer/src/pauli_math.hpp"
#endif

#include "private/gate_factory/sv_gate.hpp"
#include <stdexcept> // For std::runtime_error
#include <vector>
#include <string>
#include <unordered_map>

namespace NWQSim
{
    struct ObservableList
    {
        IdxType *zmasks;
        ValType exp_output;
        ValType *coeffs;
        IdxType numterms;
    };
    // Base class
    class QuantumState
    {
    public:
        QuantumState(SimType _sim_type) : sim_type(_sim_type)
        {
            registerGates();
        } // constructor
        virtual ~QuantumState() {} // virtual destructor

        virtual void print_config(std::string sim_backend)
        {
            Config::printConfig(i_proc, sim_backend);
        };

        virtual void reset_state() = 0;
        virtual void set_seed(IdxType seed) = 0;

        virtual void sim(std::shared_ptr<Circuit> circuit, double& time) = 0;
        virtual void sim2D(std::shared_ptr<Circuit> circuit, std::vector<int> gate_chunks, double &sim_time)
        {
            throw std::runtime_error("2D Sim Not Implemented");
        }
        virtual void sim_batch(std::shared_ptr<NWQSim::Circuit> circuit, IdxType shots, std::vector<std::vector<int32_t>>& all_results, double &sim_time)
        {
            throw std::runtime_error("Sim Batch Not Implemented");
        }
        virtual IdxType *get_results() = 0;
        virtual IdxType measure(IdxType qubit) = 0;
        virtual IdxType *measure_all(IdxType repetition)
        {
            throw std::runtime_error("Measure all not implemented.");
        }
        virtual IdxType **measure_all_long(IdxType shots = 2048)
        {
            throw std::runtime_error("Measure all for results >64 bits not implemented.");
        }
        virtual int *getSingleResult()
        {
            throw std::runtime_error("Single results for >64 qubits not supported on this backend.");
        }

        virtual void set_initial(std::string fpath, std::string format) = 0;
        virtual ValType *get_real() const = 0;
        virtual ValType *get_imag() const = 0;

        virtual ValType get_exp_z() = 0;
        virtual ValType get_exp_z(const std::vector<size_t> &in_bits) = 0;

        virtual ValType fidelity(std::shared_ptr<QuantumState> other)
        {
            throw std::runtime_error("Fidelity computation not Implemented");
        };
        virtual void print_res_state() = 0;
        virtual void dump_res_state(std::string outfile) = 0;

        virtual void save_state()
        {
            throw std::runtime_error("Save State Not Implemented");
        }

        virtual void load_state()
        {
            throw std::runtime_error("Load State Not Implemented");
        }

        virtual void clear_state()
        {
            throw std::runtime_error("Clear Buffer Not Implemented");
        }
#ifdef EIGEN
        virtual Eigen::MatrixXcd get_density_matrix()
        {
            throw std::runtime_error("Density Matrix Return Not Implemented");
        }
#endif
        virtual std::vector<std::vector<int>> get_graph_matrix()
        {
            throw std::runtime_error("Graph Matrix Return Not Implemented");
        }
        virtual std::vector<std::string> get_stabilizers()
        {
            throw std::runtime_error("Get Stabilizers Not Implemented");
        }
        virtual void replace_stabilizer(std::string pauliString, int stabPos)
        {
            throw std::runtime_error("Replace Stabilizers Not Implemented");
        }
        virtual void remove_destabilizers()
        {
            throw std::runtime_error("Remove Destabilizers Not Implemented");
        }
        virtual void delete_all_rows()
        {
            throw std::runtime_error("Delete All Rows Not Implemented");
        }
        virtual int get_num_rows()
        {
            throw std::runtime_error("Get Number of Tableau Rows Not Implemented");
        }
        virtual void stabilizer_count(std::unordered_map<std::string, std::pair<int, int>>& map)
        {
            throw std::runtime_error("Stabilizer Count Not Implemented");
        }
        virtual void add_stabilizer_bits(std::vector<int> new_x, std::vector<int> new_z, int phase_bit = 0)
        {
            throw std::runtime_error("Add Stabilizer Not Implemented");
        }
        virtual bool check_commutation(std::string& pauliString)
        {
            throw std::runtime_error("Check Commutation Not Implemented");
        }
        virtual bool check_row_commutation(std::string pauliString, int row)
        {
            throw std::runtime_error("Check Row Commutation Not Implemented");
        }
        virtual std::pair<std::string,int> get_stabilizer_line(int row)
        {
            throw std::runtime_error("Check Row Commutation Not Implemented");
        }
        virtual void add_stabilizer(std::string pauliString, int phase_bit = 0)
        {
            throw std::runtime_error("Add Stabilizer Not Implemented");
        }
        virtual IdxType get_qubits()
        {
            throw std::runtime_error("Get Qubits Not Implemented");
        }
        virtual void i_rowsum(int h, int i)
        {
            throw std::runtime_error("i_Rowsum Not Implemented");
        }
        virtual void remove_stabilizer(int row_index)
        {
            throw std::runtime_error("Remove Stabilizer Not Implemented");
        }
        virtual void apply_gate(std::string gate, int a, int b = -1)
        {
            throw std::runtime_error("Apply Gate Not Implemented");
        }
        virtual void allocate_measurement_buffers(int max_measurements)
        {
            throw std::runtime_error("Allocate Measurement Buffers Not Implemented");
        }
        virtual void free_measurement_buffers()
        {
            throw std::runtime_error("Free Measurement Buffers Not Implemented");
        }

        IdxType i_proc = 0; // process id
        ValType *buffer_state = nullptr;
        SimType sim_type;

        virtual std::vector<int32_t> get_measurement_results() {
            throw std::runtime_error("Get Measurement Results Not Implemented");
        }

        virtual void clear_measurement_results() {
            throw std::runtime_error("Clear Measurement Results Not Implemented");
        }
    };

} // namespace NWQSim
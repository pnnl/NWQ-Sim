#pragma once

#include "nwq_util.hpp"
#include "circuit.hpp"
#include "config.hpp"

#include "private/gate_factory/sv_gate.hpp"
#include <stdexcept> // For std::runtime_error
#include <vector>
#include <string>

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

        virtual void sim(std::shared_ptr<NWQSim::Circuit> circuit, double& time) = 0;

        virtual IdxType *get_results() = 0;
        virtual IdxType measure(IdxType qubit) = 0;
        virtual IdxType *measure_all(IdxType repetition) = 0;
        virtual void set_initial(std::string fpath, std::string format) = 0;
        virtual ValType *get_real() const = 0;
        virtual ValType *get_imag() const = 0;

        virtual ValType get_exp_z() = 0;
        virtual ValType get_exp_z(const std::vector<size_t> &in_bits) = 0;

        virtual ValType fidelity(std::shared_ptr<QuantumState> other)
        {
            throw std::runtime_error("Fidelity computation not implemented");
        };
        virtual void print_res_state() = 0;
        virtual void dump_res_state(std::string outfile) = 0;

        virtual void save_state()
        {
            throw std::runtime_error("Save State Not implemented");
        }

        virtual void load_state()
        {
            throw std::runtime_error("Load State Not implemented");
        }

        virtual void clear_state()
        {
            throw std::runtime_error("Clear Buffer Not implemented");
        }

        virtual std::vector<std::vector<double>> get_density_matrix()
        {
            throw std::runtime_error("Density Matrix Return Not implemented");
        }
        virtual std::vector<std::vector<int>> get_graph_matrix()
        {
            throw std::runtime_error("Graph Matrix Return Not implemented");
        }
        virtual std::vector<std::string> get_stabilizers()
        {
            throw std::runtime_error("Get Stabilizers Not implemented");
        }


        IdxType i_proc = 0; // process id
        ValType *buffer_state = nullptr;
        SimType sim_type;
    };

} // namespace NWQSim
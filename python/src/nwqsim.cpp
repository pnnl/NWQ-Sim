#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>

#include <algorithm>
#include <cctype>
#include <complex>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "backendManager.hpp"
#include "circuit.hpp"
#include "config.hpp"
#include "nwq_util.hpp"
#include "state.hpp"

namespace py = pybind11;

namespace
{
    using CircuitPtr = std::shared_ptr<NWQSim::Circuit>;
    using StatePtr = std::shared_ptr<NWQSim::QuantumState>;
    using IdxType = NWQSim::IdxType;
    using ValType = NWQSim::ValType;

    std::string to_lower(std::string value)
    {
        std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c)
                       { return static_cast<char>(std::tolower(c)); });
        return value;
    }

    template <typename Method>
    void bind_single_qubit_gate(py::class_<NWQSim::Circuit, CircuitPtr> &cls,
                                const char *name,
                                Method method)
    {
        auto binder = [method](NWQSim::Circuit &self, IdxType qubit) -> NWQSim::Circuit &
        {
            (self.*method)(qubit);
            return self;
        };
        std::string lower_name = to_lower(name);
        cls.def(lower_name.c_str(), binder, py::arg("qubit"), py::return_value_policy::reference_internal);
        cls.def(name, binder, py::arg("qubit"), py::return_value_policy::reference_internal);
    }

    template <typename Method>
    void bind_control_gate(py::class_<NWQSim::Circuit, CircuitPtr> &cls,
                           const char *name,
                           Method method)
    {
        auto binder = [method](NWQSim::Circuit &self, IdxType control, IdxType target) -> NWQSim::Circuit &
        {
            (self.*method)(control, target);
            return self;
        };
        std::string lower_name = to_lower(name);
        cls.def(lower_name.c_str(), binder, py::arg("control"), py::arg("target"), py::return_value_policy::reference_internal);
        cls.def(name, binder, py::arg("control"), py::arg("target"), py::return_value_policy::reference_internal);
    }

    template <typename Method>
    void bind_angle_gate(py::class_<NWQSim::Circuit, CircuitPtr> &cls,
                          const char *name,
                          Method method)
    {
        auto binder = [method](NWQSim::Circuit &self, ValType theta, IdxType qubit) -> NWQSim::Circuit &
        {
            (self.*method)(theta, qubit);
            return self;
        };
        std::string lower_name = to_lower(name);
        cls.def(lower_name.c_str(), binder, py::arg("theta"), py::arg("qubit"), py::return_value_policy::reference_internal);
        cls.def(name, binder, py::arg("theta"), py::arg("qubit"), py::return_value_policy::reference_internal);
    }

    template <typename Method>
    void bind_three_angle_gate(py::class_<NWQSim::Circuit, CircuitPtr> &cls,
                               const char *name,
                               Method method)
    {
        auto binder = [method](NWQSim::Circuit &self, ValType theta, ValType phi, ValType lam, IdxType qubit) -> NWQSim::Circuit &
        {
            (self.*method)(theta, phi, lam, qubit);
            return self;
        };
        std::string lower_name = to_lower(name);
        cls.def(lower_name.c_str(), binder, py::arg("theta"), py::arg("phi"), py::arg("lam"), py::arg("qubit"), py::return_value_policy::reference_internal);
        cls.def(name, binder, py::arg("theta"), py::arg("phi"), py::arg("lam"), py::arg("qubit"), py::return_value_policy::reference_internal);
    }

    template <typename Method>
    void bind_two_angle_gate(py::class_<NWQSim::Circuit, CircuitPtr> &cls,
                             const char *name,
                             Method method)
    {
        auto binder = [method](NWQSim::Circuit &self, ValType first, ValType second, IdxType qubit) -> NWQSim::Circuit &
        {
            (self.*method)(first, second, qubit);
            return self;
        };
        std::string lower_name = to_lower(name);
        cls.def(lower_name.c_str(), binder, py::arg("phi"), py::arg("lam"), py::arg("qubit"), py::return_value_policy::reference_internal);
        cls.def(name, binder, py::arg("phi"), py::arg("lam"), py::arg("qubit"), py::return_value_policy::reference_internal);
    }

    template <typename Method>
    void bind_control_angle_gate(py::class_<NWQSim::Circuit, CircuitPtr> &cls,
                                  const char *name,
                                  Method method)
    {
        auto binder = [method](NWQSim::Circuit &self, ValType theta, IdxType control, IdxType target) -> NWQSim::Circuit &
        {
            (self.*method)(theta, control, target);
            return self;
        };
        std::string lower_name = to_lower(name);
        cls.def(lower_name.c_str(), binder, py::arg("theta"), py::arg("control"), py::arg("target"), py::return_value_policy::reference_internal);
        cls.def(name, binder, py::arg("theta"), py::arg("control"), py::arg("target"), py::return_value_policy::reference_internal);
    }

    template <typename Method>
    void bind_control_phase_gate(py::class_<NWQSim::Circuit, CircuitPtr> &cls,
                                  const char *name,
                                  Method method)
    {
        auto binder = [method](NWQSim::Circuit &self, ValType theta, IdxType control, IdxType target) -> NWQSim::Circuit &
        {
            (self.*method)(theta, control, target);
            return self;
        };
        std::string lower_name = to_lower(name);
        cls.def(lower_name.c_str(), binder, py::arg("theta"), py::arg("control"), py::arg("target"), py::return_value_policy::reference_internal);
        cls.def(name, binder, py::arg("theta"), py::arg("control"), py::arg("target"), py::return_value_policy::reference_internal);
    }

    template <typename Method>
    void bind_two_qubit_angle_gate(py::class_<NWQSim::Circuit, CircuitPtr> &cls,
                                   const char *name,
                                   Method method)
    {
        auto binder = [method](NWQSim::Circuit &self, ValType theta, IdxType q0, IdxType q1) -> NWQSim::Circuit &
        {
            (self.*method)(theta, q0, q1);
            return self;
        };
        std::string lower_name = to_lower(name);
        cls.def(lower_name.c_str(), binder, py::arg("theta"), py::arg("qubit_a"), py::arg("qubit_b"), py::return_value_policy::reference_internal);
        cls.def(name, binder, py::arg("theta"), py::arg("qubit_a"), py::arg("qubit_b"), py::return_value_policy::reference_internal);
    }

    template <typename Method>
    void bind_noise_gate(py::class_<NWQSim::Circuit, CircuitPtr> &cls,
                         const char *name,
                         Method method)
    {
        auto binder = [method](NWQSim::Circuit &self, IdxType qubit, ValType param_a, ValType param_b) -> NWQSim::Circuit &
        {
            (self.*method)(qubit, param_a, param_b);
            return self;
        };
        std::string lower_name = to_lower(name);
        cls.def(lower_name.c_str(), binder, py::arg("qubit"), py::arg("param_a"), py::arg("param_b"), py::return_value_policy::reference_internal);
        cls.def(name, binder, py::arg("qubit"), py::arg("param_a"), py::arg("param_b"), py::return_value_policy::reference_internal);
    }

    template <typename Method>
    void bind_single_param_noise(py::class_<NWQSim::Circuit, CircuitPtr> &cls,
                                 const char *name,
                                 Method method)
    {
        auto binder = [method](NWQSim::Circuit &self, IdxType qubit, ValType param) -> NWQSim::Circuit &
        {
            (self.*method)(qubit, param);
            return self;
        };
        std::string lower_name = to_lower(name);
        cls.def(lower_name.c_str(), binder, py::arg("qubit"), py::arg("param"), py::return_value_policy::reference_internal);
        cls.def(name, binder, py::arg("qubit"), py::arg("param"), py::return_value_policy::reference_internal);
    }

    template <typename Method>
    void bind_control_noise_gate(py::class_<NWQSim::Circuit, CircuitPtr> &cls,
                                 const char *name,
                                 Method method)
    {
        auto binder = [method](NWQSim::Circuit &self, IdxType target, IdxType control, ValType param) -> NWQSim::Circuit &
        {
            (self.*method)(target, control, param);
            return self;
        };
        std::string lower_name = to_lower(name);
        cls.def(lower_name.c_str(), binder, py::arg("target"), py::arg("control"), py::arg("param"), py::return_value_policy::reference_internal);
        cls.def(name, binder, py::arg("target"), py::arg("control"), py::arg("param"), py::return_value_policy::reference_internal);
    }

    template <typename Method>
    void bind_probability_gate(py::class_<NWQSim::Circuit, CircuitPtr> &cls,
                                const char *name,
                                Method method)
    {
        auto binder_single = [method](NWQSim::Circuit &self, IdxType qubit, const std::vector<double> &probabilities) -> NWQSim::Circuit &
        {
            (self.*method)(qubit, probabilities);
            return self;
        };
        std::string lower_name = to_lower(name);
        cls.def(lower_name.c_str(), binder_single, py::arg("qubit"), py::arg("probabilities"), py::return_value_policy::reference_internal);
        cls.def(name, binder_single, py::arg("qubit"), py::arg("probabilities"), py::return_value_policy::reference_internal);
    }

    template <typename Method>
    void bind_control_probability_gate(py::class_<NWQSim::Circuit, CircuitPtr> &cls,
                                        const char *name,
                                        Method method)
    {
        auto binder = [method](NWQSim::Circuit &self, IdxType target, IdxType control, const std::vector<double> &probabilities) -> NWQSim::Circuit &
        {
            (self.*method)(target, control, probabilities);
            return self;
        };
        std::string lower_name = to_lower(name);
        cls.def(lower_name.c_str(), binder, py::arg("target"), py::arg("control"), py::arg("probabilities"), py::return_value_policy::reference_internal);
        cls.def(name, binder, py::arg("target"), py::arg("control"), py::arg("probabilities"), py::return_value_policy::reference_internal);
    }

    void bind_mod_noise(py::class_<NWQSim::Circuit, CircuitPtr> &cls)
    {
        auto binder = [](NWQSim::Circuit &self,
                         const std::string &op,
                         const std::string &noise,
                         ValType value,
                         const std::vector<IdxType> &qubits) -> NWQSim::Circuit &
        {
            self.MOD_NOISE(op, noise, value, qubits);
            return self;
        };
        cls.def("mod_noise", binder,
                py::arg("operation"),
                py::arg("noise_model"),
                py::arg("value"),
                py::arg("qubits"),
                py::return_value_policy::reference_internal);
        cls.def("MOD_NOISE", binder,
                py::arg("operation"),
                py::arg("noise_model"),
                py::arg("value"),
                py::arg("qubits"),
                py::return_value_policy::reference_internal);
    }

    void bind_comb_gate(py::class_<NWQSim::Circuit, CircuitPtr> &cls)
    {
        auto binder = [](NWQSim::Circuit &self, IdxType qubit, ValType phi, ValType lambda, ValType gamma) -> NWQSim::Circuit &
        {
            self.COMB(qubit, phi, lambda, gamma);
            return self;
        };
        cls.def("comb", binder, py::arg("qubit"), py::arg("phi"), py::arg("lambda"), py::arg("gamma"), py::return_value_policy::reference_internal);
        cls.def("COMB", binder, py::arg("qubit"), py::arg("phi"), py::arg("lambda"), py::arg("gamma"), py::return_value_policy::reference_internal);
    }

    void bind_cu_gate(py::class_<NWQSim::Circuit, CircuitPtr> &cls)
    {
        auto binder = [](NWQSim::Circuit &self, ValType theta, ValType phi, ValType lam, ValType gamma, IdxType control, IdxType target) -> NWQSim::Circuit &
        {
            self.CU(theta, phi, lam, gamma, control, target);
            return self;
        };
        cls.def("cu", binder,
                py::arg("theta"), py::arg("phi"), py::arg("lam"), py::arg("gamma"), py::arg("control"), py::arg("target"),
                py::return_value_policy::reference_internal);
        cls.def("CU", binder,
                py::arg("theta"), py::arg("phi"), py::arg("lam"), py::arg("gamma"), py::arg("control"), py::arg("target"),
                py::return_value_policy::reference_internal);
    }

    void bind_ccx_gate(py::class_<NWQSim::Circuit, CircuitPtr> &cls)
    {
        auto binder = [](NWQSim::Circuit &self, IdxType a, IdxType b, IdxType c) -> NWQSim::Circuit &
        {
            self.CCX(a, b, c);
            return self;
        };
        cls.def("ccx", binder, py::arg("control_1"), py::arg("control_2"), py::arg("target"), py::return_value_policy::reference_internal);
        cls.def("CCX", binder, py::arg("control_1"), py::arg("control_2"), py::arg("target"), py::return_value_policy::reference_internal);
    }

    void bind_cswap_gate(py::class_<NWQSim::Circuit, CircuitPtr> &cls)
    {
        auto binder = [](NWQSim::Circuit &self, IdxType control, IdxType a, IdxType b) -> NWQSim::Circuit &
        {
            self.CSWAP(control, a, b);
            return self;
        };
        cls.def("cswap", binder, py::arg("control"), py::arg("qubit_a"), py::arg("qubit_b"), py::return_value_policy::reference_internal);
        cls.def("CSWAP", binder, py::arg("control"), py::arg("qubit_a"), py::arg("qubit_b"), py::return_value_policy::reference_internal);
    }

    void bind_measure_all(py::class_<NWQSim::Circuit, CircuitPtr> &cls)
    {
        auto binder = [](NWQSim::Circuit &self, IdxType shots) -> NWQSim::Circuit &
        {
            self.MA(shots);
            return self;
        };
        cls.def("measure_all", binder, py::arg("shots"), py::return_value_policy::reference_internal);
        cls.def("MA", binder, py::arg("shots"), py::return_value_policy::reference_internal);
    }
}

PYBIND11_MODULE(nwqsim, m)
{
    m.doc() = "Python bindings for NWQSim";

    py::enum_<NWQSim::SimType>(m, "SimType")
        .value("SV", NWQSim::SimType::SV)
        .value("DM", NWQSim::SimType::DM)
        .value("STAB", NWQSim::SimType::STAB)
        .export_values();

    py::class_<NWQSim::CircuitMetrics>(m, "CircuitMetrics")
        .def_readonly("depth", &NWQSim::CircuitMetrics::depth)
        .def_readonly("one_q_gates", &NWQSim::CircuitMetrics::one_q_gates)
        .def_readonly("two_q_gates", &NWQSim::CircuitMetrics::two_q_gates)
        .def_readonly("gate_density", &NWQSim::CircuitMetrics::gate_density)
        .def_readonly("retention_lifespan", &NWQSim::CircuitMetrics::retention_lifespan)
        .def_readonly("measurement_density", &NWQSim::CircuitMetrics::measurement_density)
        .def_readonly("entanglement_variance", &NWQSim::CircuitMetrics::entanglement_variance);

    auto circuit = py::class_<NWQSim::Circuit, CircuitPtr>(m, "Circuit")
                        .def(py::init<IdxType>(), py::arg("num_qubits"))
                        .def("num_qubits", &NWQSim::Circuit::num_qubits)
                        .def("num_gates", &NWQSim::Circuit::num_gates)
                        .def("is_empty", &NWQSim::Circuit::is_empty)
                        .def("reset", &NWQSim::Circuit::reset, py::return_value_policy::reference_internal)
                        .def("clear", &NWQSim::Circuit::clear, py::return_value_policy::reference_internal)
                        .def("to_string", &NWQSim::Circuit::to_string)
                        .def("metrics", &NWQSim::Circuit::circuit_metrics)
                        .def("print_metrics", &NWQSim::Circuit::print_metrics);

    bind_mod_noise(circuit);
    bind_single_qubit_gate(circuit, "X", &NWQSim::Circuit::X);
    bind_single_qubit_gate(circuit, "Y", &NWQSim::Circuit::Y);
    bind_single_qubit_gate(circuit, "Z", &NWQSim::Circuit::Z);
    bind_single_qubit_gate(circuit, "H", &NWQSim::Circuit::H);
    bind_single_qubit_gate(circuit, "S", &NWQSim::Circuit::S);
    bind_single_qubit_gate(circuit, "SDG", &NWQSim::Circuit::SDG);
    bind_single_qubit_gate(circuit, "T", &NWQSim::Circuit::T);
    bind_single_qubit_gate(circuit, "TDG", &NWQSim::Circuit::TDG);
    bind_single_qubit_gate(circuit, "SX", &NWQSim::Circuit::SX);
    bind_single_qubit_gate(circuit, "ID", &NWQSim::Circuit::ID);
    bind_single_qubit_gate(circuit, "M", &NWQSim::Circuit::M);
    bind_single_qubit_gate(circuit, "RESET", &NWQSim::Circuit::RESET);
    bind_angle_gate(circuit, "RI", &NWQSim::Circuit::RI);
    bind_angle_gate(circuit, "RX", &NWQSim::Circuit::RX);
    bind_angle_gate(circuit, "RY", &NWQSim::Circuit::RY);
    bind_angle_gate(circuit, "RZ", &NWQSim::Circuit::RZ);
    bind_angle_gate(circuit, "P", &NWQSim::Circuit::P);
    bind_angle_gate(circuit, "DELAY", &NWQSim::Circuit::DELAY);
    bind_three_angle_gate(circuit, "U", &NWQSim::Circuit::U);
    bind_three_angle_gate(circuit, "U3", &NWQSim::Circuit::U3);
    bind_two_angle_gate(circuit, "U2", &NWQSim::Circuit::U2);
    bind_angle_gate(circuit, "U1", &NWQSim::Circuit::U1);

    bind_control_gate(circuit, "CX", &NWQSim::Circuit::CX);
    bind_control_gate(circuit, "CY", &NWQSim::Circuit::CY);
    bind_control_gate(circuit, "CZ", &NWQSim::Circuit::CZ);
    bind_control_gate(circuit, "CH", &NWQSim::Circuit::CH);
    bind_control_gate(circuit, "CS", &NWQSim::Circuit::CS);
    bind_control_gate(circuit, "CSDG", &NWQSim::Circuit::CSDG);
    bind_control_gate(circuit, "CT", &NWQSim::Circuit::CT);
    bind_control_gate(circuit, "CTDG", &NWQSim::Circuit::CTDG);
    bind_control_gate(circuit, "CSX", &NWQSim::Circuit::CSX);
    bind_control_gate(circuit, "ECR", &NWQSim::Circuit::ECR);
    bind_control_gate(circuit, "SWAP", &NWQSim::Circuit::SWAP);

    bind_control_angle_gate(circuit, "CRX", &NWQSim::Circuit::CRX);
    bind_control_angle_gate(circuit, "CRY", &NWQSim::Circuit::CRY);
    bind_control_angle_gate(circuit, "CRZ", &NWQSim::Circuit::CRZ);
    bind_control_phase_gate(circuit, "CP", &NWQSim::Circuit::CP);
    bind_cu_gate(circuit);
    bind_two_qubit_angle_gate(circuit, "RXX", &NWQSim::Circuit::RXX);
    bind_two_qubit_angle_gate(circuit, "RYY", &NWQSim::Circuit::RYY);
    bind_two_qubit_angle_gate(circuit, "RZZ", &NWQSim::Circuit::RZZ);

    bind_probability_gate(circuit, "CHAN1", &NWQSim::Circuit::CHAN1);
    bind_control_probability_gate(circuit, "CHAN2", &NWQSim::Circuit::CHAN2);
    bind_single_param_noise(circuit, "DEP1", &NWQSim::Circuit::DEP1);
    bind_control_noise_gate(circuit, "DEP2", &NWQSim::Circuit::DEP2);
    bind_noise_gate(circuit, "DAMP", &NWQSim::Circuit::DAMP);
    bind_single_param_noise(circuit, "T1", &NWQSim::Circuit::T1);
    bind_single_param_noise(circuit, "T2", &NWQSim::Circuit::T2);
    bind_single_param_noise(circuit, "EXC", &NWQSim::Circuit::EXC);
    bind_comb_gate(circuit);
    bind_measure_all(circuit);
    bind_ccx_gate(circuit);
    bind_cswap_gate(circuit);

    py::class_<NWQSim::QuantumState, StatePtr>(m, "QuantumState")
        .def("reset", &NWQSim::QuantumState::reset_state)
        .def("set_seed", &NWQSim::QuantumState::set_seed, py::arg("seed"))
        .def("simulate",
             [](StatePtr &state, CircuitPtr &circuit)
             {
                 double sim_time = 0.0;
                 {
                     py::gil_scoped_release release;
                     state->sim(circuit, sim_time);
                 }
                 return sim_time;
             },
             py::arg("circuit"))
        .def("simulate_batch",
             [](StatePtr &state, CircuitPtr &circuit, IdxType shots)
             {
                 std::vector<std::vector<int32_t>> results;
                 double sim_time = 0.0;
                 {
                     py::gil_scoped_release release;
                     state->sim_batch(circuit, shots, results, sim_time);
                 }
                 return py::make_tuple(sim_time, results);
             },
             py::arg("circuit"), py::arg("shots"))
        .def("measure",
             [](StatePtr &state, IdxType qubit)
             {
                 IdxType result = 0;
                 {
                     py::gil_scoped_release release;
                     result = state->measure(qubit);
                 }
                 return result;
             },
             py::arg("qubit"))
        .def("measure_all",
             [](StatePtr &state, IdxType shots)
             {
                 IdxType *raw = nullptr;
                 {
                     py::gil_scoped_release release;
                     raw = state->measure_all(shots);
                 }
                 std::vector<IdxType> values;
                 if (raw != nullptr && shots > 0)
                 {
                     values.assign(raw, raw + shots);
                 }
                 return values;
             },
             py::arg("shots"))
        .def("measurement_results",
             [](StatePtr &state)
             {
                 return state->get_measurement_results();
             })
        .def("print_config",
             [](StatePtr &state, const std::string &backend)
             {
                 state->print_config(backend);
             },
             py::arg("backend"))
        .def("dump_state",
             [](StatePtr &state, const std::string &path)
             {
                 state->dump_res_state(path);
             },
             py::arg("path"))
        .def("statevector",
             [](StatePtr &state)
             {
                 if (state->sim_type != NWQSim::SimType::SV)
                 {
                     throw std::runtime_error("statevector() is only available for statevector backends");
                 }
                 IdxType num_qubits = state->get_qubits();
                 if (num_qubits < 0)
                 {
                     throw std::runtime_error("backend does not expose qubit count");
                 }
                 const auto *real = state->get_real();
                 const auto *imag = state->get_imag();
                 if (!real || !imag)
                 {
                     throw std::runtime_error("backend does not expose amplitude buffers");
                 }
                 const size_t dim = static_cast<size_t>(1ULL << num_qubits);
                 py::array_t<std::complex<double>> out(dim);
                 auto *data = reinterpret_cast<std::complex<double> *>(out.mutable_data());
                 for (size_t i = 0; i < dim; ++i)
                 {
                     data[i] = std::complex<double>(real[i], imag[i]);
                 }
                 return out;
             })
        .def("num_qubits",
             [](StatePtr &state)
             {
                 return state->get_qubits();
             })
        .def("expectation_z",
             [](StatePtr &state, const std::vector<size_t> &bits)
             {
                 return state->get_exp_z(bits);
             },
             py::arg("bits"))
        .def("expectation_z",
             [](StatePtr &state)
             {
                 return state->get_exp_z();
             })
        .def_readonly("sim_type", &NWQSim::QuantumState::sim_type);

    m.def("create_state",
          [](const std::string &backend, IdxType num_qubits, const std::string &method)
          {
              return BackendManager::create_state(backend, num_qubits, method);
          },
          py::arg("backend"), py::arg("num_qubits"), py::arg("method"));

    m.def("available_backends", []()
          {
              std::vector<std::string> backends = {"CPU"};
#ifdef OMP_ENABLED
              backends.push_back("OPENMP");
#endif
#ifdef MPI_ENABLED
              backends.push_back("MPI");
#endif
#ifdef CUDA_ENABLED
              backends.push_back("NVGPU");
#endif
#ifdef CUDA_MPI_ENABLED
              backends.push_back("NVGPU_MPI");
#endif
#ifdef HIP_ENABLED
              backends.push_back("AMDGPU");
#endif
              return backends;
          });

    auto config = m.def_submodule("config");
    config.def("set_print_sim_trace", [](bool value)
               { NWQSim::Config::PRINT_SIM_TRACE = value; });
    config.def("set_enable_noise", [](bool value)
               { NWQSim::Config::ENABLE_NOISE = value; });
    config.def("set_enable_fusion", [](bool value)
               { NWQSim::Config::ENABLE_FUSION = value; });
    config.def("set_random_seed", [](int value)
               { NWQSim::Config::RANDOM_SEED = value; });
    config.def("set_tensor_core", [](bool value)
               { NWQSim::Config::ENABLE_TENSOR_CORE = value; });
    config.def("device_noise_file", [](const std::string &path)
               { NWQSim::Config::device_noise_file = path; });
    config.def("device_layout_file", [](const std::string &path)
               { NWQSim::Config::device_layout_file = path; });

    m.attr("__version__") = "0.1.0";
}

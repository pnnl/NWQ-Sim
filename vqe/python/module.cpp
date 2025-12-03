#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>

#include "ansatz/ansatz.hpp"
#include "execution/vqe_runner.hpp"
#include "execution/adapt_runner.hpp"
#include "hamiltonian_parser.hpp"
#include "jw_transform.hpp"
#include "backend/statevector_cpu.hpp"
#if defined(VQE_ENABLE_CUDA) || defined(VQE_ENABLE_HIP)
#include "backend/statevector_gpu.hpp"
#endif

#include <nlopt.hpp>

#include <algorithm>
#include <cctype>
#include <stdexcept>
#include <sstream>
#include <iomanip>

namespace py = pybind11;

namespace
{

  void assign_parameters(vqe::uccsd_ansatz &ansatz, const std::vector<double> &values)
  {
    auto &circuit = ansatz.mutable_circuit();
    const auto expected = circuit.parameters().size();
    if (values.size() != expected)
    {
      throw std::invalid_argument("parameter vector has incorrect length");
    }
    for (std::size_t i = 0; i < values.size(); ++i)
    {
      circuit.set_parameter(i, values[i]);
    }
  }

  double evaluate_energy_internal(vqe::uccsd_ansatz &ansatz,
                                  const std::vector<vqe::pauli_term> &terms,
                                  const std::vector<double> &parameters,
                                  bool use_gpu)
  {
    assign_parameters(ansatz, parameters);
    const auto &circuit = ansatz.get_circuit();
    const auto qubits = circuit.num_qubits();

#if defined(VQE_ENABLE_CUDA) || defined(VQE_ENABLE_HIP)
    if (use_gpu)
    {
      vqe::backend::statevector_gpu backend(qubits);
      backend.reset();
      backend.apply(circuit);
      return backend.expectation(terms).real();
    }
#else
    if (use_gpu)
    {
      throw std::runtime_error("VQE built without CUDA/HIP support; GPU backend unavailable");
    }
#endif

    vqe::backend::statevector_cpu backend(qubits);
    backend.reset();
    backend.apply(circuit);
    return backend.expectation(terms).real();
  }

  class EnergyEvaluator
  {
  public:
    EnergyEvaluator(const vqe::uccsd_ansatz &ansatz,
                    const std::vector<vqe::pauli_term> &pauli_terms,
                    bool use_gpu)
        : ansatz_(ansatz),
          pauli_terms_(pauli_terms),
          use_gpu_(use_gpu)
    {
    }

    std::size_t parameter_count() const
    {
      return ansatz_.get_circuit().parameters().size();
    }

    std::vector<double> parameters() const
    {
      return ansatz_.get_circuit().parameters();
    }

    void set_parameters(const std::vector<double> &values)
    {
      assign_parameters(ansatz_, values);
    }

    double energy(const std::vector<double> &values)
    {
      py::gil_scoped_release release;
      return evaluate_energy_internal(ansatz_, pauli_terms_, values, use_gpu_);
    }

  private:
    vqe::uccsd_ansatz ansatz_;
    std::vector<vqe::pauli_term> pauli_terms_;
    bool use_gpu_ = false;
  };

  std::string algorithm_to_string(nlopt::algorithm algo)
  {
    const char *name = nlopt::algorithm_name(algo);
    if (name)
    {
      return std::string(name);
    }
    return "";
  }

  nlopt::algorithm algorithm_from_string(const std::string &name)
  {
    const auto value = nlopt_algorithm_from_string(name.c_str());
    if (value < 0 || value >= nlopt::NUM_ALGORITHMS)
    {
      throw std::invalid_argument("Unknown NLopt algorithm: " + name);
    }
    return static_cast<nlopt::algorithm>(value);
  }

  std::string normalize_mode_string(const std::string &name)
  {
    std::string value = name;
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
      return static_cast<char>(std::tolower(c));
    });
    return value;
  }

  std::string mode_to_string(vqe::run_mode mode)
  {
    switch (mode)
    {
    case vqe::run_mode::standard:
      return "vqe";
    case vqe::run_mode::adapt:
      return "adapt";
    }
    return "vqe";
  }

  vqe::run_mode mode_from_string(const std::string &name)
  {
    const auto value = normalize_mode_string(name);
    if (value == "vqe" || value == "standard")
    {
      return vqe::run_mode::standard;
    }
    if (value == "adapt" || value == "adapt-vqe" || value == "adaptvqe")
    {
      return vqe::run_mode::adapt;
    }
    throw std::invalid_argument("Unknown mode: " + name + ". Expected 'vqe' or 'adapt'.");
  }

} // namespace

PYBIND11_MODULE(_core, m)
{
  m.doc() = "Python bindings for the NWQ-Sim VQE module";

  // ============================================================================
  // Enumerations with string representations
  // ============================================================================

  py::enum_<vqe::fermion_op_kind>(m, "FermionOpKind", "Kind of fermion operator")
      .value("CREATION", vqe::fermion_op_kind::creation, "Creation operator")
      .value("ANNIHILATION", vqe::fermion_op_kind::annihilation, "Annihilation operator")
      .def("__str__", [](vqe::fermion_op_kind kind)
           { return kind == vqe::fermion_op_kind::creation ? "creation" : "annihilation"; })
      .def("__repr__", [](vqe::fermion_op_kind kind)
           { return kind == vqe::fermion_op_kind::creation ? "FermionOpKind.CREATION" : "FermionOpKind.ANNIHILATION"; });

  py::enum_<vqe::spin>(m, "Spin", "Electron spin state")
      .value("UP", vqe::spin::up, "Spin-up electron")
      .value("DOWN", vqe::spin::down, "Spin-down electron")
      .def("__str__", [](vqe::spin s)
           { return s == vqe::spin::up ? "up" : "down"; })
      .def("__repr__", [](vqe::spin s)
           { return s == vqe::spin::up ? "Spin.UP" : "Spin.DOWN"; });

  py::enum_<vqe::orbital_kind>(m, "OrbitalKind", "Type of molecular orbital")
      .value("OCCUPIED", vqe::orbital_kind::occupied, "Occupied orbital")
      .value("VIRTUAL", vqe::orbital_kind::virtual_orbital, "Virtual (unoccupied) orbital")
      .def("__str__", [](vqe::orbital_kind kind)
           { return kind == vqe::orbital_kind::occupied ? "occupied" : "virtual"; })
      .def("__repr__", [](vqe::orbital_kind kind)
           { return kind == vqe::orbital_kind::occupied ? "OrbitalKind.OCCUPIED" : "OrbitalKind.VIRTUAL"; });

  py::enum_<vqe::operator_kind>(m, "OperatorKind", "Type of fermion operator")
      .value("CREATION", vqe::operator_kind::creation, "Creation operator")
      .value("ANNIHILATION", vqe::operator_kind::annihilation, "Annihilation operator")
      .def("__str__", [](vqe::operator_kind kind)
           { return kind == vqe::operator_kind::creation ? "creation" : "annihilation"; })
      .def("__repr__", [](vqe::operator_kind kind)
           { return kind == vqe::operator_kind::creation ? "OperatorKind.CREATION" : "OperatorKind.ANNIHILATION"; });

  // ============================================================================
  // Core data structures with Python-friendly interfaces
  // ============================================================================

  py::class_<vqe::fermion_op>(m, "BasicFermionOp", "Basic fermion operator with index and kind")
      .def(py::init<>(), "Create default fermion operator")
      .def(py::init<std::size_t, vqe::fermion_op_kind>(),
           py::arg("index"), py::arg("kind"), "Create fermion operator")
      .def_property("index", [](const vqe::fermion_op &op)
                    { return op.index; }, [](vqe::fermion_op &op, std::size_t idx)
                    { op.index = idx; }, "Orbital index")
      .def_property("kind", [](const vqe::fermion_op &op)
                    { return op.kind; }, [](vqe::fermion_op &op, vqe::fermion_op_kind k)
                    { op.kind = k; }, "Operator kind (creation/annihilation)")
      .def("__repr__", [](const vqe::fermion_op &op)
           {
        std::ostringstream oss;
        oss << "BasicFermionOp(index=" << op.index << ", kind=" 
            << (op.kind == vqe::fermion_op_kind::creation ? "CREATION" : "ANNIHILATION") << ")";
        return oss.str(); })
      .def("__str__", [](const vqe::fermion_op &op)
           {
        std::ostringstream oss;
        oss << (op.kind == vqe::fermion_op_kind::creation ? "a†" : "a") << "_" << op.index;
        return oss.str(); });

  py::class_<vqe::fermion_term>(m, "FermionTerm", "Term in a fermionic Hamiltonian")
      .def(py::init<>(), "Create default fermion term")
      .def_property("coefficient", [](const vqe::fermion_term &term)
                    { return term.coefficient; }, [](vqe::fermion_term &term, std::complex<double> coeff)
                    { term.coefficient = coeff; }, "Complex coefficient")
      .def_property("operators", [](const vqe::fermion_term &term)
                    { return term.operators; }, [](vqe::fermion_term &term, const std::vector<vqe::fermion_op> &ops)
                    { term.operators = ops; }, "List of fermion operators")
      .def("__repr__", [](const vqe::fermion_term &term)
           {
        std::ostringstream oss;
        oss << "FermionTerm(coeff=" << std::fixed << std::setprecision(6) << term.coefficient
            << ", ops=" << term.operators.size() << ")";
        return oss.str(); })
      .def("__str__", [](const vqe::fermion_term &term)
           {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(6) << term.coefficient;
        for (const auto& op : term.operators) {
          oss << " " << (op.kind == vqe::fermion_op_kind::creation ? "a†" : "a") << "_" << op.index;
        }
        return oss.str(); });

  py::class_<vqe::fermion_operator>(m, "FermionOperator", "Extended fermion operator with orbital information")
      .def(py::init<>(), "Create default fermion operator")
      .def_property("orbital_index", [](const vqe::fermion_operator &op)
                    { return op.orbital_index; }, [](vqe::fermion_operator &op, std::size_t idx)
                    { op.orbital_index = idx; }, "Orbital index")
      .def_property("orbital", [](const vqe::fermion_operator &op)
                    { return op.orbital; }, [](vqe::fermion_operator &op, vqe::orbital_kind kind)
                    { op.orbital = kind; }, "Orbital type")
      .def_property("spin_state", [](const vqe::fermion_operator &op)
                    { return op.spin_state; }, [](vqe::fermion_operator &op, vqe::spin s)
                    { op.spin_state = s; }, "Spin state")
      .def_property("op", [](const vqe::fermion_operator &op)
                    { return op.op; }, [](vqe::fermion_operator &op, vqe::operator_kind kind)
                    { op.op = kind; }, "Operator type")
      .def("__repr__", [](const vqe::fermion_operator &op)
           {
        std::ostringstream oss;
        oss << "FermionOperator(orbital=" << op.orbital_index
            << ", " << (op.orbital == vqe::orbital_kind::occupied ? "occupied" : "virtual")
            << ", " << (op.spin_state == vqe::spin::up ? "up" : "down")
            << ", " << (op.op == vqe::operator_kind::creation ? "creation" : "annihilation") << ")";
        return oss.str(); });

  py::class_<vqe::hamiltonian_data>(m, "HamiltonianData", "Fermionic Hamiltonian data structure")
      .def(py::init<>(), "Create empty Hamiltonian")
      .def_property("constant", [](const vqe::hamiltonian_data &data)
                    { return data.constant; }, [](vqe::hamiltonian_data &data, std::complex<double> c)
                    { data.constant = c; }, "Constant energy offset")
      .def_property("terms", [](const vqe::hamiltonian_data &data)
                    { return data.terms; }, [](vqe::hamiltonian_data &data, const std::vector<vqe::fermion_term> &terms)
                    { data.terms = terms; }, "List of fermion terms")
      .def_property("max_index", [](const vqe::hamiltonian_data &data)
                    { return data.max_index; }, [](vqe::hamiltonian_data &data, std::size_t idx)
                    { data.max_index = idx; }, "Maximum orbital index")
      .def_property_readonly("num_qubits", &vqe::hamiltonian_data::num_qubits, "Number of qubits needed for Jordan-Wigner encoding")
      .def("__len__", [](const vqe::hamiltonian_data &data)
           { return data.terms.size(); })
      .def("__repr__", [](const vqe::hamiltonian_data &data)
           {
        std::ostringstream oss;
        oss << "HamiltonianData(terms=" << data.terms.size() 
            << ", qubits=" << data.num_qubits() 
            << ", constant=" << std::fixed << std::setprecision(6) << data.constant << ")";
        return oss.str(); });

  py::class_<vqe::pauli_term>(m, "PauliTerm", "Pauli string term in qubit Hamiltonian")
      .def(py::init<>(), "Create identity Pauli term")
      .def(py::init<std::uint64_t, std::uint64_t, std::complex<double>>(),
           py::arg("x_mask"), py::arg("z_mask"), py::arg("coefficient") = 1.0,
           "Create Pauli term with masks")
      .def_property("x_mask", [](const vqe::pauli_term &term)
                    { return term.x_mask; }, [](vqe::pauli_term &term, std::uint64_t mask)
                    { term.x_mask = mask; }, "Bit mask for X operators")
      .def_property("z_mask", [](const vqe::pauli_term &term)
                    { return term.z_mask; }, [](vqe::pauli_term &term, std::uint64_t mask)
                    { term.z_mask = mask; }, "Bit mask for Z operators")
      .def_property("coefficient", [](const vqe::pauli_term &term)
                    { return term.coefficient; }, [](vqe::pauli_term &term, std::complex<double> coeff)
                    { term.coefficient = coeff; }, "Complex coefficient")
      .def("weight", [](const vqe::pauli_term &term)
           { return __builtin_popcountll(term.x_mask) + __builtin_popcountll(term.z_mask); }, "Number of non-identity Pauli operators")
      .def("to_string", [](const vqe::pauli_term &term, std::size_t n_qubits = 0)
           {
        if (n_qubits == 0) {
          n_qubits = std::max(64 - __builtin_clzll(term.x_mask | term.z_mask | 1), 1);
        }
        return vqe::pauli_to_string(term.x_mask, term.z_mask, n_qubits); }, py::arg("n_qubits") = 0, "Convert to Pauli string representation")
      .def("__repr__", [](const vqe::pauli_term &term)
           {
        std::ostringstream oss;
        oss << "PauliTerm(coeff=" << std::fixed << std::setprecision(6) << term.coefficient
            << ", weight=" << (__builtin_popcountll(term.x_mask) + __builtin_popcountll(term.z_mask)) << ")";
        return oss.str(); })
      .def("__str__", [](const vqe::pauli_term &term)
           {
        std::size_t n_qubits = std::max(64 - __builtin_clzll(term.x_mask | term.z_mask | 1), 1);
        return std::to_string(term.coefficient.real()) + " " + vqe::pauli_to_string(term.x_mask, term.z_mask, n_qubits); });

  py::class_<vqe::molecular_environment>(m, "System", "Molecular system specification")
      .def(py::init<>(), "Create default molecular environment")
      .def(py::init<std::size_t, std::size_t, bool, double>(),
           py::arg("n_spatial"), py::arg("n_electrons"), py::arg("xacc_indexing") = true,
           py::arg("constant_energy") = 0.0, "Create molecular environment")
      .def_property("n_spatial", [](const vqe::molecular_environment &env)
                    { return env.n_spatial; }, [](vqe::molecular_environment &env, std::size_t n)
                    { env.n_spatial = n; }, "Number of spatial orbitals")
      .def_property("n_electrons", [](const vqe::molecular_environment &env)
                    { return env.n_electrons; }, [](vqe::molecular_environment &env, std::size_t n)
                    { env.n_electrons = n; }, "Number of electrons")
      .def_property("xacc_indexing", [](const vqe::molecular_environment &env)
                    { return env.xacc_indexing; }, [](vqe::molecular_environment &env, bool use_xacc)
                    { env.xacc_indexing = use_xacc; }, "Use XACC ordering convention")
      .def_property("constant_energy", [](const vqe::molecular_environment &env)
                    { return env.constant_energy; }, [](vqe::molecular_environment &env, double energy)
                    { env.constant_energy = energy; }, "Nuclear repulsion energy")
      .def("occupied_orbitals", &vqe::molecular_environment::occupied_orbitals, "Get list of occupied orbital indices")
      .def("virtual_orbitals", &vqe::molecular_environment::virtual_orbitals, "Get list of virtual orbital indices")
      .def("total_qubits", &vqe::molecular_environment::total_qubits, "Total number of qubits needed")
      .def("qubit_index", &vqe::molecular_environment::qubit_index, py::arg("orbital_index"), py::arg("kind"), py::arg("spin"), "Get qubit index for given orbital")
      .def("__repr__", [](const vqe::molecular_environment &env)
           {
        std::ostringstream oss;
        oss << "System(orbitals=" << env.n_spatial 
            << ", electrons=" << env.n_electrons
            << ", qubits=" << env.total_qubits() << ")";
        return oss.str(); });

  py::class_<vqe::uccsd_ansatz>(m, "UCCSD", "UCCSD ansatz for VQE calculations")
      .def(py::init<const vqe::molecular_environment &, std::size_t, std::size_t>(),
           py::arg("environment"), py::arg("trotter_steps") = 1, py::arg("symmetry_level") = 3,
           "Create UCCSD ansatz")
      .def("build", &vqe::uccsd_ansatz::build, "Build the ansatz circuit")
      .def("unique_parameter_count", &vqe::uccsd_ansatz::unique_parameter_count,
           "Number of unique variational parameters")
      .def_property_readonly("parameters", [](const vqe::uccsd_ansatz &ansatz)
                             { return ansatz.get_circuit().parameters(); }, "Current parameter values")
      .def("set_parameters", [](vqe::uccsd_ansatz &ansatz, const std::vector<double> &values)
           {
            if (values.size() != ansatz.get_circuit().parameters().size()) {
              throw std::invalid_argument("Parameter vector length mismatch: expected " + 
                std::to_string(ansatz.get_circuit().parameters().size()) + 
                ", got " + std::to_string(values.size()));
            }
            assign_parameters(ansatz, values); }, py::arg("parameters"), "Set variational parameters")
      .def_property_readonly("num_qubits", [](const vqe::uccsd_ansatz &ansatz)
                             { return ansatz.get_circuit().num_qubits(); }, "Number of qubits in the circuit")
      .def_property_readonly("num_parameters", [](const vqe::uccsd_ansatz &ansatz)
                             { return ansatz.get_circuit().parameters().size(); }, "Total number of parameters")
      .def("copy", [](const vqe::uccsd_ansatz &ansatz)
           {
             return ansatz; // Creates a copy
           },
           "Create a copy of the ansatz")
      .def("__repr__", [](const vqe::uccsd_ansatz &ansatz)
           {
        std::ostringstream oss;
        oss << "UCCSD(qubits=" << ansatz.get_circuit().num_qubits()
            << ", parameters=" << ansatz.get_circuit().parameters().size() << ")";
        return oss.str(); });

  // ============================================================================
  // Configuration classes with validation and default values
  // ============================================================================

  py::class_<vqe::vqe_options>(m, "Options", "Unified configuration for VQE and ADAPT-VQE")
      .def(py::init<>(), "Create default options")
      .def_property("mode", [](const vqe::vqe_options &opts)
                    { return mode_to_string(opts.mode); }, [](vqe::vqe_options &opts, const std::string &value)
                    { opts.mode = mode_from_string(value); }, "Solver mode: 'vqe' or 'adapt'")
      .def_property("verbose", [](const vqe::vqe_options &opts)
                    { return opts.verbose; }, [](vqe::vqe_options &opts, bool verbose)
                    { opts.verbose = verbose; }, "Enable verbose logging")
      .def_property("use_gpu", [](const vqe::vqe_options &opts)
                    { return opts.use_gpu; }, [](vqe::vqe_options &opts, bool use_gpu)
                    { opts.use_gpu = use_gpu; }, "Use GPU acceleration when available")
      .def_property("use_xacc_indexing", [](const vqe::vqe_options &opts)
                    { return opts.use_xacc_indexing; }, [](vqe::vqe_options &opts, bool value)
                    { opts.use_xacc_indexing = value; }, "Use XACC orbital ordering convention")
      .def_property("random_seed", [](const vqe::vqe_options &opts) -> std::optional<unsigned>
                    { return opts.random_seed; }, [](vqe::vqe_options &opts, std::optional<unsigned> seed)
                    { opts.random_seed = seed; }, "Optional random seed for deterministic behaviours")
      .def_property("trotter_steps", [](const vqe::vqe_options &opts)
                    { return opts.trotter_steps; }, [](vqe::vqe_options &opts, std::size_t steps)
                    {
          if (steps == 0)
          {
            throw std::invalid_argument("trotter_steps must be positive");
          }
          opts.trotter_steps = steps;
        }, "Number of Trotter steps for UCCSD")
      .def_property("symmetry_level", [](const vqe::vqe_options &opts)
                    { return opts.symmetry_level; }, [](vqe::vqe_options &opts, std::size_t level)
                    { opts.symmetry_level = level; }, "Symmetry reduction level")
      .def_property("lower_bound", [](const vqe::vqe_options &opts)
                    { return opts.lower_bound; }, [](vqe::vqe_options &opts, double bound)
                    {
          opts.lower_bound = bound;
          opts.adapt_lower_bound = bound;
        }, "Lower bound for optimizer parameters")
      .def_property("upper_bound", [](const vqe::vqe_options &opts)
                    { return opts.upper_bound; }, [](vqe::vqe_options &opts, double bound)
                    {
          opts.upper_bound = bound;
          opts.adapt_upper_bound = bound;
        }, "Upper bound for optimizer parameters")
      .def_property("max_evaluations", [](const vqe::vqe_options &opts)
                    { return opts.max_evaluations; }, [](vqe::vqe_options &opts, std::size_t value)
                    {
          opts.max_evaluations = value;
          opts.adapt_max_evaluations = value;
        }, "Maximum objective evaluations")
      .def_property("relative_tolerance", [](const vqe::vqe_options &opts)
                    { return opts.relative_tolerance; }, [](vqe::vqe_options &opts, double tol)
                    {
          if (tol < 0)
          {
            throw std::invalid_argument("relative_tolerance must be non-negative");
          }
          opts.relative_tolerance = tol;
          opts.adapt_relative_tolerance = tol;
        }, "Relative convergence tolerance")
      .def_property("absolute_tolerance", [](const vqe::vqe_options &opts)
                    { return opts.absolute_tolerance; }, [](vqe::vqe_options &opts, double tol)
                    {
          if (tol < 0)
          {
            throw std::invalid_argument("absolute_tolerance must be non-negative");
          }
          opts.absolute_tolerance = tol;
          opts.adapt_absolute_tolerance = tol;
        }, "Absolute convergence tolerance")
      .def_property("stop_value", [](const vqe::vqe_options &opts)
                    { return opts.stop_value; }, [](vqe::vqe_options &opts, double value)
                    {
          opts.stop_value = value;
          opts.adapt_stop_value = value;
        }, "Target objective value for early stopping")
      .def_property("max_time", [](const vqe::vqe_options &opts)
                    { return opts.max_time; }, [](vqe::vqe_options &opts, double seconds)
                    {
          if (seconds <= 0)
          {
            throw std::invalid_argument("max_time must be positive");
          }
          opts.max_time = seconds;
          opts.adapt_max_time = seconds;
        }, "Maximum wall-clock time (seconds)")
      .def_property("optimizer", [](const vqe::vqe_options &opts)
                    { return algorithm_to_string(opts.optimizer); }, [](vqe::vqe_options &opts, const std::string &name)
                    {
          auto algo = algorithm_from_string(name);
          opts.optimizer = algo;
          opts.adapt_optimizer = algo;
        }, "Primary NLopt optimizer name")
      .def_property("algorithm_parameters", [](const vqe::vqe_options &opts)
                    { return opts.algorithm_parameters; }, [](vqe::vqe_options &opts, const std::unordered_map<std::string, double> &params)
                    {
          opts.algorithm_parameters = params;
          opts.adapt_algorithm_parameters = params;
        }, "Algorithm-specific NLopt parameters")
      .def_property("initial_parameters", [](const vqe::vqe_options &opts)
                    { return opts.initial_parameters; }, [](vqe::vqe_options &opts, const std::vector<double> &values)
                    { opts.initial_parameters = values; }, "Initial parameter guess for VQE")
      .def_property("iteration_improvement_tolerance", [](const vqe::vqe_options &opts)
                    { return opts.iteration_improvement_tolerance; }, [](vqe::vqe_options &opts, double tol)
                    {
          if (tol < 0)
          {
            throw std::invalid_argument("iteration_improvement_tolerance must be non-negative");
          }
          opts.iteration_improvement_tolerance = tol;
        }, "Minimum improvement required to log iteration progress")
      .def_property("gradient_step", [](const vqe::vqe_options &opts)
                    { return opts.gradient_step; }, [](vqe::vqe_options &opts, double step)
                    {
          if (step <= 0)
          {
            throw std::invalid_argument("gradient_step must be positive");
          }
          opts.gradient_step = step;
        }, "Forward-difference step size for optimizer gradients")
      .def_property("status_interval", [](const vqe::vqe_options &opts)
                    { return opts.status_interval; }, [](vqe::vqe_options &opts, std::size_t interval)
                    {
          if (interval == 0)
          {
            throw std::invalid_argument("status_interval must be positive");
          }
          opts.status_interval = interval;
        }, "Iteration interval for status updates")
      .def_property("adapt_max_iterations", [](const vqe::vqe_options &opts)
                    { return opts.adapt_max_iterations; }, [](vqe::vqe_options &opts, std::size_t value)
                    {
          if (value == 0)
          {
            throw std::invalid_argument("adapt_max_iterations must be positive");
          }
          opts.adapt_max_iterations = value;
        }, "Maximum number of ADAPT iterations")
      .def_property("adapt_gradient_step", [](const vqe::vqe_options &opts)
                    { return opts.adapt_gradient_step; }, [](vqe::vqe_options &opts, double step)
                    {
          if (step <= 0)
          {
            throw std::invalid_argument("adapt_gradient_step must be positive");
          }
          opts.adapt_gradient_step = step;
        }, "Finite-difference step for ADAPT gradients")
      .def_property("adapt_gradient_tolerance", [](const vqe::vqe_options &opts)
                    { return opts.adapt_gradient_tolerance; }, [](vqe::vqe_options &opts, double tol)
                    {
          if (tol < 0)
          {
            throw std::invalid_argument("adapt_gradient_tolerance must be non-negative");
          }
          opts.adapt_gradient_tolerance = tol;
        }, "Gradient norm tolerance for ADAPT termination")
      .def_property("adapt_energy_tolerance", [](const vqe::vqe_options &opts)
                    { return opts.adapt_energy_tolerance; }, [](vqe::vqe_options &opts, double tol)
                    {
          if (tol < 0)
          {
            throw std::invalid_argument("adapt_energy_tolerance must be non-negative");
          }
          opts.adapt_energy_tolerance = tol;
        }, "Energy tolerance for ADAPT convergence")
      .def_property("adapt_log_memory", [](const vqe::vqe_options &opts)
                    { return opts.adapt_log_memory; }, [](vqe::vqe_options &opts, bool enabled)
                    { opts.adapt_log_memory = enabled; }, "Log memory usage during ADAPT execution")
      .def_property("adapt_optimizer", [](const vqe::vqe_options &opts)
                    { return algorithm_to_string(opts.adapt_optimizer); }, [](vqe::vqe_options &opts, const std::string &name)
                    { opts.adapt_optimizer = algorithm_from_string(name); }, "NLopt optimizer used within ADAPT iterations")
      .def_property("adapt_lower_bound", [](const vqe::vqe_options &opts)
                    { return opts.adapt_lower_bound; }, [](vqe::vqe_options &opts, double bound)
                    { opts.adapt_lower_bound = bound; }, "Lower bound for ADAPT optimizer")
      .def_property("adapt_upper_bound", [](const vqe::vqe_options &opts)
                    { return opts.adapt_upper_bound; }, [](vqe::vqe_options &opts, double bound)
                    { opts.adapt_upper_bound = bound; }, "Upper bound for ADAPT optimizer")
      .def_property("adapt_max_evaluations", [](const vqe::vqe_options &opts)
                    { return opts.adapt_max_evaluations; }, [](vqe::vqe_options &opts, std::size_t value)
                    { opts.adapt_max_evaluations = value; }, "Maximum evaluations per ADAPT micro-optimization")
      .def_property("adapt_relative_tolerance", [](const vqe::vqe_options &opts)
                    { return opts.adapt_relative_tolerance; }, [](vqe::vqe_options &opts, double tol)
                    {
          if (tol < 0)
          {
            throw std::invalid_argument("adapt_relative_tolerance must be non-negative");
          }
          opts.adapt_relative_tolerance = tol;
        }, "Relative tolerance for ADAPT optimizer")
      .def_property("adapt_absolute_tolerance", [](const vqe::vqe_options &opts)
                    { return opts.adapt_absolute_tolerance; }, [](vqe::vqe_options &opts, double tol)
                    {
          if (tol < 0)
          {
            throw std::invalid_argument("adapt_absolute_tolerance must be non-negative");
          }
          opts.adapt_absolute_tolerance = tol;
        }, "Absolute tolerance for ADAPT optimizer")
      .def_property("adapt_stop_value", [](const vqe::vqe_options &opts)
                    { return opts.adapt_stop_value; }, [](vqe::vqe_options &opts, double value)
                    { opts.adapt_stop_value = value; }, "Stop value for ADAPT optimizer")
      .def_property("adapt_max_time", [](const vqe::vqe_options &opts)
                    { return opts.adapt_max_time; }, [](vqe::vqe_options &opts, double seconds)
                    {
          if (seconds <= 0)
          {
            throw std::invalid_argument("adapt_max_time must be positive");
          }
          opts.adapt_max_time = seconds;
        }, "Maximum time per ADAPT optimizer call")
      .def_property("adapt_iteration_improvement_tolerance", [](const vqe::vqe_options &opts)
                    { return opts.adapt_iteration_improvement_tolerance; }, [](vqe::vqe_options &opts, double tol)
                    {
          if (tol < 0)
          {
            throw std::invalid_argument("adapt_iteration_improvement_tolerance must be non-negative");
          }
          opts.adapt_iteration_improvement_tolerance = tol;
        }, "Iteration improvement tolerance for ADAPT logging")
      .def_property("adapt_status_interval", [](const vqe::vqe_options &opts)
                    { return opts.adapt_status_interval; }, [](vqe::vqe_options &opts, std::size_t interval)
                    {
          if (interval == 0)
          {
            throw std::invalid_argument("adapt_status_interval must be positive");
          }
          opts.adapt_status_interval = interval;
        }, "Status print interval for ADAPT optimizer")
      .def_property("adapt_algorithm_parameters", [](const vqe::vqe_options &opts)
                    { return opts.adapt_algorithm_parameters; }, [](vqe::vqe_options &opts, const std::unordered_map<std::string, double> &params)
                    { opts.adapt_algorithm_parameters = params; }, "Algorithm parameters for ADAPT optimizer")
      .def_property("adapt_initial_parameters", [](const vqe::vqe_options &opts)
                    { return opts.adapt_initial_parameters; }, [](vqe::vqe_options &opts, const std::vector<double> &values)
                    { opts.adapt_initial_parameters = values; }, "Initial parameter guess for ADAPT micro-optimization")
      .def("__repr__", [](const vqe::vqe_options &opts)
           {
        std::ostringstream oss;
        oss << "Options(mode='" << mode_to_string(opts.mode)
            << "', optimizer='" << algorithm_to_string(opts.optimizer)
            << "', max_eval=" << opts.max_evaluations
            << ", gpu=" << (opts.use_gpu ? "True" : "False") << ")";
        return oss.str(); });

  // ============================================================================
  // Result classes with rich information
  // ============================================================================

  py::class_<vqe::vqe_result>(m, "VQEResult", "Results from VQE optimization")
      .def(py::init<>(), "Create empty VQE result")
      .def_property_readonly("energy", [](const vqe::vqe_result &result)
                             { return result.energy; }, "Optimized ground state energy")
      .def_property_readonly("parameters", [](const vqe::vqe_result &result)
                             { return result.parameters; }, "Optimal variational parameters")
      .def_property_readonly("evaluations", [](const vqe::vqe_result &result)
                             { return result.evaluations; }, "Number of energy evaluations performed")
      .def_property_readonly("converged", [](const vqe::vqe_result &result)
                             { return result.converged; }, "Whether optimization converged")
      .def_property_readonly("success", [](const vqe::vqe_result &result)
                             { return result.converged; }, "Alias for converged")
      .def("__repr__", [](const vqe::vqe_result &result)
           {
        std::ostringstream oss;
        oss << "VQEResult(energy=" << std::fixed << std::setprecision(8) << result.energy
            << ", evaluations=" << result.evaluations 
            << ", converged=" << (result.converged ? "True" : "False") << ")";
        return oss.str(); });

  py::class_<vqe::adapt_result>(m, "AdaptResult", "Results from ADAPT-VQE optimization")
      .def(py::init<>(), "Create empty ADAPT result")
      .def_property_readonly("energy", [](const vqe::adapt_result &result)
                             { return result.energy; }, "Final ground state energy")
      .def_property_readonly("parameters", [](const vqe::adapt_result &result)
                             { return result.parameters; }, "Final variational parameters")
      .def_property_readonly("selected_indices", [](const vqe::adapt_result &result)
                             { return result.selected_indices; }, "Indices of selected operators")
      .def_property_readonly("selected_labels", [](const vqe::adapt_result &result)
                             { return result.selected_labels; }, "Labels of selected operators")
      .def_property_readonly("iterations", [](const vqe::adapt_result &result)
                             { return result.iterations; }, "Number of ADAPT iterations performed")
      .def_property_readonly("energy_evaluations", [](const vqe::adapt_result &result)
                             { return result.energy_evaluations; }, "Total number of energy evaluations")
      .def_property_readonly("converged", [](const vqe::adapt_result &result)
                             { return result.converged; }, "Whether ADAPT procedure converged")
      .def_property_readonly("success", [](const vqe::adapt_result &result)
                             { return result.converged; }, "Alias for converged")
      .def_property_readonly("circuit_depth", [](const vqe::adapt_result &result)
                             { return result.selected_indices.size(); }, "Final circuit depth (number of selected operators)")
      .def("__repr__", [](const vqe::adapt_result &result)
           {
        std::ostringstream oss;
        oss << "AdaptResult(energy=" << std::fixed << std::setprecision(8) << result.energy
            << ", iterations=" << result.iterations
            << ", depth=" << result.selected_indices.size()
            << ", converged=" << (result.converged ? "True" : "False") << ")";
        return oss.str(); });

  // ============================================================================
  // Energy evaluator with thread-safe interface
  // ============================================================================

  py::class_<EnergyEvaluator>(m, "EnergyEvaluator", "Thread-safe energy evaluator for VQE")
      .def(py::init<const vqe::uccsd_ansatz &, const std::vector<vqe::pauli_term> &, bool>(),
           py::arg("ansatz"), py::arg("pauli_terms"), py::arg("use_gpu") = false,
           "Create energy evaluator")
      .def_property_readonly("parameter_count", &EnergyEvaluator::parameter_count,
                             "Number of variational parameters")
      .def_property_readonly("num_parameters", &EnergyEvaluator::parameter_count,
                             "Alias for parameter_count")
      .def_property("parameters",
                    &EnergyEvaluator::parameters,
                    &EnergyEvaluator::set_parameters,
                    "Current parameter values")
      .def("energy", &EnergyEvaluator::energy, py::arg("parameters"),
           "Evaluate energy for given parameters")
      .def("__call__", &EnergyEvaluator::energy, py::arg("parameters"),
           "Make evaluator callable")
      .def("__repr__", [](const EnergyEvaluator &eval)
           {
        std::ostringstream oss;
        oss << "EnergyEvaluator(parameters=" << eval.parameter_count() << ")";
        return oss.str(); });

  // ============================================================================
  // Module-level functions with improved error handling and documentation
  // ============================================================================

  m.def("load", [](const std::string &path)
        {
        if (path.empty()) {
          throw std::invalid_argument("Hamiltonian file path cannot be empty");
        }
        py::gil_scoped_release release;
        try {
          return vqe::read_hamiltonian_file(path);
        } catch (const std::exception& e) {
          throw std::runtime_error("Failed to read Hamiltonian from '" + path + "': " + e.what());
        } }, py::arg("filename"), "Read fermionic Hamiltonian from file\n\n"
                                  "Parameters\n"
                                  "----------\n"
                                  "filename : str\n"
                                  "    Path to Hamiltonian file\n\n"
                                  "Returns\n"
                                  "-------\n"
                                  "HamiltonianData\n"
                                  "    Parsed Hamiltonian data structure");

  m.def("transform", [](const vqe::hamiltonian_data &data)
        {
        py::gil_scoped_release release;
        try {
          return vqe::jordan_wigner_transform(data);
        } catch (const std::exception& e) {
          throw std::runtime_error("Jordan-Wigner transformation failed: " + std::string(e.what()));
        } }, py::arg("hamiltonian"), "Apply Jordan-Wigner transformation to fermionic Hamiltonian\n\n"
                                     "Parameters\n"
                                     "----------\n"
                                     "hamiltonian : HamiltonianData\n"
                                     "    Fermionic Hamiltonian to transform\n\n"
                                     "Returns\n"
                                     "-------\n"
                                     "list[PauliTerm]\n"
                                     "    Qubit Hamiltonian as list of Pauli terms");

  m.def("pauli_to_string", &vqe::pauli_to_string,
        py::arg("x_mask"), py::arg("z_mask"), py::arg("n_qubits"),
        "Convert Pauli operator masks to string representation\n\n"
        "Parameters\n"
        "----------\n"
        "x_mask : int\n"
        "    Bit mask for X operators\n"
        "z_mask : int\n"
        "    Bit mask for Z operators\n"
        "n_qubits : int\n"
        "    Number of qubits\n\n"
        "Returns\n"
        "-------\n"
        "str\n"
        "    Pauli string (e.g., 'XIZY')");

  m.def("run_from_file", [](const std::string &path, std::size_t n_particles, const vqe::vqe_options &options)
        {
        if (path.empty()) {
          throw std::invalid_argument("Hamiltonian file path cannot be empty");
        }
        if (n_particles == 0) {
          throw std::invalid_argument("Number of particles must be positive");
        }
        py::gil_scoped_release release;
        try {
          return vqe::run_default_vqe(path, n_particles, options);
        } catch (const std::exception& e) {
          throw std::runtime_error("VQE calculation failed: " + std::string(e.what()));
        } }, py::arg("hamiltonian_file"), py::arg("n_electrons"), py::arg("options") = vqe::vqe_options{}, "Run VQE calculation from Hamiltonian file\n\n"
                                                                                                                                                "Parameters\n"
                                                                                                                                                "----------\n"
                                                                                                                                                "hamiltonian_file : str\n"
                                                                                                                                                "    Path to fermionic Hamiltonian file\n"
                                                                                                                                                "n_electrons : int\n"
                                                                                                                                                "    Number of electrons in the system\n"
                                                                                                                                                "options : Options, optional\n"
                                                                                                                                                "    VQE configuration options\n\n"
                                                                                                                                                "Returns\n"
                                                                                                                                                "-------\n"
                                                                                                                                                "VQEResult\n"
                                                                                                                                                "    Optimization results");

  m.def("run", [](vqe::uccsd_ansatz ansatz, const std::vector<vqe::pauli_term> &pauli_terms, const vqe::vqe_options &options)
        {
        if (pauli_terms.empty()) {
          throw std::invalid_argument("Pauli terms list cannot be empty");
        }
        py::gil_scoped_release release;
        try {
          auto opts = options;
          opts.mode = vqe::run_mode::standard;
          return vqe::run_default_vqe_with_ansatz(ansatz, pauli_terms, opts);
        } catch (const std::exception& e) {
          throw std::runtime_error("VQE calculation failed: " + std::string(e.what()));
        } }, py::arg("ansatz"), py::arg("hamiltonian"), py::arg("options") = vqe::vqe_options{}, "Run VQE with pre-constructed ansatz and Hamiltonian\n\n"
                                                                                                 "Parameters\n"
                                                                                                 "----------\n"
                                                                                                 "ansatz : UCCSD\n"
                                                                                                 "    Prepared UCCSD ansatz\n"
                                                                                                 "hamiltonian : list[PauliTerm]\n"
                                                                                                 "    Qubit Hamiltonian as Pauli terms\n"
                                                                                                 "options : Options, optional\n"
                                                                                                 "    VQE configuration options\n\n"
                                                                                                 "Returns\n"
                                                                                                 "-------\n"
                                                                                                 "VQEResult\n"
                                                                                                 "    Optimization results");

  m.def("adapt", [](const std::string &path, std::size_t n_particles, const vqe::vqe_options &options)
        {
        if (path.empty()) {
          throw std::invalid_argument("Hamiltonian file path cannot be empty");
        }
        if (n_particles == 0) {
          throw std::invalid_argument("Number of particles must be positive");
        }
        py::gil_scoped_release release;
        try {
          auto opts = options;
          opts.mode = vqe::run_mode::adapt;
          return vqe::run_adapt_vqe(path, n_particles, opts);
        } catch (const std::exception& e) {
          throw std::runtime_error("ADAPT-VQE calculation failed: " + std::string(e.what()));
        } }, py::arg("hamiltonian_file"), py::arg("n_electrons"), py::arg("options") = vqe::vqe_options{}, "Run ADAPT-VQE calculation from Hamiltonian file\n\n"
                                                                                                                                                  "Parameters\n"
                                                                                                                                                  "----------\n"
                                                                                                                                                  "hamiltonian_file : str\n"
                                                                                                                                                  "    Path to fermionic Hamiltonian file\n"
                                                                                                                                                  "n_electrons : int\n"
                                                                                                                                                  "    Number of electrons in the system\n"
                                                                                                                                                  "options : Options, optional\n"
                                                                                                                                                  "    Configuration options (mode will be forced to 'adapt')\n\n"
                                                                                                                                                  "Returns\n"
                                                                                                                                                  "-------\n"
                                                                                                                                                  "AdaptResult\n"
                                                                                                                                                  "    ADAPT optimization results");

  m.def("energy", [](vqe::uccsd_ansatz ansatz, const std::vector<vqe::pauli_term> &terms, const std::vector<double> &parameters, bool use_gpu)
        {
        if (terms.empty()) {
          throw std::invalid_argument("Pauli terms list cannot be empty");
        }
        if (parameters.size() != ansatz.get_circuit().parameters().size()) {
          throw std::invalid_argument("Parameter vector length mismatch");
        }
        py::gil_scoped_release release;
        try {
          return evaluate_energy_internal(ansatz, terms, parameters, use_gpu);
        } catch (const std::exception& e) {
          throw std::runtime_error("Energy evaluation failed: " + std::string(e.what()));
        } }, py::arg("ansatz"), py::arg("hamiltonian"), py::arg("parameters"), py::arg("use_gpu") = false, "Evaluate energy for given ansatz and parameters\n\n"
                                                                                                           "Parameters\n"
                                                                                                           "----------\n"
                                                                                                           "ansatz : UCCSD\n"
                                                                                                           "    Prepared UCCSD ansatz\n"
                                                                                                           "hamiltonian : list[PauliTerm]\n"
                                                                                                           "    Qubit Hamiltonian as Pauli terms\n"
                                                                                                           "parameters : list[float]\n"
                                                                                                           "    Variational parameters\n"
                                                                                                           "use_gpu : bool, optional\n"
                                                                                                           "    Use GPU acceleration if available (default: False)\n\n"
                                                                                                           "Returns\n"
                                                                                                           "-------\n"
                                                                                                           "float\n"
                                                                                                           "    Expectation value of the Hamiltonian");

  // ============================================================================
  // Utility functions
  // ============================================================================

  m.def("list_optimizers", []()
        {
        std::vector<std::string> optimizers;
        for (int i = 0; i < nlopt::NUM_ALGORITHMS; ++i) {
          const char* name = nlopt::algorithm_name(static_cast<nlopt::algorithm>(i));
          if (name && strlen(name) > 0) {
            optimizers.emplace_back(name);
          }
        }
        return optimizers; }, "Get list of available optimization algorithms\n\n"
             "Returns\n"
             "-------\n"
             "list[str]\n"
             "    Names of available NLopt algorithms");

  m.def("check_gpu_support", []()
        {
#if defined(VQE_ENABLE_CUDA) || defined(VQE_ENABLE_HIP)
          return true;
#else
          return false;
#endif
        },
        "Check if GPU support is available\n\n"
        "Returns\n"
        "-------\n"
        "bool\n"
        "    True if compiled with CUDA or HIP support");

  // Version and build information
  m.attr("__version__") = "1.0.0";
  m.attr("__build_info__") = py::dict(
      py::arg("cuda_support") =
#if defined(VQE_ENABLE_CUDA) || defined(VQE_ENABLE_HIP)
          true
#else
          false
#endif
  );
}

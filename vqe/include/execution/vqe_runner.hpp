#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "vqe_options.hpp"

namespace vqe
{

  class uccsd_ansatz;
  struct pauli_term;

  struct vqe_result
  {
    double energy = 0.0;
    std::vector<double> parameters;
    std::size_t evaluations = 0;
    bool converged = false;
    std::size_t iterations = 0;
    double initial_energy = 0.0;
    double energy_delta = 0.0;
    double apply_time = 0.0;
    double expectation_time = 0.0;
    nlopt::result status = nlopt::FAILURE;
  };

  vqe_result run_default_vqe(const std::string &hamiltonian_path,
                             std::size_t n_particles,
                             const vqe_options &options = {});

  vqe_result run_default_vqe_with_ansatz(uccsd_ansatz &ansatz,
                                         const std::vector<pauli_term> &pauli_terms,
                                         const vqe_options &options);

} // namespace vqe

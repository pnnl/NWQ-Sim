#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include "execution/vqe_runner.hpp"

namespace vqe
{

  struct adapt_result
  {
    double energy = 0.0;
    std::vector<double> parameters;
    std::vector<std::size_t> selected_indices;
    std::vector<std::string> selected_labels;
    std::size_t iterations = 0;
    std::size_t energy_evaluations = 0;
    bool converged = false;
    double initial_energy = 0.0;
  };

  adapt_result run_adapt_vqe(const std::string &hamiltonian_path,
                             std::size_t n_particles,
                             const vqe_options &options = {});

} // namespace vqe

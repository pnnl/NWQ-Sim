#pragma once

#include <complex>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "execution/vqe_runner.hpp"

std::vector<std::pair<std::string, std::complex<double>>> parseHamiltonianFile(const std::string &filename);

std::pair<double, std::vector<std::pair<std::vector<int>, double>>> qflow_nwqsim(
    const std::vector<std::pair<std::string, std::complex<double>>> &hamiltonian_ops,
    int n_part,
    std::string backend = "CPU",
    std::optional<vqe::vqe_options> options = std::nullopt);

// std::string get_termination_reason_local(int result); // MZ: "local gradient" (SPSA gradient) is not implemnted after re-organized anyway


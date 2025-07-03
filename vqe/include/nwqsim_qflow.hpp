#ifndef NWQSIM_QFLOW_HPP
#define NWQSIM_QFLOW_HPP

#include <vector>
#include <string>
#include <complex>

/**
 * @brief Parse a Hamiltonian file and extract operator-coefficient pairs
 *
 * @param filename Path to the Hamiltonian file
 * @return Vector of pairs containing operator strings and their complex coefficients
 */
std::vector<std::pair<std::string, std::complex<double>>> parseHamiltonianFile(const std::string &filename);

/**
 * @brief Run VQE optimization using NWQSim backend
 *
 * @param hamiltonian_ops Vector of Hamiltonian operators and coefficients
 * @param backend Backend type (default: "CPU")
 * @return Optimized ground state energy
 */
double qflow_nwqsim(const std::vector<std::pair<std::string, std::complex<double>>> &hamiltonian_ops, int n_part, std::string backend = "CPU");

/**
 * @brief Get termination reason for optimization result
 *
 * @param result Optimization result code
 * @return Human-readable termination reason
 */
std::string get_termination_reason_local(int result);

#endif // EXACHEM_QFLOW_HPP

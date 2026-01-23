#include "hamiltonian_parser.hpp"

#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace vqe {
namespace {

std::string trim(const std::string& input) {
  std::size_t first = 0;
  while (first < input.size() && std::isspace(static_cast<unsigned char>(input[first]))) {
    ++first;
  }
  std::size_t last = input.size();
  while (last > first && std::isspace(static_cast<unsigned char>(input[last - 1]))) {
    --last;
  }
  return input.substr(first, last - first);
}

}  // namespace

hamiltonian_data read_hamiltonian_file(const std::string& path) {
  std::ifstream file(path);
  if (!file.is_open()) {
    throw std::runtime_error("failed to open hamiltonian file: " + path);
  }

  hamiltonian_data data{};
  std::string line;
  std::size_t line_number = 0;

  while (std::getline(file, line)) {
    ++line_number;
    const std::string trimmed = trim(line);
    if (trimmed.empty() || trimmed[0] == 'c' || trimmed[0] == '#') {
      continue;
    }

    const auto closing = trimmed.find(')');
    if (trimmed.size() < 2 || trimmed.front() != '(' || closing == std::string::npos) {
      throw std::runtime_error("malformed coefficient on line " + std::to_string(line_number));
    }

    const std::string coeff_text = trimmed.substr(1, closing - 1);
    std::istringstream coeff_stream(coeff_text);
    std::string real_str;
    std::string imag_str;
    if (!std::getline(coeff_stream, real_str, ',')) {
      throw std::runtime_error("missing real coefficient on line " + std::to_string(line_number));
    }
    if (!std::getline(coeff_stream, imag_str)) {
      throw std::runtime_error("missing imaginary coefficient on line " + std::to_string(line_number));
    }

    const double real = std::stod(trim(real_str));
    const double imag = std::stod(trim(imag_str));
    const std::complex<double> coeff(real, imag);

    const std::string ops_part = trim(trimmed.substr(closing + 1));
    if (ops_part.empty()) {
      data.constant += coeff;
      continue;
    }

    fermion_term term;
    term.coefficient = coeff;
    std::istringstream ops_stream(ops_part);
    std::string token;
    while (ops_stream >> token) {
      if (token == "+") {
        continue;
      }
      const bool dagger = !token.empty() && token.back() == '^';
      const std::string index_part = dagger ? token.substr(0, token.size() - 1) : token;
      if (index_part.empty()) {
        throw std::runtime_error("missing orbital index on line " + std::to_string(line_number));
      }
      const std::size_t index = static_cast<std::size_t>(std::stoll(index_part));
      data.max_index = std::max(data.max_index, index);
      term.operators.push_back({index, dagger ? fermion_op_kind::creation : fermion_op_kind::annihilation});
    }

    if (term.operators.empty()) {
      data.constant += coeff;
    } else {
      data.terms.push_back(std::move(term));
    }
  }

  if (data.max_index == 0) {
    bool has_non_identity = false;
    for (const auto& term : data.terms) {
      for (const auto& op : term.operators) {
        has_non_identity = true;
        data.max_index = std::max(data.max_index, op.index);
      }
    }
    if (!has_non_identity) {
      data.max_index = 0;
    }
  }

  return data;
}

}  // namespace vqe

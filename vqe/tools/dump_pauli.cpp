// dump_pauli: write NWQ-Sim's Jordan-Wigner Pauli representation of an XACC
// fermionic Hamiltonian to disk.
//
// All scientific computation (file parsing, Jordan-Wigner transform, term
// accumulation/normalization) is performed by the NWQ-Sim vqe frontend sources
// this tool links against (hamiltonian_parser.cpp, jw_transform.cpp,
// pauli_term.cpp). This file only formats the resulting pauli_term list and,
// optionally, evaluates the diagonal expectation value <ref|H|ref> of a
// computational-basis reference state as a validation sum.
//
// Usage:
//   dump_pauli <input-xacc> <output-basename> <cutoff|none> [occ_qubits_csv]
//
// Outputs (data only, no headers; conventions are documented separately):
//   <output-basename>.txt   sparse text, one term per line: "coeff_re coeff_im [X0 Z3 ...]"
//   <output-basename>.json  {"n_qubits": N, "terms": [[label, [re, im]], ...]}
//
// cutoff "none" keeps every term whose coefficient is not exactly zero
// (cutoff = DBL_TRUE_MIN); a numeric cutoff is forwarded verbatim to
// NWQ-Sim's normalize_terms via jordan_wigner_transform.
//
// Build (not wired into CMake; -ffp-contract=off keeps results bit-identical
// across compilers that would otherwise fuse multiply-adds):
//   clang++ -std=c++17 -O2 -ffp-contract=off -Ivqe/include \
//     vqe/src/hamiltonian_parser.cpp vqe/src/jw_transform.cpp \
//     vqe/src/pauli_term.cpp vqe/tools/dump_pauli.cpp -o dump_pauli

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include "hamiltonian_parser.hpp"
#include "jw_transform.hpp"
#include "pauli_term.hpp"

namespace {

std::string fmt17(double v) {
  char buf[64];
  std::snprintf(buf, sizeof(buf), "%.17g", v);
  return buf;
}

// Sparse operator list for one term, e.g. "X0 Z3 Y12". Symbols are taken from
// NWQ-Sim's pauli_to_string so this tool never re-derives the X/Y/Z encoding:
// in that label, label[n_qubits-1-q] is the symbol acting on qubit q.
std::string sparse_ops(const std::string& label, std::size_t n_qubits) {
  std::string out;
  for (std::size_t q = 0; q < n_qubits; ++q) {
    const char sym = label[n_qubits - 1 - q];
    if (sym == 'I') {
      continue;
    }
    if (!out.empty()) {
      out.push_back(' ');
    }
    out.push_back(sym);
    out += std::to_string(q);
  }
  return out;
}

std::string basename_of(const std::string& path) {
  const auto pos = path.find_last_of('/');
  return pos == std::string::npos ? path : path.substr(pos + 1);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 4) {
    std::cerr << "usage: dump_pauli <input-xacc> <output-basename> <cutoff|none> [occ_qubits_csv]"
              << std::endl;
    return 1;
  }
  const std::string input_path = argv[1];
  const std::string out_base = argv[2];
  const std::string cutoff_arg = argv[3];

  double cutoff;
  std::string cutoff_desc;
  if (cutoff_arg == "none") {
    cutoff = std::numeric_limits<double>::denorm_min();
    cutoff_desc = "none (only exactly-zero coefficients dropped)";
  } else {
    cutoff = std::stod(cutoff_arg);
    cutoff_desc = fmt17(cutoff);
  }

  std::vector<std::size_t> occ_qubits;
  if (argc > 4) {
    std::stringstream ss(argv[4]);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
      if (!tok.empty()) {
        occ_qubits.push_back(static_cast<std::size_t>(std::stoull(tok)));
      }
    }
  }

  // --- Scientific pipeline: NWQ-Sim code only ---
  const auto data = vqe::read_hamiltonian_file(input_path);
  const std::size_t n_qubits = data.num_qubits();
  const auto terms = vqe::jordan_wigner_transform(data, cutoff);
  // ---------------------------------------------

  double max_abs_imag = 0.0;
  std::complex<double> identity_coeff{0.0, 0.0};
  for (const auto& t : terms) {
    max_abs_imag = std::max(max_abs_imag, std::abs(t.coefficient.imag()));
    if (t.x_mask == 0 && t.z_mask == 0) {
      identity_coeff = t.coefficient;
    }
  }

  // Diagonal validation sum <ref|H|ref> for a computational-basis reference
  // state with the given qubits occupied: only x_mask == 0 terms contribute,
  // each weighted by the parity of z_mask restricted to occupied qubits.
  std::complex<double> ref_energy{0.0, 0.0};
  std::uint64_t occ_mask = 0;
  if (!occ_qubits.empty()) {
    for (const auto q : occ_qubits) {
      occ_mask |= (static_cast<std::uint64_t>(1) << q);
    }
    for (const auto& t : terms) {
      if (t.x_mask != 0) {
        continue;
      }
      const int parity = __builtin_popcountll(t.z_mask & occ_mask) % 2;
      ref_energy += (parity ? -1.0 : 1.0) * t.coefficient;
    }
  }

  const std::string src_name = basename_of(input_path);

  // ---- sparse text output (data only) ----
  {
    std::ofstream txt(out_base + ".txt");
    if (!txt.is_open()) {
      std::cerr << "failed to open " << out_base << ".txt" << std::endl;
      return 1;
    }
    for (const auto& t : terms) {
      const std::string label = vqe::pauli_to_string(t.x_mask, t.z_mask, n_qubits);
      txt << fmt17(t.coefficient.real()) << ' ' << fmt17(t.coefficient.imag()) << " ["
          << sparse_ops(label, n_qubits) << "]\n";
    }
  }

  // ---- Qiskit-compatible JSON output ----
  {
    std::ofstream js(out_base + ".json");
    if (!js.is_open()) {
      std::cerr << "failed to open " << out_base << ".json" << std::endl;
      return 1;
    }
    js << "{\n";
    js << "  \"n_qubits\": " << n_qubits << ",\n";
    js << "  \"terms\": [\n";
    for (std::size_t i = 0; i < terms.size(); ++i) {
      const auto& t = terms[i];
      js << "    [\"" << vqe::pauli_to_string(t.x_mask, t.z_mask, n_qubits) << "\", ["
         << fmt17(t.coefficient.real()) << ", " << fmt17(t.coefficient.imag()) << "]]"
         << (i + 1 < terms.size() ? ",\n" : "\n");
    }
    js << "  ]\n";
    js << "}\n";
  }

  // ---- stats to stdout ----
  std::cout << "input            : " << src_name << "\n";
  std::cout << "n_qubits         : " << n_qubits << "\n";
  std::cout << "fermionic_terms  : " << data.terms.size() << "\n";
  std::cout << "constant_in_file : " << fmt17(data.constant.real()) << " " << fmt17(data.constant.imag()) << "\n";
  std::cout << "cutoff           : " << cutoff_desc << "\n";
  std::cout << "pauli_terms      : " << terms.size() << "\n";
  std::cout << "identity_coeff   : " << fmt17(identity_coeff.real()) << " " << fmt17(identity_coeff.imag()) << "\n";
  std::cout << "max_abs_imag     : " << fmt17(max_abs_imag) << "\n";
  if (!occ_qubits.empty()) {
    std::cout << "ref_occ_qubits   : " << argv[4] << "\n";
    std::cout << "ref_energy       : " << fmt17(ref_energy.real()) << " " << fmt17(ref_energy.imag()) << "\n";
  }
  return 0;
}

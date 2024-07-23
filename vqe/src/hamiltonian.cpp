#include "observable/hamiltonian.hpp"
#include <fstream>
#include <sstream>
#include <regex>
#include <queue>
#include <list>
namespace NWQSim{
  namespace VQE {
    struct FermiArgs {
      IdxType qubit;
      std::complex<ValType> coeff;
      FermionOpType typeval;
      IdxType term_index;
    };
    const std::regex pattern("\\(\\s*([\\d\\.e\\+-]+),\\s*([\\d\\.]+)\\)([\\d^\\s]+)");
    Hamiltonian::Hamiltonian(std::string input_path, IdxType n_particles, 
                    Transformer transform) {
      std::ifstream input_file(input_path);
      if (!input_file.is_open()) {
        throw std::runtime_error("Could not open file");
      }
      std::string line;
      std::istringstream linestream;

      std::string throwaway;
      std::complex<ValType> coeff;
      std::smatch match;
      IdxType max_index = 0;
        
      std::queue<FermiArgs> arglist;
      IdxType term_counter = 0;
      ValType energy_constant = 0;
      while (std::getline(input_file, line)) {
        if (line.length() == 0 || line[0] == 'c') {
          continue;
        }
        
        if (std::regex_search(line, match, pattern)) {
          if (match.size() != 4) {
            continue;
          }
          coeff = {std::stod(match.str(1)), std::stod(match.str(2))};
          if (match.str(3) == "") {
            energy_constant = coeff.real();
          }
          std::string fermi_ops = match.str(3);
          linestream = std::istringstream(fermi_ops);
          std::vector<std::string> term_ops;
          while(std::getline(linestream, throwaway, ' ')) {
            FermionOpType optype;
            IdxType qubit_index;
            if (throwaway == "" || throwaway == "+")  {
              continue;
            }
            if (throwaway.back() == '^') {
              optype = Creation;
              qubit_index = std::stoll(throwaway.substr(0, throwaway.length() - 1));
            } else {
              optype = Annihilation;
              qubit_index = std::stoll(throwaway);
            }
            max_index = std::max(max_index, qubit_index);

            term_ops.push_back(throwaway);
            FermiArgs args;
            args.coeff = coeff;
            args.qubit = qubit_index;
            args.term_index = term_counter;
            args.typeval = optype;
            arglist.push(args);
          }
          term_counter ++;
        } else {

        std::cout << "No match " <<  line << std::endl;
        }
      }
      if ((max_index + 1) % 2) {
        throw std::runtime_error("Odd number of spin orbitals");  
      }
        
      fermi_operators.resize(term_counter);
      IdxType n_spatial_orbitals = (max_index + 1) / 2;
      IdxType n_occ = n_particles / 2;
      IdxType n_virt = n_spatial_orbitals - n_occ;
      env = MolecularEnvironment(n_spatial_orbitals, 
                                 n_particles,
                                 xacc_scheme,
                                 energy_constant);
      IdxType coeff_index = 0;
      IdxType otypeval, spinval, orbital_index;
      Spin spin;
      OrbitalType orbitaltype;
      
      while(!arglist.empty()) {
        auto args = arglist.front();
        arglist.pop();
        
        getFermiInfoFromQubit(args.qubit, orbital_index, spinval, otypeval, n_occ, n_virt, xacc_scheme);
        spin = spinval ? Down : Up;
        orbitaltype = otypeval ? Virtual : Occupied;
        fermi_operators[args.term_index].push_back(
          FermionOperator(orbital_index, 
                          orbitaltype, 
                          spin, 
                          args.typeval, 
                          xacc_scheme,
                          args.coeff));
      };
      transform(env, fermi_operators, pauli_operators, false);
      n_ops = 0;
      for (auto& i : pauli_operators) {
        n_ops += i.size();
      }
// #ifndef NDEBUG
//       std::cout << "Hamiltonian Pauli Strings:" << std::endl;
//       for (auto& plist: pauli_operators) {
//         for (auto& p: plist) {
//           std::cout << p << "\n";
//         }
//       }
// #endif
    };
    Hamiltonian::Hamiltonian(const std::vector<std::pair<std::string, std::complex<double>>>& input_ops, 
                             IdxType n_particles, 
                             bool xacc_scheme,
                             Transformer transform) {


      std::string throwaway;
      std::complex<ValType> coeff;
      IdxType max_index = 0;
      std::queue<FermiArgs> arglist;
      std::istringstream linestream;
      IdxType term_counter = 0;
      ValType energy_constant = 0;
      for (auto& pair: input_ops) {
        

          coeff = pair.second;
          if (pair.first == "") {
            energy_constant = coeff.real();
          }
          std::string fermi_ops = pair.first;
          linestream = std::istringstream(fermi_ops);
          std::vector<std::string> term_ops;
          while(std::getline(linestream, throwaway, ' ')) {
            FermionOpType optype;
            IdxType qubit_index;
            if (throwaway == "" || throwaway == "+")  {
              continue;
            }
            if (throwaway.back() == '^') {
              optype = Creation;
              qubit_index = std::stoll(throwaway.substr(0, throwaway.length() - 1));
            } else {
              optype = Annihilation;
              qubit_index = std::stoll(throwaway);
            }
            max_index = std::max(max_index, qubit_index);

            term_ops.push_back(throwaway);
            FermiArgs args;
            args.coeff = coeff;
            args.qubit = qubit_index;
            args.term_index = term_counter;
            args.typeval = optype;
            arglist.push(args);
          }
          term_counter ++;
      }
      if ((max_index + 1) % 2) {
        throw std::runtime_error("Odd number of spin orbitals");  
      }
      fermi_operators.resize(term_counter);
      IdxType n_spatial_orbitals = (max_index + 1) / 2;
      IdxType n_occ = n_particles / 2;
      IdxType n_virt = n_spatial_orbitals - n_occ;
      env = MolecularEnvironment(n_spatial_orbitals, 
                                 n_particles,
                                 xacc_scheme,
                                 energy_constant);
      IdxType coeff_index = 0;
      IdxType otypeval, spinval, orbital_index;
      Spin spin;
      OrbitalType orbitaltype;
      
      while(!arglist.empty()) {
        auto args = arglist.front();
        arglist.pop();
        
        getFermiInfoFromQubit(args.qubit, orbital_index, spinval, otypeval, n_occ, n_virt, xacc_scheme);
        spin = spinval ? Down : Up;
        orbitaltype = otypeval ? Virtual : Occupied;
        fermi_operators[args.term_index].push_back(
          FermionOperator(orbital_index, 
                          orbitaltype, 
                          spin, 
                          args.typeval, 
                          xacc_scheme,
                          args.coeff));
      };
      transform(env, fermi_operators, pauli_operators, false);
      n_ops = 0;
      for (auto& i : pauli_operators) {
        n_ops += i.size();
      }
// #ifndef NDEBUG
//       std::cout << "Hamiltonian Pauli Strings:" << std::endl;
//       for (auto& plist: pauli_operators) {
//         for (auto& p: plist) {
//           std::cout << p << "\n";
//         }
//       }
// #endif
    };
  }; // namespace VQE
}; // namespace NWQSim
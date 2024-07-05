#include "circuit/ansatz.hpp"
#ifndef DYNAMIC_ANSATZ_HH
#define DYNAMIC_ANSATZ_HH
namespace NWQSim { 
  namespace VQE {
    using PauliOplist = std::vector<FermionOperator> >;
    using SymmetryGroup = std::vector<FermionOplist>; // for now, just assume we're preserving spin-conjugated symmetries

    class DynamicAnsatz: public Ansatz {
      protected:
        MolecularEnvironment& env;
        std::vector<SymmetryGroup> operator_pool; // list of Fermion operators, grouped according to spin symmetry
        std::vector<std::pair<IdxType, ValType> > symmetry_terms; // list of Fermion operators, grouped according to spin symmetry
        Transformer transformer;
        IdxType last_ptr;
        virtual void generate_operator_pool();
        virtual void add_operator(std::vector<std::vector<PauliOperator> > ops, ValType param) {
          theta->push_back(param);
          last_ptr = gates->size();
          for (std::vector<PauliOperator> oplist: ops) {
            std::vector<std::pair<IdxType, ValType > > idxvals = {{theta->size() - 1, 1.0}};
            for (PauliOperator pauli: oplist) {
              if (pauli.isNonTrivial() && abs(pauli.getCoeff().real()) > 1e-10) {
                ExponentialGate(pauli, OP::RZ, idxvals, 2 * pauli.getCoeff().real());
              }
            }
          }
        }; 
        virtual void remove_last_operator() {
          theta->resize(theta->size() - 1);
          last_ptr = gates->resize(last_ptr);
        }; 

    };
  };
};

#endif
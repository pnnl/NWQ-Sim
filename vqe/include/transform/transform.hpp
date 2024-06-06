#ifndef __JW_HH__
#define __JW_HH__
#include "observable/fermionic_operator.hpp"
#include "observable/pauli_operator.hpp"
#include "environment.hpp"


namespace NWQSim {
  namespace VQE {
    typedef void (*Transformer) (const MolecularEnvironment&, const std::vector<std::vector<FermionOperator> >&, std::vector<std::vector<PauliOperator> >& ,bool);
    /**
     * Common header for Fermion op to qubit (Pauli string) tranformations
     * Accept as input a list of Fermion product operators (e.g. $$ a_i^\dagger a_j^\dagger a_ra_s $$ )
    */

    // Common typedef in case we want to do some sneaky function pointer Polymorphism
    // typedef std::shared_ptr<std::vector<PauliOperator> > (*Transform) (const std::vector<FermionOperator>&, double);
    void getJordanWignerTransform(const MolecularEnvironment& env, 
                                  const std::vector<std::vector<FermionOperator> >& fermi_operators,
                                  std::vector<std::vector<PauliOperator> >& pauli_operators, 
                                  bool hermitian = false);
    void getParityTransform(const MolecularEnvironment& env, 
                                  const std::vector<std::vector<FermionOperator> >& fermi_operators,
                                  std::vector<std::vector<PauliOperator> >& pauli_operators, 
                                  bool hermitian = false);
    void getBravyiKitaevTransform(const MolecularEnvironment& env, 
                                  const std::vector<std::vector<FermionOperator> >& fermi_operators,
                                  std::vector<std::vector<PauliOperator> >& pauli_operators, 
                                  bool hermitian = false);
  }; // namespace VQE
};// namespace NWQSim
#endif
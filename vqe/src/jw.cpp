#include <iostream>
#include "transform/transform.hpp"
#include "utils.hpp"

namespace NWQSim {
  namespace VQE {
    void jwFermionToPauliSingle (
    IdxType n_occ,
    IdxType n_virt,
    IdxType n_qubits,
    FermionOperator ferm_op,
    std::vector<PauliOperator>& output) {
  size_t qubit_index = ferm_op.qubitIndex(n_occ, n_virt);
  IdxType xmask = (1 << qubit_index);
  IdxType zmask1 = (1 << qubit_index) - 1;
  IdxType zmask2 = (1 << (qubit_index + 1)) - 1;
  // we also want to subtract a, adagger
  output.push_back(PauliOperator(xmask, zmask1, n_qubits, 0.5));
  int sign = (ferm_op.getType() == Annihilation) - (ferm_op.getType() == Creation);
  output.push_back(PauliOperator(xmask, zmask2, n_qubits, std::complex(0.0, 0.5 * sign)));
}

void jwFermionToPauliSinglePair (
    size_t n_occ,
    size_t n_virt,
    const FermionOperator& ap_dagger,
    const FermionOperator& aq,
    std::vector<PauliOperator>& output,
    bool hermitian) {
  // assert(ap_dagger.getType() == Creation);
  // assert(aq.getType() == Annihilation);
  std::vector<PauliOperator> ap_dag_paulis;
  jwFermionToPauliSingle(n_occ, n_virt, 2 * (n_occ + n_virt), ap_dagger, ap_dag_paulis);
  std::vector<PauliOperator> aq_paulis;
  jwFermionToPauliSingle(n_occ, n_virt, 2 * (n_occ + n_virt),  aq, aq_paulis);
  std::complex<ValType> fermicoeff = ap_dagger.getCoeff();
  for (size_t i = 0; i < 4; i ++) {
    int i1 = (i & (1 << 0)) >> 0;
    int i2 = (i & (1 << 1)) >> 1;
    PauliOperator p_op = ap_dag_paulis[i2] * 
                         aq_paulis[i1];
    std::complex<ValType> p_coeff = fermicoeff * p_op.getCoeff();
    if (hermitian) {
      p_op.setCoeff(std::complex(0.0,  2 * p_coeff.imag()) * std::complex(0.0, -1.0));
    } else {
      p_op.setCoeff(p_coeff);
    }
    output.push_back(p_op);
  }
}

void jwFermionToPauliDoublePair (
    size_t n_occ,
    size_t n_virt,
    FermionOperator ap_dagger,
    FermionOperator aq_dagger,
    FermionOperator ar,
    FermionOperator as,
    std::vector<PauliOperator>& output,
    bool hermitian = false) {
  // assert(ap_dagger.getType() == Creation);
  // assert(aq_dagger.getType() == Creation);
  // assert(ar.getType() == Annihilation);
  // assert(as.getType() == Annihilation);
  std::vector<PauliOperator> ap_dag_paulis;
  jwFermionToPauliSingle(n_occ, n_virt,2 * (n_occ + n_virt),  ap_dagger, ap_dag_paulis);
  std::vector<PauliOperator> aq_dag_paulis;
  jwFermionToPauliSingle(n_occ, n_virt,2 * (n_occ + n_virt),  aq_dagger, aq_dag_paulis);
  std::vector<PauliOperator> ar_paulis;
  jwFermionToPauliSingle(n_occ, n_virt,2 * (n_occ + n_virt),  ar, ar_paulis);
  std::vector<PauliOperator> as_paulis;
  jwFermionToPauliSingle(n_occ, n_virt,2 * (n_occ + n_virt),  as, as_paulis);
  std::complex<ValType> fermicoeff = ap_dagger.getCoeff();
  // std::cout << fermicoeff << std::endl;
  for (int i = 15; i >= 0; i --) {
    int i1 = (i & (1 << 0)) >> 0;
    int i2 = (i & (1 << 1)) >> 1;
    int i3 = (i & (1 << 2)) >> 2;
    int i4 = (i & (1 << 3)) >> 3;

    PauliOperator p_op = ap_dag_paulis[i1] *
                         aq_dag_paulis[i2] *
                         ar_paulis[i3] *
                         as_paulis[i4];
    std::complex<ValType> p_coeff = fermicoeff * p_op.getCoeff();
    if (hermitian) {
      p_op.setCoeff(std::complex(0.0, 2 * p_coeff.imag()) * std::complex(0.0, -1.0));
    } else {
      p_op.setCoeff(p_coeff);
    }
    output.push_back(p_op);
  }
}

// used for single-
void getJordanWignerTransform(
    const MolecularEnvironment& env,
    const std::vector<std::vector<FermionOperator> >& input,
    std::vector<PauliOperator >& output,
    bool hermitian
    ) {
  for (auto& fermion_product: input) {
    if (fermion_product.size() == 2) {
      // continue;
      jwFermionToPauliSinglePair(
        env.n_occ,
        env.n_virt,
        fermion_product[0], 
        fermion_product[1],
        output,
        hermitian);
    } else if (fermion_product.size() == 4) {
      jwFermionToPauliDoublePair(
        env.n_occ,
        env.n_virt,
        fermion_product[0], 
        fermion_product[1],
        fermion_product[2],
        fermion_product[3],
        output,
        hermitian);
    }
  }
}
  }; // namespace VQE 
}; // namespace NWQSim

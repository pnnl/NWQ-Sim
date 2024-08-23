# Fermion to Qubit Transformers
To perform electronic structure simulations, the Fermionic terms in the provided second-quantized Hamiltonian needs to be converted to Pauli operators. To do so, we provide a basic `Transformer` function template to allow for future mapper implementations (e.g. parity or Bravyi-Kiteav mappings). 
## Prototype
The basic Fermion->Qubit transformer function is defined as:
```c++
typedef void (*Transformer) (const MolecularEnvironment&, const std::vector<std::vector<FermionOperator> >&, std::vector<std::vector<PauliOperator> >& ,bool);
```

That is, a Transformer takes as input a `MolecularEnvironment` (for qubit indexing), a vector of `FermionicOperator` products, and a flag to indicate that the mapped operators should be Hermitian:
$$
\tau^{ij}_{rs}-\tau^{ij,\dagger}_{rs}
$$
where $\tau^{ij}_{rs}$ is an arbitrary Fermionic operator. This flag also adds an imaginary coefficient to the operators for time evolution in the later `Ansatz` object. The output is a vector of QWC-commuting `PauliOperator` groups.
## Jordan-Wigner Mapper
The JW mapper `getJordanWignerTransform` iterates over each Fermionic term and performs the JW expansion:
$$
a^\dagger_{j}=\frac{X_j-iY_j}{2}\bigotimes_{k=1}^{j-1}Z_k
$$
$$
a_{j}=\frac{X_j+iY_j}{2}\bigotimes_{k=1}^{j-1}Z_k
$$

The function `jwFermionToPauliSingle` implements the JW transformation for a single Fermionic operator.

For second-quantized Hamiltonians, we have two types of terms with the following function handlers.
1. Singles $a^\dagger_i a_j$: `jwFermionToPauliSinglePair` takes the product over the 4 Pauli terms produced and appends them to the `output` Pauli vector
2. Doubles $a^\dagger_i a^\dagger_j a_r a_s$: `jwFermionToPauliDoublePair` takes the product over the 16 Pauli terms produced and appends them to the `output` Pauli vector

At this point, the data structure `std::vector<std::vector<PauliOperator> > temp_output` contains a vector of `PauliOperators` for each Fermionic product. 

From here, `getJordanWignerTransform` behavior differs depending on the `hermitian` flag. If enabled, the function sets `output=temp_output` and returns: these operators are used to generate the UCCSD ansatz.

Otherwise, the Pauli terms are concatenated into a single vector and grouped into QWC sets using `sorted_insertion` (See Crawford et al. 2021, implementation in vqe/src/utils.cpp). The commuting groups are returned for expectation value calculation.


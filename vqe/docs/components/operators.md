# Operators
Here we outline the implementation and design details for the three main Operator classes:
1. Fermionic Operators
2. Pauli String Operators
3. Hamiltonians


## Molecular Environments
Translation between operator types requires problem-specific context. We store this information in `MolecularEnvironment` structs:
```c++
   struct MolecularEnvironment {
      /**
       * Data structure to store environment information (orbital/particle config, energy offsets)
      */
      IdxType n_spatial; // # spatial orbitals
      IdxType n_part;    // # electrons
      IdxType n_occ;     // # occupied orbitals
      IdxType n_virt;    // # virtual orbitals
      bool xacc_scheme;  // Whether to use XACC (true) or canonical (false) qubit indexing
      ValType constant;  // constant offset
```
The first four fields encode the orbital and particle information. Note that if `n_occ` and `n_virt` are known, then `n_spatial` and `n_part` are known, and vice versa.

`xacc_scheme` determines the orbital labeling scheme: `true` for XACC ordering (spin > virtual/occupied > orbital index) and `false` for canonical (virtual/occupied > spin >  orbital index).

The main function of the `MolecularEnvironment` class is data storage, however it also adds an API wrapper for obtaining the qubit index of a molecular orbital:
```c++
IdxType addr(IdxType orbital, bool is_virtual, bool spin_down) const {
  // Return the qubit index for a given orbital index, virtual/occupied status, and spin value
  return getQubitIndex(orbital, spin_down, is_virtual, n_occ, n_virt, xacc_scheme);
}
```
## Fermionic Operators
The `FermionicOperator` class (defined in vqe/include/observable/fermionic_operator.hpp) is used to store and reference data relating to a single Fermionic operator, either a creation $a_i^\dagger$ or an annihilation $a_i$ term, specified with an `enum`:
```c++

enum FermionOpType {
  Creation,
  Annihilation
};
```

The index of a Fermionic operator depends on:
1. The spatial orbital index
2. The spin of the associated orbital
3. The type orbital (virtual or occupied)

For convenience, spin and orbital type are expressed using `enum`s:
 ```c++
enum Spin {
  Up,
  Down
};
enum OrbitalType {
  Occupied,
  Virtual
};
 ```

The `FermionOperator` class stores orbital and operator attributes along with an associated complex coefficient:
```c++
class FermionOperator {
  // Basic IR before transformation operations (JW, BK, etc.)
  protected:
    IdxType orbital_index;
    OrbitalType orb_type;
    Spin spin;
    FermionOpType type;
    bool xacc_scheme;
    std::complex<ValType> coeff;
```

Like the `MolecularEnvironment` class, `FermionicOperator` objects also offer qubit indexing utilities:
```c++
IdxType qubitIndex(IdxType n_occ, IdxType n_virt) const {
  // Flattened indexing scheme
  return getQubitIndex(orbital_index, spin, orb_type, n_occ, n_virt, xacc_scheme);
}
```

To simplify object manipulation, we add in copy and assignment constuctors, as well as an overloaded function for scalar multiplication (good for switching sign/phase).

The following ADAPT-VQE functions are also forward-declared in fermionic_operator.hpp, and are covered in adapt.md
```c++
void generate_fermionic_excitations(std::vector<std::vector< std::vector<FermionOperator> > >& _fermion_operators,
                                const MolecularEnvironment& _env);
void generate_pauli_excitations(std::vector<std::vector<PauliOperator > >& _pauli_operators,
                                const MolecularEnvironment& _env,
                                IdxType subsample = -1,
                                IdxType seed = 0);
// Construct the minimal operator pool G from Tang et al. 2021 ("Qubit-ADAPT VQE")
void generate_minimal_pauli_excitations(std::vector<std::vector<PauliOperator > >& _pauli_operators,
                                const MolecularEnvironment& _env);
```

## Pauli Operators
Pauli string manipulation is a key part of the NWQ-VQE wrapper. The `PauliOperator` class (defined in vqe/include/observables/pauli_operator.hpp) stores the data of a single Pauli string of up to 64 qubits (limited by integer mask size):
```c++
    class PauliOperator {
      private:
        IdxType dim; // number of qubits
        bool non_trivial; // stores whether the Pauli is not identity
        std::complex<ValType> coeff; // complex scalar coefficient
        IdxType xmask; // X stabilizer array
        IdxType zmask; // Z stabilizer array
        std::vector<IdxType> x_indices; // vector of X/Y qubits for convenience
```
`PauliOperator` objects can be constructed in one of two ways: either via explicit `X` and `Z` binary masks (where the system size is explicitly provided):
```c++
PauliOperator(IdxType _xmask, IdxType _zmask, IdxType _dim, std::complex<ValType> _coeff = 1.0)
```
or via std::string:
```c++
PauliOperator(std::string _opstring, std::complex<ValType> _coeff = 1.0)
```
For example, the following two initializations are equivalent:
```c++
std::complex<double> coeff (1.0, 1.0); // 1 + i
PauliOperator pauli1("XYIZ", coeff);
```
```c++
std::complex<double> coeff (1.0, 1.0); // 1 + i
// Stabilizer matrices
IdxType x_stabilizer = 0b1100;
IdxType z_stabilizer = 0b0101;
IdxType n_qubits = 4;
PauliOperator pauli1(x_stabilizer, z_stabilizer, n_qubits, coeff);
```
The `PauliOperator` class also supports products between Pauli strings:
```c++
PauliOperator pauli1("XYIZ");
PauliOperator pauli2("YYIX");
PauliOperator pauli3 = pauli1 * pauli2;
std::cout << pauli3 << std::endl;
```
Result:
```
(-1 + 0i)ZIIY
```


Most of the manipulations occur via stabilizer operations. For instance, the result of operator multiplication with another `PauliOperator` is performed via logical XOR (though we still need to iterate to find the phase):
```c++
// from PauliOperator operator*(const PauliOperator& other) const
IdxType newxmask = xmask ^ other.xmask; // the X stabilizer array of the new operator
IdxType newzmask = zmask ^ other.zmask; // the Z stabilizer array of the new operator
```
or to determine the parity of two Pauli strings:
```c++
IdxType n_ones = count_ones((xmask & other.zmask) | (zmask & other.xmask)); // take the symplectic outer product
return (n_ones % 2) == 0; // n_ones is the number of anti-commuting qubits
```



## Hamiltonians
The `Hamiltonian` class acts as a container for multiple different operator representations. There are five relevant properties:

```c++
protected:
  MolecularEnvironment env;
  std::vector<std::vector<FermionOperator> > fermi_operators;
  std::vector<std::vector<PauliOperator> > pauli_operators;
  IdxType n_ops;
  Transformer qubit_transform;
```
1. `env`: A molecular environment instance used for orbital indexing
2. `fermi_operators`: The list of Fermionic operators making up the second-quantized Hamiltonian. Each `std::vector<FermionOperator>` is a single product term, e.g. $a_i^\dagger a_j^\dagger a_ra_s$
3. `pauli_operators`: The list of Pauli operators making up the JW-mapped Hamiltonian. Each `std::vector<PauliOperator>` is a QWC-commuting group, which can be simultaneously measured with single-qubit rotations
4. `n_ops`: The number of PauliOperator objects, just used for logging
5. `qubit_transform`: The Fermion to qubit mapper used to translate `fermi_operators` into `pauli_operators` (Default: Jordan-Wigner)

The Hamiltonian class supports four different constructor types:
1. Read from an XACC-formatted file. This is the most common, and is used by the command line solver.
```c++
Hamiltonian(std::string input_path, IdxType n_particles, bool xacc_scheme,
                    Transformer transform = getJordanWignerTransform);
```
2. An explicit list of `FermionOperator` objects. Note that the convention is to store the coefficient in the first `FermionOperator` in each product.
```c++
Hamiltonian(MolecularEnvironment _env,
            const std::vector<std::vector<PauliOperator> > _pauli_operators, 
            Transformer transform = getJordanWignerTransform)
```
3. An explicit list of `PauliOperator` objects. Constructors are provided for both situations with MolecularEnvironment information as well as particle/indexing scheme information. Note: we assume that the `PauliOperator` vectors are already QWC-compatible.
```c++
Hamiltonian(const std::vector<std::vector<PauliOperator> > _pauli_operators, 
            IdxType n_particles,
            bool use_xacc,
            Transformer transform = getJordanWignerTransform);
Hamiltonian(MolecularEnvironment _env,
            const std::vector<std::vector<PauliOperator> > _pauli_operators);
```

4. An explicit list of XACC-formatted Fermionic excitation strings and complex coefficients. Provided for the (currently unused) PyBind11 interface.
```c++
Hamiltonian(const std::vector<std::pair<std::string, std::complex<double>>>& input_ops, 
            IdxType n_particles, 
            bool xacc_scheme,
            Transformer transform);
``` 

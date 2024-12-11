# Ansätze
Here we describe the `Ansatz` base class and its core data structures, as well as the `UCCSD` subclass.

## Ansatz Base Class

The `Ansatz` class inherits from the generic `NWQSim::Circuit` class, adding the following data structures to permit optimization:
```c++
class Ansatz: public Circuit {
  protected:
    std::shared_ptr<std::vector<ValType> > theta;
    std::shared_ptr<std::vector<IdxType> > parameterized_gates; // Indices of parameterized gates
    std::shared_ptr<std::vector<std::vector<std::pair<IdxType, ValType> > > > gate_parameter_pointers; // Gate parameter indices
    std::shared_ptr<std::vector<ValType> > gate_coefficients; // Gate coefficients
    std::unordered_map<std::string, IdxType> excitation_index_map;
```

The Ansatz constructor only requires the number of qubits: `Ansatz(IdxType n_qubits)`. However, after construction, the Ansatz requires a call to `buildAnsatz`. For example:
```c++
UCCSD test_ansatz(10);  // 10 qubit UCCSD
test_ansatz.buildAnsatz();
```
The separate call is used by some Ansatz subclasses (namely `DynamicAnsatz`) to build the circuit after some initial step, such as setting the operator pool. 

We now discuss each in turn. 
### Data Structures
1. `std::shared_ptr<std::vector<ValType> > theta`: The core parameter vector, with one entry for each unique parameter. `theta` is assigned to by `void setParams(const std::vector<ValType>& params)` and referenced by either `std::vector<ValType>* getParams()` (non-const) or `std::vector<ValType>& getParamRef()` (const)
2. `std::shared_ptr<std::vector<IdxType> > parameterized_gates`: The list of parameterized gate addresses within the `Circuit` class `gates` vector. 
3. `std::shared_ptr<std::vector<std::vector<std::pair<IdxType, ValType> > > > gate_parameter_pointers`: A vector of linear parameter expressions for each parameterized gate. The innermost `std::vector<std::pair<IdxType, ValType> >` describes how to compute the gate angle for each parameterized gate. For instance:
```c++
gate_parameter_pointers->at(5) = {{2, -1}, {7, 0.5}};
```
indicates that the angle for gate 5 should be 0.5 times parameter 7, minus parameter 2. That is:
```c++
gates->at(5).theta = -1 * theta->at(2) + 0.5 * theta->at(7);
```
This is how symmetry groups are explicitly implemented, explained in the UCCSD section.

4. `std::shared_ptr<std::vector<ValType> > gate_coefficients`: The coefficients for each gate angle, typically either 1.0 or 0.25 depending on the size of the `FermionOperator` product used to produce the gate
5. `std::unordered_map<std::string, IdxType> excitation_index_map`: Used to map between Fermionic operator strings and parameter indices. Used for explicit parameter inputs for `QFlow`.

### Functions
The `Ansatz` class implements the following functions:

The first two implement parameterized circuit modification:
1. `void assignGateParam(IdxType gate_index, ValType param)`: Assign a single gate angle
2. `virtual void setParams(const std::vector<ValType>& params)`: Update the parameter vector and all gate angles

The next two are general circuit utilities found to be useful


3. `std::string toQASM3()` (Work in Progress): Return the parameterized circuit as a OpenQASM3 formatted string. QASM3 is required for unbound parameters, otherwise the circuit would need to fix the gate angles.
4. `void compose(const Circuit& other, std::vector<IdxType>& qubit_mapping)`: Compose the Ansatz with another circuit (e.g. for measurement). Also just a generally useful utility for circuit construction

The following two are used for ansatz construction using Pauli operators:
5. `void OneParamGate(enum OP _op_name...)`: Construct a single parameterized gate and update the associated data structures. Called by ExponentialGate, by default assumed to be an `RZ` gate. NOTE: To use with DMSim, this would have to be modified to support the gate type as an input to match arbitrary bases
6. `void ExponentialGate(const PauliOperator& pauli_op...)`: Construct an operator evolving $e^{pauli_op}$. Here we assume `pauli_op` already has the appropriate phase/coefficients. The circuit constructs a basic CNOT ladder terminating at the first non-identity qubit. NOTE: Room for optimization via Paulihedral or Tetris for 2-qubit gate cancellation

The following three functions are not implemented in the base class and require overloading:
7. `std::vector<std::string> getFermionicOperatorStrings()`: Return the string representations of the operators used to construct the Ansatz. Despite the name, these may also be Pauli operators in the case of Qubit-ADAPT
8. `std::vector<std::pair<std::string, ValType> > getFermionicOperatorParameters()`: Return the string representations of the Fermionic operators used to construct the Ansatz along with their associated excitations. Note that there are more strings returned than parameters, as we print all members of a symmetry group and their spin-reversed forms.

## UCCSD Class
The `UCCSD` class implements the functionality defined in the `Ansatz` class to construct an operator based on Fermionic single and double excitations. 

### Properties and Data Structures
The UCCSD class tracks several properties relevant to circuit construction for debugging, memory allocation, and qubit dereferencing:
```c++
    class UCCSD: public Ansatz {
      protected:
        const MolecularEnvironment& env; // reference to a MolecularEnvironment structure
        IdxType n_singles; // Number of single excitations
        IdxType n_doubles; // Number of double excitations
        IdxType trotter_n; // Number of Trotter steps
        IdxType unique_params; // Number of unique parameters
        Transformer qubit_transform; // Fermion->Qubit mapper
```
The `n_singles` and `n_doubles` fields are used for memory allocation/reservation in the following critical data structures:
```c++
std::vector<std::vector<std::pair<IdxType, ValType> > > symmetries;
```
Tracks symmetries between Fermionic operators, used to later compute the `gate_parameter_pointers` data structure. `symmetries` has one entry per Fermionic symmetry group, whereas `gate_parameter_pointers` replicates the entries multiple times for each parameterized gate sharing an angle.
```c++
std::vector<IdxType> fermion_ops_to_params;
```
Maps symmetry groups to parameter indices, also used to construct `gate_parameter_pointers`. Note to self: `fermion_ops_to_parameters` is probably somewhat redundant in the face of `symmetries`, will have to look at that later on.
```c++
std::vector<std::vector<FermionOperator> > fermion_operators;
```
List of `FermionOperator` product terms used to generate the ansatz

### Excitation Generation
The number of single excitations can easily be computed as $N_{Occupied}\cdot N_{Virtual}$. As an example:
```c++
IdxType orbital_index = 0;
IdxType unique_parameters = 0;
// make down spin operators
FermionOperator occ_down (0, Occupied, Down, Annihilation, true);
FermionOperator virt_down (0, Virtual, Down, Creation, true);
// make spin up operators
FermionOperator occ_up (0, Occupied, Up, Annihilation, true);
FermionOperator virt_up (0, Virtual, Up, Creation, true);
// Record the operator products
fermion_operators.push_back({occ_down, virt_down});
fermion_operators.push_back({occ_down, virt_down});
// Record the symmetry (both point to )
symmetries[0] = {{0, 1.0}};
symmetries[1] = {{0, 1.0}};
fermion_ops_to_params[0] = 0;
```

A symmetry is further used to reduce the number of parameters required for both single and double excitations, in general there are two cases:
1. When two operators only differ in the spin but not spin-orbital index: $t_{I(\alpha)}^{A(\alpha)} = t_{I(\beta)}^{A(\beta)}$ and $t_{I(\alpha)J(\alpha)}^{A(\alpha)B(\alpha)} = t_{I(\beta)J(\beta)}^{A(\beta)B(\beta)}$
2. For double exicitations, when the order of indices flipped for creation and annihilation operators: $t_{I(\alpha)J(\beta)}^{A(\alpha)B(\beta)} = t_{J(\alpha)I(\beta)}^{B(\alpha)A(\beta)}$

<--
However, the number of double excitations requires a bit more combinatorics. We can delineate *mixed terms* (where all spatial orbitals are unique) from *degenerate terms* (where either the annihilation or creation operators share a spatial orbital). These two groups have different symmetry expressions:

Using bars to denote spin orientation ($\alpha_i$ for a spin up annihilation, $\beta_i$ for spin down), we can express the symmetries for mixed excitations as follows:

```math
\alpha_i^\dagger\alpha_j^\dagger\alpha_r\alpha_s=\beta_i^\dagger\beta_j^\dagger\beta_r\beta_s=\alpha_i^\dagger\beta_j^\dagger\beta_r\alpha_s+\alpha_i^\dagger\beta_j^\dagger\alpha_r\beta_s=\beta_i^\dagger\alpha_j^\dagger\alpha_r\beta_s+\beta_i^\dagger\alpha_j^\dagger\beta_r\alpha_s
```

whereas degenerate excitations have the symmetry:
```math
\beta_i^\dagger\alpha_j^\dagger\alpha_r\beta_r=\beta_j^\dagger\alpha_i^\dagger\alpha_r\beta_r
```
with a similar form for degenerate virtual orbitals.)
-->


### Construction
The UCCSD ansatz is built during the call to `buildAnsatz`. 

The circuit starts with the Hartree-Fock state, with the orbital locations depending on the indexing scheme used:
```c++
if (env.xacc_scheme) {
for (IdxType i = 0; i < env.n_occ; i++) {
  X(i);
  X(i+env.n_spatial);
}
} else {
  for (IdxType i = 0; i < 2 * env.n_occ; i++) {
    X(i);
  }
}
```
The constructor then calls the Fermion->Qubit mapper and adds an EvolutionGate for each non-trivial string (non-identity with finite coefficient):
```c++
qubit_transform(env, fermion_operators, pauli_oplist, true);  
IdxType index = 0; // parameter index, shares parameters for Pauli evolution gates corresponding to the same Fermionic operator within the same Trotter step
for (auto& fermionic_group: pauli_oplist) {
  bool wasused = 0;
  for (auto& pauli: fermionic_group) {
    assert (pauli.getCoeff().imag() == 0.0);
    double coeff = pauli.getCoeff().real(); 
    
    std::vector<std::pair<IdxType, ValType> > idxvals(symmetries[index].size());
    // construct the index map for the Ansatz gate_parameter_pointers object
    std::transform(symmetries[index].begin(), symmetries[index].end(),idxvals.begin(), 
      [&] (std::pair<IdxType, ValType> val) {
        return std::make_pair(fermion_ops_to_params[val.first], val.second);
      } );
    // If non-zero coefficient, add to the Ansatz
    if (pauli.isNonTrivial() && abs(coeff) > 1e-10) { 
      ExponentialGate(pauli, OP::RZ, idxvals, 2 * coeff);
      wasused = 1;
    }
  }
  // Sanity check
  if(!wasused) {
    printf("%lld parameter not used for operator\n", index);
  }
  index++;
}
```
The same step is repeated if multiple Trotter steps are requested. The UCCSD Ansatz is now ready for simulation.

# Expectation Values
Here we discuss the method of computing expectation values in NWQ-Sim, starting with the relevant data structures and then discussing simulation details.

For an operator in the diagonal basis of the form:
```math
P_i=c_i\bigotimes_{j}Z_j
```
where identity operators are implicit between the qubits $j$. The expectation value with respect to a statevector $\vert \psi\rangle$ can be computed as:
```math
\langle \psi\vert P_i\vert \psi\rangle=\sum_{j=1}^{2^N}c_jParity(j, P_i)|\psi_{j}|^2
```
where $Parity(j, P_i)$ measures the bitwise parity between the Z stabilizer array of $P_i$ and the binary representation of index $j$. For example, $IZZI$ yields the stabilizer $0110$, which has an even parity with $6=0110$ (+1) and an odd parity with $5=0101$ (-1). To compute each expectation value, we need to: 
1. Diagonalize the Pauli operator with a measurement circuit
2. Compute and store the summation above
3. Reverse the diagonalization to return to the previous statevector.


## ObservableLists
The `ObservableList` structure is defined in include/state.hpp:
```c++
struct ObservableList {
    IdxType* zmasks;
    ValType exp_output;
    ValType* coeffs;
    IdxType numterms;
};
```
- `zmasks` is a pointer to an array of bitmasks, each indicating the Z stabilizer of a Pauli operator in the diagonal basis.
- `exp_output` is a field to store an accumulated expectation value sum over multiple Pauli operators
- `coeffs` stores the Pauli coefficients
- `numterms` stores the number of Pauli terms to compute for which we compute the expectation value (equivalently, the size of `zmasks` and `coeffs`)

Each `ObservableList` object and its fields are allocated either in the `fill_obslist` or ``
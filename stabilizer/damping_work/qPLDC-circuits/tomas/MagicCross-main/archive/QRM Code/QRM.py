import stim
import numpy as np
from itertools import product
import sys
from typing import List, Tuple

# Reed-Muller Generator Matrix
# This function generates the Reed-Muller generator matrix for a given order r and number of variables m.
# The 'variation' parameter can be 'None', 'punctured', or 'shortened'.
# - 'None' returns the full generator matrix.
# - 'punctured' removes the first column (corresponding to the evaluation point 0...0).
# - 'shortened' removes the first column and the first row (corresponding to the constant polynomial, which is all-0 monomial). 

def RM_generator_matrix(r: int, m: int, variation: str) -> Tuple[np.ndarray, List[Tuple[int]]]:
    
    n = 2 ** m
    monomials = []

    for deg in range(r + 1):
        for bits in product([0, 1], repeat=m):
            if sum(bits) == deg:
                monomials.append(bits)

    G = []
    for mono in monomials:
        row = []
        for i in range(n):
            x = [int(b) for b in bin(i)[2:].rjust(m, '0')]
            val = 1
            for a, xi in zip(mono, x):
                if a == 1:
                    val &= xi
            row.append(val)
        G.append(row)
    
    G = np.array(G, dtype=int)
    
    if variation == "None":
        G_final = G
        monomials_final = monomials
    
    if variation == "punctured":
        # Remove the first column (corresponding to evaluation point 0...0)
        G_final = G[:, 1:]
        monomials_final = monomials

    if variation == "shortened":
        # Remove the first column (corresponding to evaluation point 0...0)
        G_final = G[:, 1:]

        # Remove the first row (corresponding to constant polynomial = all-0 monomial)
        G_final = G_final[1:, :]
        monomials_final = monomials[1:]
        
    return G_final, monomials_final


# Recursive Reed-Muller Decoder

def repetition_ml_decoder(y: np.ndarray) -> np.ndarray:
    """
    RM(0, m): all-0 or all-1 codeword; choose the better one
    """
    factor = 0.99 # Avoid extreme values of y (+1 or -1)
    score0 = np.prod((1 + factor * y) / 2)   # likelihood of all 0 codeword (+1)
    score1 = np.prod((1 - factor * y) / 2)   # likelihood of all 1 codeword (-1)

    c0 = np.ones_like(y)
    c1 = -np.ones_like(y)

    # print(f"Repetition ML scores: {score0:.4f} (0), {score1:.4f} (1)")
    if score0 > score1:
        return c0
    elif score1 > score0:
        return c1
    else:
        # If scores are equal, return with equal probability, unbiased
        return c0 if np.random.rand() < 0.5 else c1

def rm_list_decode(y, r, m, L=4) -> np.ndarray:
    """
    Recursive list decoder Ψ_m^r(L)
    y: input vector in [-1, 1]^n
    r, m: RM(r, m)
    L: list size
    Returns: list of (codeword, path_cost)
    """
    # print("r = ", r, ", m = ", m)
    factor = 0.99  # Avoid extreme values of y (+1 or -1)
    N = len(y)

    # Base cases
    if r == 0:
        return repetition_ml_decoder(y)
    if r == m:
        # RM(m,m) is full space — just use hard-decision
        return np.sign(y + (1-factor) * np.random.uniform(-1, 1, size=N))  # Add noise to avoid extreme values

    yL, yR = y[:N//2], y[N//2:]

    # Step 1: estimate v from componentwise product y_v = yL * yR
    y_v = yL * yR
    v_hat = rm_list_decode(y_v, r-1, m-1, L)

    # Step 2: for each v_hat, estimate u_hat
    y_hat = yR * v_hat  # Compute ŷ = yR * v_hat
    # Compute y_u = (yL + ŷ) / (1 + yL * ŷ)
    denom = 1 + factor * yL * y_hat
    y_u = (yL + y_hat) / denom
    u_hat = rm_list_decode(y_u, r, m-1, L)

    cw = np.concatenate([u_hat, u_hat * v_hat])

    return cw

# Utility functions for bit manipulation
bin_wt = lambda i: bin(i)[2:].count('1')
bit_rev = lambda t, m: int(bin(t)[2:].rjust(m, '0')[::-1], 2)
int2bin = lambda i, m: [int(c) for c in bin(i)[2:].rjust(m, '0')]
bin2int = lambda l: int(''.join(map(str, l)), 2)

# Non-Fault-Tolerant QRM Hypercube Encoding Circuit

def QRM_encode(rx: int, rz: int, m: int, noise_profile: list, observable_type: str, punctured: bool = True, final_meas: str = "Z") -> stim.Circuit:
    p1, p2, p_M, p_R = noise_profile
    N = 2**m
    data_coords = [(0,0)]
    for r in range(m):
        if r % 2 == 0:
            coords_shifted = [(x + 2**(r//2), y) for (x, y) in data_coords]
        else:
            coords_shifted = [(x, y + 2**(r//2)) for (x, y) in data_coords]
        data_coords += coords_shifted

    coord_to_index = {coord: i for i, coord in enumerate(data_coords)}

    if punctured:
        del coord_to_index[(0,0)] # Remove the punctured qubit at (0,0)

    qubit_indices = list(coord_to_index.values())

    # non-punctured or punctured code
    full_circuit = stim.Circuit()
    # repeat_circuit = stim.Circuit()

    # QUBIT_COORDS annotations
    for coord, index in coord_to_index.items():
        full_circuit.append_operation("QUBIT_COORDS", [index], list(coord))

    # Initialize the qubits in |0> or |+> state
    X_state_indices = []
    Z_state_indices = []
    prepared_indices = []
        
    if observable_type == "Z":
        for i in qubit_indices:
            if bin_wt(i) <= rx:
                X_state_indices.append(i)
            elif bin_wt(i) >= m - rz:
                Z_state_indices.append(i)
            else:
                prepared_indices.append(i)
        full_circuit.append_operation("RX", X_state_indices)
        full_circuit.append_operation("Z_ERROR", X_state_indices, p_R)
        full_circuit.append_operation("R", Z_state_indices + prepared_indices)
        full_circuit.append_operation("X_ERROR", Z_state_indices + prepared_indices, p_R)
        
    if observable_type == "X": # flip the state preparation
        for i in qubit_indices:
            if bin_wt(i) <= rz:
                Z_state_indices.append(i)
            elif bin_wt(i) >= m - rx:
                X_state_indices.append(i)
            else:
                prepared_indices.append(i)
        full_circuit.append_operation("RX", X_state_indices + prepared_indices)
        full_circuit.append_operation("Z_ERROR", X_state_indices + prepared_indices, p_R)
        full_circuit.append_operation("R", Z_state_indices)
        full_circuit.append_operation("X_ERROR", Z_state_indices, p_R)

    full_circuit.append("TICK")

    # non-FT hypercube encoding circuit
    for t in range(m): # rounds
        sep = 2 ** t
        CNOT_indices = []
        for i in qubit_indices:
            if int2bin(i,m)[-1-t] == 0 and i + sep < N:
                if observable_type == "Z":
                    CNOT_indices += [i, i + sep]
                else:
                    CNOT_indices += [i + sep, i]      
        full_circuit.append_operation("CNOT", CNOT_indices)
        full_circuit.append("DEPOLARIZE2", CNOT_indices, p2)
        full_circuit.append("TICK")

    # Final measurement for logical observables, should not inclue noise
    if final_meas == "Z":
        full_circuit.append_operation("M", qubit_indices)
    if final_meas == "X":
        full_circuit.append_operation("MX", qubit_indices)
    if final_meas == "Y":
        full_circuit.append_operation("MY", qubit_indices)

    return full_circuit

# An optimization: if the control qubit is zero, omit the CNOT gate at the first time step
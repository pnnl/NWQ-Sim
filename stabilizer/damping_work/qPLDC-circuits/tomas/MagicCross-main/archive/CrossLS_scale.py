# %%
import stim
from planar_PQRM import PlanarSurfaceCode, MultiPatchSurfaceCode
import pymatching
import sinter
from typing import Set, Tuple, List
from itertools import product
import copy
import pandas as pd
import seaborn as sns
import numpy as np
from numpy.typing import NDArray
import matplotlib.pyplot as plt
from stimbposd import BPOSD
from stimbposd import SinterDecoder_BPOSD, sinter_decoders

# %%
# Compute the logical error rate (LER) per round
def convert_to_per_round(total_LER: float, rounds: int) -> float:
    per_round_LER = 1/2 * (1 - np.abs(1-2*total_LER)**(1/rounds))
    return per_round_LER

# Utility functions for bit manipulation
bin_wt = lambda i: bin(i)[2:].count('1')
bit_rev = lambda t, m: int(bin(t)[2:].rjust(m, '0')[::-1], 2)
int2bin = lambda i, m: [int(c) for c in bin(i)[2:].rjust(m, '0')]
bin2int = lambda l: int(''.join(map(str, l)), 2)

# %%
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

# %%
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

def RM_decoder(r, m, noisy_cws): # punctured RM code decoder
    decoded_cws = []
    # Insert the 0-th bit (punctured)
    for noisy_cw in noisy_cws:
        y = (-1) ** noisy_cw
        y0 = np.insert(y, 0, 1)  # Insert a 1 at the 0-th bit
        y1 = np.insert(y, 0, -1)  # Insert a -1 at the 0-th bit
        decoded_cw0 = rm_list_decode(y0, r, m, L=4)
        decoded_cw1 = rm_list_decode(y1, r, m, L=4)
        overlap0 = np.sum(y0 * decoded_cw0)
        overlap1 = np.sum(y1 * decoded_cw1)
        if overlap0 > overlap1:
            decoded_cw = decoded_cw0
        else:
            decoded_cw = decoded_cw1
        decoded_cw = (1 - decoded_cw[1:]) // 2  # Convert from ±1 to {0, 1}
        decoded_cws.append(decoded_cw)
    return np.array(decoded_cws, dtype=int)

# %%
# Lattice surgery between planar surface codes and PQRM code
# The surface code will be initialized in logical |+> state, the PQRM code will be initialized in |0> or |+> state for tomography
def CrossLS(PQRM_para: List[int], d_surf: int, rounds: int, PQRM_state: str, noise_profile: List[float], punctured: bool = True) -> stim.Circuit:
    """
    Constructs a lattice surgery circuit between a PQRM code and a surface code.
    
    Parameters:
    - PQRM_para: List of parameters for the PQRM code [rx, rz, m].
    - d_surf: The distance of the surface code.
    - punctured: Whether to use a punctured PQRM code.
    
    Returns:
    - A stim.Circuit object representing the lattice surgery circuit.

    For now, we only consider PQRM with several parameters: 
    [[15,1,3]]: PQRM(1,2,4), [[31,1,3]]: PQRM(1,3,5), 
    [[63,1,3]]: PQRM(1,4,6), [[128,1,3]]: PQRM(1,5,7)
    They have the advantage of weight-4 stabilizers
    Potentially incorporate: [[128,1,7]]: PQRM(2,4,7)
    """
    
    rx, rz, m = PQRM_para
    N = 2 ** m
    log_PQRM_len_dict = {(1,2,4): 3, (1,3,5): 7, (1,4,6): 7, (1,5,7): 15, (2,4,7): 15}
    log_PQRM_len = log_PQRM_len_dict[(rx, rz, m)] # We require log_PQRM_len to be <= d_surf for ZZ interaction
    
    if log_PQRM_len and log_PQRM_len > d_surf:
        raise ValueError("We require log_PQRM_len <= d_surf for ZZ interaction.")

    # Set noise parameters
    p1, p2, p_M, p_R = noise_profile

    # Construct PQRM code
    data_coords_PQRM = [(0,0)]
    for r in range(m):
        if r % 2 == 0: # shift down first
            coords_shifted = [(x, y + 2**(r//2)) for (x, y) in data_coords_PQRM]
        else: # then shift to right
            coords_shifted = [(x + 2**(r//2), y) for (x, y) in data_coords_PQRM]
        data_coords_PQRM += coords_shifted

    # Layout
    x_range_PQRM = 2 ** (m // 2)
    y_range_PQRM = 2 ** ((m+1) // 2)

    # Scale the PQRM coordinates to match the surface code's scale
    data_coords_PQRM = [(2 * x, 2 * y) for (x, y) in data_coords_PQRM]
    coord_to_index_PQRM = {coord: i for i, coord in enumerate(data_coords_PQRM)}
    index_to_coord_PQRM = {i: coord for coord, i in coord_to_index_PQRM.items()}

    if punctured:
        data_coords_PQRM.remove((0,0))
        del coord_to_index_PQRM[(0,0)] # Remove the punctured qubit at (0,0)

    qubit_indices_PQRM = list(coord_to_index_PQRM.values())
    PQRM_logical_Z = [(0, 2 * (y + 1)) for y in range(log_PQRM_len)]


    # Construct surface code, shifted to the left
    d = d_surf
    sf_code = PlanarSurfaceCode(d, shift=(-2*d, 2))
    # Update logical Z of code 1 from left side to right side for merging
    sf_code.logical_Z[0] = PlanarSurfaceCode.shift_coords(sf_code.logical_Z[0], 2 * (d - 1), 0)

    # Construct the ancilla system for lattice surgery
    sydrome_coords_middle = [(-1, y) for y in range(2, 2 * d + 1, 2)]
    syndrome_coords_PQRM_bulk = [(2*x+1, 2*y+1) for x in range(x_range_PQRM-1) for y in range(y_range_PQRM-1)]
    syndrome_coords_PQRM_bulk.remove((1,1))
    syndrome_coords_PQRM_bound = [(8*x+3, 2*y_range_PQRM-1) for x in range(x_range_PQRM//4)] + [(2*x_range_PQRM-1, 8*y+3) for y in range(y_range_PQRM//4)]

    # Notice that these syndromes don't give the full set of Z-stabilizers
    sydrome_coords_ancsys_z = sydrome_coords_middle + syndrome_coords_PQRM_bulk + syndrome_coords_PQRM_bound
    sydrome_coords_ancsys = sydrome_coords_ancsys_z # potentially define X-stabilizers
    data_coords_ancsys = [(-1, y) for y in range(3, 2 * d, 2)]
    coord_to_index_ancsys = {coord: i for i, coord in enumerate(sydrome_coords_ancsys + data_coords_ancsys)}

    # Combine the PQRM code, surface code, ancilla system, set QUBIT_COORDS
    # qubit_coords = sf_code.qubit_coords + data_coords_PQRM
    coord_to_index = {}
    coord_to_index.update(coord_to_index_PQRM)
    offset1 = max(coord_to_index.values()) + 1
    coord_to_index_surf = sf_code.coord_to_index.copy()
    coord_to_index_shifted = {coord: i + offset1 for coord, i in coord_to_index_surf.items()}
    sf_code.coord_to_index = coord_to_index_shifted
    coord_to_index.update(coord_to_index_shifted)
    offset2 = max(coord_to_index.values()) + 1
    coord_to_index_ancsys_shifted = {coord: i + offset2 for coord, i in coord_to_index_ancsys.items()}
    coord_to_index.update(coord_to_index_ancsys_shifted)

    # Build the full circuit

    full_circuit = stim.Circuit()

    for (x, y), index in coord_to_index.items():
        full_circuit.append('QUBIT_COORDS', [index], [x, y])

    # Initialize the PQRM code in |0> state noiselessly
    X_state_indices = []
    Z_state_indices = []
    prepared_indices = []
        
    if PQRM_state == "Z":
        for i in qubit_indices_PQRM:
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
    if PQRM_state == "X": # flip the state preparation
        for i in qubit_indices_PQRM:
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
        for i in qubit_indices_PQRM:
            if int2bin(i,m)[-1-t] == 0 and i + sep < N:
                if PQRM_state == "Z":
                    CNOT_indices += [i, i + sep]
                else:
                    CNOT_indices += [i + sep, i]      
        full_circuit.append_operation("CNOT", CNOT_indices)
        full_circuit.append("DEPOLARIZE2", CNOT_indices, p2)
        full_circuit.append("TICK")

    # Initialize the surface code in |+> state, one round of noiseless SM
    data_indices_surf = [sf_code.coord_to_index[c] for c in sf_code.data_coords]
    ancilla_indices_surf = [sf_code.coord_to_index[c] for c in sf_code.ancilla_coords]
    full_circuit.append_operation("RX", data_indices_surf)
    full_circuit.append_operation("R", ancilla_indices_surf)
    full_circuit.append("TICK")
    full_circuit += sf_code.build_standard_sm_round(noise_profile)
    # the first round, detectors are the stabilizers already created by initialization
    for k, ancilla in enumerate(reversed(sf_code.ancilla_coords)):
        if ancilla in sf_code.ancilla_coords_x:
            full_circuit.append_operation("DETECTOR", [stim.target_rec(-k - 1)], list(ancilla) + [0])
    full_circuit.append("TICK")
    full_circuit = full_circuit.without_noise() # Remove the noise up to now!

    # Lattice surgery with PQRM code and surface code
    # Initialize the ancilla system for lattice surgery with noise
    sydrome_coords_ancsys_indices = [coord_to_index[c] for c in sydrome_coords_ancsys]
    data_coords_ancsys_indices = [coord_to_index[c] for c in data_coords_ancsys]
    full_circuit.append_operation("R", sydrome_coords_ancsys_indices)
    full_circuit.append_operation("X_ERROR", sydrome_coords_ancsys_indices, p_R)
    full_circuit.append_operation("RX", data_coords_ancsys_indices)
    full_circuit.append_operation("Z_ERROR", data_coords_ancsys_indices, p_R)
    full_circuit.append("TICK")

    # Incorporate the ancilla system into the surface code, build unified SM circuits
    sf_code_extended = copy.deepcopy(sf_code)
    # PQRM_data_extra = [(2*x,2*y+2) for x in range(1,4) for y in range(d_PQRM)] # only for extended surface code
    sf_code_extended.data_coords = sf_code_extended.data_coords + data_coords_ancsys + data_coords_PQRM
    sf_code_extended.ancilla_coords_z += sydrome_coords_ancsys_z
    # sf_code_extended.ancilla_coords_x += sydrome_coords_ancsys_x
    sf_code_extended.ancilla_coords += sydrome_coords_ancsys
    sf_code_extended.coord_to_index = coord_to_index
    # Update stabilizers
    stabilizers = sf_code_extended.stabilizers

    # Update the boundary X-stabilizers from weight-3 to weight-4
    for stab in stabilizers:
        anc = stab['ancilla']
        if anc[0] == -2:
            stab['data'] = [(anc[0]-1, anc[1]), (anc[0]+1, anc[1]), (anc[0], anc[1]-1), (anc[0], anc[1]+1)]

    # Create the new boundary Z/X-stabilizers 
    for syn in sydrome_coords_ancsys:
        if syn[0] == -1:
            if syn[1] == 2: # top weight-3 stabilizer
                stab = {'type': 'Z', 'ancilla': syn, 'data': [(syn[0] - 1, syn[1]), (syn[0] + 1, syn[1]), (syn[0], syn[1] + 1)]}
            elif syn[1] <= 2*log_PQRM_len-2: # middle weight-4 stabilizers
                stab = {'type': 'Z', 'ancilla': syn, 'data': [(syn[0]-1, syn[1]), (syn[0]+1, syn[1]), (syn[0], syn[1]-1), (syn[0], syn[1]+1)]}
            elif syn[1] == 2*log_PQRM_len:
                if d == log_PQRM_len: # weight-3
                    stab = {'type': 'Z', 'ancilla': syn, 'data': [(syn[0]-1, syn[1]), (syn[0]+1, syn[1]), (syn[0], syn[1]-1)]}
                else: # weight-4
                    stab = {'type': 'Z', 'ancilla': syn, 'data': [(syn[0]-1, syn[1]), (syn[0]+1, syn[1]), (syn[0], syn[1]-1), (syn[0], syn[1]+1)]}
            elif syn[1] < 2*d:
                stab = {'type': 'Z', 'ancilla': syn, 'data': [(syn[0]-1, syn[1]), (syn[0], syn[1]-1), (syn[0], syn[1]+1)]}
            else: # syn[1] == 2 * d:
                if d > log_PQRM_len: # bottom weight-2 stabilizer
                    stab = {'type': 'Z', 'ancilla': syn, 'data': [(syn[0] - 1, syn[1]), (syn[0], syn[1] - 1)]}
                else: # d = log_PQRM_len, bottom weight-3 stabilizers
                    stab = {'type': 'Z', 'ancilla': syn, 'data': [(syn[0] - 1, syn[1]), (syn[0], syn[1] - 1), (syn[0] + 1, syn[1])]}
        # if syn[0] == 0: # weight-3 X-stabilizers
        #     stab = {'type': 'X', 'ancilla': syn, 'data': [(syn[0] - 1, syn[1]), (syn[0], syn[1] + 1), (syn[0], syn[1] - 1)]}
        if syn[0] in range(1, 2*x_range_PQRM-2) and syn[1] in range(1,2*y_range_PQRM-2): # weight-4 Z-stabilizers in PQRM
            stab = {'type': 'Z', 'ancilla': syn, 'data': [(syn[0]-1, syn[1]-1), (syn[0]+1, syn[1]-1), (syn[0]-1, syn[1]+1), (syn[0]+1, syn[1]+1)]}
        if syn[0] >= 2*x_range_PQRM-1:
            stab = {'type': 'Z', 'ancilla': syn, 'data': [(syn[0]-1, syn[1]-3), (syn[0]-1, syn[1]-1), (syn[0]-1, syn[1]+1), (syn[0]-1, syn[1]+3)]}
        if syn[1] >= 2*y_range_PQRM-1:
            stab = {'type': 'Z', 'ancilla': syn, 'data': [(syn[0]-3, syn[1]-1), (syn[0]-1, syn[1]-1), (syn[0]+1, syn[1]-1), (syn[0]+3, syn[1]-1)]}
        stabilizers.append(stab)
    sf_code_extended.stabilizers = stabilizers

    # Insert one round noisy SM for the extended surface code, fix the guage
    full_circuit += sf_code_extended.build_standard_sm_round(noise_profile)
    # Detectors:
    for k, ancilla in enumerate(reversed(sf_code_extended.ancilla_coords)):
        prev = -k - 1 - len(sf_code.ancilla_coords) # not the extended ancilla coords, but the original ones!
        if ancilla[0] < -1:
            full_circuit.append_operation("DETECTOR", [stim.target_rec(-k - 1), stim.target_rec(prev)], list(ancilla) + [0])
        if ancilla[0] > 0:
            full_circuit.append_operation("DETECTOR", [stim.target_rec(-k - 1)], list(ancilla) + [0])

    # Later rounds, repeated meausrements, all stabilizers gives detectors
    repeat_circuit = stim.Circuit()
    repeat_circuit.append("TICK")
    repeat_circuit += sf_code_extended.build_standard_sm_round(noise_profile)
    repeat_circuit.append_operation("SHIFT_COORDS", [], [0,0,1])

    for k, ancilla in enumerate(reversed(sf_code_extended.ancilla_coords)):
        prev = -k - 1 - len(sf_code_extended.ancilla_coords) # now it's the extended ancilla coords
        repeat_circuit.append_operation("DETECTOR", [stim.target_rec(-k - 1), stim.target_rec(prev)], list(ancilla) + [0])

    full_circuit += repeat_circuit * (rounds-1)

    # Final measurements 
    # Measure PQRM code in X basis along with the middle data qubits in the ancilla system
    full_circuit.append("TICK")
    final_data1 = data_coords_PQRM + data_coords_ancsys
    final_data1_indices = qubit_indices_PQRM + data_coords_ancsys_indices
    full_circuit.append_operation("Z_ERROR", final_data1_indices, p_M)
    full_circuit.append_operation("MX", final_data1_indices)


    # Measure the surface code patch in Z or X basis depending on the the initial state of PQRM (the state to be teleported)
    final_data2 = sf_code.data_coords
    final_data2_indices = [sf_code.coord_to_index[coord] for coord in final_data2]
    if PQRM_state == "Z":
        full_circuit.append_operation("X_ERROR", final_data2_indices, p_M)
        full_circuit.append_operation("M", final_data2_indices)
    if PQRM_state == "X":
        full_circuit.append_operation("Z_ERROR", final_data2_indices, p_M)
        full_circuit.append_operation("MX", final_data2_indices)

    # Final detectors for the surface code patch
    for stab in sf_code_extended.stabilizers:
        if stab['type'] == PQRM_state and stab['ancilla'][0] <= -2:
            ancilla_index = sf_code_extended.ancilla_coords[::-1].index(stab['ancilla'])
            anc_rec = stim.target_rec(-ancilla_index - 1 - len(final_data1+final_data2))
            data_rec_targets = [stim.target_rec(-(final_data1+final_data2)[::-1].index(q) - 1) for q in stab['data']]
            full_circuit.append_operation("DETECTOR", [anc_rec] + data_rec_targets, list(stab['ancilla']) + [1])
    
    # Constructing detectors for PQRM
    # All the X stabilizers extended by the data in the ancsys
    # X stabilizers in PQRM span on rows or columns; for row-stabilizers, combine the ancsys data in between to form stabilizers
    # X stabilizers will be postselected in the end
    PQRM_stabilizer_x, _ = RM_generator_matrix(rx, m, variation="shortened") # noitce that they are already punctured
    stabilizers_x_indices = []
    logical_Z_indices = [coord_to_index_PQRM[coord] for coord in PQRM_logical_Z]
    binary_logical_Z = np.zeros(N, dtype=int)
    binary_logical_Z[logical_Z_indices] = 1
    binary_logical_Z = binary_logical_Z[1:]

    for stab_x in PQRM_stabilizer_x:
        bin_prod = stab_x * binary_logical_Z
        if np.any(bin_prod != 0):
            indices = np.where(bin_prod != 0)[0] + 1
            intersect_qubits_coord = [index_to_coord_PQRM[ind] for ind in indices]
            y_coords = [coord[1] for coord in intersect_qubits_coord]
            data_added = []
            for i in range(0, len(y_coords), 2):
                start, end = y_coords[i], y_coords[i+1]
                data_added.extend([(-1,y) for y in range(start + 1, end) if y % 2 == 1])
            data_added_indices = [coord_to_index[coord] for coord in data_added]
            stab_x_indices = list(np.where(stab_x != 0)[0] + 1) + data_added_indices
        else:
            stab_x_indices = list(np.where(stab_x != 0)[0] + 1)
        stabilizers_x_indices.append(stab_x_indices)

    for stab_x_indices in stabilizers_x_indices:
        data_rec_targets = [stim.target_rec(-(final_data1_indices+final_data2_indices)[::-1].index(q) - 1) for q in stab_x_indices]
        full_circuit.append_operation("DETECTOR", data_rec_targets, list(index_to_coord_PQRM[stab_x_indices[0]]) + [1])
        # Add the 4th coordinate for postselection
    
    # Logical observables
    if PQRM_state == "Z": # Z gauge + logical Z in surface code
        log_ob_gauge_indices = [sf_code_extended.ancilla_coords[::-1].index(anc) for anc in sydrome_coords_ancsys_z if anc[0] < 0]
        log_ob_gauge_records = [stim.target_rec(-k-1-len(final_data1+final_data2)) for k in log_ob_gauge_indices]
        log_ob_data = [(-2,2*y+2) for y in range(d)]
        log_ob_data_indices = [final_data2[::-1].index(data) for data in log_ob_data]
        log_ob_data_records = [stim.target_rec(-k-1) for k in log_ob_data_indices]
        full_circuit.append_operation("OBSERVABLE_INCLUDE", log_ob_gauge_records + log_ob_data_records, 0)
    
    if PQRM_state == "X":
        # Logical observale: merged logical X of patch 1 and 2
        log_ob_x_sf = [(-2*x,2) for x in range(1, d+1)]
        log_ob_x_PQRM = [(2*x,2*y) for x in range(x_range_PQRM//2) for y in range(y_range_PQRM)]
        log_ob_x_PQRM.remove((0,0))
        log_ob_x_middle = [(-1,4*y+5) for y in range((log_PQRM_len + 1)//2-1)]
        log_ob_data_indices_sf = [(final_data1+final_data2)[::-1].index(data) for data in log_ob_x_sf]
        log_ob_data_indices_PQRM = [(final_data1+final_data2)[::-1].index(data) for data in log_ob_x_PQRM]
        log_ob_middle_indices = [(final_data1+final_data2)[::-1].index(data) for data in log_ob_x_middle]
        log_ob_data_records = [stim.target_rec(-k-1) for k in log_ob_data_indices_sf + log_ob_middle_indices + log_ob_data_indices_PQRM]
        full_circuit.append_operation("OBSERVABLE_INCLUDE", log_ob_data_records, 0)

    # probably pass the detector ranges that will be modified..., the measurement results corresponding to the PQRM
    
    return full_circuit

# %%
if __name__ == "__main__":
    tasks = []
    rx = 1
    rz = 3
    m = 5
    d_surf = 7
    rounds = d_surf 
    noise = 0.001
    noise_profile = [1e-6, noise, noise, noise]
    for PQRM_state in ["Z", "X"]:
        circuit = CrossLS([rx, rz, m], d_surf, rounds, PQRM_state, noise_profile, punctured=True)
        tasks.append(sinter.Task(circuit=circuit, json_metadata={'d_surf': d_surf, 'PQRM_para': [rx,rz,m], 'rounds': rounds, 'PQRM_state': PQRM_state, 'p1': noise_profile[0], 'p2': noise_profile[1], 'pM': noise_profile[2], 'pR': noise_profile[3]}))

    collected_stats: List[sinter.TaskStats] = sinter.collect(
        num_workers=6,
        tasks=tasks,
        decoders=['bposd'],
        max_shots=1000, # change to large number for full scale experiment
        max_errors=10000, # Whichever first hits the threshold stops the simulation
        custom_decoders=sinter_decoders(),
        print_progress=True,
    )

    # Convert collected_stats into Dataframe
    records = []
    for stat in collected_stats:
        record = {
            'decoder': stat.decoder,
            'shots': stat.shots,
            'errors': stat.errors,
            'seconds': stat.seconds,
            'LER': convert_to_per_round(stat.errors / stat.shots, 1) if stat.shots > 0 else None, # don't average
            'rounds': rounds,
            'd_surf': stat.json_metadata['d_surf'],
            'PQRM_state': stat.json_metadata['PQRM_state'],
            'p1': stat.json_metadata['p1'], 'p2': stat.json_metadata['p2'], 'pM': stat.json_metadata['pM'], 'pR': stat.json_metadata['pR']
        }
        records.append(record)
    df_sf = pd.DataFrame(records)
    df_sf.to_csv("results_729.csv", index=False) # If you wanna save the experiment data to csv
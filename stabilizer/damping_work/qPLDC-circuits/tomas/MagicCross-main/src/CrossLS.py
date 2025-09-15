# No initialization noise on surface code and PQRM code

import stim
from planar_PQRM import PlanarSurfaceCode
from typing import Tuple, List
from itertools import product
import copy
import numpy as np

# Utility functions for bit manipulation
bin_wt = lambda i: bin(i)[2:].count('1')
bit_rev = lambda t, m: int(bin(t)[2:].rjust(m, '0')[::-1], 2)
int2bin = lambda i, m: [int(c) for c in bin(i)[2:].rjust(m, '0')]
bin2int = lambda l: int(''.join(map(str, l)), 2)

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
    syndrome_coords_PQRM_bound_dict = {(4,4): [(3,7),(7,3)], (4,8):[(3,15),(7,3),(7,5),(7,7),(7,11)],(8,8):[(3,15),(5,15),(7,15),(11,15),(15,3),(15,5),(15,7),(15,11)]}
    syndrome_coords_PQRM_bound = syndrome_coords_PQRM_bound_dict[(x_range_PQRM, y_range_PQRM)]

    # Notice that these syndromes give the full set of Z-stabilizers
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


    # Initialize the PQRM code in |0> state with noise
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
    
    full_circuit = full_circuit.without_noise() # if you wanna remove the noise on PQRM initialization as well

    # Lattice surgery with PQRM code and surface code
    # Initialize the ancilla system for lattice surgery with noise
    sydrome_coords_ancsys_indices = [coord_to_index[c] for c in sydrome_coords_ancsys]
    data_coords_ancsys_indices = [coord_to_index[c] for c in data_coords_ancsys]
    full_circuit.append_operation("R", sydrome_coords_ancsys_indices)
    full_circuit.append_operation("X_ERROR", sydrome_coords_ancsys_indices, p_R)
    full_circuit.append_operation("RX", data_coords_ancsys_indices)
    full_circuit.append_operation("Z_ERROR", data_coords_ancsys_indices, p_R)
    full_circuit.append("TICK")

    # Incorporate the ancilla system and PQRM into the surface code, build unified SM circuits
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
    if PQRM_state == "Z": # Z gauge (+ some Z stabilizer of PQRM) + logical Z in surface code
        log_ob_gauge_indices = [sf_code_extended.ancilla_coords[::-1].index(anc) for anc in sydrome_coords_ancsys_z if anc[0] < 0]
        log_ob_gauge_records = [stim.target_rec(-k-1-len(final_data1+final_data2)) for k in log_ob_gauge_indices]
        log_ob_data = [(-2,2*y+2) for y in range(d)]
        log_ob_data_indices = [final_data2[::-1].index(data) for data in log_ob_data]
        log_ob_data_records = [stim.target_rec(-k-1) for k in log_ob_data_indices]
        if y_range_PQRM == 8:
            log_ob_stab_z = [(2*x+1,y) for x in range(x_range_PQRM-1) for y in [9,13]] + [(2*x_range_PQRM-1,11)]
            log_ob_stab_z_indices = [sf_code_extended.ancilla_coords[::-1].index(anc) for anc in log_ob_stab_z]
            log_ob_stab_z_records = [stim.target_rec(-k-1-len(final_data1+final_data2)) for k in log_ob_stab_z_indices]
            log_ob_gauge_records += log_ob_stab_z_records
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
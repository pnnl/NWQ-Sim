from typing import Tuple, Set, Dict, List
import stim

# classes of surface codes

class PlanarSurfaceCode:
    def __init__(self, d: int, shift: list[int] = [0,0]):
        
        self.d = d

        self.data_coords = set()
        self.ancilla_coords_x = set()
        self.ancilla_coords_z = set()
        self.stabilizers = []
        self.logical_X = []
        self.logical_Z = []

        self._generate_planar_surface_code()
        if shift != [0,0]:
            self.patch_shift(shift)
    

    def _generate_planar_surface_code(self):
        # maybe can extend it to m*n surface code patch, default 1*1
        
        d = self.d

        for y in range(2 * d - 1):
            if y % 2 == 0: # Even rows Z-stabilizers
                for x in range(d-1):
                    anc = (2 * x + 1, y)
                    self.ancilla_coords_z.add(anc)
                    if y == 0: # top row Z-stabilizers
                        self.stabilizers.append({'type': 'Z', 'ancilla': anc, 'data': [(2 * x, y), (2 * x + 2, y), (2 * x + 1, y + 1)]})
                    elif y == 2 * d - 2: # bottom row Z-stabilizers
                        self.stabilizers.append({'type': 'Z', 'ancilla': anc, 'data': [(2 * x, y), (2 * x + 2, y), (2 * x + 1, y - 1)]})
                    else: # Other Z-stabilizers
                        self.stabilizers.append({'type': 'Z', 'ancilla': anc, 'data': [(2 * x, y), (2 * x + 2, y), (2 * x + 1, y - 1), (2 * x + 1, y + 1)]})
                [self.data_coords.add((2 * x, y)) for x in range(d)]
            else: # Odd rows X-stabilizers
                for x in range(d):
                    anc = (2 * x, y)
                    self.ancilla_coords_x.add(anc)
                    if x == 0:
                        self.stabilizers.append({'type': 'X', 'ancilla': anc, 'data': [(2 * x + 1, y), (2 * x, y + 1), (2 * x, y - 1)]})
                    elif x == d - 1:
                        self.stabilizers.append({'type': 'X', 'ancilla': anc, 'data': [(2 * x - 1, y), (2 * x, y + 1), (2 * x, y - 1)]})
                    else:
                        self.stabilizers.append({'type': 'X', 'ancilla': anc, 'data': [(2 * x - 1, y), (2 * x + 1, y), (2 * x, y + 1), (2 * x, y - 1)]})
                [self.data_coords.add((2 * x + 1, y)) for x in range(d - 1)]

        self.qubit_coords = PlanarSurfaceCode.sort_coords(self.data_coords | self.ancilla_coords_x | self.ancilla_coords_z)
        self.data_coords = PlanarSurfaceCode.sort_coords(self.data_coords)
        self.ancilla_coords = PlanarSurfaceCode.sort_coords(self.ancilla_coords_x | self.ancilla_coords_z)
        self.ancilla_coords_x = PlanarSurfaceCode.sort_coords(self.ancilla_coords_x)
        self.ancilla_coords_z = PlanarSurfaceCode.sort_coords(self.ancilla_coords_z)
        
        self.coord_to_index = {coord: i for i, coord in enumerate(self.qubit_coords)}
        self.index_to_coord = {i: coord for coord, i in self.coord_to_index.items()}

        self.logical_Z.append([(0, 2 * y) for y in range(d)])
        self.logical_X.append([(2 * x, 0) for x in range(d)])

    def get_info(self):
        # The returned coordinates are all sorted
        return {
            'qubit_coords': self.qubit_coords,
            'data_coords': self.data_coords,
            'ancilla_coords': self.ancilla_coords,
            'ancilla_coords_x': self.ancilla_coords_x,
            'ancilla_coords_z': self.ancilla_coords_z,
            'stabilizers': self.stabilizers,
            'coord_to_index': self.coord_to_index,
            'index_to_coord': self.index_to_coord,
            'logical_Z': self.logical_Z,
            'logical_X': self.logical_X
        }
    
    def patch_shift(self, coords_diff: list[int]):
        for key in self.get_info():
            if key in ['qubit_coords', 'data_coords', 'ancilla_coords', 'ancilla_coords_x', 'ancilla_coords_z']:
                new_value = PlanarSurfaceCode.shift_coords(getattr(self, key), coords_diff[0], coords_diff[1])
                setattr(self, key, new_value)
            elif key == 'stabilizers':
                new_value = [{'type': stab['type'], 'ancilla': PlanarSurfaceCode.shift_coord(stab['ancilla'], coords_diff[0], coords_diff[1]), 
                                       'data': PlanarSurfaceCode.shift_coords(stab['data'], coords_diff[0], coords_diff[1])} for stab in getattr(self, key)]
                setattr(self, key, new_value)
            elif key in ['logical_Z', 'logical_X']:
                new_value = [PlanarSurfaceCode.shift_coords(logical_op, coords_diff[0], coords_diff[1]) for logical_op in getattr(self, key)]
                setattr(self, key, new_value)
            else:
                continue
        
        self.coord_to_index = {coord: i for i, coord in enumerate(self.qubit_coords)}
        self.index_to_coord = {i: coord for coord, i in self.coord_to_index.items()}

    def get_parity_check_matrix(self):
        pass
    
    def build_standard_sm_round(self, noise_profile: list, data_error = True, code_capacity = False) -> stim.Circuit:
        
        p1, p2, p_M, p_R = noise_profile

        # This one-round circuit only contains operations and noise, without detectors
        circuit = stim.Circuit()

        # Before a sm round, insert data qubit errors
        if data_error:
            data_indices = [self.coord_to_index[c] for c in self.data_coords]
            circuit.append_operation("X_ERROR", data_indices, p1)

        # TICK 1: Hadamard on X ancillas
        H_targets = [self.coord_to_index[coord] for coord in self.ancilla_coords_x]
        circuit.append_operation("H", H_targets)
        if not code_capacity:
            circuit.append_operation("X_ERROR", H_targets, p1)
        circuit.append("TICK")

        # TICKs 2–7: CNOT gate scheduling, following Li's paper
        # for different gate schedulings, change the tick_deltas
        tick_deltas = [((-1, 0), (0, 0)), ((+1, 0), (0, 0)), ((0, +1), (0, +1)), ((0, -1), (0, -1)), ((0, 0), (-1, 0)), ((0, 0), (+1, 0))]

        for dx_z, dx_x in tick_deltas:
            CNOT_pairs = []
            for ancilla_x in self.ancilla_coords_x:
                anc = self.coord_to_index[ancilla_x]
                data_coord = (ancilla_x[0] + dx_x[0], ancilla_x[1] + dx_x[1])
                if data_coord in self.data_coords:
                    data = self.coord_to_index[data_coord]
                    circuit.append_operation("CNOT", [anc, data])
                    CNOT_pairs += [anc, data]
            for ancilla_z in self.ancilla_coords_z:
                anc = self.coord_to_index[ancilla_z]
                data_coord = (ancilla_z[0] + dx_z[0], ancilla_z[1] + dx_z[1])
                if data_coord in self.data_coords:
                    data = self.coord_to_index[data_coord]
                    circuit.append_operation("CNOT", [data, anc])
                    CNOT_pairs += [data, anc]
            if not code_capacity:
                circuit.append_operation("X_ERROR", CNOT_pairs, p2)
            circuit.append("TICK")

        # TICK 8: Hadamard on X ancillas
        circuit.append_operation("H", H_targets)
        if not code_capacity:
            circuit.append_operation("X_ERROR", H_targets, p1)
        circuit.append("TICK")

        # TICK 9: MR with SPAM noise
        measure_targets = [self.coord_to_index[coord] for coord in self.ancilla_coords]
        if not code_capacity:
            circuit.append_operation("X_ERROR", measure_targets, p_M)
        circuit.append_operation("MR", measure_targets)
        if not code_capacity:
            circuit.append_operation("X_ERROR", measure_targets, p_R)
        # Don't insert TICK here, insert TICK after specifying the detectors in the construction of the full circuit

        return circuit
    
    def build_full_surface_code_circuit(self, rounds: int, noise_profile: list, observable_type: str, data_error = True, code_capacity = False) -> stim.Circuit:
    
        p1, p2, p_M, p_R = noise_profile

        full_circuit = stim.Circuit()
        repeat_circuit = stim.Circuit()

        # QUBIT_COORDS annotations
        for coord, index in self.coord_to_index.items():
            full_circuit.append_operation("QUBIT_COORDS", [index], list(coord))

        # Initialization
        data_indices = [self.coord_to_index[c] for c in self.data_coords]
        ancilla_indices = [self.coord_to_index[c] for c in self.ancilla_coords]
        if observable_type == "Z":
            full_circuit.append_operation("R", data_indices)
            if not code_capacity:
                full_circuit.append_operation("X_ERROR", data_indices, p_R)
        else:
            full_circuit.append_operation("RX", data_indices)
            if not code_capacity:
                full_circuit.append_operation("Z_ERROR", data_indices, p_R)
        
        full_circuit.append_operation("R", ancilla_indices)
        if not code_capacity:
            full_circuit.append_operation("X_ERROR", ancilla_indices, p_R)
        full_circuit.append("TICK")

        # First round
        full_circuit += self.build_standard_sm_round(noise_profile, data_error)
        # the first round, detectors are the stabilizers already created by initialization
        for k, ancilla in enumerate(reversed(self.ancilla_coords)):
            if observable_type == "Z":
                if ancilla in self.ancilla_coords_z:
                    full_circuit.append_operation("DETECTOR", [stim.target_rec(-k - 1)], list(ancilla) + [0])
            else:
                if ancilla in self.ancilla_coords_x:
                    full_circuit.append_operation("DETECTOR", [stim.target_rec(-k - 1)], list(ancilla) + [0])

        # Later rounds
        repeat_circuit.append("TICK")
        repeat_circuit += self.build_standard_sm_round(noise_profile, data_error, code_capacity=code_capacity)
        repeat_circuit.append_operation("SHIFT_COORDS", [], [0,0,1])
        # Insert detectors
        for k, ancilla in enumerate(reversed(self.ancilla_coords)):
            # ancilla_index = self.coord_to_index[ancilla]
            prev = -k - 1 - len(self.ancilla_coords)
            repeat_circuit.append_operation("DETECTOR", [stim.target_rec(-k - 1), stim.target_rec(prev)], list(ancilla) + [0])

        full_circuit += repeat_circuit * (rounds-1)

        # Final measurement of data qubits in Z or X basis
        full_circuit.append("TICK")
        if observable_type == "X":
            if not code_capacity:
                full_circuit.append_operation("Z_ERROR", data_indices, p_M)
            full_circuit.append_operation("MX", data_indices),
        else: # if it's "Z"
            if not code_capacity:
                full_circuit.append_operation("X_ERROR", data_indices, p_M)
            full_circuit.append_operation("M", data_indices)

        # Final detectors
        for stab in self.stabilizers:
            if stab['type'] == observable_type:
                anc = stab['ancilla']
                try:
                    ancilla_index = self.ancilla_coords[::-1].index(anc)
                    anc_rec = stim.target_rec(-ancilla_index - 1 - len(self.data_coords))
                    data_rec_targets = [stim.target_rec(-self.data_coords[::-1].index(q) - 1) for q in stab['data']]
                    full_circuit.append_operation("DETECTOR", [anc_rec] + data_rec_targets, list(anc) + [1])
                except ValueError:
                    continue
        
        # logical observables
        logical_key = 'logical_' + observable_type
        for l, logical_data in enumerate(getattr(self, logical_key)):
            data_indices = [self.data_coords[::-1].index(data) for data in logical_data]
            full_circuit.append_operation("OBSERVABLE_INCLUDE", [stim.target_rec(-k-1) for k in data_indices], l)

        return full_circuit
    
    @staticmethod
    def sort_coords(coords: Set[Tuple[int, int]]) -> List[Tuple[int, int]]:
        return sorted(coords, key=lambda c: (c[1], c[0]))
    
    @staticmethod
    def shift_coords(coords: List[Tuple[int, int]], dx: int, dy: int) -> List[Tuple[int, int]]:
        return [(x + dx, y + dy) for (x, y) in coords]
    
    @staticmethod
    def shift_coord(coord: Tuple[int, int], dx: int, dy: int) -> Tuple[int, int]:
        x, y = coord
        return (x + dx, y + dy)


class MultiPatchSurfaceCode(PlanarSurfaceCode):
    def __init__(self, code_list: List[PlanarSurfaceCode]):
        import copy
        self.code_list = [copy.deepcopy(code) for code in code_list]

        self.data_coords = []
        self.ancilla_coords_x = []
        self.ancilla_coords_z = []
        self.stabilizers = []
        self.coord_to_index = {}
        self.logical_X = []
        self.logical_Z = []
        self.gauge_X = []
        self.gauge_Z = []

        # combine layout
        self._patch_combine()
    
    def _patch_combine(self):
        offset = 0
        for code in self.code_list:
            # combine all the info
            self.data_coords += code.data_coords
            self.ancilla_coords_x += code.ancilla_coords_x
            self.ancilla_coords_z += code.ancilla_coords_z
            self.ancilla_coords = self.ancilla_coords_x + self.ancilla_coords_z
            self.qubit_coords = self.data_coords + self.ancilla_coords
            self.stabilizers += code.stabilizers
            self.logical_X += code.logical_X
            self.logical_Z += code.logical_Z

            # update shift index
            coord_to_index_shifted = {coord: i + offset for coord, i in code.coord_to_index.items()}
            self.coord_to_index.update(coord_to_index_shifted)
            self.index_to_coord = {i: coord for coord, i in self.coord_to_index.items()}
            offset = max(self.coord_to_index.values()) + 1

    def build_standard_sm_round(self, noise_profile: List[float], code_capacity: bool = False) -> stim.Circuit:
        # Reuse the parent class's method
        return super().build_standard_sm_round(noise_profile=noise_profile, code_capacity=code_capacity)
    

    def build_full_surface_code_circuit(self, rounds: int, noise_profile: List[float], observable_type: str, code_capacity = False) -> stim.Circuit:
        # Reuse the parent class's method
        return super().build_full_surface_code_circuit(rounds=rounds, noise_profile=noise_profile, observable_type=observable_type, code_capacity=code_capacity)
    
    
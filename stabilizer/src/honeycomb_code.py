''' some qec circuit generation '''
# from helper.tab_state import Tableau
from typing import Union, List
import numpy as np
import stim
#from abc import abstractmethod
import math
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon

class NotQubitLocationError(Exception):
    """Raised when the location is not a qubit"""
    pass

class NoQubitActiveError(Exception):
    '''Raised when there is no qubits being activated in the code'''
    pass

class ConsistencyError(Exception):
    pass

class Qubit:
    def __init__(self):
        self.location = (0, 0, 0)
        self.index = 0
        self.active = False
    
    def __str__(self):
        if self.active:
            output = f'active qubit{self.location} {self.index}'
        else:
            output = f'inactivate qubit{self.location} {self.index}'
        
        return output



def check_qubits(location, ts, q_type='x'):
    """
    Return the location of the ts-th neighboring data qubit
    around the ancilla located at `location`.

    Parameters
    ----------
    location : tuple
        Coordinates of the ancilla qubit (x, y).
    ts : int
        Index from 0 to 5 indicating which of the 6 neighboring data qubits.
    q_type : str
        'x' or 'z' — could affect direction symmetry (optional, here unused).

    Returns
    -------
    tuple
        (x, y) location of the neighboring data qubit.
    """

    #Directions around a regular hexagon (pointy-top orientation)
    #These offsets assume the hex center is at (x, y) and each side is unit length.
    
    if q_type == 'x':
        directions = [
        (+2.2, 0),                   # Right
        (+1.2, +2),         # Top-right
        (-.8, +2),         # Top-left
        (-1.8, 0),                   # Left
        (-.8, -2),         # Bottom-left
        (+1.2, -2),         # Bottom-right
        ]
    elif q_type == 'z':
        directions = [
        (+1.8, 0),                   # Right
        (+.8, +2),         # Top-right
        (-1.2, +2),         # Top-left
        (-2.2, 0),                   # Left
        (-1.2, -2),         # Bottom-left
        (+.8, -2),         # Bottom-right
        ]
    else:
        raise ValueError(f"Unsupported q_type: {q_type}")
    
    offset = directions[ts % 6]
    x, y = location
    dx, dy = offset

    x = round(x + dx,1) #integer-float issue
        
    return (x, y + dy)


class HoneycombQubits:
    def __init__(self):
        self.qubit_index = {} #indices mapped by location for all qubits {(x,y)->i, (x,y)->i...}
        self.qubit_location = {} #locations mapped by index for all qubits {i->(x,y), i->(x,y)...}

        self.data_qubits = [] #locations of only data qubits [[x,y],[x,y]]
        self.x_qubits = [] #locations+colors of only x check ancillas [((x,y),color), ((x,y),color)...]
        self.z_qubits = [] #locaitons+colors of only z check ancillas [((x,y),color), ((x,y),color)...]

        self.n_x_checks = 0
        self.n_z_checks = 0
        self.n_data = 0
        

        self.x_indices = [] #[1,2,3,4...]
        self.z_indices = [] #[1,2,3,4...]
        self.data_indices = [] #[1,2,3,4...]

        self.logical_locations = [] #[(i, (x,y)), (i, (x,y))...]
        self.x_neighbours = {} #x_indices: [r0_data_index, r1_data_index, r2_data_index, r3_data_index]
        self.z_neighbours = {} #z_indices: [r0_data_index, r1_data_index, r2_data_index, r3_data_index]

        self.bell_pairs = []

        self.hex_side_length = math.sqrt(5)

        
        self.distance = None
        self.off = ''

        return 
    
    def _reset(self):
        self.qubit_index = {} #indices mapped by location for all qubits {(x,y)->i, (x,y)->i...}
        self.qubit_location = {} #locations mapped by index for all qubits {i->(x,y), i->(x,y)...}

        self.data_qubits = [] #locations of only data qubits [[x,y],[x,y]]
        self.x_qubits = [] #locations+colors of only x check ancillas [((x,y),color), ((x,y),color)...]
        self.z_qubits = [] #locaitons+colors of only z check ancillas [((x,y),color), ((x,y),color)...]

        self.n_x_checks = 0
        self.n_z_checks = 0
        self.n_data = 0
        

        self.x_indices = [] #[1,2,3,4...]
        self.z_indices = [] #[1,2,3,4...]
        self.data_indices = [] #[1,2,3,4...]

        self.logical_locations = [] #[(i, (x,y)), (i, (x,y))...]
        self.x_neighbours = {} #x_indices: [r0_data_index, r1_data_index, r2_data_index, r3_data_index]
        self.z_neighbours = {} #z_indices: [r0_data_index, r1_data_index, r2_data_index, r3_data_index]

        self.bell_pairs = []

        self.hex_side_length = math.sqrt(5)

        
        self.distance = None
        self.off = ''

        return 
    
    # def get_observables(self, shots, circuit, postselect=False):
    #     sampler = circuit.compile_detector_sampler()
    #     results = sampler.sample(shots)
    #     observables = [False] * shots

    #     for shot in range(shots):
    #         # print(results[shot])
    #         for i in range(len(self.logical_locations)):
    #             # print([-1*i-1])
    #             # print(results[shot][-1*i-1])
    #             observables[shot] ^= results[shot][-1*i-1]

    #     # if postselect:
            

    #     return observables
    
    def inject_and_measure(self, hook_opt=False, injection='0', rounds=1, method='hex', logical_basis='z', detector_reset=True):
        """
        Create the full circuit for injection, syndrome extraction, and logical observable basis.

        Parameters
        ----------
        injection : str
            Injected state.
        rounds : int
            Number of rounds to repeat the syndrome measurement.
        method : int
            Syndrome extraction method.
            - 'hex', 'willow'
        logical_basis: str
            Basis to measure the observable in.
        output_format : str
            'str', 'list', or 'stim' to return a stim.Circuit object.

        Returns
        -------
        str, list, or stim.Circuit
            The syndrome check block in the specified format.
        """

        circuit = stim.Circuit()

        if method.lower() == 'hex':
            injection_circuit = self.arb_state_injection_bell_hex(output_format='str', inj=injection, detector_reset=detector_reset)
            circuit.append_from_stim_program_text(injection_circuit)

            syndrome_circuit = self.syndrome_check_block_hex(round=rounds, hook_opt=hook_opt, output_format='str', detector_reset=detector_reset)
            circuit.append_from_stim_program_text(syndrome_circuit)

            readout_circuit = self.readout_block_hex(output_format='str', basis=logical_basis, detector_reset=detector_reset)
            circuit.append_from_stim_program_text(readout_circuit)
        
        return circuit
        

    def select_active(self, distance=3):
        """
        Set up a 6.6.6 triangular honeycomb color code layout.

        - Qubits are arranged in honeycomb hexagons bounded by a triangle, which creates some trapezoids on the boundary
        - Each hexagon (trapezoid) contains 6(4) data qubits 
        - One X and one Z ancilla qubit sit inside each honeycomb (trapezoid)
        - ** Numbering prioritizes left first, then prioritizes lowest second**

        Parameters
        ----------
        distance : int
            Number of stabilizer faces along one axis.
        """
        if self.distance is not None:
            self._reset()
        self.distance = distance

        self.data_qubits = set()
        self.logical_locations = set()

        color = 'red'
        
        for i in range(0, distance+(distance//2), 1): #base of the triangle
            for j in range(i+1): #grows the pyramid ontop

                if (i+j)%3 == 1:
                    if 2*j%6 == 0:
                        color = 'red'
                    if 2*j%6 == 2: 
                        color = 'green'
                    if 2*j%6 == 4: 
                        color = 'blue'

                    self.x_qubits.append(((round(2*(i-.1-j*.5),1), 2*j), color))#integer-float issue
                    self.z_qubits.append(((round(2*(i+.1-j*.5),1), 2*j), color))#integer-float issue
                else:
                    data_qubit = (2*(i-j*.5), 2*j)
                    self.data_qubits.add(data_qubit)
                    if j == i:
                        self.logical_locations.add(data_qubit)
                j-=1

        self.data_qubits = np.array(sorted(self.data_qubits))
        self.n_data = len(self.data_qubits)
        self.n_x_checks = len(self.x_qubits)
        self.n_z_checks = len(self.z_qubits)

        self._process()

        return

    def _process(self):
        """
        Handles the numbering of qubits once they've been selected

        Parameters
        ----------
        N/A
        """
        if not (self.n_x_checks or self.n_z_checks or self.n_data):
            raise NoQubitActiveError(
                f'_process called without any active qubits: '
                f'nx: {self.n_x_checks}, nz: {self.n_z_checks}, nd: {self.n_data}'
            )

        self.qubit_index = {}       # Map: (x, y) → index
        self.qubit_location = {}    # Map: index → (x, y)
        self.data_indices = []      # List of data qubit indices
        self.x_indices = []         # List of x ancilla indices
        self.z_indices = []         # List of z ancilla indices
        self.x_neighbours = {}      # Map: x_index → [data indices]
        self.z_neighbours = {}      # Map: z_index → [data indices]
        logical_indices = []

        ct = 0

        self.logical_locations = sorted(self.logical_locations, key=lambda x: x[0])

        # --- Index Data Qubits ---
        for loc in self.data_qubits:
            loc_t = tuple(loc)
            self.qubit_index[loc_t] = ct
            if loc_t in self.logical_locations:
                logical_indices.append(ct)
            self.qubit_location[ct] = loc_t
            self.data_indices.append(ct)
            ct += 1

        self.logical_locations = list(zip(logical_indices, self.logical_locations))

        # --- Index X Ancilla Qubits ---
        for loc, color in self.x_qubits:
            loc_t = tuple(loc)
            self.qubit_index[loc_t] = ct
            self.qubit_location[ct] = loc_t
            self.x_indices.append(ct)
            ct += 1
            

        # --- Index Z Ancilla Qubits ---
        for loc, color in self.z_qubits:
            loc_t = tuple(loc)
            self.qubit_index[loc_t] = ct
            self.qubit_location[ct] = loc_t
            self.z_indices.append(ct)
            ct += 1

        # --- Find Neighbors for X Checks ---
        for anc_idx, (loc, color) in zip(self.x_indices, self.x_qubits):
            # print(loc)
            # print(anc_idx)
            nearby_data = []
            for ts in range(6):  #6 directions around the hex
                # print(loc)
                loc_data = check_qubits(loc, ts, q_type='x')
                # print(loc_data)
                if tuple(loc_data) in self.qubit_index:
                    nearby_data.append(self.qubit_index[tuple(loc_data)])
                else:
                    nearby_data.append(-1)  #Placeholder for boundary qubit (nothing there)
            self.x_neighbours[anc_idx] = tuple(nearby_data)

        # --- Find Neighbors for Z Checks ---
        for anc_idx, (loc, color) in zip(self.z_indices, self.z_qubits):
            # print(loc)
            # print(anc_idx)
            nearby_data = []
            for ts in range(6):
                # print(loc)
                loc_data = check_qubits(loc, ts, q_type='z')
                # print(loc_data)
                if tuple(loc_data) in self.qubit_index:
                    nearby_data.append(self.qubit_index[tuple(loc_data)])
                else:
                    nearby_data.append(-1)
            self.z_neighbours[anc_idx] = tuple(nearby_data)

        '''
        
        Bell pair setup
        
        '''

        #Find the bell pairs 
        self._gen_bell_pairs()

        return
    
    
    
    def _gen_bell_pairs(self):
        '''Generate the Bell pairs for state injection.'''

        #Red hex stores all of the qubits around each red hexagon, grouped by the specific hexagon. Since the code is 3-colorable, red hexagons never share qubits.
        #This does the same neighbor find that filled x_neighbours and z_neighbours
        red_hex = []        
        for z_coord, color in self.z_qubits:
            if color == 'red':
                hex_verts = []
                for i in range(6):
                    loc = check_qubits(z_coord, i, 'z')
                    if tuple(loc) in self.qubit_index:
                        qubit = self.qubit_index[tuple(loc)]
                        hex_verts.append(qubit)
                red_hex.append(hex_verts)
        
        #The red boundary should needs to be included too, since there aren't proper red hexagons on one side of the triangle
        flat_red_hex = {x for sublist in red_hex for x in sublist}
        right_boundary = [x for x in self.data_indices if x not in flat_red_hex]
        red_hex.append(right_boundary)
        
        
        #Map the logical operator and remove it from the red verticies. 
        #The logical operator is a special case handled independently, since the bell pairs cross the middle of a hexagon.
        for i in range(1, len(self.logical_locations), 2):
            first, _ = self.logical_locations[i]
            second, _ = self.logical_locations[i+1]
            self.bell_pairs.append((first, second))
            for sublist in red_hex:
                if 0 in sublist:
                    sublist.remove(0)
                if first in sublist:
                    sublist.remove(first)
                if second in sublist:
                    sublist.remove(second)
                    

        def distance(a, b):
            x1, y1 = self.qubit_location[a]
            x2, y2 = self.qubit_location[b]
            return math.sqrt(abs(x1 - x2)**2 + abs(y1 - y2)**2)
        
        #Map each value to its group index
        value_to_group = {}
        for i, group in enumerate(red_hex):
            for val in group:
                value_to_group[val] = i

        all_values = [val for group in red_hex for val in group]
        unpaired = set(all_values)

        while unpaired:
            paired_this_round = set()
            for a in list(unpaired):
                if a in paired_this_round:
                    continue
                group_a = value_to_group[a]

                #Find the unique b qubit to pair from a different red group. It is always 1 hexagon side length away.
                found_b = None
                for b in unpaired:
                    if b == a or b in paired_this_round:
                        continue
                    if value_to_group[b] == group_a:
                        continue
                    if distance(a,b) <= self.hex_side_length:
                        found_b = b
                        break

                if found_b is not None:
                    self.bell_pairs.append((a, found_b))
                    paired_this_round.add(a)
                    paired_this_round.add(found_b)

            #Remove paired points from unpaired
            unpaired -= paired_this_round

            #If no pairs were formed this round, break to avoid infinite loop
            if not paired_this_round:
                break

        # print("Pairs:", self.bell_pairs)
        

    def draw(self, fontsize: int = 12, show: bool = True, off = ''):
        """
        Visualise the currently‑selected patch.

        Parameters
        ----------
        fontsize : int
            Font size used for the 'D', 'X', 'Z' labels.
        show : bool
            If True (default) call plt.show() at the end.
        off : str
            Used to turn off a color (i.e. red) or type of check (X or Z)
        """

        if not (self.n_x_checks or self.n_z_checks or self.n_data):
            raise NoQubitActiveError("No qubits selected: call select_active() first")
        
        coords_D = [tuple(q) for q in self.data_qubits]
        coords_X = [tuple(q) for q, _ in self.x_qubits]
        coords_Z = [tuple(q) for q, _ in self.z_qubits]

        avg_coords = [
            ((x + z) / 2, (y + w) / 2)
            for (x, y), (z, w) in zip(coords_X, coords_Z)
        ]
        Z_colors = [color for _, color in self.z_qubits]

        def _plot(coords, label_char, scatter_kwargs, legend_label=None):
            if not coords:
                return
            xs, ys = zip(*coords)

            if legend_label:
                plt.scatter(xs, ys, **scatter_kwargs, label=legend_label)
            else:
                plt.scatter(xs, ys, **scatter_kwargs)
            for x, y in coords:
                plt.text(x, y, f"{self.qubit_index[(x, y)]}",
                        ha="center", va="center", fontsize=fontsize)
                
        for (x, y), color in zip(avg_coords, Z_colors):
            if self.off.lower() != color:
                hex = RegularPolygon(
                    (x, y),              #center
                    numVertices=6,
                    radius=2,   #controls size
                    orientation=np.radians(30),  #flat
                    facecolor=color,
                    edgecolor=None,
                    alpha=0.4
                )
                plt.gca().add_patch(hex)
            

        _plot(coords_D, "D", dict(marker="s", s=120, color="white", edgecolors="k"), legend_label="Data")
        if (self.off != 'x') and (self.off != 'X'):
            _plot(coords_X, "X", dict(marker="o", s=120, color="tab:red", edgecolors="k"), legend_label="X Ancilla")
        if (self.off != 'z') and (self.off != 'Z'):
            _plot(coords_Z, "Z", dict(marker="^", s=120, color="tab:blue", edgecolors="k"), legend_label="Z Ancilla")

        
        first = True
        #Plot bell pairs
        if self.bell_pairs:
            for a, b in self.bell_pairs:
                x1, y1 = self.qubit_location[a]
                x2, y2 = self.qubit_location[b]
                plt.plot([x1, x2], [y1, y2], '-', color='gold', linewidth=3, markersize=8, label='Bell States' if first else "")
                first = False

        first = True
        #Plot logical op

        if self.logical_locations:
            prev_x, prev_y = None, None
            for qubit, _ in self.logical_locations:
                x, y = self.qubit_location[qubit]
                if prev_x is not None:
                    plt.plot([prev_x, x], [prev_y, y], 'k:', linewidth=2, label='Logical X, Z' if first else "")
                    first = False
                prev_x, prev_y = x, y
        

        plt.gca().set_aspect("equal", adjustable="box")
        plt.grid(True, linestyle="--", linewidth=0.5)
        plt.xlabel("x  (rightwards)")
        plt.ylabel("y  (downwards)")
        plt.title("6.6.6 Triangular Honeycomb Code")
        plt.legend()

        if show:
            plt.show()

    
    def syndrome_check_block_hex(self, round=3, hook_opt=False, output_format='stim', detector_reset=True):
        """
        Generate the syndrome check blocks for a non-local connection. Turn off one type of ancilla qubit. 
        Every remaining ancilla qubit connects to every data qubit for syndrome measurment.

        Order of CNOT extraction cycle to minimize hook errors after post-selection GIVEN our injection site 
        is in the bottom left and the triangle is left handed: (perform the first syndrome extraction cnot at the 1 position, f in the paper)
            USTC: North vertex and left handed
                - Clockwise, starting from north-east, for X-check: 5-6-3-2-1-4
                - Clockwise, starting from north-east, for Z-check: 4-7-6-5-2-3 (1 is H on the ancilla qubit to init)
            Us (120 degree rotation counterclockwise vs USTC): Southwest vertex and left handed (j in the paper)
                - Clockwise, starting from north-east, for X-check

        Parameters
        ----------
        round : int
            Number of rounds to repeat the syndrome measurement.
        output_format : str
            'str', 'list', or 'stim' to return a stim.Circuit object.

        Returns
        -------
        str, list, or stim.Circuit
            The syndrome check block in the specified format.
        """
        if len(self.data_qubits) == 0:
            raise NoQubitActiveError("Code patch is not initialized. Run select_active before this call.")
        if round < 1:
            raise ValueError("round must be >= 1.")
        
        #Turn off X ancillas
        self.off = 'X'

        gate_sequence = ''

        x_indices_str = " ".join(map(str, self.x_indices))


        if hook_opt:
            '''
            
            Z Basis round
            
            '''
            if detector_reset:
                gate_sequence += ("\nR " + x_indices_str)
            gate_sequence += ("\nCX")


            for x_index in self.x_neighbours:
                data_qubits = []
                for data_index in self.x_neighbours[x_index]:
                    if data_index != -1:
                        data_qubits.append(data_index)
                data_qubits.sort()
                for qubit in data_qubits:
                    gate_sequence += (f" {qubit} {x_index}")

            gate_sequence += ("\nM " + x_indices_str)
            annotation = -1
            j = self.n_x_checks
            for loc, color in self.x_qubits:
                if color == 'red': annotation = 0
                elif color == 'green': annotation = 1
                elif color == 'blue': annotation = 2
                else: ValueError("Invalid color. How?")
                gate_sequence += (f"\nDETECTOR{loc[0], loc[1], 0, annotation} rec[-{j}] rec[-{j + 2*self.n_x_checks}]")
                j-=1

            '''
            
            X Basis round
            
            '''

            gate_sequence += "\nTICK "
            if detector_reset:
                gate_sequence += ("\nR " + x_indices_str)
            gate_sequence += ("\nH " + x_indices_str)
            gate_sequence += "\nCX"

            for x_index in self.x_neighbours:
                data_qubits = []
                for data_index in self.x_neighbours[x_index]:
                    if data_index != -1:
                        data_qubits.append(data_index)
                data_qubits.sort()
                for qubit in data_qubits:
                    gate_sequence+=(f" {x_index} {qubit}")

            gate_sequence += ("\nH " + x_indices_str)
            gate_sequence += ("\nM " + x_indices_str)
            
            j = self.n_x_checks
            for loc, color in self.x_qubits:
                if color == 'red': annotation = 0
                elif color == 'green': annotation = 1
                elif color == 'blue': annotation = 2
                else: ValueError("Invalid color. How?")
                if detector_reset:
                    gate_sequence += (f"\nDETECTOR{loc[0], loc[1], 0, annotation} rec[-{j}] rec[-{j + 2*self.n_x_checks}]")
                j-=1
            gate_sequence += "\nTICK "
        
        else:
            '''
            
            Z Basis round
            
            '''
            if detector_reset:
                gate_sequence += ("\nR " + x_indices_str)
            gate_sequence += ("\nCX")


            for x_index in self.x_neighbours:
                data_qubits = []
                for data_index in self.x_neighbours[x_index]:
                    if data_index != -1:
                        data_qubits.append(data_index)
                data_qubits.sort()
                for qubit in data_qubits:
                    gate_sequence += (f" {qubit} {x_index}")

            gate_sequence += ("\nM " + x_indices_str)
            annotation = -1
            j = self.n_x_checks
            for loc, color in self.x_qubits:
                if color == 'red': annotation = 0
                elif color == 'green': annotation = 1
                elif color == 'blue': annotation = 2
                else: ValueError("Invalid color. How?")
                if detector_reset:
                    gate_sequence += (f"\nDETECTOR{loc[0], loc[1], 0, annotation} rec[-{j}] rec[-{j + 2*self.n_x_checks}]")
                j-=1

            '''
            
            X Basis round
            
            '''

            gate_sequence += "\nTICK "
            if detector_reset:
                gate_sequence += ("\nR " + x_indices_str)
            gate_sequence += ("\nH " + x_indices_str)
            gate_sequence += "\nCX"

            for x_index in self.x_neighbours:
                data_qubits = []
                for data_index in self.x_neighbours[x_index]:
                    if data_index != -1:
                        data_qubits.append(data_index)
                data_qubits.sort()
                for qubit in data_qubits:
                    gate_sequence+=(f" {x_index} {qubit}")

            gate_sequence += ("\nH " + x_indices_str)
            gate_sequence += ("\nM " + x_indices_str)
            
            j = self.n_x_checks
            for loc, color in self.x_qubits:
                if color == 'red': annotation = 0
                elif color == 'green': annotation = 1
                elif color == 'blue': annotation = 2
                else: ValueError("Invalid color. How?")
                if detector_reset:
                    gate_sequence += (f"\nDETECTOR{loc[0], loc[1], 0, annotation} rec[-{j}] rec[-{j + 2*self.n_x_checks}]")
                j-=1
            gate_sequence += "\nTICK "

        if output_format == 'stim':
            if round > 1:
                block_str = gate_sequence
                repeat_block = f"REPEAT {round} {{\n{block_str}\n}}"
                circuit = stim.Circuit(repeat_block)
            else:
                circuit = stim.Circuit(gate_sequence)
            return circuit

        elif output_format == 'list':
            if round > 1:
                return [f"REPEAT {round} {{"] + gate_sequence + ["}}"]
            else:
                return gate_sequence

        elif output_format == 'str':
            if round > 1:
                return f"REPEAT {round} {{\n" + (gate_sequence) + "\n}"
            else:
                return gate_sequence

        else:
            raise ValueError(f"Unsupported output_format: {output_format}")

    def syndrome_check_block_willow(self, round=3, output_format='stim', detector_reset=True):
        """
        Generate the syndrome check blocks (excluding initialization & measurement at the end and beginning).
        
        This block follows the dense syndrome extract method outlined by Google in https://arxiv.org/pdf/2412.14256v1.
        This method permits a local square lattice topology by warping the shape of the honeycombs and including dedicated X and Z ancilla qubits.

        The generated circuit will be in Stim format or as a list or string of gate strings.

        Parameters
        ----------
        round : int
            Number of rounds to repeat the syndrome measurement.
        output_format : str
            'str', 'list', or 'stim' to return a stim.Circuit object.

        Returns
        -------
        str, list, or stim.Circuit
            The syndrome check block in the specified format.
        """
        if len(self.data_qubits) == 0:
            raise NoQubitActiveError("Code patch is not initialized. Run select_active before this call.")
        if round < 1:
            raise ValueError("round must be >= 1.")
        
        #Turn on all stabilizer checks again
        self.off = ''

        gate_sequence = []

        x_indices_str = " ".join(map(str, self.x_indices))
        z_indices_str = " ".join(map(str, self.z_indices))
        if detector_reset:
            gate_sequence.append("R " + x_indices_str + " " + z_indices_str)
        gate_sequence.append("H " + x_indices_str)

        for control, target in zip(self.z_indices, self.x_indices):
            # print(f"CX {control} {target}")
            gate_sequence.append(f"CX {control} {target}")
        
        for x_index in self.x_neighbours:
            data_qubits = []
            for data_index in self.x_neighbours[x_index]:
                if data_index != -1:
                    data_qubits.append(data_index)
            data_qubits.sort()
            for qubit in data_qubits[:len(data_qubits)//2]:
                
                gate_sequence.append(f"CX {x_index} {qubit}")
                gate_sequence.append(f"CX {qubit} {x_index}")

        for z_index in self.z_neighbours:
            data_qubits = []
            for data_index in self.z_neighbours[z_index]:
                if data_index != -1:
                    data_qubits.append(data_index)
            data_qubits.sort()
            for qubit in data_qubits[len(data_qubits)//2:]:
                gate_sequence.append(f"CX {qubit} {z_index}")
                gate_sequence.append(f"CX {z_index} {qubit}")
        
        for control, target in zip(self.x_indices, self.z_indices):
            gate_sequence.append(f"CX {control} {target}")

        gate_sequence.append("H " + x_indices_str)
        gate_sequence.append("M " + x_indices_str + " " + z_indices_str)
        gate_sequence.append("TICK")

        total_ancillas = len(self.x_indices) + len(self.z_indices)
        if detector_reset:
            for j in range(total_ancillas):
                gate_sequence.append(f"DETECTOR rec[-{j+1}] rec[-{j+1 + total_ancillas}]")

        if output_format == 'stim':
            if round > 1:
                block_str = '\n'.join(gate_sequence)
                repeat_block = f"REPEAT {round} {{\n{block_str}\n}}"
                circuit = stim.Circuit(repeat_block)
            else:
                circuit = stim.Circuit('\n'.join(gate_sequence))
            return circuit

        elif output_format == 'list':
            if round > 1:
                return [f"REPEAT {round} {{"] + gate_sequence + ["}}"]
            else:
                return gate_sequence

        elif output_format == 'str':
            if round > 1:
                return f"REPEAT {round} {{\n" + '\n'.join(gate_sequence) + "\n}"
            else:
                return '\n'.join(gate_sequence)

        else:
            raise ValueError(f"Unsupported output_format: {output_format}")
        

    
    def readout_block_willow(self, output_format='stim'):
        """
        Generate the readout block.

        The generated circuit will be in Stim format or as a list or string of gate strings.

        Parameters
        ----------
        output_format : str
            'str', 'list', or 'stim' to return a stim.Circuit object.

        Returns
        -------
        str, list, or stim.Circuit
            The syndrome check block in the specified format.
        """

        #Generate Bell states w/ circuit
        gate_sequence = 'H '
        gate_sequence += ''.join(str(x)+ " " for x in self.x_indices)
        gate_sequence += '\nM '
        gate_sequence += ''.join(str(x)+ " " for x in self.x_indices)
        gate_sequence += ''.join(str(x)+ " " for x in self.z_indices)
        gate_sequence += ''.join(str(d)+ " " for d in self.data_indices)
        
        # for d in self.data_indices:
        #     if d in [index[0] for _, index in self.logical_locations]:
        #         gate_sequence += f'\nOBSERVABLE_INCLUDE(0) rec[{d-len(self.data_indices)}] '


        if output_format == 'stim':
            circuit = stim.Circuit((gate_sequence))
            return circuit

        elif output_format == 'list':
            return gate_sequence

        elif output_format == 'str':
            return ''.join(gate_sequence)

        else:
            raise ValueError(f"Unsupported output_format: {output_format}")
        
    def readout_block_hex(self, output_format='stim', basis='z', detector_reset=True):
        """
        Generate the readout block.

        The generated circuit will be in Stim format or as a list or string of gate strings.

        Parameters
        ----------
        output_format : str
            'str', 'list', or 'stim' to return a stim.Circuit object.

        Returns
        -------
        str, list, or stim.Circuit
            The syndrome check block in the specified format.
        """

        #Generate Bell states w/ circuit
        gate_sequence = ''
        

        if basis.lower() == 'x':
            gate_sequence += '\nH '
            gate_sequence += ''.join(str(d)+ " " for d in self.data_indices)
        elif basis.lower() == 'y':
            gate_sequence = '\nS '
            gate_sequence += ''.join(str(d)+ " " for d in self.data_indices)
            gate_sequence += '\nH '
            gate_sequence += ''.join(str(d)+ " " for d in self.data_indices)

        gate_sequence += '\nM '
        gate_sequence += ''.join(str(d)+ " " for d in self.data_indices)
        # print(gate_sequence)

        if detector_reset:
            for index, loc in self.logical_locations:
                gate_sequence += f"\nOBSERVABLE_INCLUDE(0) rec[{index-len(self.data_indices)}]"


        if output_format == 'stim':
            circuit = stim.Circuit((gate_sequence))
            return circuit

        elif output_format == 'list':
            return gate_sequence

        elif output_format == 'str':
            return ''.join(gate_sequence)

        else:
            raise ValueError(f"Unsupported output_format: {output_format}")
        
    
    def arb_state_injection_single(self, basis='z', output_format="stim"):
        '''Prepares all of the qubits for state injection cycle via single qubit stabilizers
        Data qubit 0 is taken to be the injection site.'''

        data_qubits = ""
        anc_qubits = ""
        gates = "\nR "

        if basis.lower() == 'z':
            self.off = 'x'
        elif basis.lower() == 'x':
            self.off = 'z'
        else:
            ValueError("Invalid basis choice for state injection! Basis {basis}")
            
        
        for loc in self.data_qubits:
            if(self.qubit_index[tuple(loc)] > 0):
                data_qubits += str(self.qubit_index[tuple(loc)]) + " "
        if self.off == 'x':
            for loc, _ in self.z_qubits:
                anc_qubits += str(self.qubit_index[tuple(loc)]) + " "
        if self.off == 'z':
            for loc, _ in self.x_qubits:
                anc_qubits += str(self.qubit_index[tuple(loc)]) + " "
        

        gates+= data_qubits
        gates+= anc_qubits

        if basis.lower() == 'x':
            gates += (f"\nH {anc_qubits}")   

        if output_format == 'str':
            return ''.join(gates)
        elif output_format == 'list':
            return gates
        else:
            return gates
        
    def arb_state_injection_bell_hex(self, output_format="stim", inj='0', I_error=True, detector_reset=True):
        '''Prepares all of the qubits for state injection cycle via two-qubit stabilizers, then a syndrome check skipping Red ancillas.
        Builds out the rest of the stabilizers by preparing Bell states. This results in 66% stabilizer coverage for the first round.
        We take Red to be the color that data qubit 0 is apart of. For the injection round, Red stabilizers cannot be determined,
        but Blue and Green are still active.
        Data qubit 0 is taken to be the injection site.'''


        gate_sequence = ''
        for loc, index in zip(self.qubit_index, self.qubit_location):
            if loc not in [pos for pos, color in self.z_qubits]:
                gate_sequence += f"QUBIT_COORDS({loc[0]}, {loc[1]}) {index}\n"
            # gate_sequence += "\n"
        
        # print(gate_sequence)
        # if I_error:
            # gate_sequence += "I "

            # for loc, index in zip(self.qubit_index, self.qubit_location):
            #     if loc not in [pos for pos, color in self.z_qubits]:
            #         gate_sequence += f"{index} "
        # for index, loc in self.logical_locations:
        #     gate_sequence += f'\nOBSERVABLE_INCLUDE(0) rec[{index-len(self.data_qubits)-len(self.x_qubits)}]'

        #Inject a state
        if inj == '1':
            gate_sequence += "\nX 0"
        elif inj == '+':
            gate_sequence += "\nH 0"
        elif inj == '-':
            gate_sequence += "\nH 0\nZ 0"
        elif inj == 'i' or inj.lower() == '+i':
            gate_sequence += "\nH 0\nS 0"
        elif inj == '-i':
            gate_sequence += "\nH 0\nS_DAG 0"
        else:
            ValueError(f"Invalid injection state: {inj}")

        gate_sequence += '\nH '
        gate_sequence += ' '.join(str(i) + "" for i, j in self.bell_pairs)
        gate_sequence += '\nCNOT '
        gate_sequence += ' '.join(str(i)+ " " +str(j) for i, j in self.bell_pairs)

        colors = dict(self.x_qubits)
        x_indices_str = " ".join(
            str(x_index)
            for x_index in self.x_indices
            if self.off != colors[self.qubit_location[x_index]]
        )

        gate_sequence += "\nTICK "
        '''
        
        Z Basis round
        
        '''

        # gate_sequence += "\nR " + x_indices_str
        gate_sequence += '\nCX' 
        for x_index in self.x_neighbours:
            data_qubits = []
            for data_index in self.x_neighbours[x_index]:
                if data_index != -1:
                    data_qubits.append(data_index)
            data_qubits.sort()
            for qubit in data_qubits:
                gate_sequence += (f" {qubit} {x_index}")

        gate_sequence += "\nM " + x_indices_str
        annotation = -1
        j = len(self.x_qubits)
        for loc, color in self.x_qubits:
            if color != 'red':
                if color == 'green':
                    annotation = 4
                elif color == 'blue':
                    annotation = 5
                else:
                    ValueError("Invalid color. How?")

                if detector_reset:
                    gate_sequence += f"\nDETECTOR{loc[0], loc[1], 0, annotation} rec[-{j}]"
            j-=1

        gate_sequence += "\nTICK "
        '''
        
        X Basis round
        
        '''
        if detector_reset:
            gate_sequence += ("\nR " + x_indices_str)
        gate_sequence += ("\nH " + x_indices_str)
        gate_sequence += "\nCX"
        for x_index in self.x_neighbours:
            data_qubits = []
            for data_index in self.x_neighbours[x_index]:
                if data_index != -1:
                    data_qubits.append(data_index)
            data_qubits.sort()
            for qubit in data_qubits:
                gate_sequence+=(f" {x_index} {qubit}")

        gate_sequence += ("\nH " + x_indices_str)
        gate_sequence += ("\nM " + x_indices_str)
        
        j = len(self.x_qubits)
        for loc, color in self.x_qubits:
            if color != 'red':
                if color == 'green':
                    annotation = 1
                elif color == 'blue':
                    annotation = 2
                else:
                    ValueError("Invalid color. How?")
                if detector_reset:
                    gate_sequence += f"\nDETECTOR{loc[0], loc[1], 0, annotation} rec[-{j}]"
            j-=1

        gate_sequence += "\nTICK "

        # print(gate_sequence)

        if output_format == 'stim':
            circuit = stim.Circuit((gate_sequence))
            return circuit

        elif output_format == 'list':
            return gate_sequence

        elif output_format == 'str':
            return ''.join(gate_sequence)

        else:
            raise ValueError(f"Unsupported output_format: {output_format}")
        
    def arb_state_injection_bell_willow(self, output_format="stim"):
        '''Prepares all of the qubits for state injection cycle via two-qubit stabilizers, then a syndrome check skipping Red ancillas.
        Builds out the rest of the stabilizers by preparing Bell states. This results in 66% stabilizer coverage for the first round.
        We take Red to be the color that data qubit 0 is apart of. For the injection round, Red stabilizers cannot be determined,
        but Blue and Green are still active.
        Data qubit 0 is taken to be the injection site.'''

        #Find the bell pairs
        self._gen_bell_pairs()

        #Turn off red stabilizers
        self.off = 'red'

        #Generate Bell states w/ circuit
        gate_sequence = 'R '
        gate_sequence += ' '.join(str(i) + "" for i  in self.data_indices if i > 0)
        gate_sequence += '\nH '
        gate_sequence += ' '.join(str(i) + "" for i, j in self.bell_pairs)
        gate_sequence += '\nCNOT '
        gate_sequence += ' '.join(str(i)+ " " +str(j) for i, j in self.bell_pairs)

        x_indices_str = " ".join(
            str(x_index)
            for x_index in self.x_indices
        )
        z_indices_str = " ".join(
            str(z_index)
            for z_index in self.z_indices
        )

        gate_sequence += ("\nR " + x_indices_str + " " + z_indices_str)
        gate_sequence += ("\nH " + x_indices_str)

        for control, target in zip(self.z_indices, self.x_indices):
            gate_sequence += (f"\nCX {control} {target}")
        
        for x_index in self.x_neighbours:
            # print(colors[self.qubit_location[x_index]])
            data_qubits = []
            for data_index in self.x_neighbours[x_index]:
                if data_index != -1:
                    data_qubits.append(data_index)
            data_qubits.sort()
            for qubit in data_qubits[:len(data_qubits)//2]:
                gate_sequence += (f"\nCX {x_index} {qubit}")
                gate_sequence += (f"\nCX {qubit} {x_index}")

        
        for z_index in self.z_neighbours:
            # print(colors[self.qubit_location[z_index]])
            data_qubits = []
            for data_index in self.z_neighbours[z_index]:
                if data_index != -1:
                    data_qubits.append(data_index)
            data_qubits.sort()
            for qubit in data_qubits[len(data_qubits)//2:]:
                gate_sequence += (f"\nCX {qubit} {z_index}")
                gate_sequence += (f"\nCX {z_index} {qubit}")
            
        for control, target in zip(self.x_indices, self.z_indices):
            gate_sequence += (f"\nCX {control} {target}")

        gate_sequence += ("\nH " + x_indices_str)
        gate_sequence += "\nM " +  " ".join(str(i) + "" for i in self.x_indices) + " " + " ".join(str(j) + "" for j in self.z_indices)
        # print(gate_sequence)
        

        
        colors = dict(self.x_qubits)
        j = len(colors)
        for loc in colors:
            if self.off != colors[loc]:
                gate_sequence += (f"\nDETECTOR({loc[0]}, {loc[1]}) rec[-{j+1}]")
                j-=1
        colors = dict(self.z_qubits)
        for loc in colors:
            if self.off != colors[loc]:
                gate_sequence += (f"\nDETECTOR({loc[0]}, {loc[1]}) rec[-{j+1}]")
                j-=1

        gate_sequence += ("\nTICK")
        

        if output_format == 'stim':
            circuit = stim.Circuit((gate_sequence))
            return circuit

        elif output_format == 'list':
            return gate_sequence

        elif output_format == 'str':
            return ''.join(gate_sequence)

        else:
            raise ValueError(f"Unsupported output_format: {output_format}")

        
    # def code_stabilizer(self, output_format='stim') -> list:
    #     '''
    #     return the stabilizers of the code

    #     Parameters:
    #     -----------
    #     output_format: str
    #         'stim': for list of stim.PauliString obj (default)
    #         'tab': for the Tableau class instance
    #         'list': for the list of strings.
    #     '''
    #     stabilizers = []
    #     for key, val in self.x_neighbours.items():
    #         stabilizers.append(pauli_string(self.n_data, val, 'X'))
        
    #     for key, val in self.z_neighbours.items():
    #         stabilizers.append(pauli_string(self.n_data, val, 'Z'))

    #     if output_format == 'stim':
    #         stim_stab = []
    #         for ps in stabilizers:
    #             stim_stab.append(stim.PauliString(ps))
    #         return stim_stab
        
    #     elif output_format == 'tab':
    #         stab_tab = []
    #         for ps in stabilizers:
    #             stab_tab.append(Tableau.convert_back(ps.replace('_', "I")))
    #         temp_stab = stab_tab[0]
    #         for tab in stab_tab[1:]:
    #             temp_stab.append(tab)
    #         return temp_stab
        
    #     else:
    #         return stabilizers
    
    # def logical_ops(self, basis, output_format='stim'):
    #     '''return the logical operator'''
    #     loc = np.array([0,0], dtype=int)
        
    #     if basis == 'z' or basis == 'Z':
    #         indices = []
    #         for shift in range(self.distance):
    #             loc_op = loc + np.array([2*shift, 0], dtype=int)
    #             indices.append(self.qubit_index[tuple(loc_op)])
    #         ops = ['Z'] * self.distance

    #     elif basis == 'x' or basis == 'X':
    #         indices = []
    #         for shift in range(self.distance):
    #             loc_op = loc + np.array([0, 2*shift], dtype=int)
    #             indices.append(self.qubit_index[tuple(loc_op)])
    #         ops = ["X"] * self.distance
        
    #     else:
    #         indices = []
    #         for shift in range(1, self.distance):
    #             loc_op = loc + np.array([0, 2*shift], dtype=int)
    #             indices.append(self.qubit_index[tuple(loc_op)])
    #         ops = ["X"] * (self.distance - 1)

    #         for shift in range(1, self.distance):
    #             loc_op = loc + np.array([2*shift, 0], dtype=int)
    #             indices.append(self.qubit_index[tuple(loc_op)])
    #         ops += ['Z'] * (self.distance - 1)

    #         indices.append(self.qubit_index[tuple(loc)])
    #         ops.append('Y')

    #     ps = pauli_string(self.n_data, indices, ops)
    #     if output_format == 'stim':
    #         return stim.PauliString(ps)
    #     elif output_format == 'tab':
    #         return Tableau.convert_back(ps)
    #     elif output_format == 'index':
    #         return indices
    #     else:
    #         return ps
        

    # def code_measure_out(self, basis, output_format='list', full_checks=True):
    #     '''
    #     the measurement along a given basis

    #     Parameters:
    #     -----------
    #     basis: str
    #         'x', 'y', or 'z': y-basis has not been implemented yet.
    #     output_format: str
    #         'list': list of gates
    #         'str': a single string
    #         whether there are full syndrome checks before the measurement circuits
    #         when round=1, in stim circuit generation, there is no full round checks.
    #         the number of measurements can be different.
    #     '''
    #     #full_check_bit = int(full_checks)
    #     gates = []
    #     ## Z or X measurements: measure all qubits into the corresponding basis
    #     ## then apply a syndrome check using the measured outcomes
    #     data_indices_str = " ".join(map(str, self.data_indices))

    #     if basis == 'Y' or basis == 'y':
    #         ## this needs to be implemented later
    #         #print('this has not been implemented yet, using Z ')
    #         #basis = 'z'

    #         ## the strategy here is:
    #         ## using the non-fault-tolerant gates to implement a S^dagger gate, 
    #         ## then the Y is converted into X basis measurements
    #         logical_z_indices = self.logical_ops(basis='z', output_format='index')
    #         for j in range(len(logical_z_indices)-1):
    #             gates.append(f'CX {logical_z_indices[j]} {logical_z_indices[j+1]}')
            
    #         gates.append("TICK")
    #         gates.append(f'S_DAG {logical_z_indices[-1]}')
    #         gates.append("TICK")
    #         for j in range(len(logical_z_indices)-1, 0, -1):
    #             gates.append(f'CX {logical_z_indices[j-1]} {logical_z_indices[j]}')
            
    #         gates.append("TICK")
    #         #basis = 'x'

    #     if basis == 'x' or basis == 'X' or basis == 'y' or basis == 'Y':
    #         gates.append('H ' + data_indices_str)
        
    #     gates.append('M ' + data_indices_str)

    #     ## now consider the syndrome checks at the final measurement results.
    #     loc = np.array([0,0], dtype=int) ## this is for tracking the locations of the data qubits for the logical operators

    #     # if basis == 'Y' or basis == 'y':
    #     #     ## need to be implemented
    #     #     pass
    #     if basis == 'z' or basis == 'Z':
    #         for j, z_index in enumerate(self.z_indices):
    #             data_indices = self.z_neighbours[z_index]
    #             rec_indices = [k - self.n_data for k in data_indices if k != -1]
    #             gate_str = f'DETECTOR ' + ' '.join(map(lambda x: f'rec[{x}]', rec_indices)) \
    #                 + f' rec[{j - (self.n_data + self.n_z_checks)}]'
    #             gates.append(gate_str)
                
    #         ## generate the obsrvable
    #         indices = []
    #         for shift in range(self.distance):
    #             loc_op = loc + np.array([2*shift, 0], dtype=int)
    #             indices.append(self.qubit_index[tuple(loc_op)])
                
    #     elif basis == 'x' or basis == 'X' or basis == 'y' or basis == 'Y':
    #     #elif basis == 'x' or basis == 'X':
    #         for j, x_index in enumerate(self.x_indices):
    #             data_indices = self.x_neighbours[x_index]
    #             rec_indices = [k - self.n_data for k in data_indices if k != -1]
    #             if basis == 'x' or basis == 'X':
    #                 gate_str = f'DETECTOR ' + ' '.join(map(lambda x: f'rec[{x}]', rec_indices)) \
    #                     + f' rec[{j - (self.n_data + self.n_z_checks + self.n_x_checks)}]'
    #                 gates.append(gate_str)

    #         indices = []
    #         for shift in range(self.distance):
    #             loc_op = loc + np.array([0, 2*shift], dtype=int)
    #             indices.append(self.qubit_index[tuple(loc_op)])
        

    #     logical_rec_indices = [k - self.n_data for k in indices if k != -1]
    #     ## add the final observable
    #     gates.append(f'OBSERVABLE_INCLUDE(0) ' + " ".join(map(lambda x: f'rec[{x}]', logical_rec_indices)))

    #     if output_format == 'str':
    #         return '\n'.join(gates)
    #     else:
    #         return gates

    # #def encoding_circuit(self, output_format='list'):
    # #    pass

    # def decoding_circuit(self, output_format='list'):
    #     pass
    
    # def stim_circ_init(self, output_format='list'):
    #     '''generate the informaiton about the qubit location, etc., for stim circuits'''
    #     gate_sequence = []
    #     for loc in self.qubit_index:
    #         index = self.qubit_index[loc]
    #         gate_sequence.append(f'QUBIT_COORDS({loc[0]}, {loc[1]}) {index}')
        
    #     if output_format == 'list':
    #         return gate_sequence
    #     else:
    #         return '\n'.join(gate_sequence) + '\n'

    # def logical_state_vector(self, basis = 'z'):
    #     '''generate the state vector for the code basis state |+sigma>'''
    #     code_stab = self.code_stabilizer(output_format='stim')
    #     logical_op = self.logical_ops(basis=basis, output_format='stim')
    #     code_stab.append(logical_op)

    #     stab_state = stim.Tableau.from_stabilizers(code_stab)

    #     return stab_state.to_state_vector(endian='little')
    
    # def stim_circ_aux_info(self, output_format='list'):
    #     '''include the coordinate and index information for stim circuits'''
    #     output = []
    #     for index, coord in self.qubit_location.items():
    #         output.append(f'QUBIT_COORDS{tuple(coord)} {index}')
        
    #     if output_format == 'str':
    #         return '\n'.join(output)
    #     else:
    #         return output
 
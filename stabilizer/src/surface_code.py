''' some qec circuit generation '''
from helper.tab_state import Tableau
from typing import Union, List
import numpy as np
import stim
#from abc import abstractmethod
import math
import matplotlib.pyplot as plt

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
    

def qubit_type(location):
    loc = np.array(location[:2])

    for j,val in enumerate(loc):
        if not math.isclose(val, round(val), abs_tol=1e-8):
            raise NotQubitLocationError(f"location {location} is not a valid qubit location.")
        else:
            loc[j] = round(val)
    
    loc = loc.astype(int)

    ## check if the given location is a data qubit
    if not ((loc[0]&1) | (loc[1]&1)):
        return 'data', loc
    
    loc = loc % 4
    if sum(loc) % 4 == 0:
        ## only (1,3) and (3,1), which are assigned as X
        return 'x', loc
    else:
        return 'z', loc

def is_qubit(location):
    try:
        _ =qubit_type(location)
    except NotQubitLocationError:
        return False
    
    return True

def is_data(location):
    '''check if the given location is a data qubit'''
    try:
        q_type, _ = qubit_type(location)
    except NotQubitLocationError:
        return False
    
    if q_type == 'data':
        return True
    else:
        return False
def is_x_syndrome(location):
    '''check if the given location corresponding to a X-stabilizer check qubits'''
    try:
        q_type, _ = qubit_type(location)
    except NotQubitLocationError:
        return False
    
    if q_type == 'x':
        return True
    else:
        return False

def is_z_syndrome(location):
    '''check if the given location corresponding to a Z-stabilizer check qubits'''
    try:
        q_type, _ = qubit_type(location)
    except NotQubitLocationError:
        return False
    
    if q_type == 'z':
        return True
    else:
        return False

def check_qubits(location, syn_round, q_type='x', consistency_check=False):
    '''
    find the data qubits for the syndrome check gate operations 
    Assuming vertical edges are logical X operator (Z edge), and horizontal edges are Z.

    Parameters:
    -----------
    syn_round: int, (0,1,2,3):
        if syn_round == 0: left-up
        if syn_round == 3: right-down
        if syn_round == 1: if Z check: right-up, if X check: left-down
        if syn_round == 2: if Z check: left-down, if X check: right-up
    '''
    if consistency_check:
        q_type_checking, loc = qubit_type(location)
        if q_type != q_type_checking:
            raise ConsistencyError(f'Inconsistency detected! loc {location} is a {q_type_checking}-type syndrome qubits. \n \
                                   The given q_type is {q_type}-type.')
    loc = np.array(location).astype('int')
    b0 = syn_round & 1
    b1 = (syn_round >> 1)&1
    
    if q_type == 'z':
        loc[0] -= (-1) ** b1
        loc[1] -= (-1) ** b0
    else:
        loc[0] -= (-1) ** b0
        loc[1] -= (-1) ** b1

    return loc


class SurfaceCodeQubits:
    def __init__(self):
        self.qubit_index = {}
        self.qubit_location = {}
    
        self.data_qubits = []
        self.x_qubits = []
        self.z_qubits = []

        self.n_x_checks = 0
        self.n_z_checks = 0
        self.n_data = 0

        self.x_indices = []
        self.z_indices = []
        self.data_indices = []
        self.x_neighbours = {} # x_indices: [r0_data_index, r1_data_index, r2_data_index, r3_data_index]
        self.z_neighbours = {} # z_indices: [r0_data_index, r1_data_index, r2_data_index, r3_data_index]

        self.distance = None

        return 
    
    def _reset(self):

        self.qubit_index = {}
        self.qubit_location = {}
    
        self.data_qubits = []
        self.x_qubits = []
        self.z_qubits = []

        self.n_x_checks = 0
        self.n_z_checks = 0
        self.n_data = 0

        self.x_indices = []
        self.z_indices = []
        self.data_indices = []
        self.x_neighbours = {} # x_indices: [r0_data_index, r1_data_index, r2_data_index, r3_data_index]
        self.z_neighbours = {} # z_indices: [r0_data_index, r1_data_index, r2_data_index, r3_data_index]

        self.distance = None
        return 
 

    def select_active(self, distance = 3):
        """
        Build the lists of active qubits for a rotated surface code patch
        of distance d.

        Data qubits: (even even) sites: (0,0)to (2d-2,2d-2)
        Horizontal edges: X edges (logical Z operator)
        Vertical edges: Z edges (logical X operator)

        Parameters
        ----------
        distance : int, default 3
            code distance for the rotated surface code patch
        """
        if self.distance is not None:
            self._reset()
        self.distance = distance
        
        coord0 = np.array([0, 0], dtype=int)
        for j in range(distance):
            for k in range(distance):
                self.data_qubits.append(coord0 + np.array([j,k]) * 2)
        
        ## now consider the syndrome check qubits inside the qubit region
        coord1 = np.array([1,1], dtype=int)
        for j in range(distance - 1):
            for k in range(distance - 1):
                syn_coord = coord1 + np.array([j, k]) * 2
                if is_x_syndrome(syn_coord):
                    self.x_qubits.append(syn_coord)
                else:
                    self.z_qubits.append(syn_coord)
        
        ## now consider the x_edges:
        ## horizontal edges: need X
        coord1 = np.array([1, -1], dtype=int) ## upper edge
        coord2 = coord1 + np.array([0, distance * 2], dtype = int)
        for j in range(0, distance-1, 2):
            self.x_qubits.append(coord1 + np.array([2*j, 0]))
            self.x_qubits.append(coord2 + np.array([2*(j+1), 0]))
        
        ## vertical edge: need Z
        coord1 =  np.array([-1, 1], dtype=int)
        coord2 = coord1 + np.array([2 * distance, 0], dtype=int)
        for j in range(0, distance-1, 2):
            self.z_qubits.append(coord2 + np.array([0, 2*j]))
            self.z_qubits.append(coord1 + np.array([0, 2*(j+1)]))
        
        ## indexing the qubits:
        self.x_qubits = np.array(self.x_qubits).astype('int')
        self.z_qubits = np.array(self.z_qubits).astype('int')
        self.data_qubits = np.array(self.data_qubits).astype('int')

        self.n_x_checks = len(self.x_qubits)
        self.n_z_checks = len(self.z_qubits)
        self.n_data = len(self.data_qubits)

        self._process()
        return

    def _process(self):
        if not (self.n_x_checks or self.n_z_checks or self.n_data):
            print(f'nx: {self.n_x_checks}, nz: {self.n_z_checks}, nd: {self.n_data}')
            raise NoQubitActiveError('the process is incorrectly applied when there are no active qubits.')

        self.qubit_index = {}
        self.qubit_location = {}
        self.x_indices = []
        self.z_indices = []
        self.data_indices = []

        ct = 0
        for loc in self.data_qubits:
            self.qubit_index[tuple(loc)] = ct
            self.qubit_location[ct] = loc
            self.data_indices.append(ct)
            ct += 1
        for loc in self.x_qubits:
            self.qubit_index[tuple(loc)] = ct
            self.qubit_location[ct] = loc
            self.x_indices.append(ct)
            ct += 1
        for loc in self.z_qubits:
            self.qubit_index[tuple(loc)] = ct
            self.qubit_location[ct] = loc
            self.z_indices.append(ct)
            ct += 1

        for loc in self.x_qubits:
            nearby_data = []
            x_index = self.qubit_index[tuple(loc)]
            for ts in range(4):
                ## solve the nearby data qubits' indices
                ## if not active data qubits found, use -1 for the index
                loc_data = check_qubits(loc, ts, 'x')
                if tuple(loc_data) in self.qubit_index:
                    nearby_data.append(self.qubit_index[tuple(loc_data)])
                else:
                    nearby_data.append(-1)
            self.x_neighbours[x_index] = tuple(nearby_data)
        
        for loc in self.z_qubits:
            nearby_data = []
            z_index = self.qubit_index[tuple(loc)]
            for ts in range(4):
                ## solve the nearby data qubits' indices
                ## if not active data qubits found, use -1 for the index
                loc_data = check_qubits(loc, ts, 'z')
                if tuple(loc_data) in self.qubit_index:
                    nearby_data.append(self.qubit_index[tuple(loc_data)])
                else:
                    nearby_data.append(-1)
            self.z_neighbours[z_index] = tuple(nearby_data)
        return
    

    def draw(self, fontsize: int = 14, show: bool = True):
        """
        Visualise the currently‑selected patch.

        Parameters
        ----------
        fontsize : int
            Font size used for the 'D', 'X', 'Z' labels.
        show : bool
            If True (default) call plt.show() at the end.  Pass False when
            composing figures manually.
        """
        if not (self.n_x_checks or self.n_z_checks or self.n_data):
            raise NoQubitActiveError("No qubits selected: call select_active() first")

        # collect all coordinates; convert numpy arrays to tuples
        coords_D = [tuple(q) for q in self.data_qubits]
        coords_X = [tuple(q) for q in self.x_qubits]
        coords_Z = [tuple(q) for q in self.z_qubits]

        # helper to plot a set of (x,y) points with a text label
        def _plot(coords, label, **kw):
            if not coords:
                return
            xs, ys = zip(*coords)
            plt.scatter(xs, ys, **kw)
            for x, y in coords:
                plt.text(x, y, label+f':{self.qubit_index[(x,y)]}', ha="center", va="center", fontsize=fontsize)

        # Data qubits as squares
        _plot(coords_D, "D", marker="s", s=120, color="white", edgecolors="k")
        # X checks as circles
        _plot(coords_X, "X", marker="o", s=120, color="tab:red", edgecolors="k")
        # Z checks as triangles
        _plot(coords_Z, "Z", marker="^", s=120, color="tab:blue", edgecolors="k")

        plt.gca().set_aspect("equal", adjustable="box")
        plt.grid(True, linestyle="--", linewidth=0.5)
        plt.xlabel("x  (rightwards)")
        plt.ylabel("y  (downwards)")
        plt.title("Surface‑code patch")

        # Invert the y‑axis so that *down* is positive
        plt.gca().invert_yaxis()

        if show:
            plt.show()


    def syndrome_check_block(self, round=1, output_format='str'):
        '''
        generate the syndrome check blocks. No initialization and final measurements.
        the generated circuit will be in the stim format, or QASM.
        '''
        if len(self.data_qubits) == 0:
            raise NoQubitActiveError("Code patch is not initialized. Run select_active before this call.")

        if round < 1:
            raise ValueError("round can only be larger than 1.")
        
        n_total_checks = self.n_x_checks + self.n_z_checks

        ## generating the single round of syndrome check
        gate_sequence = []
        x_indices_str = " ".join(map(str, self.x_indices))
        z_indices_str = " ".join(map(str, self.z_indices))

        ## resettting all the syndrome qubits
        gate_sequence.append("R "+x_indices_str+ " "+z_indices_str)
            
        ## reset the X syndrome into |+> state by applying H gates
        gate_sequence.append("H " + x_indices_str)

        gate_sequence.append("TICK")

        ## now consider performing the syndrome checks
        for ts in range(4):
            gate_str = 'CX '
            for x_index in self.x_indices:
                data_index = self.x_neighbours[x_index][ts]
                if data_index != -1:
                    gate_str += f'{x_index} {data_index} '
            
            for z_index in self.z_indices:
                data_index = self.z_neighbours[z_index][ts]
                if data_index != -1:
                    gate_str += f'{data_index} {z_index} '
            gate_sequence.append(gate_str)
            gate_sequence.append('TICK')
        
        ## now consider applying syndrome measurements
        gate_sequence.append('H '+x_indices_str)

        ## perform measurements
        gate_sequence.append('M '+ x_indices_str + " " + z_indices_str)

        gate_sequence.append('TICK')

        ## now consider comparing the syndrome check results from the last checks
        for j in range(n_total_checks):
            gate_sequence.append(f'DETECTOR rec[-{j+1}] rec[-{j+1+n_total_checks}]')

        if round > 1:
            #gate_str = f"REPEAT {round} {{\n {'\n'.join(gate_sequence)} \n}}\n"
            #gate_str = f"REPEAT {round} " + "{ \n" + f"{'\n'.join(gate_sequence)}" + "\n}\n"
            if output_format == 'str':
                gate_str = f"REPEAT {round} " + \
                        "{\n" + \
                        f"{'\n'.join(gate_sequence)}" + \
                        "\n}\n"
                return gate_str
            else:
                gate_sequence = [f"REPEAT {round} " + "{"] + gate_sequence
                gate_sequence.append('}')
                return gate_sequence
        else:
            if output_format == 'list':
                return gate_sequence
            else:
                gate_str = '\n'.join(gate_sequence) + '\n'

                return gate_str
    
    def syndrome_check_block_compressed(self, round=1, output_format='list'):
        '''generate the same syndrome check block, but in qasm format, and reuse the measurement qubits so we only use one aux qubit.'''
        if len(self.data_qubits) == 0:
            raise NoQubitActiveError("Code patch is not initialized. Run select_active before this call.")

        if round < 1:
            raise ValueError("round can only be larger than 1.")
        
        n_total_checks = self.n_x_checks + self.n_z_checks
        #classical_rec_size = n_total_checks * round
        
        aux_index = self.n_data
        #single_round_gates = [f'creg rec[{classical_rec_size}]']
        single_round_gates = []
        ct = 0
        for x_index in self.x_indices:
            single_round_gates.append(f'reset q[{aux_index}];')
            single_round_gates.append(f'h q[{aux_index}];')
            data_indices = self.x_neighbours[x_index]
            for data_index in data_indices:
                if data_index != -1:
                    single_round_gates.append(f'cx q[{aux_index}] q[{data_index}];')
            single_round_gates.append(f'h q[{aux_index}];')
            single_round_gates.append(f'measure q[{aux_index}] -> rec[{ct}]')
            ct += 1
        
        for z_index in self.z_indices:
            single_round_gates.append(f'reset q[{aux_index}];')
            data_indices = self.z_neighbours[z_index]
            for data_index in data_indices:
                if data_index != -1:
                    single_round_gates.append(f'cx q[{data_index}] q[{aux_index}];')
            single_round_gates.append(f'measure q[{aux_index}] -> rec[{ct}]')
            ct += 1
        
        if round == 1:
            total_gates = single_round_gates
        else:
            total_gates = list(single_round_gates)
            for _ in range(1, round):
                for temp in total_gates:
                    if temp.startswith('measure'):
                        total_gates.append(f'measure q[{aux_index}] -> rec[{ct}]')
                        ct +=1
                    else:
                        total_gates.append(temp)

        if output_format == 'list':
            return total_gates
        elif output_format == 'str':
            return '\n'.join(total_gates)
        else:
            print('unknown output format, using list')
            return total_gates
    
    def state_injection(self, basis, output_format="list"):
        '''non-fault-tolerant state injection'''
        x_checks = self.x_indices
        z_checks = self.z_indices

        #satisfied_checks = set({})
        satisfied_checks_order = [] # record the position in the total checks, used later in qec
        #satisfied_z_checks = []

        # loc_inj = np.array([0,0], dtype=int)

        # ## disable the injection site's nearby syndrome checks
        # for i in (-1, 1):
        #     for j in (-1, 1):
        #         check_loc = tuple(loc_inj + np.array([i,j], dtype=int))
        #         if check_loc in self.qubit_index:
        #             ## this is within the active code patch
        #             ## this syndrome check should be either x or z checks.
        #             index = self.qubit_index
        #             if is_x_syndrome(check_loc):
        #                 x_checks.remove(index)
        #             else:
        #                 z_checks.remove(index)
        #         else:
        #             continue
        
        ## x-logical op edge should be in x_init to start with
        x_init = set(self.logical_ops(basis='x', output_format='index'))
        x_init.remove(0) ## this is the injection site

        ## similar for z
        z_init = set(self.logical_ops(basis='z', output_format='index'))
        z_init.remove(0)

        ## find all un-assinged data qubit indices
        unassigned_data = []
        for x in self.data_indices:
            if (x!=0 and x not in x_init and x not in z_init):
                unassigned_data.append(x)
        unassigned_data = set(unassigned_data)

        ## check all the syndrome checks, if it is possible to make it determined, then do so
        ## check the x_checks
        for j, x_check_index in enumerate(x_checks):
            determined = True
            x_nearby = set(self.x_neighbours[x_check_index])
            if -1 in x_nearby:
                x_nearby.remove(-1)
            for x_nearby_data in x_nearby:
                if x_nearby_data in z_init or x_nearby_data == 0:
                    determined = False
                    break
                elif x_nearby_data in x_init:
                    continue
                else:
                    ## this means the data qubit is unassigned.
                    ## assign it to x
                    x_init.add(x_nearby_data)
                    unassigned_data.remove(x_nearby_data)
            
            if determined:
                #satisfied_checks.add(x_check_index)
                satisfied_checks_order.append(j)
            #x_checks.remove(x_check_index)

        for j, z_check_index in enumerate(z_checks):
            determined = True
            z_nearby = set(self.z_neighbours[z_check_index])
            if -1 in z_nearby:
                z_nearby.remove(-1)
            for z_nearby_data in z_nearby:
                if z_nearby_data in x_init or z_nearby_data == 0:
                    determined = False
                    break
                elif z_nearby_data in z_init:
                    continue
                else:
                    z_init.add(z_nearby_data)
                    unassigned_data.remove(z_nearby_data)
            
            if determined:
                #satisfied_checks.add(z_check_index)
                satisfied_checks_order.append(j + len(x_checks))
            #z_checks.remove()

        ## if there are any data qubits left, we init them in Z
        if len(unassigned_data) != 0:
            z_init = z_init + list(unassigned_data)
        
        ## now construct the initialization circuit
        gates = []
        x_data_str = " ".join(map(str, x_init))
        z_data_str = " ".join(map(str, z_init))

        x_indices_str = " ".join(map(str, self.x_indices)) # x-syndrome-checks
        z_indices_str = " ".join(map(str, self.z_indices)) # z-syndrome-checks
 
        gates.append(f'R 0 {x_data_str} {z_data_str}')
        gates.append(f'H {x_data_str}')
        if basis == 'x' or basis == 'X':
            gates.append(f'H 0')
        elif basis == 'y' or basis == 'Y':
            gates.append(f'H 0')
            gates.append(f'S 0')
        gates.append('TICK')

        gates.append(f'R {x_indices_str} {z_indices_str}')
        gates.append(f'H {x_indices_str}')
        gates.append('TICK')
        ## perform one-rounds of syndrome checks.
        for ts in range(4):
            gate_str = 'CX '
            for x_index in self.x_indices:
                data_index = self.x_neighbours[x_index][ts]
                if data_index != -1:
                    gate_str += f'{x_index} {data_index} '
            
            for z_index in self.z_indices:
                data_index = self.z_neighbours[z_index][ts]
                if data_index != -1:
                    gate_str += f'{data_index} {z_index} '
            gates.append(gate_str)
            gates.append('TICK')
        
        ## now consider applying syndrome measurements
        gates.append('H '+x_indices_str)

        ## perform measurements
        gates.append('M '+ x_indices_str + " " + z_indices_str) 
        
        ## error correction only require the determined ones
        for det_order in satisfied_checks_order:
            gates.append(f'DETECTOR rec[{det_order - (len(x_checks) + len(z_checks))}]')
        
        if output_format == 'str':
            return '\n'.join(gates)
        elif output_format == 'list':
            return gates
        else:
            return gates

    def code_initialization(self, basis, output_format='str'):
        ## currently only supporting X and Z basis w/o problems.
        gate_sequence = []
        x_indices_str = " ".join(map(str, self.x_indices))
        z_indices_str = " ".join(map(str, self.z_indices))
        data_indices_str = " ".join(map(str, self.data_indices))

        ## initialize all qubits into |0> state, both data qubits and syndrome qubits.
        gate_sequence.append('R '+x_indices_str + " "+z_indices_str)
        gate_sequence.append('R '+data_indices_str)
        if basis == 'x' or basis == 'X':
            ## initialize the data qubits into |+> state
            gate_sequence.append('H '+data_indices_str)
        elif basis == 'z' or basis == 'Z':
            pass
        elif basis == 'y' or basis == 'Y':
            print('Y basis is not fault tolerant, using state injection instead.')
            return self.state_injection(basis='y', output_format=output_format)
        # elif basis == 'y' or basis == 'Y':
        #     ## initialize all the data qubits into |+Y> state
        #     #gate_sequence.append('SQRT_X_DAG '+data_indices_str)
        #     x_init_indices = []
        #     y_init_indices = []
        #     for loc in self.data_qubits:
        #         if loc[0] == loc[1]:
        #             ## diagonal qubits are intiilized into physical |+Y> basis
        #             index = self.qubit_index[tuple(loc)]
        #             y_init_indices.append(index)
        #         elif loc[0] < loc[1]:
        #             ## this is the X region, for example (0, 1) will be on X edge
        #             index = self.qubit_index[tuple(loc)]
        #             x_init_indices.append(index)
        #     gate_sequence.append('SQRT_X_DAG '+" ".join(map(str, y_init_indices)))
        #     gate_sequence.append('H '+" ".join(map(str, x_init_indices)))
        else:
            print('basis is unknown, only can take x, z.')
            print('using z instead')
            basis = 'z'

        ## initialize the X checks to be |+> state
        gate_sequence.append("H "+x_indices_str)

        ## perform one-rounds of syndrome checks.
        for ts in range(4):
            gate_str = 'CX '
            for x_index in self.x_indices:
                data_index = self.x_neighbours[x_index][ts]
                if data_index != -1:
                    gate_str += f'{x_index} {data_index} '
            
            for z_index in self.z_indices:
                data_index = self.z_neighbours[z_index][ts]
                if data_index != -1:
                    gate_str += f'{data_index} {z_index} '
            gate_sequence.append(gate_str)
            gate_sequence.append('TICK')
        
        ## now consider applying syndrome measurements
        gate_sequence.append('H '+x_indices_str)

        ## perform measurements
        gate_sequence.append('M '+ x_indices_str + " " + z_indices_str) 
        
        ## now consider comparing the syndrome check results from the last checks
        ## note that for X and Z initialization, we focus on X and Z checks only
        if basis == 'x' or basis == 'X':
            for j in range(self.n_x_checks):
                gate_sequence.append(f'DETECTOR rec[{j - (self.n_x_checks + self.n_z_checks)}]') 
        
        elif basis == 'z' or basis == 'Z':
            for j in range(self.n_z_checks):
                gate_sequence.append(f'DETECTOR rec[{j - self.n_z_checks}]') 
        
        else:
            raise ValueError('FT state preparation can only support x or z.')
        
        if output_format == "list":
            return gate_sequence
        else:
            return "\n".join(gate_sequence)

    def code_initialization_compressed(self, basis, output_format='list'):
        """use only one extra qubits for syndrome measurements, output qasm circuit"""
        gates = []
       
        if basis == 'x' or basis == 'X':
            for data_index in self.data_indices:
                gates.append(f'h q[{data_index}];')
        
        elif basis == 'y' or basis == 'Y':
            for loc in self.data_qubits:
                if loc[0] == loc[1]:
                    ## diagonal qubits are initialized into physical |+Y> basis
                    index = self.qubit_index[tuple(loc)]
                    gates.append(f'sxdg q[{index}];')
                elif loc[0] < loc[1]:
                    ## this is the x region, for example (0, 1) will be in the X operator
                    index = self.qubit_index[tuple(loc)]
                    gates.append(f'h q[{index}];')
                else:
                    ## this is the Z region, no need to apply extra gates.
                    pass
        elif basis == 'z' or basis == 'Z':
            pass
        else:
            print('basis not known. only accept x, y, z.')
            print('using z basis instead.')
            
        syn_gates = self.syndrome_check_block_compressed(round=1, output_format=output_format)
        if output_format == 'list':
            gates += syn_gates
            return gates
        elif output_format == 'str':
            gates = '\n'.join(gates) + '\n' + syn_gates
            return gates
        else:
            raise ValueError('output format not implemented')

    def code_initialization_stab_check(self, basis):
        '''check the stabilizer generator to verify the scheme under ideal cases.'''
        mtx_eye = np.eye(self.n_data, dtype=int)
        mtx_zero = np.zeros((self.n_data, self.n_data), dtype=int)
        mtx_sign = np.zeros((self.n_data, 1), dtype=int)

        if basis == 'x' or basis == 'X':
            ## all data qubits are initialized into |+> states.
            init_mtx = np.bmat([mtx_eye, mtx_zero, mtx_sign])
            init_tab = Tableau(init_mtx)
        
        elif basis == 'y' or basis == 'Y':
            init_mtx = np.zeros((self.n_data, 2*self.n_data+1), dtype=int)
            ct = 0
            for loc in self.data_qubits:
                index = self.qubit_index[tuple(loc)]
                if loc[0] == loc[1]:
                    ## diagonal qubits are intiilized into physical |+Y> basis
                    init_mtx[ct, index] = init_mtx[ct, self.n_data+index] = 1
                elif loc[0] < loc[1]:
                    ## this is the X region, for example (0, 1) will be on X Op
                    init_mtx[ct, index] = 1
                else:
                    ## otherwise, it is in the Z region, for example (1,0) will be on the Z logical op
                    init_mtx[ct, self.n_data + index] = 1
                ct += 1
            init_tab = Tableau(init_mtx)
        else:
            ## all data qubits are initilized into |0> states.
            init_mtx = np.bmat([mtx_zero, mtx_eye, mtx_sign])
            init_tab = Tableau(init_mtx)

            if basis != 'z' or basis != 'Z':
                print('The basis in unknown. Use Z basis instead.')
        
        ## generate the measurement tableau
        n_syn = self.n_x_checks + self.n_z_checks
        measure_mtx = np.zeros((n_syn, 2*self.n_data+1), dtype=int)

        ct = 0
        for x_index in self.x_indices:
            data_indices = self.x_neighbours[x_index]
            for data_index in data_indices:
                if data_index != -1:
                    measure_mtx[ct, data_index] = 1
            #print(ct)
            #print(Tableau(measure_mtx[ct]))
            ct += 1
        
        for z_index in self.z_indices:
            data_indices = self.z_neighbours[z_index]
            for data_index in data_indices:
                if data_index != -1:
                    measure_mtx[ct, self.n_data + data_index] = 1
            #print(ct)
            #print(Tableau(measure_mtx[ct]))
            ct += 1
        measure_tab = Tableau(measure_mtx)

        return init_tab, measure_tab 
        
    def code_stabilizer(self, output_format='stim') -> list:
        '''
        return the stabilizers of the code

        Parameters:
        -----------
        output_format: str
            'stim': for list of stim.PauliString obj (default)
            'tab': for the Tableau class instance
            'list': for the list of strings.
        '''
        stabilizers = []
        for key, val in self.x_neighbours.items():
            stabilizers.append(pauli_string(self.n_data, val, 'X'))
        
        for key, val in self.z_neighbours.items():
            stabilizers.append(pauli_string(self.n_data, val, 'Z'))

        if output_format == 'stim':
            stim_stab = []
            for ps in stabilizers:
                stim_stab.append(stim.PauliString(ps))
            return stim_stab
        
        elif output_format == 'tab':
            stab_tab = []
            for ps in stabilizers:
                stab_tab.append(Tableau.convert_back(ps.replace('_', "I")))
            temp_stab = stab_tab[0]
            for tab in stab_tab[1:]:
                temp_stab.append(tab)
            return temp_stab
        
        else:
            return stabilizers
    
    def logical_ops(self, basis, output_format='stim'):
        '''return the logical operator'''
        loc = np.array([0,0], dtype=int)
        
        if basis == 'z' or basis == 'Z':
            indices = []
            for shift in range(self.distance):
                loc_op = loc + np.array([2*shift, 0], dtype=int)
                indices.append(self.qubit_index[tuple(loc_op)])
            ops = ['Z'] * self.distance

        elif basis == 'x' or basis == 'X':
            indices = []
            for shift in range(self.distance):
                loc_op = loc + np.array([0, 2*shift], dtype=int)
                indices.append(self.qubit_index[tuple(loc_op)])
            ops = ["X"] * self.distance
        
        else:
            indices = []
            for shift in range(1, self.distance):
                loc_op = loc + np.array([0, 2*shift], dtype=int)
                indices.append(self.qubit_index[tuple(loc_op)])
            ops = ["X"] * (self.distance - 1)

            for shift in range(1, self.distance):
                loc_op = loc + np.array([2*shift, 0], dtype=int)
                indices.append(self.qubit_index[tuple(loc_op)])
            ops += ['Z'] * (self.distance - 1)

            indices.append(self.qubit_index[tuple(loc)])
            ops.append('Y')

        ps = pauli_string(self.n_data, indices, ops)
        if output_format == 'stim':
            return stim.PauliString(ps)
        elif output_format == 'tab':
            return Tableau.convert_back(ps)
        elif output_format == 'index':
            return indices
        else:
            return ps
        

    def code_measure_out(self, basis, output_format='list', full_checks=True):
        '''
        the measurement along a given basis

        Parameters:
        -----------
        basis: str
            'x', 'y', or 'z': y-basis has not been implemented yet.
        output_format: str
            'list': list of gates
            'str': a single string
            whether there are full syndrome checks before the measurement circuits
            when round=1, in stim circuit generation, there is no full round checks.
            the number of measurements can be different.
        '''
        #full_check_bit = int(full_checks)
        gates = []
        ## Z or X measurements: measure all qubits into the corresponding basis
        ## then apply a syndrome check using the measured outcomes
        data_indices_str = " ".join(map(str, self.data_indices))

        if basis == 'Y' or basis == 'y':
            ## this needs to be implemented later
            #print('this has not been implemented yet, using Z ')
            #basis = 'z'

            ## the strategy here is:
            ## using the non-fault-tolerant gates to implement a S^dagger gate, 
            ## then the Y is converted into X basis measurements
            logical_z_indices = self.logical_ops(basis='z', output_format='index')
            for j in range(len(logical_z_indices)-1):
                gates.append(f'CX {logical_z_indices[j]} {logical_z_indices[j+1]}')
            
            gates.append("TICK")
            gates.append(f'S_DAG {logical_z_indices[-1]}')
            gates.append("TICK")
            for j in range(len(logical_z_indices)-1, 0, -1):
                gates.append(f'CX {logical_z_indices[j-1]} {logical_z_indices[j]}')
            
            gates.append("TICK")
            #basis = 'x'

        if basis == 'x' or basis == 'X' or basis == 'y' or basis == 'Y':
            gates.append('H ' + data_indices_str)
        
        gates.append('M ' + data_indices_str)

        ## now consider the syndrome checks at the final measurement results.
        loc = np.array([0,0], dtype=int) ## this is for tracking the locations of the data qubits for the logical operators

        # if basis == 'Y' or basis == 'y':
        #     ## need to be implemented
        #     pass
        if basis == 'z' or basis == 'Z':
            for j, z_index in enumerate(self.z_indices):
                data_indices = self.z_neighbours[z_index]
                rec_indices = [k - self.n_data for k in data_indices if k != -1]
                gate_str = f'DETECTOR ' + ' '.join(map(lambda x: f'rec[{x}]', rec_indices)) \
                    + f' rec[{j - (self.n_data + self.n_z_checks)}]'
                gates.append(gate_str)
                
            ## generate the obsrvable
            indices = []
            for shift in range(self.distance):
                loc_op = loc + np.array([2*shift, 0], dtype=int)
                indices.append(self.qubit_index[tuple(loc_op)])
                
        elif basis == 'x' or basis == 'X' or basis == 'y' or basis == 'Y':
        #elif basis == 'x' or basis == 'X':
            for j, x_index in enumerate(self.x_indices):
                data_indices = self.x_neighbours[x_index]
                rec_indices = [k - self.n_data for k in data_indices if k != -1]
                if basis == 'x' or basis == 'X':
                    gate_str = f'DETECTOR ' + ' '.join(map(lambda x: f'rec[{x}]', rec_indices)) \
                        + f' rec[{j - (self.n_data + self.n_z_checks + self.n_x_checks)}]'
                    gates.append(gate_str)

            indices = []
            for shift in range(self.distance):
                loc_op = loc + np.array([0, 2*shift], dtype=int)
                indices.append(self.qubit_index[tuple(loc_op)])
        

        logical_rec_indices = [k - self.n_data for k in indices if k != -1]
        ## add the final observable
        gates.append(f'OBSERVABLE_INCLUDE(0) ' + " ".join(map(lambda x: f'rec[{x}]', logical_rec_indices)))

        if output_format == 'str':
            return '\n'.join(gates)
        else:
            return gates

    #def encoding_circuit(self, output_format='list'):
    #    pass

    def decoding_circuit(self, output_format='list'):
        pass
    
    def stim_circ_init(self, output_format='list'):
        '''generate the informaiton about the qubit location, etc., for stim circuits'''
        gate_sequence = []
        for loc in self.qubit_index:
            index = self.qubit_index[loc]
            gate_sequence.append(f'QUBIT_COORDS({loc[0]}, {loc[1]}) {index}')
        
        if output_format == 'list':
            return gate_sequence
        else:
            return '\n'.join(gate_sequence) + '\n'

    def logical_state_vector(self, basis = 'z'):
        '''generate the state vector for the code basis state |+sigma>'''
        code_stab = self.code_stabilizer(output_format='stim')
        logical_op = self.logical_ops(basis=basis, output_format='stim')
        code_stab.append(logical_op)

        stab_state = stim.Tableau.from_stabilizers(code_stab)

        return stab_state.to_state_vector(endian='little')
    
    def stim_circ_aux_info(self, output_format='list'):
        '''include the coordinate and index information for stim circuits'''
        output = []
        for index, coord in self.qubit_location.items():
            output.append(f'QUBIT_COORDS{tuple(coord)} {index}')
        
        if output_format == 'str':
            return '\n'.join(output)
        else:
            return output
        ''' some qec circuit generation '''
from helper.tab_state import Tableau

import numpy as np
import stim
#from abc import abstractmethod
import math
import matplotlib.pyplot as plt

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
    

def qubit_type(location):
    loc = np.array(location[:2])

    for j,val in enumerate(loc):
        if not math.isclose(val, round(val), abs_tol=1e-8):
            raise NotQubitLocationError(f"location {location} is not a valid qubit location.")
        else:
            loc[j] = round(val)
    
    loc = loc.astype(int)

    ## check if the given location is a data qubit
    if not ((loc[0]&1) | (loc[1]&1)):
        return 'data', loc
    
    loc = loc % 4
    if sum(loc) % 4 == 0:
        ## only (1,3) and (3,1), which are assigned as X
        return 'x', loc
    else:
        return 'z', loc

def is_qubit(location):
    try:
        _ =qubit_type(location)
    except NotQubitLocationError:
        return False
    
    return True

def is_data(location):
    '''check if the given location is a data qubit'''
    try:
        q_type, _ = qubit_type(location)
    except NotQubitLocationError:
        return False
    
    if q_type == 'data':
        return True
    else:
        return False
def is_x_syndrome(location):
    '''check if the given location corresponding to a X-stabilizer check qubits'''
    try:
        q_type, _ = qubit_type(location)
    except NotQubitLocationError:
        return False
    
    if q_type == 'x':
        return True
    else:
        return False

def is_z_syndrome(location):
    '''check if the given location corresponding to a Z-stabilizer check qubits'''
    try:
        q_type, _ = qubit_type(location)
    except NotQubitLocationError:
        return False
    
    if q_type == 'z':
        return True
    else:
        return False

def check_qubits(location, syn_round, q_type='x', consistency_check=False):
    '''
    find the data qubits for the syndrome check gate operations 
    Assuming vertical edges are logical X operator (Z edge), and horizontal edges are Z.

    Parameters:
    -----------
    syn_round: int, (0,1,2,3):
        if syn_round == 0: left-up
        if syn_round == 3: right-down
        if syn_round == 1: if Z check: right-up, if X check: left-down
        if syn_round == 2: if Z check: left-down, if X check: right-up
    '''
    if consistency_check:
        q_type_checking, loc = qubit_type(location)
        if q_type != q_type_checking:
            raise ConsistencyError(f'Inconsistency detected! loc {location} is a {q_type_checking}-type syndrome qubits. \n \
                                   The given q_type is {q_type}-type.')
    loc = np.array(location).astype('int')
    b0 = syn_round & 1
    b1 = (syn_round >> 1)&1
    
    if q_type == 'z':
        loc[0] -= (-1) ** b1
        loc[1] -= (-1) ** b0
    else:
        loc[0] -= (-1) ** b0
        loc[1] -= (-1) ** b1

    return loc


class SurfaceCodeQubits:
    def __init__(self):
        self.qubit_index = {}
        self.qubit_location = {}
    
        self.data_qubits = []
        self.x_qubits = []
        self.z_qubits = []

        self.n_x_checks = 0
        self.n_z_checks = 0
        self.n_data = 0

        self.x_indices = []
        self.z_indices = []
        self.data_indices = []
        self.x_neighbours = {} # x_indices: [r0_data_index, r1_data_index, r2_data_index, r3_data_index]
        self.z_neighbours = {} # z_indices: [r0_data_index, r1_data_index, r2_data_index, r3_data_index]

        self.distance = None

        return 
    
    def _reset(self):

        self.qubit_index = {}
        self.qubit_location = {}
    
        self.data_qubits = []
        self.x_qubits = []
        self.z_qubits = []

        self.n_x_checks = 0
        self.n_z_checks = 0
        self.n_data = 0

        self.x_indices = []
        self.z_indices = []
        self.data_indices = []
        self.x_neighbours = {} # x_indices: [r0_data_index, r1_data_index, r2_data_index, r3_data_index]
        self.z_neighbours = {} # z_indices: [r0_data_index, r1_data_index, r2_data_index, r3_data_index]

        self.distance = None
        return 
 

    def select_active(self, distance = 3):
        """
        Build the lists of active qubits for a rotated surface code patch
        of distance d.

        Data qubits: (even even) sites: (0,0)to (2d-2,2d-2)
        Horizontal edges: X edges (logical Z operator)
        Vertical edges: Z edges (logical X operator)

        Parameters
        ----------
        distance : int, default 3
            code distance for the rotated surface code patch
        """
        if self.distance is not None:
            self._reset()
        self.distance = distance
        
        coord0 = np.array([0, 0], dtype=int)
        for j in range(distance):
            for k in range(distance):
                self.data_qubits.append(coord0 + np.array([j,k]) * 2)
        
        ## now consider the syndrome check qubits inside the qubit region
        coord1 = np.array([1,1], dtype=int)
        for j in range(distance - 1):
            for k in range(distance - 1):
                syn_coord = coord1 + np.array([j, k]) * 2
                if is_x_syndrome(syn_coord):
                    self.x_qubits.append(syn_coord)
                else:
                    self.z_qubits.append(syn_coord)
        
        ## now consider the x_edges:
        ## horizontal edges: need X
        coord1 = np.array([1, -1], dtype=int) ## upper edge
        coord2 = coord1 + np.array([0, distance * 2], dtype = int)
        for j in range(0, distance-1, 2):
            self.x_qubits.append(coord1 + np.array([2*j, 0]))
            self.x_qubits.append(coord2 + np.array([2*(j+1), 0]))
        
        ## vertical edge: need Z
        coord1 =  np.array([-1, 1], dtype=int)
        coord2 = coord1 + np.array([2 * distance, 0], dtype=int)
        for j in range(0, distance-1, 2):
            self.z_qubits.append(coord2 + np.array([0, 2*j]))
            self.z_qubits.append(coord1 + np.array([0, 2*(j+1)]))
        
        ## indexing the qubits:
        self.x_qubits = np.array(self.x_qubits).astype('int')
        self.z_qubits = np.array(self.z_qubits).astype('int')
        self.data_qubits = np.array(self.data_qubits).astype('int')

        self.n_x_checks = len(self.x_qubits)
        self.n_z_checks = len(self.z_qubits)
        self.n_data = len(self.data_qubits)

        self._process()
        return

    def _process(self):
        if not (self.n_x_checks or self.n_z_checks or self.n_data):
            print(f'nx: {self.n_x_checks}, nz: {self.n_z_checks}, nd: {self.n_data}')
            raise NoQubitActiveError('the process is incorrectly applied when there are no active qubits.')

        self.qubit_index = {}
        self.qubit_location = {}
        self.x_indices = []
        self.z_indices = []
        self.data_indices = []

        ct = 0
        for loc in self.data_qubits:
            self.qubit_index[tuple(loc)] = ct
            self.qubit_location[ct] = loc
            self.data_indices.append(ct)
            ct += 1
        for loc in self.x_qubits:
            self.qubit_index[tuple(loc)] = ct
            self.qubit_location[ct] = loc
            self.x_indices.append(ct)
            ct += 1
        for loc in self.z_qubits:
            self.qubit_index[tuple(loc)] = ct
            self.qubit_location[ct] = loc
            self.z_indices.append(ct)
            ct += 1

        for loc in self.x_qubits:
            nearby_data = []
            x_index = self.qubit_index[tuple(loc)]
            for ts in range(4):
                ## solve the nearby data qubits' indices
                ## if not active data qubits found, use -1 for the index
                loc_data = check_qubits(loc, ts, 'x')
                if tuple(loc_data) in self.qubit_index:
                    nearby_data.append(self.qubit_index[tuple(loc_data)])
                else:
                    nearby_data.append(-1)
            self.x_neighbours[x_index] = tuple(nearby_data)
        
        for loc in self.z_qubits:
            nearby_data = []
            z_index = self.qubit_index[tuple(loc)]
            for ts in range(4):
                ## solve the nearby data qubits' indices
                ## if not active data qubits found, use -1 for the index
                loc_data = check_qubits(loc, ts, 'z')
                if tuple(loc_data) in self.qubit_index:
                    nearby_data.append(self.qubit_index[tuple(loc_data)])
                else:
                    nearby_data.append(-1)
            self.z_neighbours[z_index] = tuple(nearby_data)
        return
    

    def draw(self, fontsize: int = 14, show: bool = True):
        """
        Visualise the currently‑selected patch.

        Parameters
        ----------
        fontsize : int
            Font size used for the 'D', 'X', 'Z' labels.
        show : bool
            If True (default) call plt.show() at the end.  Pass False when
            composing figures manually.
        """
        if not (self.n_x_checks or self.n_z_checks or self.n_data):
            raise NoQubitActiveError("No qubits selected: call select_active() first")

        # collect all coordinates; convert numpy arrays to tuples
        coords_D = [tuple(q) for q in self.data_qubits]
        coords_X = [tuple(q) for q in self.x_qubits]
        coords_Z = [tuple(q) for q in self.z_qubits]

        # helper to plot a set of (x,y) points with a text label
        def _plot(coords, label, **kw):
            if not coords:
                return
            xs, ys = zip(*coords)
            plt.scatter(xs, ys, **kw)
            for x, y in coords:
                plt.text(x, y, label+f':{self.qubit_index[(x,y)]}', ha="center", va="center", fontsize=fontsize)

        # Data qubits as squares
        _plot(coords_D, "D", marker="s", s=120, color="white", edgecolors="k")
        # X checks as circles
        _plot(coords_X, "X", marker="o", s=120, color="tab:red", edgecolors="k")
        # Z checks as triangles
        _plot(coords_Z, "Z", marker="^", s=120, color="tab:blue", edgecolors="k")

        plt.gca().set_aspect("equal", adjustable="box")
        plt.grid(True, linestyle="--", linewidth=0.5)
        plt.xlabel("x  (rightwards)")
        plt.ylabel("y  (downwards)")
        plt.title("Surface‑code patch")

        # Invert the y‑axis so that *down* is positive
        plt.gca().invert_yaxis()

        if show:
            plt.show()


    def syndrome_check_block(self, round=3, output_format='str'):
        '''
        generate the syndrome check blocks. No initialization and final measurements.
        the generated circuit will be in the stim format, or QASM.
        '''
        if len(self.data_qubits) == 0:
            raise NoQubitActiveError("Code patch is not initialized. Run select_active before this call.")

        if round < 1:
            raise ValueError("round can only be larger than 1.")
        
        n_total_checks = self.n_x_checks + self.n_z_checks

        ## generating the single round of syndrome check
        gate_sequence = []
        x_indices_str = " ".join(map(str, self.x_indices))
        z_indices_str = " ".join(map(str, self.z_indices))

        ## resettting all the syndrome qubits
        gate_sequence.append("R "+x_indices_str+ " "+z_indices_str)
            
        ## reset the X syndrome into |+> state by applying H gates
        gate_sequence.append("H " + x_indices_str)

        gate_sequence.append("TICK")

        ## now consider performing the syndrome checks
        for ts in range(4):
            gate_str = 'CX '
            for x_index in self.x_indices:
                data_index = self.x_neighbours[x_index][ts]
                if data_index != -1:
                    gate_str += f'{x_index} {data_index} '
            
            for z_index in self.z_indices:
                data_index = self.z_neighbours[z_index][ts]
                if data_index != -1:
                    gate_str += f'{data_index} {z_index} '
            gate_sequence.append(gate_str)
            gate_sequence.append('TICK')
        
        ## now consider applying syndrome measurements
        gate_sequence.append('H '+x_indices_str)

        ## perform measurements
        gate_sequence.append('M '+ x_indices_str + " " + z_indices_str)

        gate_sequence.append('TICK')

        ## now consider comparing the syndrome check results from the last checks
        for j in range(n_total_checks):
            gate_sequence.append(f'DETECTOR rec[-{j+1}] rec[-{j+1+n_total_checks}]')

        if round > 1:
            #gate_str = f"REPEAT {round} {{\n {'\n'.join(gate_sequence)} \n}}\n"
            #gate_str = f"REPEAT {round} " + "{ \n" + f"{'\n'.join(gate_sequence)}" + "\n}\n"
            if output_format == 'str':
                gate_str = f"REPEAT {round} " + \
                        "{\n" + \
                        f"{'\n'.join(gate_sequence)}" + \
                        "\n}\n"
                return gate_str
            else:
                gate_sequence = [f"REPEAT {round} " + "{"] + gate_sequence
                gate_sequence.append('}')
                return gate_sequence
        else:
            if output_format == 'list':
                return gate_sequence
            else:
                gate_str = '\n'.join(gate_sequence) + '\n'

                return gate_str
    
    def syndrome_check_block_compressed(self, round=1, output_format='list'):
        '''generate the same syndrome check block, but in qasm format, and reuse the measurement qubits so we only use one aux qubit.'''
        if len(self.data_qubits) == 0:
            raise NoQubitActiveError("Code patch is not initialized. Run select_active before this call.")

        if round < 1:
            raise ValueError("round can only be larger than 1.")
        
        n_total_checks = self.n_x_checks + self.n_z_checks
        #classical_rec_size = n_total_checks * round
        
        aux_index = self.n_data
        #single_round_gates = [f'creg rec[{classical_rec_size}]']
        single_round_gates = []
        ct = 0
        for x_index in self.x_indices:
            single_round_gates.append(f'reset q[{aux_index}];')
            single_round_gates.append(f'h q[{aux_index}];')
            data_indices = self.x_neighbours[x_index]
            for data_index in data_indices:
                if data_index != -1:
                    single_round_gates.append(f'cx q[{aux_index}] q[{data_index}];')
            single_round_gates.append(f'h q[{aux_index}];')
            single_round_gates.append(f'measure q[{aux_index}] -> rec[{ct}]')
            ct += 1
        
        for z_index in self.z_indices:
            single_round_gates.append(f'reset q[{aux_index}];')
            data_indices = self.z_neighbours[z_index]
            for data_index in data_indices:
                if data_index != -1:
                    single_round_gates.append(f'cx q[{data_index}] q[{aux_index}];')
            single_round_gates.append(f'measure q[{aux_index}] -> rec[{ct}]')
            ct += 1
        
        if round == 1:
            total_gates = single_round_gates
        else:
            total_gates = list(single_round_gates)
            for _ in range(1, round):
                for temp in total_gates:
                    if temp.startswith('measure'):
                        total_gates.append(f'measure q[{aux_index}] -> rec[{ct}]')
                        ct +=1
                    else:
                        total_gates.append(temp)

        if output_format == 'list':
            return total_gates
        elif output_format == 'str':
            return '\n'.join(total_gates)
        else:
            print('unknown output format, using list')
            return total_gates
    
    def state_injection(self, basis, output_format="list"):
        '''non-fault-tolerant state injection'''
        x_checks = self.x_indices
        z_checks = self.z_indices

        #satisfied_checks = set({})
        satisfied_checks_order = [] # record the position in the total checks, used later in qec
        #satisfied_z_checks = []

        # loc_inj = np.array([0,0], dtype=int)

        # ## disable the injection site's nearby syndrome checks
        # for i in (-1, 1):
        #     for j in (-1, 1):
        #         check_loc = tuple(loc_inj + np.array([i,j], dtype=int))
        #         if check_loc in self.qubit_index:
        #             ## this is within the active code patch
        #             ## this syndrome check should be either x or z checks.
        #             index = self.qubit_index
        #             if is_x_syndrome(check_loc):
        #                 x_checks.remove(index)
        #             else:
        #                 z_checks.remove(index)
        #         else:
        #             continue
        
        ## x-logical op edge should be in x_init to start with
        x_init = set(self.logical_ops(basis='x', output_format='index'))
        x_init.remove(0) ## this is the injection site

        ## similar for z
        z_init = set(self.logical_ops(basis='z', output_format='index'))
        z_init.remove(0)

        ## find all un-assinged data qubit indices
        unassigned_data = []
        for x in self.data_indices:
            if (x!=0 and x not in x_init and x not in z_init):
                unassigned_data.append(x)
        unassigned_data = set(unassigned_data)

        ## check all the syndrome checks, if it is possible to make it determined, then do so
        ## check the x_checks
        for j, x_check_index in enumerate(x_checks):
            determined = True
            x_nearby = set(self.x_neighbours[x_check_index])
            if -1 in x_nearby:
                x_nearby.remove(-1)
            for x_nearby_data in x_nearby:
                if x_nearby_data in z_init or x_nearby_data == 0:
                    determined = False
                    break
                elif x_nearby_data in x_init:
                    continue
                else:
                    ## this means the data qubit is unassigned.
                    ## assign it to x
                    x_init.add(x_nearby_data)
                    unassigned_data.remove(x_nearby_data)
            
            if determined:
                #satisfied_checks.add(x_check_index)
                satisfied_checks_order.append(j)
            #x_checks.remove(x_check_index)

        for j, z_check_index in enumerate(z_checks):
            determined = True
            z_nearby = set(self.z_neighbours[z_check_index])
            if -1 in z_nearby:
                z_nearby.remove(-1)
            for z_nearby_data in z_nearby:
                if z_nearby_data in x_init or z_nearby_data == 0:
                    determined = False
                    break
                elif z_nearby_data in z_init:
                    continue
                else:
                    z_init.add(z_nearby_data)
                    unassigned_data.remove(z_nearby_data)
            
            if determined:
                #satisfied_checks.add(z_check_index)
                satisfied_checks_order.append(j + len(x_checks))
            #z_checks.remove()

        ## if there are any data qubits left, we init them in Z
        if len(unassigned_data) != 0:
            z_init = z_init + list(unassigned_data)
        
        ## now construct the initialization circuit
        gates = []
        x_data_str = " ".join(map(str, x_init))
        z_data_str = " ".join(map(str, z_init))

        x_indices_str = " ".join(map(str, self.x_indices)) # x-syndrome-checks
        z_indices_str = " ".join(map(str, self.z_indices)) # z-syndrome-checks
 
        gates.append(f'R 0 {x_data_str} {z_data_str}')
        gates.append(f'H {x_data_str}')
        if basis == 'x' or basis == 'X':
            gates.append(f'H 0')
        elif basis == 'y' or basis == 'Y':
            gates.append(f'H 0')
            gates.append(f'S 0')
        gates.append('TICK')

        gates.append(f'R {x_indices_str} {z_indices_str}')
        gates.append(f'H {x_indices_str}')
        gates.append('TICK')
        ## perform one-rounds of syndrome checks.
        for ts in range(4):
            gate_str = 'CX '
            for x_index in self.x_indices:
                data_index = self.x_neighbours[x_index][ts]
                if data_index != -1:
                    gate_str += f'{x_index} {data_index} '
            
            for z_index in self.z_indices:
                data_index = self.z_neighbours[z_index][ts]
                if data_index != -1:
                    gate_str += f'{data_index} {z_index} '
            gates.append(gate_str)
            gates.append('TICK')
        
        ## now consider applying syndrome measurements
        gates.append('H '+x_indices_str)

        ## perform measurements
        gates.append('M '+ x_indices_str + " " + z_indices_str) 
        
        ## error correction only require the determined ones
        for det_order in satisfied_checks_order:
            gates.append(f'DETECTOR rec[{det_order - (len(x_checks) + len(z_checks))}]')
        
        if output_format == 'str':
            return '\n'.join(gates)
        elif output_format == 'list':
            return gates
        else:
            return gates

    def code_initialization(self, basis, output_format='str'):
        ## currently only supporting X and Z basis w/o problems.
        gate_sequence = []
        x_indices_str = " ".join(map(str, self.x_indices))
        z_indices_str = " ".join(map(str, self.z_indices))
        data_indices_str = " ".join(map(str, self.data_indices))

        ## initialize all qubits into |0> state, both data qubits and syndrome qubits.
        gate_sequence.append('R '+x_indices_str + " "+z_indices_str)
        gate_sequence.append('R '+data_indices_str)
        if basis == 'x' or basis == 'X':
            ## initialize the data qubits into |+> state
            gate_sequence.append('H '+data_indices_str)
        elif basis == 'z' or basis == 'Z':
            pass
        elif basis == 'y' or basis == 'Y':
            print('Y basis is not fault tolerant, using state injection instead.')
            return self.state_injection(basis='y', output_format=output_format)
        # elif basis == 'y' or basis == 'Y':
        #     ## initialize all the data qubits into |+Y> state
        #     #gate_sequence.append('SQRT_X_DAG '+data_indices_str)
        #     x_init_indices = []
        #     y_init_indices = []
        #     for loc in self.data_qubits:
        #         if loc[0] == loc[1]:
        #             ## diagonal qubits are intiilized into physical |+Y> basis
        #             index = self.qubit_index[tuple(loc)]
        #             y_init_indices.append(index)
        #         elif loc[0] < loc[1]:
        #             ## this is the X region, for example (0, 1) will be on X edge
        #             index = self.qubit_index[tuple(loc)]
        #             x_init_indices.append(index)
        #     gate_sequence.append('SQRT_X_DAG '+" ".join(map(str, y_init_indices)))
        #     gate_sequence.append('H '+" ".join(map(str, x_init_indices)))
        else:
            print('basis is unknown, only can take x, z.')
            print('using z instead')
            basis = 'z'

        ## initialize the X checks to be |+> state
        gate_sequence.append("H "+x_indices_str)

        ## perform one-rounds of syndrome checks.
        for ts in range(4):
            gate_str = 'CX '
            for x_index in self.x_indices:
                data_index = self.x_neighbours[x_index][ts]
                if data_index != -1:
                    gate_str += f'{x_index} {data_index} '
            
            for z_index in self.z_indices:
                data_index = self.z_neighbours[z_index][ts]
                if data_index != -1:
                    gate_str += f'{data_index} {z_index} '
            gate_sequence.append(gate_str)
            gate_sequence.append('TICK')
        
        ## now consider applying syndrome measurements
        gate_sequence.append('H '+x_indices_str)

        ## perform measurements
        gate_sequence.append('M '+ x_indices_str + " " + z_indices_str) 
        
        ## now consider comparing the syndrome check results from the last checks
        ## note that for X and Z initialization, we focus on X and Z checks only
        if basis == 'x' or basis == 'X':
            for j in range(self.n_x_checks):
                gate_sequence.append(f'DETECTOR rec[{j - (self.n_x_checks + self.n_z_checks)}]') 
        
        elif basis == 'z' or basis == 'Z':
            for j in range(self.n_z_checks):
                gate_sequence.append(f'DETECTOR rec[{j - self.n_z_checks}]') 
        
        else:
            raise ValueError('FT state preparation can only support x or z.')
        
        if output_format == "list":
            return gate_sequence
        else:
            return "\n".join(gate_sequence)

    def code_initialization_compressed(self, basis, output_format='list'):
        """use only one extra qubits for syndrome measurements, output qasm circuit"""
        gates = []
       
        if basis == 'x' or basis == 'X':
            for data_index in self.data_indices:
                gates.append(f'h q[{data_index}];')
        
        elif basis == 'y' or basis == 'Y':
            for loc in self.data_qubits:
                if loc[0] == loc[1]:
                    ## diagonal qubits are initialized into physical |+Y> basis
                    index = self.qubit_index[tuple(loc)]
                    gates.append(f'sxdg q[{index}];')
                elif loc[0] < loc[1]:
                    ## this is the x region, for example (0, 1) will be in the X operator
                    index = self.qubit_index[tuple(loc)]
                    gates.append(f'h q[{index}];')
                else:
                    ## this is the Z region, no need to apply extra gates.
                    pass
        elif basis == 'z' or basis == 'Z':
            pass
        else:
            print('basis not known. only accept x, y, z.')
            print('using z basis instead.')
            
        syn_gates = self.syndrome_check_block_compressed(round=1, output_format=output_format)
        if output_format == 'list':
            gates += syn_gates
            return gates
        elif output_format == 'str':
            gates = '\n'.join(gates) + '\n' + syn_gates
            return gates
        else:
            raise ValueError('output format not implemented')

    def code_initialization_stab_check(self, basis):
        '''check the stabilizer generator to verify the scheme under ideal cases.'''
        mtx_eye = np.eye(self.n_data, dtype=int)
        mtx_zero = np.zeros((self.n_data, self.n_data), dtype=int)
        mtx_sign = np.zeros((self.n_data, 1), dtype=int)

        if basis == 'x' or basis == 'X':
            ## all data qubits are initialized into |+> states.
            init_mtx = np.bmat([mtx_eye, mtx_zero, mtx_sign])
            init_tab = Tableau(init_mtx)
        
        elif basis == 'y' or basis == 'Y':
            init_mtx = np.zeros((self.n_data, 2*self.n_data+1), dtype=int)
            ct = 0
            for loc in self.data_qubits:
                index = self.qubit_index[tuple(loc)]
                if loc[0] == loc[1]:
                    ## diagonal qubits are intiilized into physical |+Y> basis
                    init_mtx[ct, index] = init_mtx[ct, self.n_data+index] = 1
                elif loc[0] < loc[1]:
                    ## this is the X region, for example (0, 1) will be on X Op
                    init_mtx[ct, index] = 1
                else:
                    ## otherwise, it is in the Z region, for example (1,0) will be on the Z logical op
                    init_mtx[ct, self.n_data + index] = 1
                ct += 1
            init_tab = Tableau(init_mtx)
        else:
            ## all data qubits are initilized into |0> states.
            init_mtx = np.bmat([mtx_zero, mtx_eye, mtx_sign])
            init_tab = Tableau(init_mtx)

            if basis != 'z' or basis != 'Z':
                print('The basis in unknown. Use Z basis instead.')
        
        ## generate the measurement tableau
        n_syn = self.n_x_checks + self.n_z_checks
        measure_mtx = np.zeros((n_syn, 2*self.n_data+1), dtype=int)

        ct = 0
        for x_index in self.x_indices:
            data_indices = self.x_neighbours[x_index]
            for data_index in data_indices:
                if data_index != -1:
                    measure_mtx[ct, data_index] = 1
            #print(ct)
            #print(Tableau(measure_mtx[ct]))
            ct += 1
        
        for z_index in self.z_indices:
            data_indices = self.z_neighbours[z_index]
            for data_index in data_indices:
                if data_index != -1:
                    measure_mtx[ct, self.n_data + data_index] = 1
            #print(ct)
            #print(Tableau(measure_mtx[ct]))
            ct += 1
        measure_tab = Tableau(measure_mtx)

        return init_tab, measure_tab 
        
    def code_stabilizer(self, output_format='stim') -> list:
        '''
        return the stabilizers of the code

        Parameters:
        -----------
        output_format: str
            'stim': for list of stim.PauliString obj (default)
            'tab': for the Tableau class instance
            'list': for the list of strings.
        '''
        stabilizers = []
        for key, val in self.x_neighbours.items():
            stabilizers.append(pauli_string(self.n_data, val, 'X'))
        
        for key, val in self.z_neighbours.items():
            stabilizers.append(pauli_string(self.n_data, val, 'Z'))

        if output_format == 'stim':
            stim_stab = []
            for ps in stabilizers:
                stim_stab.append(stim.PauliString(ps))
            return stim_stab
        
        elif output_format == 'tab':
            stab_tab = []
            for ps in stabilizers:
                stab_tab.append(Tableau.convert_back(ps.replace('_', "I")))
            temp_stab = stab_tab[0]
            for tab in stab_tab[1:]:
                temp_stab.append(tab)
            return temp_stab
        
        else:
            return stabilizers
    
    def logical_ops(self, basis, output_format='stim'):
        '''return the logical operator'''
        loc = np.array([0,0], dtype=int)
        
        if basis == 'z' or basis == 'Z':
            indices = []
            for shift in range(self.distance):
                loc_op = loc + np.array([2*shift, 0], dtype=int)
                indices.append(self.qubit_index[tuple(loc_op)])
            ops = ['Z'] * self.distance

        elif basis == 'x' or basis == 'X':
            indices = []
            for shift in range(self.distance):
                loc_op = loc + np.array([0, 2*shift], dtype=int)
                indices.append(self.qubit_index[tuple(loc_op)])
            ops = ["X"] * self.distance
        
        else:
            indices = []
            for shift in range(1, self.distance):
                loc_op = loc + np.array([0, 2*shift], dtype=int)
                indices.append(self.qubit_index[tuple(loc_op)])
            ops = ["X"] * (self.distance - 1)

            for shift in range(1, self.distance):
                loc_op = loc + np.array([2*shift, 0], dtype=int)
                indices.append(self.qubit_index[tuple(loc_op)])
            ops += ['Z'] * (self.distance - 1)

            indices.append(self.qubit_index[tuple(loc)])
            ops.append('Y')

        ps = pauli_string(self.n_data, indices, ops)
        if output_format == 'stim':
            return stim.PauliString(ps)
        elif output_format == 'tab':
            return Tableau.convert_back(ps)
        elif output_format == 'index':
            return indices
        else:
            return ps
        

    def code_measure_out(self, basis, output_format='list', full_checks=True):
        '''
        the measurement along a given basis

        Parameters:
        -----------
        basis: str
            'x', 'y', or 'z': y-basis has not been implemented yet.
        output_format: str
            'list': list of gates
            'str': a single string
            whether there are full syndrome checks before the measurement circuits
            when round=1, in stim circuit generation, there is no full round checks.
            the number of measurements can be different.
        '''
        #full_check_bit = int(full_checks)
        gates = []
        ## Z or X measurements: measure all qubits into the corresponding basis
        ## then apply a syndrome check using the measured outcomes
        data_indices_str = " ".join(map(str, self.data_indices))

        if basis == 'Y' or basis == 'y':
            ## this needs to be implemented later
            #print('this has not been implemented yet, using Z ')
            #basis = 'z'

            ## the strategy here is:
            ## using the non-fault-tolerant gates to implement a S^dagger gate, 
            ## then the Y is converted into X basis measurements
            logical_z_indices = self.logical_ops(basis='z', output_format='index')
            for j in range(len(logical_z_indices)-1):
                gates.append(f'CX {logical_z_indices[j]} {logical_z_indices[j+1]}')
            
            gates.append("TICK")
            gates.append(f'S_DAG {logical_z_indices[-1]}')
            gates.append("TICK")
            for j in range(len(logical_z_indices)-1, 0, -1):
                gates.append(f'CX {logical_z_indices[j-1]} {logical_z_indices[j]}')
            
            gates.append("TICK")
            #basis = 'x'

        if basis == 'x' or basis == 'X' or basis == 'y' or basis == 'Y':
            gates.append('H ' + data_indices_str)
        
        gates.append('M ' + data_indices_str)

        ## now consider the syndrome checks at the final measurement results.
        loc = np.array([0,0], dtype=int) ## this is for tracking the locations of the data qubits for the logical operators

        # if basis == 'Y' or basis == 'y':
        #     ## need to be implemented
        #     pass
        if basis == 'z' or basis == 'Z':
            for j, z_index in enumerate(self.z_indices):
                data_indices = self.z_neighbours[z_index]
                rec_indices = [k - self.n_data for k in data_indices if k != -1]
                gate_str = f'DETECTOR ' + ' '.join(map(lambda x: f'rec[{x}]', rec_indices)) \
                    + f' rec[{j - (self.n_data + self.n_z_checks)}]'
                gates.append(gate_str)
                
            ## generate the obsrvable
            indices = []
            for shift in range(self.distance):
                loc_op = loc + np.array([2*shift, 0], dtype=int)
                indices.append(self.qubit_index[tuple(loc_op)])
                
        elif basis == 'x' or basis == 'X' or basis == 'y' or basis == 'Y':
        #elif basis == 'x' or basis == 'X':
            for j, x_index in enumerate(self.x_indices):
                data_indices = self.x_neighbours[x_index]
                rec_indices = [k - self.n_data for k in data_indices if k != -1]
                if basis == 'x' or basis == 'X':
                    gate_str = f'DETECTOR ' + ' '.join(map(lambda x: f'rec[{x}]', rec_indices)) \
                        + f' rec[{j - (self.n_data + self.n_z_checks + self.n_x_checks)}]'
                    gates.append(gate_str)

            indices = []
            for shift in range(self.distance):
                loc_op = loc + np.array([0, 2*shift], dtype=int)
                indices.append(self.qubit_index[tuple(loc_op)])
        

        logical_rec_indices = [k - self.n_data for k in indices if k != -1]
        ## add the final observable
        gates.append(f'OBSERVABLE_INCLUDE(0) ' + " ".join(map(lambda x: f'rec[{x}]', logical_rec_indices)))

        if output_format == 'str':
            return '\n'.join(gates)
        else:
            return gates

    #def encoding_circuit(self, output_format='list'):
    #    pass

    def decoding_circuit(self, output_format='list'):
        pass
    
    def stim_circ_init(self, output_format='list'):
        '''generate the informaiton about the qubit location, etc., for stim circuits'''
        gate_sequence = []
        for loc in self.qubit_index:
            index = self.qubit_index[loc]
            gate_sequence.append(f'QUBIT_COORDS({loc[0]}, {loc[1]}) {index}')
        
        if output_format == 'list':
            return gate_sequence
        else:
            return '\n'.join(gate_sequence) + '\n'

    def logical_state_vector(self, basis = 'z'):
        '''generate the state vector for the code basis state |+sigma>'''
        code_stab = self.code_stabilizer(output_format='stim')
        logical_op = self.logical_ops(basis=basis, output_format='stim')
        code_stab.append(logical_op)

        stab_state = stim.Tableau.from_stabilizers(code_stab)

        return stab_state.to_state_vector(endian='little')
    
    def stim_circ_aux_info(self, output_format='list'):
        '''include the coordinate and index information for stim circuits'''
        output = []
        for index, coord in self.qubit_location.items():
            output.append(f'QUBIT_COORDS{tuple(coord)} {index}')
        
        if output_format == 'str':
            return '\n'.join(output)
        else:
            return output
        

    def state_injection_part(
        self,
        basis: str,
        injection_method: str = "flexion",
        output_format: str = "list",
    ) -> Union[List[str], str]:
        """
        Returns the Stim lines doing a single injection‐check for the given basis.
        
        Parameters
        ----------
        basis : {'x','y','z'}
            Which logical basis to inject.
        injection_method : str
            Flavor of injection‐check.  Currently only 'flexion' is implemented.
        output_format : {'list','str'}
            Whether to return a Python list of strings, or a single '\n'.join‐ed string.
        """
        method = injection_method.lower()
        if method != "flexion":
            raise NotImplementedError(f"Injection method '{injection_method}' not supported")

        lines: List[str] = self.stim_circ_init(output_format="list")

        lines += self.code_initialization(basis=basis, output_format="list")

        lines += self.syndrome_check_block(round=1, output_format="list")

        if output_format == "list":
            return lines
        else:
            return "\n".join(lines) + "\n"
    
    def build_flexion_rounds(
        self,
        basis: str,
        injection_method: str = "flexion",
        rounds: int = 1,
        output_format: str = "list",
    ) -> Union[List[str], stim.Circuit]:
        """
        1) Do a single FT injection-check for `basis` (via state_injection_part).
        2) Append `rounds` copies of one syndrome_check_block.
        3) Append the final code_measure_out for the same basis.
        """
        # 1) injection + check
        lines = self.state_injection_part(
            basis=basis,
            injection_method=injection_method,
            output_format="list",
        )

        # 2) main QEC rounds
        synd_part = self.syndrome_check_block(round=1, output_format="list")
        for _ in range(rounds):
            lines += synd_part

        # 3) final readout
        m_part = self.code_measure_out(basis=basis, output_format="list")
        lines += m_part

        if output_format == "list":
            return lines

        # otherwise wrap into a stim.Circuit
        return stim.Circuit("\n".join(lines))
    

    def encoding_circuit(
        self,
        inj_basis: str,
        meas_basis: str,
        injection_method: str = "flexion",
        rounds: int = 1,
        output_format: str = "list",
    ) -> Union[List[str], stim.Circuit]:
        """
        1) inject logical |inj_basis>_L  (via state_injection_part)
        2) do `rounds` syndrome rounds
        3) measure in `meas_basis`  (via code_measure_out)

        Parameters
        ----------
        inj_basis : {'x','y','z'}
            which logical basis state to prepare
        meas_basis : {'x','y','z'}
            which basis to rotate into before the final measurement
        injection_method : str
        rounds : int
        output_format : {'list','stim'}
        """

        lines = self.state_injection_part(
            basis=inj_basis,
            injection_method=injection_method,
            output_format="list",
        )

        # 2) repeating syndrome rounds
        synd = self.syndrome_check_block(round=rounds, output_format="list")
        lines += synd

        # 3) final measurement in meas_basis
        m_part = self.code_measure_out(basis=meas_basis, output_format="list")
        lines += m_part

        if output_format == "list":
            return lines
        else:
            return stim.Circuit("\n".join(lines))



    def noisy_encoding_circuit(
        self,
        inj_basis: str,
        meas_basis: str,
        injection_method: str,
        rounds: int,
        p1: float,
        p2: float,
        p_meas: float,
        *,
        skip_prep_noise: bool = False,
        skip_meas_noise: bool = False,
        output_format: str = "list",
    ) -> Union[List[str], stim.Circuit]:
        prep_block = self.state_injection_part(
            basis=inj_basis,
            injection_method=injection_method,
            output_format="list",
        )

        qec_block = self.syndrome_check_block(round=rounds, output_format="list")

        meas_block = self.code_measure_out(
            basis=meas_basis, output_format="list"
        )

        one_q = {"R", "H", "S", "S_DAG", "SQRT_X", "SQRT_X_DAG"}
        two_q = {"CX", "CNOT"}

        noisy: List[str] = []

        for line in prep_block:
            noisy.append(line)
            if not skip_prep_noise:
                op, *targs = line.split()
                if op in one_q:
                    noisy.append(f"DEPOLARIZE1({p1}) " + " ".join(targs))
                elif op in two_q:
                    noisy.append(f"DEPOLARIZE2({p2}) " + " ".join(targs))

        for line in qec_block:
            op, *targs = line.split()
            if op == "M":
                noisy.append(f"X_ERROR({p_meas}) " + " ".join(targs))
                noisy.append(line)
            else:
                noisy.append(line)
                if op in one_q:
                    noisy.append(f"DEPOLARIZE1({p1}) " + " ".join(targs))
                elif op in two_q:
                    noisy.append(f"DEPOLARIZE2({p2}) " + " ".join(targs))

        for line in meas_block:
            noisy.append(line)
            if not skip_meas_noise:
                op, *targs = line.split()
                if op == "M":
                    noisy.append(f"X_ERROR({p_meas}) " + " ".join(targs))
                elif op in one_q:
                    noisy.append(f"DEPOLARIZE1({p1}) " + " ".join(targs))
                elif op in two_q:
                    noisy.append(f"DEPOLARIZE2({p2}) " + " ".join(targs))

        if output_format == "list":
             return noisy
        else:
             return stim.Circuit("\n".join(noisy))
        

def pauli_string(nq:int, bits:list, ps, sign:bool=True):
    '''
    generate a pauli-string 

    Parameters:
    -----------
    nq: int
        total number of qubits
    bits: list
        the list of indices for nontrivial pauli operators
    ps: list or str
        if str: all nontrivial pauli oeprators will be the given ps
        if list: should matching the length of bits
    sign: bool:
        true: + sign, false: - sign
    
    Returns:
    -------
    p_string: str
        pauli string with sign
    '''
    temp = ['_'] * nq
    if isinstance(ps, str):
        ps = [ps] * len(bits)

    for j, s in enumerate(bits):
        if s >= 0:
            temp[s] = ps[j]
    
    if sign:
        return f'+{''.join(temp)}'
    else:
        return f'-{''.join(temp)}'



def surface_code_stabs(sc_qubits:SurfaceCodeQubits):

    stabilizers = []
    for x_index in sc_qubits.x_indices:
        x_nn = sc_qubits.x_neighbours[x_index]
        stab = ['_'] * sc_qubits.n_data
        for j in x_nn:
            if j == -1:
                continue
            else:
                stab[j] = "X"
        stabilizers.append('+'+''.join(stab))

    for z_index in sc_qubits.z_indices:
        z_nn = sc_qubits.z_neighbours[z_index]
        stab = ['_'] * sc_qubits.n_data
        for j in z_nn:
            if j == -1:
                continue
            else:
                stab[j] = "Z"
        stabilizers.append('+'+''.join(stab))
    
    return [stim.PauliString(p) for p in stabilizers]
    
def surface_code_logical(sc_qubits:SurfaceCodeQubits, log_op = 'z'):
    loc0 = np.array([0,0], dtype=int)
    if log_op == 'z' or log_op == 'Z' or log_op == 'y' or log_op == 'Y':
        z_index = []
        for j in range(sc_qubits.distance):
            loc = loc0 + np.array([2*j, 0], dtype=int)
            z_index.append(sc_qubits.qubit_index[tuple(loc)])
        
        op = ['_'] * sc_qubits.n_data
        for k in z_index:
            op[k] = 'Z'
        
        z_op = stim.PauliString("+" + ''.join(op))

        if log_op == 'z' or log_op == 'Z':
            return z_op

    if log_op == 'x' or log_op == 'X' or log_op == 'y' or log_op =='Y':
        x_index = []
        for j in range(sc_qubits.distance):
            loc = loc0 + np.array([0, 2*j], dtype=int)
            x_index.append(sc_qubits.qubit_index[tuple(loc)])
        
        op = ['_'] * sc_qubits.n_data

        for k in x_index:
            op[k] = "X"
        
        x_op = stim.PauliString("+"+"".join(op))

    if log_op == "x" or log_op == 'X':
        return x_op

    else:
        return 1j * x_op * z_op
        

    

    
        

def pauli_string(nq:int, bits:list, ps, sign:bool=True):
    '''
    generate a pauli-string 

    Parameters:
    -----------
    nq: int
        total number of qubits
    bits: list
        the list of indices for nontrivial pauli operators
    ps: list or str
        if str: all nontrivial pauli oeprators will be the given ps
        if list: should matching the length of bits
    sign: bool:
        true: + sign, false: - sign
    
    Returns:
    -------
    p_string: str
        pauli string with sign
    '''
    temp = ['_'] * nq
    if isinstance(ps, str):
        ps = [ps] * len(bits)

    for j, s in enumerate(bits):
        if s >= 0:
            temp[s] = ps[j]
    
    if sign:
        return f'+{''.join(temp)}'
    else:
        return f'-{''.join(temp)}'



def surface_code_stabs(sc_qubits:SurfaceCodeQubits):

    stabilizers = []
    for x_index in sc_qubits.x_indices:
        x_nn = sc_qubits.x_neighbours[x_index]
        stab = ['_'] * sc_qubits.n_data
        for j in x_nn:
            if j == -1:
                continue
            else:
                stab[j] = "X"
        stabilizers.append('+'+''.join(stab))

    for z_index in sc_qubits.z_indices:
        z_nn = sc_qubits.z_neighbours[z_index]
        stab = ['_'] * sc_qubits.n_data
        for j in z_nn:
            if j == -1:
                continue
            else:
                stab[j] = "Z"
        stabilizers.append('+'+''.join(stab))
    
    return [stim.PauliString(p) for p in stabilizers]
    
def surface_code_logical(sc_qubits:SurfaceCodeQubits, log_op = 'z'):
    loc0 = np.array([0,0], dtype=int)
    if log_op == 'z' or log_op == 'Z' or log_op == 'y' or log_op == 'Y':
        z_index = []
        for j in range(sc_qubits.distance):
            loc = loc0 + np.array([2*j, 0], dtype=int)
            z_index.append(sc_qubits.qubit_index[tuple(loc)])
        
        op = ['_'] * sc_qubits.n_data
        for k in z_index:
            op[k] = 'Z'
        
        z_op = stim.PauliString("+" + ''.join(op))

        if log_op == 'z' or log_op == 'Z':
            return z_op

    if log_op == 'x' or log_op == 'X' or log_op == 'y' or log_op =='Y':
        x_index = []
        for j in range(sc_qubits.distance):
            loc = loc0 + np.array([0, 2*j], dtype=int)
            x_index.append(sc_qubits.qubit_index[tuple(loc)])
        
        op = ['_'] * sc_qubits.n_data

        for k in x_index:
            op[k] = "X"
        
        x_op = stim.PauliString("+"+"".join(op))

    if log_op == "x" or log_op == 'X':
        return x_op

    else:
        return 1j * x_op * z_op
        
import numpy as np
import stim

def stim_state_tomog(
    distance: int,
    p1: float,
    p2: float,
    p_meas: float,
    inj_basis: str,
    shots: int = 500_000,
    skip_prep_noise: bool = False,
    skip_meas_noise: bool = False,
    injection_method: str = "flexion",
    draw_timelines: bool = False,
):
    """
    Prepare and tomograph a logical Surface-Code qubit with depolarizing + measurement noise.

    Parameters
    ----------
    distance : int
        Surface-code distance (also number of rounds).
    p1 : float
        Single-qubit depolarizing probability.
    p2 : float
        Two-qubit depolarizing probability.
    p_meas : float
        Measurement bit-flip probability.
    inj_basis : {'x','y','z'}
        Which logical basis to inject (|+>_L, |i>_L or |0>_L).
    shots : int
        Number of samples per circuit.
    injection_method : str
        Passed through to noisy_encoding_circuit (e.g. 'flexion').
    draw_timelines : bool
        If True, call .diagram('timeline-svg') on each circuit.

    Returns
    -------
    rho : 2×2 complex ndarray
        Reconstructed single-qubit density matrix.
    bloch : dict
        Bloch components {'X':…, 'Y':…, 'Z':…}.
    circuits : dict
        The three Stim circuits {'z':…, 'x':…, 'y':…}.
    """
    # 1) build circuits
    sc = SurfaceCodeQubits()
    sc.select_active(distance=distance)

    circuits = {}
    for meas_basis in ("z","x","y"):
        c = sc.noisy_encoding_circuit(
            inj_basis=inj_basis,
            meas_basis=meas_basis,
            injection_method=injection_method,
            rounds=distance,
            p1=p1, p2=p2, p_meas=p_meas,
            skip_prep_noise=skip_prep_noise,
            skip_meas_noise=skip_meas_noise,
            output_format="stim",
        )

        if draw_timelines:
            c.diagram("timeline-svg")
        circuits[meas_basis] = c

    # 2) tomography
    def estimate_expectation(logical_bits: np.ndarray) -> float:
        return np.mean(1 - 2*logical_bits)

    def tomo_rho(
        circ_z: stim.Circuit,
        circ_x: stim.Circuit,
        circ_y: stim.Circuit,
        shots: int
    ):
        sampler_z = circ_z.compile_detector_sampler()
        sampler_x = circ_x.compile_detector_sampler()
        sampler_y = circ_y.compile_detector_sampler()

        samps_z = sampler_z.sample(shots, append_observables=True)
        samps_x = sampler_x.sample(shots, append_observables=True)
        samps_y = sampler_y.sample(shots, append_observables=True)

        bits_z = samps_z[:, -1]
        bits_x = samps_x[:, -1]
        bits_y = samps_y[:, -1]

        ez = estimate_expectation(bits_z)
        ex = estimate_expectation(bits_x)
        ey = estimate_expectation(bits_y)

        rho = 0.5 * np.array([
            [1 + ez,     ex - 1j*ey],
            [ex + 1j*ey, 1 - ez    ],
        ])
        return rho, {"X": ex, "Y": ey, "Z": ez}

    rho, bloch = tomo_rho(
        circuits["z"],
        circuits["x"],
        circuits["y"],
        shots
    )

    return rho, bloch, circuits

    
def stim_process_tomog(
    distance: int,
    p1: float,
    p2: float,
    p_meas: float,
    shots: int = 500_000,
    skip_prep_noise: bool = False,
    skip_meas_noise: bool = False,
    injection_method: str = "flexion",
    draw_timelines: bool = False,
):
    """
    Perform process tomography of the logical channel implemented by the
    noisy surface‐code memory.  Returns:

      M : 3×3 real ndarray
          Pauli‐transfer matrix so that r_out = M r_in on Bloch vectors.
      blochs : dict[str→dict]
          For each inj_basis in {'x','y','z'}, the reconstructed Bloch
          components {'X':…, 'Y':…, 'Z':…}.
      circuits : dict[str→dict]
          The Stim circuits used for each measurement basis.
    """
    import numpy as np


    basis_order = ["x", "y", "z"]
    blochs = {}
    circuits = {}

    M = np.zeros((3, 3), dtype=float)

    for j, inj in enumerate(basis_order):
        rho, bloch, circs = stim_state_tomog(
            distance=distance,
            p1=p1,
            p2=p2,
            p_meas=p_meas,
            inj_basis=inj,
            shots=shots,
            skip_prep_noise=skip_prep_noise,
            skip_meas_noise=skip_meas_noise,
            injection_method=injection_method,
            draw_timelines=draw_timelines,
        )
        blochs[inj] = bloch
        circuits[inj] = circs

        M[:, j] = [bloch["X"], bloch["Y"], bloch["Z"]]

    return M, blochs, circuits

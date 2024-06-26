from mpi4py import MPI
import sys
sys.path.insert(0, 'build/vqe')
import nwqflow
from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.units import DistanceUnit
from qiskit_nature.second_q.transformers import ActiveSpaceTransformer
from qiskit_nature.second_q.mappers import JordanWignerMapper
import numpy as np

comm = MPI.COMM_WORLD
# rank = MPI
# Example Hamiltonian from PySCF
driver = PySCFDriver(
    atom='Li 0.0, 0.0, 0.0; H 0.0, 0.0 1.2',
    basis='sto-3g',
    charge=0,
    spin=0,
    unit=DistanceUnit.ANGSTROM
)
problem = driver.run()
problem = ActiveSpaceTransformer(problem.num_particles, 5, range(5)).transform(problem)

# Convert the PySCF format to XACC (e.g. "+_3 +_2 -_0 -_1" -> 3^ 2^ 0 1)
oplist = []
for op, coeff in problem.second_q_ops()[0].items():
    new_op = []
    for o in op.split():
        ostr = o[2:]
        if o[0:2] == '+_':
            ostr += "^"
        new_op.append(ostr)
    oplist.append((' '.join(new_op), coeff))

n_params = nwqflow.get_param_count(num_particles=4, num_spatial_orbitals=5)  # number of optimization parameters
initial_point = np.random.rand(n_params)  # initial parameter point
nwqflow.optimize_effective_hamiltonian(oplist, 4, initial_point, backend='MPI', xacc=True, comm=comm)
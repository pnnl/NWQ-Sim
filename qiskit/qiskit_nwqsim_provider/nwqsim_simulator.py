"""Backend for NWQSim simulation environment."""

import json
import time
import uuid

from qiskit.providers import BackendV2 as Backend
from qiskit.transpiler import Target
from qiskit.providers import Options
from qiskit.result import Result
from qiskit_nwqsim_provider.nwqsim_job import NWQSimJob
from qiskit import qobj as qobj_mod
from qiskit.compiler import assemble
import subprocess
import importlib_resources # for finding the nwq_qasm 
import os

#======================================== NWQSim Simulator ================================
class NWQSimSimulator(Backend):
    def __init__(self):
        super().__init__()
        # Create Target
        self._target = Target("NWQSim Target")
        # Set option validators
        self.options.set_validator("shots", (1, 65536))

    configuration = {
        'backend_name': 'nwqsim',
        'backend_version': '2.0',
        'n_qubits': 30,
        'basis_gates': ['x', 'y', 'z', 'h', 's', 'sdg', 't', 'tdg', 'ri', 'rx', 'ry', 'rz', 'sx', 'p', 'u', 'cx', 'cy', 'cz', 'ch', 'cs', 'csdg', 'ct', 'ctdg', 'crx', 'cry', 'crz', 'csx', 'cp', 'cu', 'id', 'swap', 'm', 'ma', 'reset', 'u1', 'u2', 'u3', 'ccx', 'cswap', 'rxx', 'ryy', 'rzz'],
        'coupling_map': None,
        'url': 'https://github.com/pnnl/nwq-sim',
        'simulator': True,
        'local': True,
        'conditional': True,
        'open_pulse': False,
        'memory': True,
        'max_shots': 65536,
        'description': 'NWQSim simulation environment (including DM-Sim and SV-Sim)',
        'gates': []
    }

    @property
    def target(self):
        return self._target

    @property
    def max_circuits(self):
        return 1024

    @classmethod
    def _default_options(cls):
        return Options(shots=1024, seed=334)

    def run_experiment(self, experiment, options):
        if os.path.isfile("./build/qasm/nwq_qasm"): # meaning we are in the repo/qiskit folder
            nwq_qasm_path = "./build/qasm/nwq_qasm"
        else:
            nwq_qasm_path = str(importlib_resources.files("qiskit_nwqsim_provider").parent.joinpath("qasm/nwq_qasm"))
        # print(" (DEBUG)", nwq_qasm_path)
        # if not os.path.isfile(nwq_qasm_path):
        #     Exception("Please use the absolute path for nwq_qis_path to the cloned repo.")
        qasmjson = json.dumps(experiment.to_dict())
        if len(qasmjson) > 1000000:
            qasmjsonbin = qasmjson.encode('utf-8')
            with open("expriment.json","wb") as outfile:
                outfile.write(qasmjsonbin)
            # cmd = str( "./build/qasm/nwq_qasm -j expriment.json") # shots parameter is not included
            cmd = str( nwq_qasm_path+" -shots {:d} -j expriment.json").format(options['shots'])
        else:
            # cmd = str( "./build/qasm/nwq_qasm -js '") + qasmjson + "'"  # shots parameter is not included
            cmd = str( nwq_qasm_path+" -shots {:d} -js '").format(options['shots']) + qasmjson + "'" 
        #print (cmd)
        output = subprocess.getoutput(cmd)
        #print(output)
        pos = output.find('nwq_sim_counts=')
        res_counts = json.loads(output[output.find('{',pos):output.find('}',pos)+1])
        pos = output.find('state_vector=')
        vals = output[output.find('[',pos)+1:output.find(']',pos)-2] #last is "x+yj, ]" that's why -2
        state_vector = []
        for val in vals.split(','):
            components = val.split('+')
            state_vector.append(complex(float(components[0]),float(components[1][:-1]))) #-1 is to take out "j"
        #print (output)
        self.start_time = time.time()
        self.end_time = time.time()
        result_dict = {'header': {'name': experiment.header.name,
                                  'memory_slots': experiment.config.memory_slots,
                                  'creg_sizes': experiment.header.creg_sizes
                                  },
                       'status': 'DONE', 'time_taken': self.end_time - self.start_time,
                       'seed': options['seed'], 
                       'shots': options['shots'],
                       'data': {'counts': res_counts, 'statevector':state_vector},
                       'success': True
                       }
        return result_dict

    def run(self, circuits, **kwargs):
        # serialize circuits submit to backend and create a job
        for kwarg in kwargs:
            # if not hasattr(kwarg, self.options): # Reversed, see https://docs.python.org/3.8/library/functions.html#hasattr
            if not hasattr(self.options, kwarg):
                warnings.warn(
                    "Option %s is not used by this backend" % kwarg, UserWarning, stacklevel=2)
        options = {
            'shots': kwargs.get('shots', self.options.shots),
            'seed': kwargs.get('seed', self.options.seed),
        }
        job_id = str(uuid.uuid4())
        local_job = NWQSimJob(self, job_id, self._run_job, circuits, options)
        local_job.submit()
        return local_job

    def _run_job(self, job_id, circuits, options):
        """Run circuits in q_job"""
        start = time.time()
        if isinstance(circuits, qobj_mod.QasmQobj):
            qobj = circuits
        else:
            qobj = assemble(circuits)
        end = time.time()
        #print("assemble time is:" + str(end-start))
        result_list = []

        start = time.time()
        for experiment in qobj.experiments:
            result_list.append(self.run_experiment(experiment, options))
        end = time.time()
        #print("Time is:" + str(end-start))

        result = {'backend_name': self.configuration["backend_name"],
                  'backend_version': self.configuration["backend_version"],
                  'qobj_id': qobj.qobj_id,
                  'job_id': job_id,
                  'results': result_list,
                  'status': 'COMPLETED',
                  'success': True,
                  'time_taken': (end - start)}
        return Result.from_dict(result)


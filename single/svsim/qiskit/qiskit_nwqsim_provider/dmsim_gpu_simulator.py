# -*- coding: utf-8 -*-

# Copyright 2018, IBM.
#
# This source code is licensed under the Apache License, Version 2.0 found in
# the LICENSE.txt file in the root directory of this source tree.

"""Backend for the DM-Sim GPU simulation environment."""

import itertools
import json
import logging
import operator
import os
import platform
import random
import subprocess
import time
import uuid
import numpy as np

from math import log2
from qiskit.util import local_hardware_info
from qiskit.providers import BackendV1
from qiskit.providers.models import BackendConfiguration, BackendStatus
from qiskit.result import Result
from qiskit.compiler import assemble
from qiskit import qobj as qobj_mod

from qiskit_nwqsim_provider.qobj_to_nwqsim import run_qobj_instructions
from qiskit_nwqsim_provider.nwqsim_job import NWQSimJob
import qiskit_nwqsim_provider.libdmsim as dmsim

RUN_MESSAGE = """NWQSim simulator by Ang Li, PNNL, USA
Developer: Ang Li
For more information, please visit http://github.com/pnnl/svsim"""

# this class handles the conversion to and from qiskit style data
class DMSimWrapper:
    """Converter to and from DMSim_gpu simulator"""

    def __init__(self, qobj, shell):
        self.shots = qobj.config.shots
        self.n_qubits = qobj.config.n_qubits
        self.qobj = qobj
        self.memory_mapping = self._build_memory_mapping()
        #self.sim = dmsim.Simulation(self.n_qubits, shell.n_gpus)
        self.sim = dmsim.Simulation()
        self.additional_output_data = []

    def final_statevector(self):
        return self.sim.get_sv()

    def _build_memory_mapping(self):
        qu2cl = {}
        if isinstance(self.qobj, qobj_mod.QasmQobj):
            for instruction in self.qobj.experiments[0].instructions:
                if instruction.name == 'measure':
                    qu2cl[instruction.qubits[0]] = instruction.memory[0]
            return qu2cl
        qubit_map = {}
        count = 0

        # If a list of quantum circuits use the first element
        # since we only can have a maximum of a single
        # circuit per job.
        if isinstance(self.qobj, list):
            self.qobj = self.qobj[0]

        for bit in self.qobj.qubits:
            qubit_map[bit] = count
            count += 1
        clbit_map = {}
        count = 0
        for bit in self.qobj.clbits:
            clbit_map[bit] = count
            count += 1
        for instruction in self.qobj.data:
            if instruction[0].name == 'measure':
                for index, qubit in enumerate(instruction[1]):
                    qu2cl[qubit_map[qubit]] = clbit_map[instruction[2][index]]
        return qu2cl

    def _rearrange_result(self, input):
        if isinstance(self.qobj, qobj_mod.QasmQobj):
            length = self.qobj.experiments[0].header.memory_slots
        else:
            length = self.qobj.num_clbits
        bin_output = list('0' * length)
        bin_input = list(bin(input)[2:].rjust(length, '0'))
        bin_input.reverse()
        for qu, cl in self.memory_mapping.items():
            bin_output[cl] = bin_input[qu]
        bin_output.reverse()
        return hex(int(''.join(bin_output), 2))

    def _format_counts(self, samples):
        counts = {}
        for result in samples:
            #h_result = self._rearrange_result(result)
            h_result = result
            if h_result not in counts:
                counts[h_result] = 1
            else:
                counts[h_result] += 1

        #print(counts)
        return counts

    def run_experiment(self, experiment):
        header = experiment.header
        measurement_data = {'mapping': {},
                            'clbits': header.creg_sizes,
                            'clbits_num': len(header.clbit_labels),
                            'qubits_num': len(header.qubit_labels)}
        self.sim.reset_sim()
        self.start_time = time.time()
        run_qobj_instructions(experiment, self.sim, self.shots, measurement_data['mapping'])
        self.end_time = time.time()

        #print("Compile time is:" + str(self.end_time - self.start_time))

        self.start_time = time.time()
        run_output = self.sim.measure_all(self.shots)
        self.end_time = time.time()
        #print("Run time is:" + str(self.end_time - self.start_time))

        result_dict = {'header': {'name': experiment.header.name,
                                  'memory_slots': experiment.config.memory_slots,
                                  'creg_sizes': experiment.header.creg_sizes
                                  },
                       'status': 'DONE', 'time_taken': self.end_time - self.start_time,
                       #'seed': self.seed, 
                       'shots': self.shots,
                       'data': {'counts': self._format_counts(run_output), 'statevector':self.final_statevector()},
                       #'data': {'counts': self._format_counts(run_output)},
                       'success': True
                       }
        return result_dict

    # parsing nwqsim output
    def parse_output(self, run_output, measurement_data):
        result = run_output
        qubits = measurement_data['qubits_num']
        if 'counts' in result:
            result['counts'] = self.convert_counts(result['counts'], measurement_data)
        return result

    def convert_snapshot(self, snapshot_data, translation_table):
        if 'statevector' in snapshot_data:
            snapshot_data['statevector'] = self.convert_statevector_data(
                snapshot_data['statevector'], translation_table)
        if 'probabilities' in snapshot_data:
            probs_data = snapshot_data.pop('probabilities')
            if 'probabilities' in self.additional_output_data:
                snapshot_data['probabilities'] = self.convert_probabilities(probs_data, translation_table)
        if 'probabilities_ket' in snapshot_data:
            probs_ket_data = snapshot_data.pop('probabilities_ket')
            if 'probabilities_ket' in self.additional_output_data:
                snapshot_data['probabilities_ket'] = self.convert_probabilities_ket(probs_ket_data)
        return snapshot_data

    def to_qiskit_complex(self, num_string):
        num = complex(num_string.replace('i', 'j'))  # first obtain an actual number
        return [num.real, num.imag]

    def convert_statevector_data(self, statevector, translation_table):
        return [self.to_qiskit_complex(statevector[translation_table[i]])
                for i in range(len(translation_table))]

    def convert_probabilities(self, probs_data, translation_table):
        return [probs_data[translation_table[i]] for i in range(len(translation_table))]

    def convert_probabilities_ket(self, probs_ket_data):
        return dict([(key[::-1], value) for key, value in probs_ket_data.items()])

    def convert_counts(self, counts, measurement_data):
        result = {}
        for qubits, count in counts.items():
            clbits = self.qubits_to_clbits(qubits, measurement_data)
            if clbits is not None:
                result[clbits] = result.get(clbits, 0) + count
        return result

    # converting the actual measurement results for all qubits to clbits the user expects to see
    def qubits_to_clbits(self, qubits, measurement_data):
        clbits = list('0' * measurement_data['clbits_num'])
        for (qubit, clbit) in measurement_data['mapping'].items():
            clbits[clbit] = qubits[qubit]
        s = "".join(clbits)[::-1]
        if s == '':
            return None
        return hex(int(s, 2))

logger = logging.getLogger(__name__)

VERSION_PATHS = [
    os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "VERSION.txt")),
    os.path.abspath(os.path.join(os.path.dirname(__file__), "VERSION.txt")),
]

for version_path in VERSION_PATHS:
    if os.path.exists(version_path):
        with open(version_path, "r") as version_file:
            VERSION = version_file.read().strip()

#======================================== DMSim Simulator ================================
class DMSimGpuSimulator(BackendV1):
    """Python interface to NWQSim's DMSim gpu simulator"""
    MAX_QUBITS_MEMORY = int(log2(local_hardware_info()['memory']*(1024**3)/16))
    DEFAULT_CONFIGURATION = {
        'backend_name': 'dmsim',
        'backend_version': VERSION,
        'n_qubits': min(16,MAX_QUBITS_MEMORY),
        'url': 'https://github.com/pnnl/dmsim',
        'simulator': True,
        'local': True,
        'conditional': True,
        'open_pulse': False,
        'memory': True,
        'max_shots': 1000000,
        'coupling_map': None,
        'description': 'NWQSim density matrix simulator with NVIDIA GPU cluster backend',
        'basis_gates': ['x', 'sx', 'rz', 'cx', 'id', 'reset', 'measure', 'ma'],
        'gates': [
            {
                'name': 'TODO',
                'parameters': [],
                'qasm_def': 'TODO'
            }
        ]
    }
    def __init__(self, configuration=None, provider=None):
        self.do_noise_1q = False
        self.do_noise_2q = False
        self.n_gpus = 1
        self.noise_1q = np.zeros(16)
        self.noise_2q = np.zeros(256)
        super().__init__(configuration=(
            configuration or BackendConfiguration.from_dict(self.DEFAULT_CONFIGURATION)),
                         provider=provider)
    def run(self, qc, **kwargs):
        if isinstance(qc, qobj_mod.QasmQobj):
            qobj = qc
        else:
            qobj = assemble(qc)
        qobj.config.shots = kwargs.get('shots',qobj.config.shots)
        job_id = str(uuid.uuid4())
        local_job = NWQSimJob(self, job_id, self._run_job, qobj)
        local_job.submit()
        return local_job

    def _run_job(self, job_id, qc):
        """Run circuits in q_job"""
        if isinstance(qc, qobj_mod.QasmQobj):
            qobj = qc
        else:
            qobj = assemble(qc)
        result_list = []
        self._validate(qobj)
        s = DMSimWrapper(qobj, self)
        start = time.time()
        for experiment in qobj.experiments:
            result_list.append(s.run_experiment(experiment))
        end = time.time()
        #print("Time is:" + str(end-start))
        result = {'backend_name': self._configuration.backend_name,
                  'backend_version': VERSION,
                  'qobj_id': qobj.qobj_id,
                  'job_id': job_id,
                  'results': result_list,
                  'status': 'COMPLETED',
                  'success': True,
                  'time_taken': (end - start)}
        return Result.from_dict(result)

    def set_noise_1q_sop(self, noise_1q):
        self.do_noise_1q = True
        self.noise_1q = noise_1q
    
    def set_noise_2q_sop(self, noise_2q):
        self.do_noise_2q = True
        self.noise_2q = noise_2q

    def set_n_gpus(self, n_gpus):
        self.n_gpus = n_gpus

    def _validate(self, qobj):
        # Do any necessary validation
        return

    def status(self):
        return BackendStatus(backend_name=self.name(),
                             backend_version=self.configuration().backend_version,
                             operational=True,
                             pending_jobs=0,
                             status_msg='')

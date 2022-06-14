# -*- coding: utf-8 -*-

# Copyright 2018, IBM.
#
# This source code is licensed under the Apache License, Version 2.0 found in
# the LICENSE.txt file in the root directory of this source tree.

import random
import unittest

import numpy
from scipy.stats import chi2_contingency

from qiskit import execute
from qiskit import QuantumCircuit
from qiskit import QuantumRegister
from qiskit import ClassicalRegister
from qiskit import BasicAer
from qiskit_nwqsim_provider import NWQSimProvider

try:
    global_pq_simulator = NWQSimProvider('DMSimSimulator')['dmsim']
except ImportError:
    _skip_class = True
else:
    _skip_class = False


from ._random_circuit_generator import RandomCircuitGenerator
from .common import QiskitTestCase


@unittest.skipIf(_skip_class, 'NWQSim simulator unavailable')
class TestNWQSimSimulator(QiskitTestCase):
    """
    Test NWQSim simulator.
    """

    # noinspection PyPep8Naming
    @classmethod
    def setUpClass(cls):
        super().setUpClass()

        # Set up random circuits
        n_circuits = 20
        min_depth = 1
        max_depth = 50
        min_qubits = 1
        max_qubits = 4
        random_circuits = RandomCircuitGenerator(min_qubits=min_qubits,
                                                 max_qubits=max_qubits,
                                                 min_depth=min_depth,
                                                 max_depth=max_depth,
                                                 seed=None)
        for _ in range(n_circuits):
            basis = list(random.sample(random_circuits.op_signature.keys(),
                                       random.randint(2, 7)))
            #if 'reset' in basis:
                #basis.remove('reset')
            #if 'u0' in basis:
                #basis.remove('u0')
            #if 'measure' in basis:
                #basis.remove('measure')
            random_circuits.add_circuits(1, basis=basis)
        cls.rqg = random_circuits

    def run_on_simulators(self, qc, pq_simulator, qk_simulator, shots, seed):
        job_pq = execute(qc, pq_simulator, shots=shots, seed_simulator=seed)
        job_qk = execute(qc, qk_simulator, shots=shots, seed_simulator=seed)
        counts_pq = job_pq.result().get_counts()
        counts_qk = job_qk.result().get_counts()
        states = counts_qk.keys() | counts_pq.keys()
        # contingency table
        ctable = numpy.array([[counts_pq.get(key, 0) for key in states],
                              [counts_qk.get(key, 0) for key in states]])
        result = chi2_contingency(ctable)
        return counts_pq, counts_qk, result

    def test_gate_x(self):
        shots = 100
        qr = QuantumRegister(1)
        cr = ClassicalRegister(1)
        qc = QuantumCircuit(qr, cr, name='test_gate_x')
        qc.x(qr[0])
        qc.measure(qr, cr)
        job = execute(qc, global_pq_simulator, shots=shots)
        result_pq = job.result(timeout=30)
        self.assertEqual(result_pq.get_counts(),
                         {'1': shots})

    def test_entangle(self):
        shots = 100
        N = 5
        qr = QuantumRegister(N)
        cr = ClassicalRegister(N)
        qc = QuantumCircuit(qr, cr, name='test_entangle')

        qc.h(qr[0])
        for i in range(1, N):
            qc.cx(qr[0], qr[i])
        qc.measure(qr, cr)
        timeout = 30
        job = execute(qc, global_pq_simulator, shots=shots)
        result = job.result(timeout=timeout)
        counts = result.get_counts()
        self.log.info(counts)
        for key, _ in counts.items():
            with self.subTest(key=key):
                self.assertTrue(key in ['0' * N, '1' * N])

    def test_output_style(self):
        qk_simulator = BasicAer.get_backend('qasm_simulator', )

        qr = QuantumRegister(2)
        cr = ClassicalRegister(2)
        qc = QuantumCircuit(qr, cr, name='test_output_order')
        qc.h(qr[0])
        qc.measure(qr[0], cr[0])
        qc.measure(qr[1], cr[1])
        shots = 100

        counts_pq, counts_qk, result = self.run_on_simulators(qc, global_pq_simulator,
                                                              qk_simulator, shots=shots, seed=1)
        self.assertGreater(result[1], 0.01)

        cr1 = ClassicalRegister(1)
        cr2 = ClassicalRegister(1)
        qc = QuantumCircuit(qr, cr1, cr2, name='test_output_separation')
        qc.h(qr[0])
        qc.measure(qr[0], cr1[0])
        qc.measure(qr[1], cr2[0])

        counts_pq, counts_qk, result = self.run_on_simulators(qc, global_pq_simulator,
                                                              qk_simulator, shots=shots, seed=1)
        self.log.info('chi2_contingency: %s', str(result))
        self.assertGreater(result[1], 0.01)

    def test_random_circuits(self):
        qk_simulator = BasicAer.get_backend('qasm_simulator', )
        for circuit in self.rqg.get_circuits(format_='QuantumCircuit'):
            self.log.info(circuit.qasm())
            shots = 100
            job_pq = execute(circuit, global_pq_simulator, shots=shots, seed_simulator=1)
            job_qk = execute(circuit, qk_simulator, shots=shots, seed_simulator=1)
            result_pq = job_pq.result()
            result_qk = job_qk.result()
            counts_pq = result_pq.get_counts()
            counts_qk = result_qk.get_counts()
            self.log.info('local_nwqsim_simulator: %s', str(counts_pq))
            self.log.info('local_qiskit_simulator: %s', str(counts_qk))
            states = counts_qk.keys() | counts_pq.keys()
            # contingency table
            ctable = numpy.array([[counts_pq.get(key, 0) for key in states],
                                  [counts_qk.get(key, 0) for key in states]])
            result = chi2_contingency(ctable)
            self.log.info('chi2_contingency: %s', str(result))
            with self.subTest(circuit=circuit):
                self.assertGreater(result[1], 0.01)


if __name__ == '__main__':
    unittest.main(verbosity=2)

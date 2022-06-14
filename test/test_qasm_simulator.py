# -*- coding: utf-8 -*-

# Copyright 2017, IBM.
#
# This source code is licensed under the Apache License, Version 2.0 found in
# the LICENSE.txt file in the root directory of this source tree.

# pylint: disable=missing-docstring,redefined-builtin

import unittest
import os
from qiskit import QuantumCircuit, transpile
from common import QiskitTestCase
from qiskit_nwqsim_provider import NWQSimProvider
from qiskit import execute


class TestNWQSimSimulatorBasic(QiskitTestCase):
    """Runs the Basic qasm_simulator tests from Terra on NWQSim."""

    def setUp(self):
        self.seed = 88
        self.backend = NWQSimProvider('DMSimSimulator').backends['dmsim']
        qasm_filename = os.path.join(os.path.dirname(__file__), 'qasms', 'example.qasm')
        compiled_circuit = QuantumCircuit.from_qasm_file(qasm_filename)
        transpiled_circuit = transpile(compiled_circuit, self.backend)
        compiled_circuit.name = 'test'
        self.circuit = transpiled_circuit

    def test_nwqsim_simulator_single_shot(self):
        """Test single shot run."""
        result = execute(self.circuit, self.backend, seed_transpiler=34342, shots=1).result()
        self.assertEqual(result.success, True)

    def test_nwqsim_simulator(self):
        """Test data counts output for single circuit run against reference."""
        shots = 1024
        result = execute(self.circuit, self.backend, seed_transpiler=34342, shots=shots).result()
        threshold = 0.04 * shots
        counts = result.get_counts(self.circuit)
        target = {'100 100': shots / 8, '011 011': shots / 8,
                  '101 101': shots / 8, '111 111': shots / 8,
                  '000 000': shots / 8, '010 010': shots / 8,
                  '110 110': shots / 8, '001 001': shots / 8}
        self.assertDictAlmostEqual(counts, target, threshold)


if __name__ == '__main__':
    unittest.main()

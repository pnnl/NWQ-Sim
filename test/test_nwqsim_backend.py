# -*- coding: utf-8 -*-

# Copyright 2019, IBM.
#
# This source code is licensed under the Apache License, Version 2.0 found in
# the LICENSE.txt file in the root directory of this source tree.

"""Test NWQSim backend."""

from qiskit.circuit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit import execute
from qiskit_nwqsim_provider import NWQSimProvider

from common import QiskitTestCase

class NWQSimBackendTestCase(QiskitTestCase):
    """Tests for the NWQSim backend."""

    def setUp(self):
        super().setUp()
        self.backend = NWQSimProvider('DMSimSimulator')['dmsim']

    def test_configuration(self):
        """Test backend.configuration()."""
        configuration = self.backend.configuration()
        return configuration

    def test_properties(self):
        """Test backend.properties()."""
        properties = self.backend.properties()
        self.assertEqual(properties, None)

    def test_status(self):
        """Test backend.status()."""
        status = self.backend.status()
        return status

    def test_run_circuit(self):
        """Test running a single circuit."""
        result = execute(bell(), self.backend, seed_transpiler=34342).result()
        self.assertEqual(result.success, True)
        return result


def bell():
    """Return a Bell circuit."""
    qr = QuantumRegister(2, name='qr')
    cr = ClassicalRegister(2, name='qc')
    qc = QuantumCircuit(qr, cr, name='bell')
    qc.h(qr[0])
    qc.cx(qr[0], qr[1])
    qc.measure(qr, cr)
    return qc

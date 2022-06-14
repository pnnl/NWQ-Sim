# -*- coding: utf-8 -*-

# Copyright 2018, IBM.
#
# This source code is licensed under the Apache License, Version 2.0 found in
# the LICENSE.txt file in the root directory of this source tree.

# pylint: disable=unused-import
# pylint: disable=redefined-builtin

"""Test Qiskit's QuantumCircuit class for multiple registers."""

from qiskit import QuantumRegister, ClassicalRegister, QuantumCircuit
from qiskit import execute
#from qiskit.quantum_info import state_fidelity, basis_state
from qiskit.test import QiskitTestCase
from qiskit_nwqsim_provider import NWQSimProvider

class TestCircuitMultiRegs(QiskitTestCase):
    """QuantumCircuit Qasm tests."""

    def test_circuit_multi(self):
        """Test circuit multi regs declared at start.
        """
        qreg0 = QuantumRegister(2, 'q0')
        creg0 = ClassicalRegister(2, 'c0')
        qreg1 = QuantumRegister(2, 'q1')
        creg1 = ClassicalRegister(2, 'c1')
        circ = QuantumCircuit(qreg0, qreg1)
        circ.x(qreg0[1])
        circ.x(qreg1[0])

        meas = QuantumCircuit(qreg0, qreg1, creg0, creg1)
        meas.measure(qreg0, creg0)
        meas.measure(qreg1, creg1)

        qc = circ + meas

        backend_sim = NWQSimProvider('DMSimSimulator')['dmsim']

        result = execute(qc, backend_sim, seed_transpiler=34342).result()
        counts = result.get_counts(qc)

        target = {'01 10': 1024}
        self.assertEqual(counts, target)

        #result = execute(circ, backend_sim, seed_transpiler=3438).result()
        #state = result.get_statevector(circ)
        #self.assertAlmostEqual(state_fidelity(basis_state('0110', 4), state), 1.0, places=7)

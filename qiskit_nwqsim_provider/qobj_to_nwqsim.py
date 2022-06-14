# -*- coding: utf-8 -*-

# This code is part of Qiskit.
#
# (C) Copyright IBM 2019.
#
# This code is licensed under the Apache License, Version 2.0. You may
# obtain a copy of this license in the LICENSE.txt file in the root directory
# of this source tree or at http://www.apache.org/licenses/LICENSE-2.0.
#
# Any modifications or derivative works of this code must retain this
# copyright notice, and modified files need to carry a notice indicating
# that they have been altered from the originals.

from numpy import pi

def run_qobj_instructions(experiment, sim, shots, mapping):
    ops = experiment.instructions
    #print(ops)
    for op in ops:
        if op.name == 'x':
            sim.X(op.qubits[0])
        elif op.name == 'y':
            sim.Y(op.qubits[0])
        elif op.name == 'z':
            sim.Z(op.qubits[0])
        elif op.name == 'h':
            sim.H(op.qubits[0])
        elif op.name == 's':
            sim.S(op.qubits[0])
        elif op.name == 'sdg':
            sim.SDG(op.qubits[0])
        elif op.name == 't':
            sim.T(op.qubits[0])
        elif op.name == 'tdg':
            sim.TDG(op.qubits[0])
        elif op.name == 'rx':
            sim.RX(op.params[0], op.qubits[0])
        elif op.name == 'ry':
            sim.RY(op.params[0], op.qubits[0])
        elif op.name == 'rz':
            sim.RZ(op.params[0], op.qubits[0])
        elif op.name == 'sx':
            sim.SX(op.qubits[0])
        elif op.name == 'p':
            sim.P(op.params[0], op.qubits[0])
        elif op.name == 'u':
            sim.U(op.params[0], op.params[1], op.params[2], op.qubits[0])
        elif op.name == 'cx':
            sim.CX(op.qubits[0], op.qubits[1])
        elif op.name == 'cy':
            sim.CY(op.qubits[0], op.qubits[1])
        elif op.name == 'cz':
            sim.CZ(op.qubits[0], op.qubits[1])
        elif op.name == 'ch':
            sim.CH(op.qubits[0], op.qubits[1])
        elif op.name == 'crx':
            sim.CRX(op.params[0], op.qubits[0], op.qubits[1])
        elif op.name == 'cry':
            sim.CRY(op.params[0], op.qubits[0], op.qubits[1])
        elif op.name == 'crz':
            sim.CRZ(op.params[0], op.qubits[0], op.qubits[1])
        elif op.name == 'csx':
            sim.CSX(op.qubits[0], op.qubits[1])
        elif op.name == 'cp':
            sim.CP(op.params[0], op.qubits[0], op.qubits[1])
        elif op.name == 'cu':
            sim.CU(op.params[0], op.params[1], op.params[2], op.params[3], op.qubits[0], op.qubits[1])
        elif op.name == 'swap':
            sim.SWAP(op.qubits[0], op.qubits[1])
        #elif op.name == 'measure':
            #sim.measureZ(op.qubits[0])
            #print('\n=============')
            #print(op.memory, op.register)
            #print('=============')
            #mapping[op.qubits[0]] = op.memory[0]



        #elif op.name == 'cx':
            #sim.CX(op.conditional, op.qubits[0])
        #elif op.name == 'measure':
            #sim.MA(shots)
            #mapping[op.qubits[0]] = op.memory[0]

    #'basis_gates': ['x', 'y', 'z', 'h', 's', 'sdg', 't', 'tdg', 'r', 'rx', 'ry', 'rz', 
        #'sx', 'p', 'u', 'cx', 'cy', 'cz', 'ch', 'cs', 'csdg', 'ct', 'ctdg', 'cr', 
        #'crx', 'cry', 'crz', 'csx', 'cp', 'cu', 'id', 'swap', 'm', 'ma', 'reset']


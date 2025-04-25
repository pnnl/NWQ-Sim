from qiskit import QuantumCircuit


def stim_to_qiskit(stim_circuit):
    n_qubits = stim_circuit.num_qubits
    qc = QuantumCircuit(n_qubits)

    def process_instruction(op):
        name = op.name
        targets = [t.value for t in op.targets_copy()]
        if name == 'H':
            qc.h(targets[0])
        elif name in ['CX', 'CNOT']:
            qc.cx(targets[0], targets[1])
        elif name == 'S':
            qc.s(targets[0])
        elif name == 'X':
            qc.x(targets[0])
        elif name == 'Z':
            qc.z(targets[0])
        elif name == 'Y':
            qc.y(targets[0])
        # Add other gate translations if needed

    for op in stim_circuit:
        if isinstance(op, stim.CircuitInstruction):
            process_instruction(op)
        elif isinstance(op, stim.CircuitRepeatBlock):
            for _ in range(op.repeat_count):
                for subop in op.body:
                    process_instruction(subop)

    return qc

# Example usage:
import stim
stim_circ = stim.Circuit.from_file("NWQ-Sim/stabilizer/stim_files/logical_cnot.stim")
qiskit_circ = stim_to_qiskit(stim_circ)
print(qiskit_circ.qasm())

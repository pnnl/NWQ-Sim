import stim

# Load circuit from file
with open("NWQ-Sim/stabilizer/stim_files/logical_cnot.stim", "r") as f:
    circuit_text = f.read()

circuit = stim.Circuit(circuit_text)

# Strip unsupported instructions for QASM export
clean_circuit = circuit.without_noise()

# Convert to QASM
qasm_str = clean_circuit.to_qasm(open_qasm_version=2, skip_dets_and_obs=True)
# Print QASM output
print(qasm_str)

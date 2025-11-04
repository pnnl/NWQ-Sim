from __future__ import annotations

from collections import Counter

from nwqsim import Circuit, create_state


def main() -> None:
    num_qubits = 16
    backend = "cpu"
    method = "stab"

    circuit = Circuit(num_qubits)
    
    for i in range(num_qubits):
        circuit.h(i-1)
        circuit.m(i-1)

    state = create_state(backend, num_qubits, method)
    sim_time_ms = state.simulate(circuit)
    # samples = state.measure_all(1024)
    samples = state.measurement_results()

    print(f"Simulated {len(samples)} shots in {sim_time_ms:.3f} ms")
    counts = Counter(samples)
    for bitstring, count in sorted(counts.items()):
        print(f"{bitstring:0{num_qubits}b}: {count}")


if __name__ == "__main__":
    main()

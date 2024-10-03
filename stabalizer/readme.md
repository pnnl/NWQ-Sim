This branch implements a tableau object to efficiently simulate the evolution of Pauli & Clifford gates. The Clifford set (S, H, CNOT) combined with T gates form a universal set for quantum computation, which incentivizes the creation of such a simulator.

*** How to Use ***
Using NWQSim normally involves the creation of a circuit object containing the qubit gates, followed by a state object to track their simulation and evolution. This tableau implementation is an augmentation of the state object. In addition to the usual state->sim function, there is a state->clifford_sim which creates a tableau and evolves it according to the same circuit object that state->sim uses.
Using the clifford sim should basically follow the typical steps in the user manual, aside from substituing out the state->sim call. Currently the tableau supports X, Y, Z, S, H, CNOT and M gates, but features and testing are currently in development.

Integrated Option:

1. Create an NWQSim 'Circuit' object and add gates.
2. Create an NWQSim 'State' object as usual for simulation.
3. Instead of 'state'->sim('circuit') use 'state'->clifford_sim('circuit') to simulate using the stabilizer simulator backend.
4. Call 'state'->get_results() as usual to see the results returned by any measurement gates in the circuit.

    Example:
        auto circuit = std::make_shared<NWQSim::Circuit>(n_qubits);
        // Add some gates to the circuit
        for(int i = 0; i < n_qubits; i++)
        {
            circuit -> H(i);
        }
        for(int i = 0; i < n_qubits-1; i++)
        {
            circuit -> CX(i,i+1);
        }
        for(int i = 0; i < n_qubits; i++)
        {
            circuit -> S(i);
        }
        for(int i = 0; i < n_qubits; i++)
        {
            circuit -> M(i);
        }

        std::string backend = "CPU";
        std::string sim_method = "sv";
        auto state = BackendManager::create_state(backend, n_qubits, sim_method);
        double state_time;
        state2->clifford_sim(circuit, state_time);
        long long int *results2 = state2->get_results();
        std::cout << "Stabalizer sim: " << results2[0] << std::endl;
        std::cout << state2_time << std::endl;

Manual Option: In this mode, you can add gates to the tableau, or construct it with a Pauli string of stabilizers.

1. Create an NWQSim::tableau object with an initial NWQSim::Circuit (this can be empty, gates can be appended later).
2. Use the functions in tableau.hpp to directly interact with the tableau.
3. Add gates with 'tableau'.add_gates(circuit). Each time gates are added, call 'tableau'.simulate() to ensure the tableau is up to date.

This branch implements a tableau object to efficiently (polynomial time as opposed to exponential) simulate the evolution of Pauli & Clifford gates. The Clifford set (S, H, CNOT) combined with T gates form a universal set for quantum computation, which incentivizes the creation of such a simulator.

*** How to Use ***
Using NWQSim normally involves the creation of a circuit object containing the qubit gates, followed by a state object to track their simulation and evolution. This tableau implementation is an augmentation of the state object. In addition to the usual state->sim function, there is a state->clifford_sim which creates a tableau and evolves it according to the same circuit object that state->sim uses.
Using the clifford sim should basically follow the typical steps in the user manual, aside from substituing out the state->sim call. Currently the tableau supports X, Y, Z, S, H, CNOT and M gates, but features and testing are currently in development.

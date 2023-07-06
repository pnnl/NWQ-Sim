import json


def tfFromUnit(unit):  # time factor from unit, convert time unis to seconds
    if unit == "ns":
        return 1e-9
    elif unit == "us":
        return 1e-6
    elif unit == "ms":
        return 1e-3
    elif unit == "s":
        return 1
    elif unit == "ks":
        return 1000


def procBackendProp(backend, path_to_file=None):
    """The function that output a file for backend noise properties in json.
    The file is named "{backend name}Config.json"
    @param backend: an Qiskit IBMQBackend object, see Validation/qiskit_backend_noise_model_validation.ipynb for example
    @param path_to_file: a string that stores the path to file, e.g., "../Data/"
    @return info_dict: the same dictionary the stores in the output "{backend name}Config.json"
    """
    backend_prop = backend.properties()
    backend_config = backend.configuration()
    backend_name = backend_prop.backend_name  # str name
    backend_version = backend_prop.backend_version  # str, actually not needed
    num_qubits = len(backend_prop.qubits)  # int
    basis_gates = backend_config.basis_gates  # list
    info_dict = {}
    info_dict["name"] = backend_name
    info_dict["version"] = backend_version
    info_dict["num_qubits"] = num_qubits
    info_dict["basis_gates"] = basis_gates  # Qubit properties
    info_dict["T1"] = {}  # unit: sec, from backend.properties().t1()
    info_dict["T2"] = {}  # unit: sec, from backend.properties().t2()
    info_dict["freq"] = {}  # unit: Hz
    info_dict["readout_length"] = {}
    info_dict["prob_meas0_prep1"] = {}
    info_dict["prob_meas1_prep0"] = {}

    # Record
    for qubit_index in range(num_qubits):
        qubit_prop_dict = backend_prop.qubit_property(qubit_index)
        info_dict["T1"][str(qubit_index)] = float(
            qubit_prop_dict["T1"][0]
        )  # the unit by default is sec
        info_dict["T2"][str(qubit_index)] = float(
            qubit_prop_dict["T2"][0]
        )  # the unit by default is sec
        info_dict["freq"][str(qubit_index)] = float(qubit_prop_dict["frequency"][0])
        info_dict["readout_length"][str(qubit_index)] = float(
            qubit_prop_dict["readout_length"][0]
        )  # the unit by default is sec
        info_dict["prob_meas0_prep1"][str(qubit_index)] = float(
            qubit_prop_dict["prob_meas0_prep1"][0]
        )
        info_dict["prob_meas1_prep0"][str(qubit_index)] = float(
            qubit_prop_dict["prob_meas1_prep0"][0]
        )

        # Gate Properties, for both 1-qubit and 2-qubit
        info_dict["gate_lens"] = {}
        info_dict["gate_errs"] = {}
        info_dict["cx_coupling"] = []

        for gate_obj in backend_prop.gates:
            gate_name = gate_obj.gate  # e.g., "id", "cx"
            gate_qubits = (
                gate_obj.qubits
            )  # e.g, [1], [4,5]         # create a key includes name and operating qubits

            if len(gate_qubits) == 1:
                gate_key = gate_name + str(gate_qubits[0])
            elif len(gate_qubits) == 2:
                coupling = "{:d}_{:d}".format(gate_qubits[0], gate_qubits[1])
                gate_key = gate_name + coupling
                info_dict["cx_coupling"].append(coupling)

            # Obtain gate error rate and gate length

            if (
                gate_obj.parameters[0].name == "gate_error"
            ):  # one is gate_error, the other is gate_length
                info_dict["gate_errs"][gate_key] = gate_obj.parameters[0].value
                try:
                    info_dict["gate_lens"][gate_key] = gate_obj.parameters[
                        1
                    ].value * tfFromUnit(gate_obj.parameters[1].unit)
                except:
                    info_dict["gate_lens"][gate_key] = 0
            else:
                info_dict["gate_lens"][gate_key] = gate_obj.parameters[
                    0
                ].value * tfFromUnit(gate_obj.parameters[0].unit)

                try:
                    info_dict["gate_errs"][gate_key] = gate_obj.parameters[1].value
                except:
                    info_dict["gate_errs"][gate_key] = 0

    if path_to_file != None:
        if path_to_file != "" and path_to_file[-1] != "/":
            path_to_file = path_to_file + "/"

        with open(path_to_file + info_dict["name"] + "Config.json", "w") as f:
            json.dump(info_dict, f)
    return info_dict

import sys
import json
import random
import collections
from collections import OrderedDict

# how to use
# python generate_dummy_topology.py num_qubits machine
# machine: ibmq, rigetti, honeywell, quantinuum
# e.g. python generate_dummy_topology.py 7 ibmq

num_qubits = sys.argv[1]
machine = sys.argv[2]

version = "0.1"

config = OrderedDict()

config['name'] = machine+"_dummy"
config['version']  = version
config['num_qubits'] = int(num_qubits)

two_qubits_gate = ""

if machine == "ibmq":
    config['basis_gates'] = ["id", "rz", "sx", "x", "cx", "reset"]
    two_qubits_gate = "cx" 

elif machine == "rigetti":
    config['basis_gates'] = ["rx","rz","xy"]
    two_qubits_gate = "xy"

elif machine == "honeywell":
    config['basis_gates'] = ["gpi","gpi2","gz","ms"]
    two_qubits_gate = "ms"

elif machine == "quantinuum":
    config['basis_gates'] = ["rx","rz","zz"]
    two_qubits_gate = "zz"

else:
    print("Architecture not supported!")
    exit()

config["T1"] = OrderedDict()
config["T2"] = OrderedDict()
config["freq"] = OrderedDict()
config["readout_length"] = OrderedDict()
config["prob_meas0_prep1"] = OrderedDict()
config["prob_meas1_prep0"] = OrderedDict()
config["gate_lens"] = OrderedDict()
config["gate_errs"] = OrderedDict()
for i in range(int(num_qubits)):
    config["T1"][str(i)] = random.uniform(0, 0.00001)
    config["T2"][str(i)] = random.uniform(0, 0.00001)
    config["freq"][str(i)] = random.uniform(5000000000, 5300000000)
    config["readout_length"][str(i)] = random.uniform(0,0.0000006)
    config["prob_meas0_prep1"][str(i)] = random.uniform(0,0.1)
    config["prob_meas1_prep0"][str(i)] = random.uniform(0,0.1)
    for gate in config['basis_gates']:
        if machine == "ibmq":
            if gate != "rz" and gate != two_qubits_gate:
                config["gate_lens"][gate+str(i)] = random.uniform(0,0.00000001)
                config["gate_errs"][gate+str(i)] = random.uniform(0,0.00000001)
            else:
                if gate == "rz":
                    config["gate_lens"][gate+str(i)] = 0
                    config["gate_errs"][gate+str(i)] = 0
        else:
            if gate != two_qubits_gate:
                config["gate_lens"][gate+str(i)] = random.uniform(0,0.00000001)
                config["gate_errs"][gate+str(i)] = random.uniform(0,0.00000001)
    config["gate_lens"]["reset"+str(i)] = random.uniform(0,0.00000001) 
    config["gate_errs"]["reset"+str(i)] = random.uniform(0,0.00000001) 

config[two_qubits_gate+"_coupling"] = []
for i in range(int(num_qubits)):
    for j in range(int(num_qubits)):
        if i != j:
            config["gate_lens"][two_qubits_gate+str(i)+"_"+str(j)] = random.uniform(0,0.00000001)
            config["gate_errs"][two_qubits_gate+str(i)+"_"+str(j)] = random.uniform(0,0.00000001)
            config[two_qubits_gate+"_coupling"].append(str(i)+"_"+str(j))

with open('dummy'+"_"+machine+str(num_qubits)+'.json', 'w') as outfile:
    json.dump(config, outfile)



import itertools
import stim
import numpy as np
import re
from collections import defaultdict
from typing import List, Tuple, Union, Dict, Optional

NoiseSpecs = List[Tuple[int, Tuple[str, ...], float]]

NOISE_GATES = {
    "X_ERROR", "Y_ERROR", "Z_ERROR",
    "DEPOLARIZE1", "DEPOLARIZE2",
    "PAULI_CHANNEL_1", "PAULI_CHANNEL_2",
    "AMPLITUDE_DAMP",  # allow in settings; we won't emit into Stim text
    "E",  # Correlated Pauli error
}

SINGLE_QUBIT_GATES = {
    "H", "X", "Y", "Z",
    "S", "SDAG", "T", "TDAG",
    "RX", "RY", "RZ",
}

# Common two‑qubit Clifford gates (extend as needed)
TWO_QUBIT_GATES = {
    "CX", "CNOT", "ZCX",
    "CZ", "ZZ", "SWAP",
}

PAULI_DOUBLE = [
    ("I", "I"),
    ("I", "X"), ("I", "Y"), ("I", "Z"),
    ("X", "I"), ("X", "X"), ("X", "Y"), ("X", "Z"),
    ("Y", "I"), ("Y", "X"), ("Y", "Y"), ("Y", "Z"),
    ("Z", "I"), ("Z", "X"), ("Z", "Y"), ("Z", "Z"),
]

# ...existing code...
# ...existing code...
def stim_to_qasm_with_depolarize_noise(circuit: stim.Circuit) -> str:
    """
    Converts a stim circuit with DEPOLARIZE1/DEPOLARIZE2/AMPLITUDE_DAMP instructions
    to a QASM 2.0 string, defining custom gates for these noise channels.

    Note: REPEAT blocks are unsupported in QASM; we unroll them by using
    circuit.flattened() so all operations are emitted explicitly.
    """
    qasm_body = ""
    # Unroll REPEAT blocks by flattening the circuit
    for instruction in circuit.flattened():
        name = instruction.name
        targets = instruction.targets_copy()
        args = instruction.gate_args_copy()

        qasm_inst = ""
        if name in ["X_ERROR", "Y_ERROR", "Z_ERROR"]:
            gate_name = name.lower()
            params = f"({args[0]})" if args else ""
            for target in targets:
                qasm_body += f"{gate_name}{params} q[{target.value}];\n"
        elif name == "DEPOLARIZE1":
            gate_name = "dep1"
            params = f"({args[0]})"
            for target in targets:
                qasm_body += f"{gate_name}{params} q[{target.value}];\n"
        elif name == "DEPOLARIZE2":
            gate_name = "dep2"
            params = f"({args[0]})"
            # DEPOLARIZE2 targets pairs of qubits
            for i in range(0, len(targets), 2):
                q1 = targets[i].value
                q2 = targets[i+1].value
                qasm_body += f"{gate_name}{params} q[{q1}],q[{q2}];\n"
        elif name == "AMPLITUDE_DAMP":  # map Stim op (if present) to QASM damp
            gate_name = "damp"
            a = list(args)
            while len(a) < 2:
                a.append(0.0)
            params = f"({a[0]}, {a[1]})"
            for target in targets:
                qasm_body += f"{gate_name}{params} q[{target.value}];\n"
        elif name == "PAULI_CHANNEL_1":
            gate_name = "chan1"
            a = list(args)
            while len(a) < 3:
                a.append(0.0)
            params = f"({a[0]},{a[1]},{a[2]})"
            for target in targets:
                qasm_body += f"{gate_name}{params} q[{target.value}];\n"
        elif name == "PAULI_CHANNEL_2":
            gate_name = "chan2"
            a = list(args)
            params = "(" + ",".join(str(x) for x in a) + ")"
            if len(targets) % 2 != 0:
                raise ValueError("PAULI_CHANNEL_2 requires an even number of targets (pairs)")
            for i in range(0, len(targets), 2):
                q1 = targets[i].value
                q2 = targets[i+1].value
                qasm_body += f"{gate_name}{params} q[{q1}],q[{q2}];\n"
        else:
            # Fallback: rely on Stim's per-instruction QASM conversion
            temp_circ = stim.Circuit()
            temp_circ.append(instruction)
            try:
                qasm_lines = temp_circ.to_qasm(open_qasm_version=2).splitlines()
                for line in qasm_lines:
                    if 'q[' in line and not line.startswith(('qreg', 'creg', 'OPENQASM', 'include')):
                         qasm_body += line + "\n"
            except ValueError:
                if name not in ["DETECTOR", "OBSERVABLE_INCLUDE", "SHIFT_COORDS", "TICK"]:
                     qasm_body += f"# Unsupported for QASM: {name}\n"

    # --- Assemble the full QASM string ---
    num_qubits = circuit.num_qubits
    qasm_header = "OPENQASM 2.0;\n"
    qasm_header += 'include "qelib1.inc";\n'
    qasm_header += f"qreg q[{num_qubits}];\n"
    # creg is not always needed and can be large.
    # qasm_header += f"creg c[{num_qubits}];\n"

    return qasm_header + qasm_body

def remove_noise_statements(circuit: stim.Circuit) -> stim.Circuit:
    """
    Remove all noise instructions from the input stim.Circuit.
    Returns a new circuit containing only the non-noise instructions.
    """
    noiseless_circ = stim.Circuit()
    circ = circuit.flattened()
    for inst in circ:
        # Remove instructions that are recognized as noise gates
        if inst.name in NOISE_GATES:
            continue
        noiseless_circ.append(inst)
    return noiseless_circ

def decompose_operations(circuit: stim.Circuit) -> stim.Circuit:
    """
    Decompose non-Z-basis measurement/reset operations into Z-basis versions with Clifford rotations.
    """
    decomposed_circ = stim.Circuit()
    circ = circuit.flattened()
    for inst in circ:
        name = inst.name.upper()
        targets = [t.value for t in inst.targets_copy()] #type: ignore
        # Decompose measurements
        if name in ("M", "MZ"):
            decomposed_circ.append("M", targets)    #type:ignore
        elif name == "MX":
            decomposed_circ.append("H", targets) #type:ignore
            decomposed_circ.append("M", targets) #type:ignore
        elif name == "MY":
            decomposed_circ.append("S_DAG", targets) #type:ignore
            decomposed_circ.append("H", targets) #type:ignore
            decomposed_circ.append("M", targets) #type:ignore
        elif name == "MR" or name == "MRZ":
            decomposed_circ.append("M", targets) #type:ignore
            decomposed_circ.append("R", targets) #type:ignore
        elif name == "MRX":
            decomposed_circ.append("H", targets) #type:ignore
            decomposed_circ.append("M", targets) #type:ignore
            decomposed_circ.append("R", targets) #type:ignore
            decomposed_circ.append("H", targets) #type:ignore
        elif name == "MRY":
            decomposed_circ.append("S_DAG", targets) #type:ignore
            decomposed_circ.append("H", targets) #type:ignore
            decomposed_circ.append("M", targets) #type:ignore
            decomposed_circ.append("R", targets) #type:ignore
            decomposed_circ.append("S", targets) #type:ignore
            decomposed_circ.append("H", targets) #type:ignore
        elif name in ("R", "RZ"):
            decomposed_circ.append("R", targets) #type:ignore
        elif name == "RX":
            decomposed_circ.append("R", targets) #type:ignore
            decomposed_circ.append("H", targets) #type:ignore
        elif name == "RY":
            decomposed_circ.append("S_DAG", targets) #type:ignore
            decomposed_circ.append("H", targets) #type:ignore
            decomposed_circ.append("R", targets) #type:ignore
            decomposed_circ.append("H", targets) #type:ignore
            decomposed_circ.append("S", targets) #type:ignore
        else:
            decomposed_circ.append(inst)
    return decomposed_circ


class ErrorModel:
    '''use this for tracking the error propagation within a stim circuit'''

    def __init__(self, stim_circ: stim.Circuit):
        self.noiseless_circ = remove_noise_statements(stim_circ)
        self.noiseless_circ_list = str(self.noiseless_circ).splitlines()
        ## this is to set weather the special type of the error are considered or not
        self._error_config = {'Two_qubit': True,
                             'Single_qubit': True,
                             'Measurement': True,
                             'Reset': True,
                             'Identity': True}

        ## this is to set the error type and the error rate
        self._error_settings: dict[str, Optional[str]] = {'Two_qubit': 'DEPOLARIZE2(0.01)',
                             'Single_qubit': 'DEPOLARIZE1(0.01)',
                             'Measurement': 'X_ERROR(0.01)',
                             'Reset': 'X_ERROR(0.01)',
                             'Identity': 'DEPOLARIZE1(0.01)'}

        self.conv = self.noiseless_circ.compile_m2d_converter()
        self.noise_specs = None
        self.total_measurements = 0
        self.total_detectors = 0
        self.refresh_circuit_statistics()
        return
    
    # ------------------------------------------------------------------
    # Helper: recompute cached measurement / detector counts
    # ------------------------------------------------------------------
    def refresh_circuit_statistics(self) -> None:
        """
        Re‑scan ``self.noiseless_circ`` and update
    
            • ``self.total_measurements`` – total number of measurement results
              (all M‑family instructions) in the *noise‑free* circuit;

            • ``self.total_detectors``    – total number of ``DETECTOR`` statements
              in the same circuit.

        Call this helper **after** you have rebuilt or modified
        ``self.noiseless_circ`` so that the cached counters remain consistent.
        """
        n_meas = 0
        n_det  = 0
        for inst in self.noiseless_circ:
            name = inst.name.upper()
            if name.startswith("M"):               # M, MZ, MX, MY, etc.
                n_meas += len(inst.targets_copy())  # one bit per target # type:ignore
            elif name == "DETECTOR":
                n_det += 1

        self.total_measurements = n_meas
        self.total_detectors    = n_det
        return
    
    
    def setting_error(self, error_place: str, make_error: bool, error_statements: str | None = None):
        """
        Configure whether and how to inject specific types of noise into the circuit.

        This function allows the user to enable or disable noise injection at different 
        locations in the circuit (e.g., single-qubit gates, two-qubit gates, measurement, 
        or reset operations), and to specify the corresponding noise model (e.g., DEPOLARIZE1, 
        X_ERROR, PAULI_CHANNEL_2, etc.).

        Parameters
        ----------
        error_place : str
            The category of circuit component to which the noise setting applies. 
            Must be one of: 'Two_qubit', 'Single_qubit', 'Measurement', 'Reset'.

        make_error : bool
            Flag indicating whether to include noise of the specified type in the circuit.
            If True, noise of the given type will be inserted (based on `error_statements`); 
            if False, no noise will be inserted for this component.

        error_statements : str
            A Stim noise instruction string specifying the noise model and its parameters.
            Examples: "DEPOLARIZE2(0.01)", "X_ERROR(0.02)", "PAULI_CHANNEL_1(0.01,0.01,0.01)".
            If set to None, no error instruction will be applied (even if `make_error` is True).

        Raises
        ------
        ValueError
            If `error_place` is not a valid key, or if `error_statements` does not match the
            expected format for the given `error_place`.

        Notes
        -----
        - Two-qubit errors must be specified with 'DEPOLARIZE2(...)' or 'PAULI_CHANNEL_2(...)'.
        - Single-qubit, measurement, and reset errors can use any single-qubit noise instruction
          supported by Stim: X/Y/Z_ERROR, DEPOLARIZE1, or PAULI_CHANNEL_1.
        - The function updates the internal configuration used later by `generate_noisy_circuit`.

        Examples
        --------
        >>> model.setting_error('two_qubit', True, 'DEPOLARIZE2(0.01)')
        >>> model.setting_error('Measurement', True, 'X_ERROR(0.005)')
        >>> model.setting_error('Reset', False, None)
        """
        allowed_keys = self._error_config.keys()
        if error_place not in allowed_keys:
            raise ValueError(f"error_place must be one of {list(allowed_keys)}; got '{error_place}'.")
        self._error_config[error_place] = make_error
        if not make_error:
            ## if the error is not active. the error statement does not needed
            return

        if error_statements is not None:
            s = error_statements.strip().upper()
            if error_place == "Two_qubit":
                if not (s.startswith("DEPOLARIZE2") or s.startswith("PAULI_CHANNEL_2")):
                    raise ValueError(
                        "For two_qubit, error_statements must be 'DEPOLARIZE2(...)' or 'PAULI_CHANNEL_2(...)'."
                    )
            else:
                # allow single-qubit Pauli or depolarize1 or pauli_channel_1 or amplitude_damp
                if not (s.startswith("X_ERROR") or s.startswith("Y_ERROR") or s.startswith("Z_ERROR") or
                        s.startswith("DEPOLARIZE1") or s.startswith("PAULI_CHANNEL_1") or
                        s.startswith("AMPLITUDE_DAMP")):
                    raise ValueError(
                        f"For {error_place}, error_statements must be single-qubit Pauli error, 'DEPOLARIZE1(...)', 'PAULI_CHANNEL_1(...)', or 'AMPLITUDE_DAMP(p_phase,p_amp)'."
                    )
            self._error_settings[error_place] = s
        else:
            self._error_settings[error_place] = None
        return

    def setting_error_from_dict(self, error_config: dict = {}):
        """
        Bulk-configure the noise model for the *initialisation circuit* (`self.init_circ`)
        using a single dictionary.

        This helper iterates over the four recognised error locations
        - ``'Two_qubit'``, ``'Single_qubit'``, ``'Measurement'``, and ``'Reset'`` - and
        calls :py:meth:`ErrorModel.setting_error` once for each entry that appears in
        ``error_config``.  It is a convenience wrapper that lets the user specify the
        complete error recipe in one shot instead of four individual calls.

        Parameters
        ----------
        error_config : dict, optional
            Mapping ``str → str | None`` that overrides the default noise settings.
            Keys must be drawn from
            ``{'Two_qubit', 'Single_qubit', 'Measurement', 'Reset'}`` - any key not
            present in the dictionary is left unchanged.

            * **Value is a *string*** - interpreted as a Stim noise instruction and
              passed to :py:meth:`setting_error` with ``make_error=True``.
            * **Value is ``None``** - disables that error channel
              (equivalent to ``make_error=False``).

            Example
            ````python
            {
                "Two_qubit"   : "DEPOLARIZE2(0.01)",
                "Single_qubit": "DEPOLARIZE1(0.001)",
                "Measurement" : None,           # disable measurement noise
                # "Reset" key omitted → keep current setting
            }
            ````

        Notes
        -----
        * The method affects **only** the *initialisation* circuit stored in
          ``self.init_circ``.  If you intend to apply the same configuration to the
          syndrome or measurement circuits, call the method on those
          :pyclass:`ErrorModel` instances as needed.
        * All sanity checks (valid key, instruction syntax, etc.) are delegated to
          :py:meth:`ErrorModel.setting_error`.

        Examples
        --------
        >>> cfg = {
        ...     "Two_qubit":    "DEPOLARIZE2(0.02)",
        ...     "Measurement":  None,
        ... }
        >>> em = ErrorModel(my_circuit)
        >>> em.setting_error_from_dict(cfg)
        # After the call:
        #   • two-qubit errors set to DEPOLARIZE2(0.02)
        #   • measurement noise disabled
        #   • single-qubit and reset noise keep their previous settings
        """
        for key in self._error_config:
            if key in error_config:
                statement = error_config[key]
                if statement is not None:
                    self.setting_error(error_place=key,
                                                 make_error=True,
                                                 error_statements=statement)
                else:
                    self.setting_error(error_place=key,
                                                 make_error=False)
        return
    
    def generate_noisy_circuit(self) -> stim.Circuit:
        """
        按照当前 _error_config / _error_settings 在 self.noiseless_circ 中插入噪声。
        返回 Stim 电路文本行列表。
        """
        circ_lines: list[str] = []
        qubits = [str(t) for t in range(self.noiseless_circ.num_qubits)] 
        qubit_str = " ".join(qubits)
        # print(qubits)

        for inst in self.noiseless_circ:
            name = inst.name.upper()
            targets = [str(t.value) for t in inst.targets_copy()]  # type: ignore
            target_str = " ".join(targets)
            # print(name)

            if name in {"DETECTOR", "OBSERVABLE_INCLUDE"}:
                circ_lines.append(str(inst))
                continue

            if name in {"QUBIT_COORDS"}:
                circ_lines.append(str(inst))
                continue

            # ────────────────────────────────────────────────────────────────
            # 1) 两比特 Clifford 门 ── 在门后对每一对目标插入 two-qubit 噪声
            # ────────────────────────────────────────────────────────────────
            if name in TWO_QUBIT_GATES and self._error_config["Two_qubit"]:
                circ_lines.append(f"{name} {target_str}")        # 先写原门
                err_stmt = self._error_settings["Two_qubit"]
                if err_stmt:
                    if len(targets) % 2 != 0:
                        raise ValueError(f"{name} target must be even: {targets}")
                    for i in range(0, len(targets), 2):
                        circ_lines.append(f"{err_stmt} {targets[i]} {targets[i+1]}")
                continue

            # ────────────────────────────────────────────────────────────────
            # 2) 单比特 Clifford 门 ── 在门后插入 single-qubit 噪声
            # ────────────────────────────────────────────────────────────────
            if name in SINGLE_QUBIT_GATES and self._error_config["Single_qubit"]:
                circ_lines.append(f"{name} {target_str}")
                err_stmt = self._error_settings["Single_qubit"]
                if err_stmt:
                    # Skip injecting AMPLITUDE_DAMP into Stim text; will inject in QASM later
                    if not err_stmt.strip().upper().startswith("AMPLITUDE_DAMP"):
                        circ_lines.append(f"{err_stmt} {target_str}")
                continue

            if name in {"I"}:
                circ_lines.append(f"{name} {target_str}")
                err_stmt = self._error_settings["Identity"]
                if err_stmt:
                    if not err_stmt.strip().upper().startswith("AMPLITUDE_DAMP"):
                        circ_lines.append(f"{err_stmt} {target_str}")
                continue

            # ────────────────────────────────────────────────────────────────
            # 3) 测量类 (M, MX, …) ── 在测量 **前** 插入 Measurement 噪声
            # ────────────────────────────────────────────────────────────────
            if name.startswith("MR") and self._error_config["Measurement"]:
                err_stmt = self._error_settings["Measurement"]
                if err_stmt and (not err_stmt.strip().upper().startswith("AMPLITUDE_DAMP")):
                    circ_lines.append(f"{err_stmt} {qubit_str}")
                circ_lines.append(f"{name} {target_str}")
                # if err_stmt and (not err_stmt.strip().upper().startswith("AMPLITUDE_DAMP")):
                #     circ_lines.append(f"{err_stmt} {target_str}")
                # if err_stmt and (not err_stmt.strip().upper().startswith("AMPLITUDE_DAMP")): #for long reset
                #     circ_lines.append(f"{err_stmt} {target_str}")
                continue

            # if name.startswith("M") and self._error_config["Measurement"]:
            #     err_stmt = self._error_settings["Measurement"]
            #     if err_stmt and (not err_stmt.strip().upper().startswith("AMPLITUDE_DAMP")):
            #         circ_lines.append(f"{err_stmt} {target_str}")
            #     circ_lines.append(f"{name} {target_str}")
            #     continue

            # ────────────────────────────────────────────────────────────────
            # 4) 复位类 (R, RX, …) ── 在复位 **后** 插入 Reset 噪声
            # ────────────────────────────────────────────────────────────────
            if name in {"R", "RX", "RY", "RZ"} and self._error_config["Reset"]:
                circ_lines.append(f"{name} {target_str}")
                err_stmt = self._error_settings["Reset"]
                if err_stmt and (not err_stmt.strip().upper().startswith("AMPLITUDE_DAMP")):
                    circ_lines.append(f"{err_stmt} {target_str}")
                continue

            # 5) 其它门：直接抄过去
            circ_lines.append(f"{name} {target_str}")
            # print(circ_lines)

        self.noisy_circ = stim.Circuit("\n".join(circ_lines))
        return self.noisy_circ
    
    def single_noise_specs(self, overwrite=False):
        '''
        generate the single-noise specs for later use based on the noisy circuit.
        The specs 
        '''
        if overwrite:
            self.noise_specs = None
        
        if self.noise_specs is not None:
            return self.noise_specs
        
        try:
            noisy_circ = self.noisy_circ
        except AttributeError:
            noisy_circ = self.generate_noisy_circuit()

        self.noise_specs = generate_single_error_specs(noisy_circ)
        return self.noise_specs
    
    def single_noise_profile(self):
        '''
        looping over all the single-noise configurations. Run a single shot and determine
        the measurement outcome. We then use the ideal case as a reference to mark the
        affected detectors. This is similar to generate a rough (un-merged) detector
        error model file. 
        '''
        ## this is the noise spects
        noise_specs = self.single_noise_specs()

        ## use the spect to generate the circuits, and compute the measurement shots
        ## and use the noiseless convert to convert to the detector events.

        measure = []
        det = []

        err_det_results = {}
        for j, spec in enumerate(noise_specs):
            line_idx, tokens, err_p = spec
            clean_circ = list(self.noiseless_circ_list)
            gate_strs = [f"{tok[0]} {tok[1:]}" for tok in tokens]
            clean_circ[line_idx:line_idx] = gate_strs

            noisy_circ = stim.Circuit('\n'.join(clean_circ))
            #print(noisy_circ)

            ## reset the simulator
            sim = stim.TableauSimulator()
            #sim.reset()
            sim.do(noisy_circ)
            measure = sim.current_measurement_record()
            #print(measure)

            ## convert to detector results
            det, _ = self.conv.convert(measurements=np.array([measure], dtype=np.bool_), separate_observables=True)
            ## check if this measurement results has been visited or not
            det_packed = tuple(np.packbits(det, bitorder='big'))
            if det_packed in err_det_results:
                val = err_det_results[det_packed]
                val.append((j, err_p))
            
            else:
                err_det_results[det_packed] = [(j, err_p)]
        
        self.total_measurements = len(measure)
        self.total_detectors = len(det)
    
        self.err_det_result = err_det_results
        return err_det_results
    
    def generate_noise_circuit_from_specs(self,
                                        noise_specs: NoiseSpecs, 
                                        return_stim: bool = False,
        ) -> Union[stim.Circuit, List[str]]:
        """
        Build a new circuit by inserting *explicit* Pauli gates according to
        ``noise_specs``.

        A ``noise_specs`` entry is a triple:

            (line_idx, pauli_tokens, probability)

        where
            • line_idx       - 0-based index in the **noise-free** circuit
                               (`self.noiseless_circ_list`) where the Pauli gates
                               should be **inserted *before*** that line;
            • pauli_tokens   - tuple like ('X9',) or ('Z1', 'Y4');
            • probability    - branch weight (ignored here but kept for completeness).

        Notes
        -----
        *   DETECTOR / OBSERVABLE_INCLUDE lines are never present in
            ``noise_profile`` - caller guarantees this.
        *   If ``line_idx`` ≥ depth of the noise-free circuit, the entry is ignored.
        *   Setting ``return_stim=True`` returns a ``stim.Circuit`` object;
            otherwise the flattened list[str] representing the circuit.

        Parameters
        ----------
        noise_specs : list[tuple[int, tuple[str, ...], float]]
            The Pauli-error insertion specification.

        return_stim : bool, default ``False``
            If True  - return a ready-to-use ``stim.Circuit``.  
            If False - return the list of operations (strings).

        Returns
        -------
        stim.Circuit | list[str]
            The newly built circuit in the requested form.

        Raises
        ------
        ValueError
            If an invalid Pauli token appears in ``noise_profile``.
        """
        ## sweeping the noise_profile, group the operations based on the inserting positions.
        noise_dict: dict[int, List[str]] = defaultdict(list)
        token_re = re.compile(r"^[XYZ]([+-]?\d+)$")        # X9, Y-1, Z12 …

        for pos, pauli_tokens, _p in noise_specs:
            for tok in pauli_tokens:
                # Basic sanity check
                if not token_re.match(tok):
                    raise ValueError(f"Invalid Pauli token in noise_profile: '{tok}'")
                gate_str = f"{tok[0]} {tok[1:]}"           # 'X9' -> "X 9"
                noise_dict[pos].append(gate_str)
        
        # ── 2) Rebuild circuit line-by-line
        new_lines: List[str] = []
        depth = len(self.noiseless_circ_list)

        for idx, line in enumerate(self.noiseless_circ_list):
            # Insert (if any) *before* current line
            if idx in noise_dict:
                new_lines.extend(noise_dict[idx])
            new_lines.append(line)

        # ── 3) Optionally append noise for line == depth  (≥ depth is ignored)
        if depth in noise_dict:
            new_lines.extend(noise_dict[depth])

        # ── 4) Return in desired form
        return stim.Circuit("\n".join(new_lines)) if return_stim else new_lines

        # noise_dict = {}
        # for position, ops, _ in noise_profile:
        #     if position not in noise_dict:
        #         noise_dict[position] = []

        #     operations = noise_dict[position]
        #     for op in ops:
        #         operations.append(f'{op[0]} {op[1:]}')
        
        # new_circuit = []
        # current_line_number = 0
        
        # for gate in self.noiseless_circ_list:
        #     if current_line_number in noise_dict:
        #         new_circuit += noise_dict[current_line_number]
            
        #     new_circuit.append(gate)
        #     current_line_number += 1
        
        # if current_line_number in noise_dict:
        #     new_circuit += noise_dict[current_line_number]
        
        # if return_stim:
        #     return stim.Circuit('\n'.join(new_circuit))
        # else:
        #     return new_circuit
                 
    def generate_error_detector_model(self):
        '''generate the error detector model'''
        try:
            err_det_result = self.err_det_result
        except AttributeError:
            err_det_result = self.single_noise_profile()
        
        err_tracking = []
        error_statements = []
        counter = 0
        for key, val in err_det_result.items():
            if np.sum(key) == 0:
                continue
            det_result = np.unpackbits(np.array(key, dtype=np.uint8), bitorder='big')
            locs = np.where(det_result == 1)[0]

            val = np.array(val)
            total_p = np.sum(val[:, 1])
            err_tracking.append(val[:,0].astype(int))
            
            es = f'error({total_p}) '
            es = es + ' '.join(map(lambda x: f'D{int(x)}', locs))
            es = es + f' L{counter:d}'
            error_statements.append(es)
            
            counter += 1

        self.dem = error_statements
        self.err_tracking = err_tracking
        return
    
    def single_error_logical_state(self, measure_labels):
        '''
        focus on all the single-error configurations.
        Check the no-detection events. find out their probability, 
        and if the final state is error-free or not.

        Parameters:
        -----------
        measure_labels: int or list
            if given as an int: this menas the last n 
                measurements are corresponding to the code stabilizer
            if given as a list: this means the measurements with these indices 
                should be counted as the final stabilizer measurement outcomes.
        '''
        pass
        


    def mc_sample_errors(self):
        '''
        using MC to sample error configurations, and accumulate statistics for logical error rates.
        '''
        pass
            


# ───────────────────────────────────────────────────────────────────────────────
# Helper: expand *one* noise instruction into alternative Pauli‑gate lists
# ───────────────────────────────────────────────────────────────────────────────
def _expand_error_instruction(inst: stim.CircuitInstruction) -> List[Tuple[List[str], float]]:
    """
    Given a Stim noise instruction, return all possible single-error alternatives as
    (gates, probability) tuples.

    Returns
    -------
    List[Tuple[List[str], float]]
        Each tuple is (gates, branch_probability).
    """
    name = inst.name.upper()
    args = inst.gate_args_copy()
    # Default: probability for each branch
    if name in {"X_ERROR", "Y_ERROR", "Z_ERROR"}:
        pauli = name[0]
        p = args[0] if args else 0.0
        targets = [t.value for t in inst.targets_copy()]
        n = len(targets)
        branch_prob = p / n if n > 0 else 0.0
        return [([f"{pauli} {q}"], branch_prob) for q in targets]

    if name == "DEPOLARIZE1":
        p = args[0] if args else 0.0
        targets = [t.value for t in inst.targets_copy()]
        # Each qubit experiences an independent depolarizing channel of strength p.
        # For that qubit the probability of each single‑Pauli error (X, Y, or Z) is p/3.
        branch_prob = p / 3
        branches = []
        for q in targets:
            for pauli in ("X", "Y", "Z"):
                branches.append(([f"{pauli} {q}"], branch_prob))
        return branches

    if name == "PAULI_CHANNEL_1":
        # Args: p_x, p_y, p_z
        args = inst.gate_args_copy()
        px, py, pz = (args + [0.0, 0.0, 0.0])[:3]
        targets = [t.value for t in inst.targets_copy()]
        branches = []
        for q in targets:
            if px > 0:
                branches.append(([f"X {q}"], px))
            if py > 0:
                branches.append(([f"Y {q}"], py))
            if pz > 0:
                branches.append(([f"Z {q}"], pz))
        return branches

    if name == "DEPOLARIZE2":
        # Stim 允许在一次指令里写多个 (q0 q1 q2 q3 …)，
        # 每两个连续索引视为一对独立的两比特通道。
        p = args[0] if args else 0.0
        qs = [t.value for t in inst.targets_copy()]

        if len(qs) % 2 != 0:
            raise ValueError("DEPOLARIZE2 target list must have EVEN length")

        branches = []
        branch_prob = p / 15 if p > 0 else 0.0

        # 逐对处理：(q0,q1), (q2,q3), ...
        for i in range(0, len(qs), 2):
            q0, q1 = qs[i], qs[i + 1]
            for p0, p1 in PAULI_DOUBLE:
                if p0 == p1 == "I":
                    continue
                gates = []
                if p0 != "I":
                    gates.append(f"{p0} {q0}")
                if p1 != "I":
                    gates.append(f"{p1} {q1}")
                branches.append((gates, branch_prob))

        return branches

    if name == "E":
        # E(p) X1 Y2 Z3 ...
        args = inst.gate_args_copy()
        p = args[0] if args else 0.0
        tokens = str(inst).split()[1:]  # skip "E(p)"
        paulis = [(tok[0], int(tok[1:])) for tok in tokens]  # ('X',3), …
        # only one branch, prob = p
        return [([f"{pauli} {q}" for pauli, q in paulis], p)]

    raise NotImplementedError(f"Unsupported noise gate: {name}")


# ───────────────────────────────────────────────────────────────────────────────
# Main routine: expand a noisy circuit into *single‑error* circuit instances
# ───────────────────────────────────────────────────────────────────────────────
def generate_single_error_circuits(noisy_circ: stim.Circuit):
    """
    Expand *one* noisy Stim circuit into a list of *deterministic* circuits,
    each containing **exactly one** explicit Pauli error *branch* corresponding to the
    single-error configurations implied by the original noise instructions.

    Workflow (requested by Chen)
    ----------------------------
    0. Convert ``noisy_circ`` into a flattened sequence of instructions.
    1. Sweep once to build
       (a)-``clean_ops``      : list[str] of *non-noise* operations;
       (b) ``error_ops``      : list[stim.CircuitInstruction] of noise terms;
       (c) ``err2idx``        : map: error-op-index → insertion position in
           ``clean_ops`` where the error originally occurred.
    2. For every entry in ``error_ops``:
       * Decompose it into explicit Pauli gate strings.
       * Insert those gates at the recorded position *inside a copy* of
         ``clean_ops``.
       * Convert the resulting list back into a ``stim.Circuit`` and append to
         the output list.

    Parameters
    ----------
    noisy_circ : stim.Circuit
        The *noisy* circuit whose single-error expansions are to be generated.

    Returns
    -------
    list[stim.Circuit]
        A list where each element is a Stim circuit encoding one concrete,
        deterministic single-error configuration.
    """
    # ── step 0: prepare containers
    clean_ops: list[str] = []
    error_ops: list[stim.CircuitInstruction] = []
    err2idx: list[int] = []

    # ── step 1: single pass over the flattened circuit
    for inst in noisy_circ.flattened():
        if inst.name in NOISE_GATES:
            # Record where this noise instruction *would have* appeared
            err_idx = len(error_ops)
            err2idx.append(len(clean_ops))
            error_ops.append(inst) #type: ignore
        else:
            clean_ops.append(str(inst))

    # ── step 2: build single‑error circuits and collect probabilities
    single_error_circs: list[stim.Circuit] = []
    single_error_probs: list[float] = []
    for err_idx, err_inst in enumerate(error_ops):
        insert_at = err2idx[err_idx]
        for gates, prob in _expand_error_instruction(err_inst):
            new_ops = clean_ops.copy()
            new_ops[insert_at:insert_at] = gates
            single_error_circs.append(stim.Circuit("\n".join(new_ops)))
            single_error_probs.append(prob)
    return single_error_circs, single_error_probs


# ───────────────────────────────────────────────────────────────────────────────
# Utility: convert gate strings like "X 0" into compact tokens "X0"
# ───────────────────────────────────────────────────────────────────────────────
def _compact_gate_tokens(gates: list[str]) -> tuple[str, ...]:
    """Convert a list of gate strings (e.g. "X 0") into a tuple like ('X0', 'Z3')."""
    tokens: list[str] = []
    for g in gates:
        parts = g.split()
        if len(parts) == 2:
            tokens.append(f"{parts[0]}{parts[1]}")
        else:
            # Fallback – keep original text if unexpected format
            tokens.append(parts[0])
    return tuple(tokens)


# ───────────────────────────────────────────────────────────────────────────────
# New routine: produce *specification list* describing single‑error branches
# ───────────────────────────────────────────────────────────────────────────────
def generate_single_error_specs(noisy_circ: stim.Circuit):
    """
    Produce a compact list of specifications for every single‑error configuration
    implied by ``noisy_circ``.

    Each element in the returned list is a tuple:

        (line_idx, gates_tuple, probability)

    where
        • line_idx      – index (0‑based) in the *noise‑free* circuit where the
                          Pauli error(s) should be inserted;
        • gates_tuple   – tuple such as ('X0',) or ('Z1', 'X4') describing the
                          explicit Pauli gates that realise this branch;
        • probability   – the branch probability associated with this error.

    This helper does **not** build new circuits; it only returns metadata that
    can be stored, analysed, or fed into later processing pipelines.
    """
    # Containers for a single sweep
    clean_ops: list[str] = []
    error_ops: list[stim.CircuitInstruction] = []
    err2idx: list[int] = []

    # 1) Separate noise vs non‑noise instructions
    for inst in noisy_circ.flattened():
        if inst.name in NOISE_GATES:
            err2idx.append(len(clean_ops))
            error_ops.append(inst) #type: ignore
        else:
            clean_ops.append(str(inst))

    # 2) Build the specification list
    specs: list[tuple[int, tuple[str, ...], float]] = []
    for err_idx, err_inst in enumerate(error_ops):
        line_idx = err2idx[err_idx]
        for gates, prob in _expand_error_instruction(err_inst):
            tokens = _compact_gate_tokens(gates)
            specs.append((line_idx, tokens, prob))

    return specs


def dem_det_coords(dem_text: str) -> Dict[int, Tuple[int,int]]:
    """
    Parse a Stim detector-error-model string and map detector indices to (x, y).

    Parameters
    ----------
    dem_text : str
        Full DEM text containing lines like
            detector(3, 5, 1) D8
            detector(-2.0, 4.0) D17

    Returns
    -------
    dict[int, (float, float)]
        Mapping: detector index → (x, y) coordinates.
    """
    # 捕获 detector 行：允许有或没有第三个坐标
    det_re = re.compile(
        r"detector\(\s*([-\d.]+)\s*,\s*([-\d.]+)"
        r"(?:\s*,\s*([-\d.]+))?\s*\)\s*D(\d+)",      # 第三项 (z) 可选
        re.IGNORECASE
    )

    det_map: Dict[int, Tuple[int,int]] = {}
    for match in det_re.finditer(dem_text):
        x_str, y_str, _z_str, d_str = match.groups()
        det_map[int(d_str)] = (int(x_str), int(y_str))

    return det_map


# --- Injection of AMPLITUDE_DAMP into QASM based on ErrorModel settings ---

def _parse_amp_damp_args(stmt: Optional[str]) -> Optional[tuple[float, float]]:
    if not stmt:
        return None
    s = stmt.strip().upper()
    if not s.startswith("AMPLITUDE_DAMP"):
        return None
    try:
        inside = s[s.find("(") + 1:s.find(")")]
        parts = [p.strip() for p in inside.split(',') if p.strip()]
        if not parts:
            return None
        p1 = float(parts[0])
        p2 = float(parts[1]) if len(parts) > 1 else 0.0
        return (p1, p2)
    except Exception:
        return None

def inject_amplitude_damp(qasm_text: str, error_model: "ErrorModel") -> str:
    """
    Insert damp(p1,p2) into a QASM string based on the noise model settings.
    This function iterates through the QASM file and injects `damp` noise
    instructions according to the policies defined in the `ErrorModel`.

    - For 'Reset' noise: `damp` is inserted *after* each `reset` operation.
    - For 'Measurement' noise: `damp` is inserted *before* each block of
      `measure` operations, once per qubit in the circuit.
    """
    p_meas_args = _parse_amp_damp_args(error_model._error_settings.get("Measurement"))
    p_reset_args = _parse_amp_damp_args(error_model._error_settings.get("Reset"))

    if not p_meas_args and not p_reset_args:
        return qasm_text

    lines = qasm_text.splitlines()
    new_lines = []
    
    # Regex to find qubit targets like q[0], q[12], etc.
    qubit_target_re = re.compile(r"q\[(\d+)\]")

    in_measurement_block = False
    for line in lines:
        stripped_line = line.strip()
        is_measurement = stripped_line.startswith("measure") or stripped_line.startswith("m ")

        # --- Handle Measurement ---
        if p_meas_args and is_measurement:
            if not in_measurement_block:
                # This is the start of a new measurement block.
                # Inject damp for all qubits.
                in_measurement_block = True
                indent = line[:line.find(stripped_line[0])] if stripped_line else ""
                
                num_qubits = error_model.noiseless_circ.num_qubits
                for i in range(num_qubits):
                    damp_inst = f"damp({p_meas_args[0]},{p_meas_args[1]}) q[{i}];"
                    new_lines.append(indent + damp_inst)
            
            new_lines.append(line)
        
        # --- Handle Reset ---
        elif p_reset_args and stripped_line.startswith("reset"):
            in_measurement_block = False # Reset ends a measurement block
            targets = qubit_target_re.findall(line)
            new_lines.append(line) # Add reset instruction first
            indent = line[:line.find(stripped_line[0])] if stripped_line else ""
            for target in targets:
                damp_inst = f"damp({p_reset_args[0]},{p_reset_args[1]}) q[{target}];"
                new_lines.append(indent + damp_inst)
        
        else:
            # Any other instruction also ends a measurement block
            in_measurement_block = False
            new_lines.append(line)
            
    return "\n".join(new_lines)

    # Track the first contiguous block of RESET statements
    in_first_reset_block = False
    first_reset_block_done = False

    for raw_line in qasm_text.splitlines():
        # Preserve completely empty lines
        if not raw_line.strip():
            out_lines.append(raw_line)
            continue

        # Split trailing comment (// ...)
        code_part = raw_line
        comment_part = ""
        cidx = raw_line.find("//")
        if cidx != -1:
            code_part = raw_line[:cidx]
            comment_part = raw_line[cidx:]

        # If the line contains only comment after stripping code_part
        if not code_part.strip():
            out_lines.append(raw_line)
            continue

        indent_match = indent_re.match(code_part)
        indent = indent_match.group(1) if indent_match else ""

        # Split into statements while preserving order; we will re-emit per statement
        stmts = [s.strip() for s in code_part.split(';') if s.strip()]

        for stmt in stmts:
            r = reset_stmt_re.match(stmt)
            if r:
                qtok = r.group(1)
                # Seeing reset: if first block not finished, we are inside it
                if not first_reset_block_done:
                    in_first_reset_block = True
                # Emit RESET
                out_lines.append(f"{indent}{stmt};")
                # If inside first contiguous reset block, inject POST-RESET damp
                if in_first_reset_block:
                    p = _pick_post_reset_params()
                    if p:
                        out_lines.append(f"{indent}damp({p[0]},{p[1]}) {qtok};")
                # Do not inject for later resets (handled before measurements instead)
                continue

            m = meas_stmt_re.match(stmt)
            if m:
                qtok = m.group(1)
                # Any non-reset breaks the initial reset block if it was in progress
                if in_first_reset_block and not first_reset_block_done:
                    in_first_reset_block = True #flip to make all resets have error after
                    first_reset_block_done = False
                # Inject PRE-MEAS damp (for all subsequent rounds)
                p = _pick_pre_meas_params()
                if p:
                    out_lines.append(f"{indent}damp({p[0]},{p[1]}) {qtok};")
                # Emit MEAS after the damp
                out_lines.append(f"{indent}{stmt};")
                continue

            if in_first_reset_block and not first_reset_block_done:
                in_first_reset_block = True #flip to make all resets have error after
                first_reset_block_done = False

            out_lines.append(f"{indent}{stmt};")

        # Attach trailing comment as its own line to avoid corrupting code layout
        if comment_part:
            out_lines.append(f"{indent}{comment_part}")

    return "\n".join(out_lines)

if __name__ == '__main__':
    # Test code that constructs an example circuit with DEPOLARIZE1 on multiple targets,
    # invokes generate_single_error_circuits, and prints each resulting circuit and probability.
    example_circuit = stim.Circuit("""
    H 0
    DEPOLARIZE1(0.1) 0 1
    CX 0 1
    E(0.1) X1
    M 0
    DETECTOR rec[-1]
    """)
    single_error_circuits, single_error_probs = generate_single_error_circuits(example_circuit)
    for i, (circ, prob) in enumerate(zip(single_error_circuits, single_error_probs)):
        print(f"Single error circuit {i+1}: (prob={prob})")
        print(circ)
        print()
    
    specs = generate_single_error_specs(example_circuit)
    print(specs)
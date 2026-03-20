"""Automatic noise labeling for Clifford circuits and stabilizer codes.

Assigns each noise source (each single-qubit Pauli channel, uniquely
indexed) one of four QStab IR error types plus descriptive tags:

**Primary label** (exactly one):

- **Type 0** (``data_qubit_error``, :math:`\\lambda_0`): Single-qubit
  Pauli error on a data qubit **before** any stabilizer CX gate starts
  or **after** the qubit's last CX — absorbed into initial/final error
  flow :math:`\\widetilde{E}`.

- **Type I** (``ghost_error``, :math:`\\lambda_1`): Error on a data
  qubit **during** the stabilizer measurement phase where the qubit
  still has future CX gates with an ancilla.  The error propagates
  through a future CX, flipping the measurement outcome while
  persisting on the data qubit.

- **Type II** (``hook_error``, :math:`\\lambda_2`): Error on the ancilla
  qubit with remaining CX gates — back-propagates through the CNOT
  chain to data qubits (:math:`e_B \\in \\mathcal{B}^1(\\mathcal{T}_s)`,
  :math:`|e_B| \\leq r`).

- **Type III** (``measurement_error``, :math:`\\lambda_3`): Error on
  the ancilla after its last CX — only flips the classical measurement
  register :math:`G[c_l]` without affecting data qubit error flow.

**Tags** (zero or more, orthogonal to primary label):

- ``data_qubit``: noise acts on a data qubit.
- ``ancilla_qubit``: noise acts on an ancilla qubit.

For DEPOLARIZE2 noise, each qubit receives its own noise index and is
classified independently.

The error propagation invariant from the fault-tolerance proof:

.. math::

    \\lambda_0 + \\lambda_1 + r\\lambda_2 + \\lambda_3 \\geq \\lambda_E

Typical usage::

    from scalerqec.QEC.autolabel import auto_label_from_circuit, auto_label_from_stabcode

    label_map = auto_label_from_circuit(circuit, num_data_qubits=7)
    label_map = auto_label_from_stabcode(code)
"""

from __future__ import annotations

from ..Clifford.clifford import (
    CliffordCircuit,
    Measurement,
    Reset,
    SingleQGate,
    TwoQGate,
    pauliNoise,
)
from ..Clifford.noiselabel import NoiseLabelMap

# The four canonical error types from the QStab IR semantics.
TYPE0 = "data_qubit_error"   # λ₀: data qubit error outside measurement phase
TYPE1 = "ghost_error"        # λ₁: data qubit error during phase, future CX
TYPE2 = "hook_error"         # λ₂: ancilla error with remaining CX
TYPE3 = "measurement_error"  # λ₃: ancilla error after last CX (flips M)

# Descriptive tags (orthogonal to primary label)
TAG_DATA = "data_qubit"
TAG_ANCILLA = "ancilla_qubit"


def auto_label_from_circuit(
    circuit: CliffordCircuit,
    num_data_qubits: int | None = None,
) -> NoiseLabelMap:
    """Label every noise source using the four QStab IR error types.

    Classification is **positional**: it depends on where the noise
    occurs relative to the CX schedule of stabilizer gadgets, not just
    the adjacent gate.

    For data qubits:

    - **Before** the first CX(data↔ancilla) in the circuit → Type 0
    - **During** the measurement phase, with future CX remaining for
      this qubit → Type I (ghost)
    - **During** the measurement phase, with no future CX → Type 0
    - **After** all ancilla measurements → Type 0

    For ancilla qubits:

    - With future CX remaining → Type II (hook, back-propagation)
    - After last CX, before/at measurement → Type III

    For DEPOLARIZE2 noise (2 consecutive pauliNoise after a CX on the
    CX's qubits), the noise is associated with the preceding CX for
    context, but classification still uses positional rules.

    Args:
        circuit: A compiled :class:`CliffordCircuit` with gates and noise.
        num_data_qubits: Number of data qubits (indices ``0..n-1``).

    Returns:
        A populated :class:`NoiseLabelMap`.
    """
    label_map = NoiseLabelMap()
    gatelist = circuit._gatelists
    n = len(gatelist)

    if num_data_qubits is None:
        # Without layout info, fall back to simple gate-based heuristic
        return _label_without_layout(circuit)

    # ------------------------------------------------------------------
    # Pre-scan: build positional information for the CX schedule
    # ------------------------------------------------------------------

    # Collect positions of CX(data↔ancilla) gates per qubit
    cx_positions: dict[int, list[int]] = {}
    first_cx_pos: int = n   # sentinel: no CX found
    last_ancilla_m_pos: int = -1  # sentinel: no ancilla M found

    for i, gate in enumerate(gatelist):
        if isinstance(gate, TwoQGate) and gate._name == "CNOT":
            ctrl_data = gate._control < num_data_qubits
            tgt_data = gate._target < num_data_qubits
            if ctrl_data != tgt_data:  # mixed data↔ancilla CX
                cx_positions.setdefault(gate._control, []).append(i)
                cx_positions.setdefault(gate._target, []).append(i)
                first_cx_pos = min(first_cx_pos, i)
        elif isinstance(gate, Measurement) and gate._qubitindex >= num_data_qubits:
            last_ancilla_m_pos = max(last_ancilla_m_pos, i)

    # ------------------------------------------------------------------
    # Pass 1: Claim DEPOLARIZE2 noise (for context enrichment only)
    # ------------------------------------------------------------------
    claimed_context: dict[int, TwoQGate] = {}

    for i, gate in enumerate(gatelist):
        if not isinstance(gate, TwoQGate):
            continue
        if (i + 2 < n
                and isinstance(gatelist[i + 1], pauliNoise)
                and isinstance(gatelist[i + 2], pauliNoise)):
            q_pair = {
                gatelist[i + 1]._qubitindex,
                gatelist[i + 2]._qubitindex,
            }
            if q_pair == {gate._control, gate._target}:
                claimed_context[i + 1] = gate
                claimed_context[i + 2] = gate

    # ------------------------------------------------------------------
    # Pass 2: Classify every noise entry by position
    # ------------------------------------------------------------------
    for i, gate in enumerate(gatelist):
        if not isinstance(gate, pauliNoise):
            continue

        qubit = gate._qubitindex
        is_data = qubit < num_data_qubits
        tags: set[str] = {TAG_DATA if is_data else TAG_ANCILLA}

        # Build context from claimed DEPOLARIZE2 or from following gate
        context: dict = {}
        if i in claimed_context:
            assoc = claimed_context[i]
            context = {
                "gate": assoc._name,
                "control": assoc._control,
                "target": assoc._target,
            }

        # Get future CX positions for this qubit
        qubit_cxs = cx_positions.get(qubit, [])
        has_future_cx = any(cx_pos > i for cx_pos in qubit_cxs)

        if is_data:
            label = _classify_data_noise(
                i, first_cx_pos, last_ancilla_m_pos, has_future_cx
            )
        else:
            label = _classify_ancilla_noise(has_future_cx)

        label_map.add(
            gate._noiseindex,
            label,
            qubit=qubit,
            tags=tags,
            context=context,
        )

    return label_map


def _classify_data_noise(
    noise_pos: int,
    first_cx_pos: int,
    last_ancilla_m_pos: int,
    has_future_cx: bool,
) -> str:
    """Classify data qubit noise by position in the CX schedule.

    - Before first CX → Type 0 (initial error flow)
    - After all ancilla M → Type 0 (final error state)
    - During measurement phase, future CX → Type I (ghost)
    - During measurement phase, no future CX → Type 0 (idle/done)
    """
    # Before any stabilizer CX has started
    if noise_pos < first_cx_pos:
        return TYPE0

    # After all ancilla measurements
    if last_ancilla_m_pos >= 0 and noise_pos > last_ancilla_m_pos:
        return TYPE0

    # During measurement phase
    if has_future_cx:
        return TYPE1  # ghost: will propagate through future CX
    return TYPE0  # idle or done with all CX


def _classify_ancilla_noise(has_future_cx: bool) -> str:
    """Classify ancilla qubit noise by remaining CX gates.

    - Future CX remaining → Type II (hook, back-propagation)
    - No future CX → Type III (measurement error)
    """
    if has_future_cx:
        return TYPE2  # hook: error will back-propagate through remaining CX
    return TYPE3  # measurement: only flips classical register


def _label_without_layout(circuit: CliffordCircuit) -> NoiseLabelMap:
    """Fallback labeling when num_data_qubits is unknown.

    Uses a simple gate-adjacency heuristic: classifies noise by
    the next non-noise gate.
    """
    label_map = NoiseLabelMap()
    gatelist = circuit._gatelists
    pending: list[pauliNoise] = []

    for gate in gatelist:
        if isinstance(gate, pauliNoise):
            pending.append(gate)
            continue

        for noise in pending:
            if isinstance(gate, Measurement):
                label = TYPE3
            elif isinstance(gate, TwoQGate):
                label = TYPE2
            elif isinstance(gate, SingleQGate):
                label = TYPE2
            else:
                label = TYPE0
            label_map.add(noise._noiseindex, label, qubit=noise._qubitindex)
        pending.clear()

    for noise in pending:
        label_map.add(noise._noiseindex, TYPE0, qubit=noise._qubitindex)

    return label_map


def _get_gate_qubits(gate: object) -> set[int]:
    """Return the set of qubits a gate acts on."""
    if isinstance(gate, TwoQGate):
        return {gate._control, gate._target}
    if isinstance(gate, (SingleQGate, Measurement, Reset)):
        return {gate._qubitindex}
    return set()


def auto_label_from_stabcode(code: object) -> NoiseLabelMap:
    """Label noise sources in a compiled :class:`StabCode` with full context.

    Each noise entry receives:

    - A primary label (Type 0/I/II/III).
    - Tags (``data_qubit`` / ``ancilla_qubit``).
    - Context enriched with ``round``, ``stabilizer_index``,
      ``stabilizer`` from the IR.

    Args:
        code: A compiled :class:`~scalerqec.QEC.qeccircuit.StabCode`.

    Returns:
        A populated :class:`NoiseLabelMap`.
    """
    from .qeccircuit import StabCode, StabPropInstruction

    assert isinstance(code, StabCode), "Expected a StabCode instance"
    assert code._circuit_compiled, (
        "StabCode must be compiled before auto-labeling. "
        "Call code.construct_circuit() first."
    )

    circuit: CliffordCircuit = code._circuit
    num_data_qubits = code._n
    num_stabs = len(code._stabs)

    # Step 1: gate-level labeling with four types + tags
    label_map = auto_label_from_circuit(circuit, num_data_qubits)

    # Step 2: enrich with round/stabilizer context from IR
    if not code._IR_compiled or not code._IRList:
        return label_map

    # Build round schedule from IR
    round_stab_info: list[tuple[int, int, str]] = []
    for irinst in code._IRList:
        if isinstance(irinst, StabPropInstruction):
            round_stab_info.append(
                (irinst.round, irinst.get_stabindex(), irinst.stab)
            )

    # Track current round per ancilla by counting resets
    ancilla_current_round: dict[int, int] = {}
    for gate in circuit._gatelists:
        if isinstance(gate, Reset):
            q = gate._qubitindex
            if num_data_qubits <= q < num_data_qubits + num_stabs:
                ancilla_current_round[q] = (
                    ancilla_current_round.get(q, -1) + 1
                )

    # Re-walk and assign round/stabilizer context
    ancilla_current_round.clear()
    for gate in circuit._gatelists:
        if isinstance(gate, Reset):
            q = gate._qubitindex
            if num_data_qubits <= q < num_data_qubits + num_stabs:
                ancilla_current_round[q] = (
                    ancilla_current_round.get(q, -1) + 1
                )

        elif isinstance(gate, pauliNoise):
            noise_idx = gate._noiseindex
            qubit = gate._qubitindex
            entry = label_map.get_entry(noise_idx)
            if entry is None:
                continue

            # Noise on an ancilla qubit — direct context
            if num_data_qubits <= qubit < num_data_qubits + num_stabs:
                ancilla_q = qubit
                rnd = ancilla_current_round.get(ancilla_q, 0)
                stab_idx = ancilla_q - num_data_qubits
                entry.context["round"] = rnd
                entry.context["stabilizer_index"] = stab_idx
                if stab_idx < len(code._stabs):
                    entry.context["stabilizer"] = code._stabs[stab_idx]

            # Noise on a data qubit during entangling — inherit
            # the associated ancilla's round context
            elif (
                qubit < num_data_qubits
                and "control" in entry.context
                and "target" in entry.context
            ):
                ctrl = entry.context["control"]
                tgt = entry.context["target"]
                ancilla_q = None
                if num_data_qubits <= ctrl < num_data_qubits + num_stabs:
                    ancilla_q = ctrl
                elif num_data_qubits <= tgt < num_data_qubits + num_stabs:
                    ancilla_q = tgt

                if ancilla_q is not None:
                    rnd = ancilla_current_round.get(ancilla_q, 0)
                    stab_idx = ancilla_q - num_data_qubits
                    entry.context["round"] = rnd
                    entry.context["stabilizer_index"] = stab_idx
                    if stab_idx < len(code._stabs):
                        entry.context["stabilizer"] = code._stabs[stab_idx]

    return label_map

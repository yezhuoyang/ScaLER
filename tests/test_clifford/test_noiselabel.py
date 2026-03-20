"""Tests for the noise labeling system (four QStab IR error types + tags)."""

import json
import pytest
import numpy as np

from scalerqec.Clifford.noiselabel import NoiseLabel, NoiseLabelMap
from scalerqec.Clifford.clifford import CliffordCircuit
from scalerqec.QEC.autolabel import (
    auto_label_from_circuit,
    TYPE0,
    TYPE1,
    TYPE2,
    TYPE3,
    TAG_DATA,
    TAG_ANCILLA,
)


class TestNoiseLabel:
    def test_basic_creation(self):
        label = NoiseLabel(noise_index=0, label="hook_error", qubit=3)
        assert label.noise_index == 0
        assert label.label == "hook_error"
        assert label.qubit == 3
        assert label.tags == set()
        assert label.context == {}

    def test_with_tags_and_context(self):
        label = NoiseLabel(
            noise_index=5,
            label="measurement_error",
            qubit=7,
            tags={"ancilla_qubit", "round_2"},
            context={"round": 2, "stabilizer_index": 1},
        )
        assert "ancilla_qubit" in label.tags
        assert label.context["round"] == 2


class TestNoiseLabelMap:
    def test_add_and_query(self):
        lm = NoiseLabelMap()
        lm.add(0, "hook_error", qubit=3, tags={"ancilla_qubit"})
        lm.add(1, "measurement_error", qubit=5, tags={"ancilla_qubit"})

        assert lm.get_label(0) == "hook_error"
        assert lm.get_label(1) == "measurement_error"
        assert lm.get_label(99) is None

    def test_get_indices(self):
        lm = NoiseLabelMap()
        lm.add(0, "hook_error", qubit=3)
        lm.add(1, "hook_error", qubit=4)
        lm.add(2, "measurement_error", qubit=5)

        assert lm.get_indices("hook_error") == [0, 1]
        assert lm.get_indices("measurement_error") == [2]
        assert lm.get_indices("nonexistent") == []

    def test_get_indices_by_tag(self):
        lm = NoiseLabelMap()
        lm.add(0, "hook_error", tags={"ancilla_qubit"})
        lm.add(1, "ghost_error", tags={"data_qubit"})
        lm.add(2, "measurement_error", tags={"ancilla_qubit"})

        assert lm.get_indices_by_tag("ancilla_qubit") == [0, 2]
        assert lm.get_indices_by_tag("data_qubit") == [1]
        assert lm.get_indices_by_tag("nonexistent") == []

    def test_add_tag(self):
        lm = NoiseLabelMap()
        lm.add(0, "hook_error", qubit=3)
        lm.add_tag(0, "ancilla_qubit")
        lm.add_tag(0, "round_0")

        assert lm.get_tags(0) == {"ancilla_qubit", "round_0"}
        assert lm.get_indices_by_tag("ancilla_qubit") == [0]

    def test_replace_existing(self):
        lm = NoiseLabelMap()
        lm.add(0, "hook_error", qubit=3, tags={"ancilla_qubit"})
        lm.add(0, "measurement_error", qubit=3, tags={"ancilla_qubit"})

        assert lm.get_label(0) == "measurement_error"
        assert len(lm) == 1
        assert lm.get_indices("hook_error") == []
        assert lm.get_indices_by_tag("ancilla_qubit") == [0]

    def test_categories(self):
        lm = NoiseLabelMap()
        lm.add(0, "hook_error")
        lm.add(1, "measurement_error")
        lm.add(2, "hook_error")

        assert lm.categories == ["hook_error", "measurement_error"]

    def test_all_tags(self):
        lm = NoiseLabelMap()
        lm.add(0, "hook_error", tags={"ancilla_qubit"})
        lm.add(1, "ghost_error", tags={"data_qubit"})

        assert set(lm.all_tags) == {"ancilla_qubit", "data_qubit"}

    def test_category_counts(self):
        lm = NoiseLabelMap()
        lm.add(0, "hook_error")
        lm.add(1, "hook_error")
        lm.add(2, "measurement_error")

        assert lm.category_counts() == {"hook_error": 2, "measurement_error": 1}

    def test_tag_counts(self):
        lm = NoiseLabelMap()
        lm.add(0, "hook_error", tags={"ancilla_qubit"})
        lm.add(1, "ghost_error", tags={"data_qubit"})
        lm.add(2, "measurement_error", tags={"ancilla_qubit"})

        assert lm.tag_counts() == {"ancilla_qubit": 2, "data_qubit": 1}

    def test_contains(self):
        lm = NoiseLabelMap()
        lm.add(0, "hook_error")
        assert 0 in lm
        assert 1 not in lm

    def test_iter(self):
        lm = NoiseLabelMap()
        lm.add(0, "hook_error", qubit=3)
        lm.add(1, "measurement_error", qubit=5)

        entries = list(lm)
        assert len(entries) == 2
        assert entries[0].noise_index == 0

    def test_serialization_roundtrip(self):
        lm = NoiseLabelMap()
        lm.add(0, "hook_error", qubit=3, tags={"ancilla_qubit"}, context={"round": 0})
        lm.add(1, "ghost_error", qubit=0, tags={"data_qubit"})

        json_str = lm.to_json()
        lm2 = NoiseLabelMap.from_json(json_str)

        assert len(lm2) == 2
        assert lm2.get_label(0) == "hook_error"
        assert lm2.get_tags(0) == {"ancilla_qubit"}
        assert lm2.get_entry(0).context == {"round": 0}
        assert lm2.get_label(1) == "ghost_error"
        assert lm2.get_tags(1) == {"data_qubit"}

    def test_repr(self):
        lm = NoiseLabelMap()
        lm.add(0, "hook_error")
        lm.add(1, "hook_error")
        lm.add(2, "measurement_error")

        r = repr(lm)
        assert "hook_error=2" in r
        assert "measurement_error=1" in r


class TestAutoLabelFourTypes:
    """Test auto-labeling with the four QStab IR error types + tags.

    Classification is POSITIONAL — depends on where the noise occurs
    relative to the CX schedule:

      - Data noise before first CX(data↔ancilla) → Type 0 (initial Ẽ)
      - Data noise during CX phase, with future CX → Type I (ghost)
      - Data noise during CX phase, no future CX → Type 0
      - Ancilla noise with future CX → Type II (hook)
      - Ancilla noise with no future CX → Type III (measurement)
    """

    def test_z_type_stabilizer_gadget(self):
        """Z-type stabilizer ZZ with ancilla=2.

        Gate list: R(2), n0(q0), n1(q2), CX(0,2), n2(q1), n3(q2), CX(1,2), n4(q2), M(2)

        Positional classification:
          n0: data q0, before first CX → Type 0
          n1: ancilla q2, has future CX(0,2) and CX(1,2) → Type II
          n2: data q1, after CX(0,2), has future CX(1,2) → Type I (ghost)
          n3: ancilla q2, has future CX(1,2) → Type II
          n4: ancilla q2, no future CX → Type III
        """
        circ = CliffordCircuit(3)
        circ.add_reset(2)
        circ._track_noise_location(0)   # n0
        circ._track_noise_location(2)   # n1
        circ.add_cnot(0, 2)
        circ._track_noise_location(1)   # n2
        circ._track_noise_location(2)   # n3
        circ.add_cnot(1, 2)
        circ._track_noise_location(2)   # n4
        circ.add_measurement(2)
        circ._qubit_num = 3

        lm = auto_label_from_circuit(circ, num_data_qubits=2)

        assert lm.get_label(0) == TYPE0  # data q0, before first CX
        assert lm.get_label(1) == TYPE2  # ancilla q2, future CX
        assert lm.get_label(2) == TYPE1  # data q1, ghost (future CX)
        assert lm.get_label(3) == TYPE2  # ancilla q2, future CX
        assert lm.get_label(4) == TYPE3  # ancilla q2, no future CX

        assert TAG_DATA in lm.get_tags(0)
        assert TAG_ANCILLA in lm.get_tags(1)
        assert TAG_DATA in lm.get_tags(2)
        assert TAG_ANCILLA in lm.get_tags(3)
        assert TAG_ANCILLA in lm.get_tags(4)

    def test_x_type_stabilizer_gadget(self):
        """X-type stabilizer XX with ancilla=2.

        R 2, H 2, CX 2 0, CX 2 1, H 2, M 2
        """
        circ = CliffordCircuit(3)
        circ.add_reset(2)
        circ._track_noise_location(2)   # n0: ancilla before H
        circ.add_hadamard(2)
        circ._track_noise_location(2)   # n1: ancilla before CX
        circ._track_noise_location(0)   # n2: data before first CX
        circ.add_cnot(2, 0)
        circ._track_noise_location(2)   # n3: ancilla between CX
        circ._track_noise_location(1)   # n4: data between CX gates
        circ.add_cnot(2, 1)
        circ._track_noise_location(2)   # n5: ancilla after last CX
        circ.add_hadamard(2)
        circ._track_noise_location(2)   # n6: ancilla before M
        circ.add_measurement(2)
        circ._qubit_num = 3

        lm = auto_label_from_circuit(circ, num_data_qubits=2)

        assert lm.get_label(0) == TYPE2  # ancilla, future CX
        assert lm.get_label(1) == TYPE2  # ancilla, future CX
        assert lm.get_label(2) == TYPE0  # data q0, before first CX
        assert lm.get_label(3) == TYPE2  # ancilla, future CX(2,1)
        assert lm.get_label(4) == TYPE1  # data q1, after CX(2,0), future CX(2,1)
        assert lm.get_label(5) == TYPE3  # ancilla, no future CX
        assert lm.get_label(6) == TYPE3  # ancilla, no future CX

        for i in [0, 1, 3, 5, 6]:
            assert TAG_ANCILLA in lm.get_tags(i)
        for i in [2, 4]:
            assert TAG_DATA in lm.get_tags(i)

    def test_data_qubit_idle_and_final_measurement(self):
        """Noise on data qubits with no CX → all Type 0."""
        circ = CliffordCircuit(3)
        circ._track_noise_location(0)   # n0: data before H
        circ.add_hadamard(0)
        circ._track_noise_location(1)   # n1: data before M
        circ.add_measurement(1)
        circ._qubit_num = 3

        lm = auto_label_from_circuit(circ, num_data_qubits=2)

        assert lm.get_label(0) == TYPE0
        assert TAG_DATA in lm.get_tags(0)
        assert lm.get_label(1) == TYPE0
        assert TAG_DATA in lm.get_tags(1)

    def test_depolarize2_creates_two_indices(self):
        """Two noise locations before CX: data → Type 0, ancilla → Type II."""
        circ = CliffordCircuit(3)
        circ.add_reset(2)
        circ._track_noise_location(0)   # n0: data
        circ._track_noise_location(2)   # n1: ancilla
        circ.add_cnot(0, 2)
        circ._qubit_num = 3

        lm = auto_label_from_circuit(circ, num_data_qubits=2)

        assert lm.total == 2
        assert lm.get_label(0) == TYPE0  # data before first CX
        assert lm.get_label(1) == TYPE2  # ancilla with future CX
        assert TAG_DATA in lm.get_tags(0)
        assert TAG_ANCILLA in lm.get_tags(1)

    def test_error_count_invariant(self):
        """Verify λ₀ + λ₁ + λ₂ + λ₃ = total noise sources."""
        circ = CliffordCircuit(3)
        circ.add_reset(2)
        circ._track_noise_location(0)
        circ._track_noise_location(2)
        circ.add_cnot(0, 2)
        circ._track_noise_location(1)
        circ._track_noise_location(2)
        circ.add_cnot(1, 2)
        circ._track_noise_location(2)
        circ.add_measurement(2)
        circ._track_noise_location(0)
        circ.add_hadamard(0)
        circ._qubit_num = 3

        lm = auto_label_from_circuit(circ, num_data_qubits=2)
        counts = lm.category_counts()

        total = sum(counts.values())
        assert total == circ._totalnoise

        all_indices = set()
        for cat in (TYPE0, TYPE1, TYPE2, TYPE3):
            for idx in lm.get_indices(cat):
                assert idx not in all_indices
                all_indices.add(idx)
        assert len(all_indices) == circ._totalnoise

    def test_every_noise_has_qubit_tag(self):
        """When num_data_qubits is given, every noise gets data/ancilla tag."""
        circ = CliffordCircuit(3)
        circ.add_reset(2)
        circ._track_noise_location(0)
        circ._track_noise_location(2)
        circ.add_cnot(0, 2)
        circ._track_noise_location(2)
        circ.add_measurement(2)
        circ._qubit_num = 3

        lm = auto_label_from_circuit(circ, num_data_qubits=2)

        for entry in lm:
            has_role = TAG_DATA in entry.tags or TAG_ANCILLA in entry.tags
            assert has_role, f"noise {entry.noise_index} missing qubit role tag"

    def test_stim_circuit_all_labeled(self):
        """Every noise in a compiled stim circuit gets a valid label."""
        stim_str = """R 0
R 1
R 2
H 0
CX 0 1
CX 0 2
H 0
M 0
M 1
M 2
DETECTOR rec[-1] rec[-2]
DETECTOR rec[-2] rec[-3]
OBSERVABLE_INCLUDE(0) rec[-1]"""

        circ = CliffordCircuit(3)
        circ.compile_from_stim_circuit_str(stim_str)

        lm = auto_label_from_circuit(circ, num_data_qubits=3)
        assert lm.total == circ._totalnoise
        for i in range(circ._totalnoise):
            label = lm.get_label(i)
            assert label in (TYPE0, TYPE1, TYPE2, TYPE3)


class TestAutoLabelNoisyCircuits:
    """Test auto-labeling on noisy circuits (noise after gates)."""

    def test_depolarize2_after_cx(self):
        """DEPOLARIZE2 after CX → classified by position in CX schedule."""
        stim_str = """R 0
R 1
R 2
CX 0 2
DEPOLARIZE2(0.01) 0 2
CX 1 2
DEPOLARIZE2(0.01) 1 2
M 2
DETECTOR rec[-1]
OBSERVABLE_INCLUDE(0) rec[-1]"""

        circ = CliffordCircuit(3)
        circ.compile_from_noisy_stim_circuit_str(stim_str)

        lm = auto_label_from_circuit(circ, num_data_qubits=2)

        # n0 (q0): after CX(0,2), no future CX → Type 0
        assert lm.get_label(0) == TYPE0
        # n1 (q2): after CX(0,2), has future CX(1,2) → Type II
        assert lm.get_label(1) == TYPE2
        # n2 (q1): after CX(1,2), no future CX → Type 0
        assert lm.get_label(2) == TYPE0
        # n3 (q2): after CX(1,2), no future CX → Type III
        assert lm.get_label(3) == TYPE3

    def test_x_error_before_measurement(self):
        """X_ERROR before M → ancilla with no future CX → Type III."""
        stim_str = """R 0
R 1
R 2
CX 0 2
CX 1 2
X_ERROR(0.01) 2
M 2
DETECTOR rec[-1]
OBSERVABLE_INCLUDE(0) rec[-1]"""

        circ = CliffordCircuit(3)
        circ.compile_from_noisy_stim_circuit_str(stim_str)

        lm = auto_label_from_circuit(circ, num_data_qubits=2)

        assert lm.get_label(0) == TYPE3
        assert TAG_ANCILLA in lm.get_tags(0)

    def test_depolarize2_and_x_error_same_qubit(self):
        """DEPOLARIZE2 after CX + X_ERROR before M on the same ancilla."""
        stim_str = """R 0
R 1
R 2
CX 0 2
DEPOLARIZE2(0.01) 0 2
X_ERROR(0.01) 2
M 2
DETECTOR rec[-1]
OBSERVABLE_INCLUDE(0) rec[-1]"""

        circ = CliffordCircuit(3)
        circ.compile_from_noisy_stim_circuit_str(stim_str)

        lm = auto_label_from_circuit(circ, num_data_qubits=2)

        assert lm.get_label(0) == TYPE0   # data q0, no future CX
        assert lm.get_label(1) == TYPE3   # ancilla q2, no future CX
        assert lm.get_label(2) == TYPE3   # X_ERROR on ancilla q2

    def test_depolarize1_after_h(self):
        """DEPOLARIZE1 after H on ancilla → with future CX → Type II."""
        stim_str = """R 0
R 1
R 2
H 2
DEPOLARIZE1(0.01) 2
CX 2 0
M 2
DETECTOR rec[-1]
OBSERVABLE_INCLUDE(0) rec[-1]"""

        circ = CliffordCircuit(3)
        circ.compile_from_noisy_stim_circuit_str(stim_str)

        lm = auto_label_from_circuit(circ, num_data_qubits=2)

        assert lm.get_label(0) == TYPE2
        assert TAG_ANCILLA in lm.get_tags(0)

    def test_full_noisy_repetition_code(self):
        """Full noisy d=3 repetition code with correct positional labels."""
        stim_str = """R 0
R 1
R 2
R 3
R 4
DEPOLARIZE1(0.001) 0
DEPOLARIZE1(0.001) 1
DEPOLARIZE1(0.001) 2
X_ERROR(0.001) 3
X_ERROR(0.001) 4
TICK
CX 0 3
DEPOLARIZE2(0.001) 0 3
CX 1 4
DEPOLARIZE2(0.001) 1 4
DEPOLARIZE1(0.001) 2
TICK
CX 1 3
DEPOLARIZE2(0.001) 1 3
CX 2 4
DEPOLARIZE2(0.001) 2 4
DEPOLARIZE1(0.001) 0
TICK
X_ERROR(0.001) 3
M 3
X_ERROR(0.001) 4
M 4
TICK
M 0
M 1
M 2
DETECTOR rec[-5] rec[-3] rec[-2]
DETECTOR rec[-4] rec[-2] rec[-1]
OBSERVABLE_INCLUDE(0) rec[-3]"""

        circ = CliffordCircuit(5)
        circ.compile_from_noisy_stim_circuit_str(stim_str)

        lm = auto_label_from_circuit(circ, num_data_qubits=3)

        assert lm.total == 17

        # Before first CX: all data noise → Type 0
        assert lm.get_label(0) == TYPE0   # q0, before CX
        assert lm.get_label(1) == TYPE0   # q1, before CX
        assert lm.get_label(2) == TYPE0   # q2, before CX

        # Ancilla after R, has future CX → Type II
        assert lm.get_label(3) == TYPE2   # q3
        assert lm.get_label(4) == TYPE2   # q4

        # DEPOLARIZE2 from CX(0,3): q0 done, q3 has future CX(1,3)
        assert lm.get_label(5) == TYPE0   # q0, no future CX
        assert lm.get_label(6) == TYPE2   # q3, has future CX(1,3)

        # DEPOLARIZE2 from CX(1,4): q1 has future CX(1,3)! → ghost
        assert lm.get_label(7) == TYPE1   # q1, future CX(1,3) → ghost!
        assert lm.get_label(8) == TYPE2   # q4, has future CX(2,4)

        # Idle q2: has future CX(2,4) → ghost
        assert lm.get_label(9) == TYPE1   # q2, future CX(2,4) → ghost!

        # DEPOLARIZE2 from CX(1,3): q1 done, q3 done
        assert lm.get_label(10) == TYPE0  # q1, no future CX
        assert lm.get_label(11) == TYPE3  # q3, no future CX → measurement

        # DEPOLARIZE2 from CX(2,4): q2 done, q4 done
        assert lm.get_label(12) == TYPE0  # q2, no future CX
        assert lm.get_label(13) == TYPE3  # q4, no future CX → measurement

        # Idle q0: no future CX → Type 0
        assert lm.get_label(14) == TYPE0  # q0, idle

        # X_ERROR before M → Type III
        assert lm.get_label(15) == TYPE3  # q3, measurement
        assert lm.get_label(16) == TYPE3  # q4, measurement

        # Verify counts
        counts = lm.category_counts()
        assert counts[TYPE0] == 7   # n0,1,2,5,10,12,14
        assert counts[TYPE1] == 2   # n7,9 (ghost)
        assert counts[TYPE2] == 4   # n3,4,6,8 (hook)
        assert counts[TYPE3] == 4   # n11,13,15,16 (measurement)

    def test_noisy_circuit_every_noise_tagged(self):
        """Every noise in a noisy circuit gets a data/ancilla tag."""
        stim_str = """R 0
R 1
R 2
CX 0 2
DEPOLARIZE2(0.01) 0 2
CX 1 2
DEPOLARIZE2(0.01) 1 2
X_ERROR(0.01) 2
M 2
DETECTOR rec[-1]
OBSERVABLE_INCLUDE(0) rec[-1]"""

        circ = CliffordCircuit(3)
        circ.compile_from_noisy_stim_circuit_str(stim_str)

        lm = auto_label_from_circuit(circ, num_data_qubits=2)

        for entry in lm:
            has_role = TAG_DATA in entry.tags or TAG_ANCILLA in entry.tags
            assert has_role, f"noise {entry.noise_index} missing qubit role tag"

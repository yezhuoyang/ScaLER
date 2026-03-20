"""Noise labeling system for circuit noise source categorization.

Each noise source (single-qubit Pauli channel) receives:

1. A **primary label** — exactly one of the four QStab IR error types:

   - ``data_qubit_error``  (Type 0, :math:`\\lambda_0`)
   - ``ghost_error``       (Type I, :math:`\\lambda_1`)
   - ``hook_error``        (Type II, :math:`\\lambda_2`)
   - ``measurement_error`` (Type III, :math:`\\lambda_3`)

2. **Tags** — additional descriptive labels such as ``"ancilla_qubit"``
   or ``"data_qubit"`` that are orthogonal to the primary type.  A noise
   source can have zero or more tags.

3. **Context** — metadata dict with keys like ``round``,
   ``stabilizer_index``, ``gate``, etc.

For DEPOLARIZE2 noise, each qubit gets its own noise index and is
classified independently (two noise indices per two-qubit noise event).

Typical usage::

    from scalerqec.Clifford.noiselabel import NoiseLabelMap

    label_map = NoiseLabelMap()
    label_map.add(0, "hook_error", qubit=3, tags={"ancilla_qubit"})
    label_map.add(1, "ghost_error", qubit=0, tags={"data_qubit"})
    label_map.add(2, "measurement_error", qubit=3, tags={"ancilla_qubit"})

    # Query by primary label
    hook_indices = label_map.get_indices("hook_error")

    # Query by tag
    ancilla_indices = label_map.get_indices_by_tag("ancilla_qubit")
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from typing import Any


@dataclass
class NoiseLabel:
    """A single labeled noise source.

    Attributes:
        noise_index: Index into the circuit's noise source list
            (``CliffordCircuit._index_to_noise``).  For DEPOLARIZE2,
            each qubit gets a separate index.
        label: Primary QStab IR error type (exactly one of
            ``data_qubit_error``, ``ghost_error``, ``hook_error``,
            ``measurement_error``).
        qubit: The qubit this noise source acts on.
        tags: Additional descriptive tags (e.g. ``"ancilla_qubit"``,
            ``"data_qubit"``).  Orthogonal to the primary label.
        context: Optional metadata such as round number, stabilizer
            index, or gate type.
    """

    noise_index: int
    label: str
    qubit: int = -1
    tags: set[str] = field(default_factory=set)
    context: dict[str, Any] = field(default_factory=dict)


class NoiseLabelMap:
    """Bidirectional mapping between noise source indices and labels.

    Maintains an ordered list of :class:`NoiseLabel` entries and
    provides efficient lookup by index, by primary label, or by tag.

    The map is constructed *after* circuit compilation and references
    noise sources by their integer index.  It adds zero overhead to
    the unlabeled simulation path.
    """

    def __init__(self) -> None:
        self._labels: list[NoiseLabel] = []
        self._index_to_label: dict[int, NoiseLabel] = {}
        self._category_to_indices: dict[str, list[int]] = {}
        self._tag_to_indices: dict[str, list[int]] = {}

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def add(
        self,
        noise_index: int,
        label: str,
        qubit: int = -1,
        tags: set[str] | None = None,
        context: dict[str, Any] | None = None,
    ) -> None:
        """Add a label for a noise source.

        If the noise index already has a label, it is replaced.

        Args:
            noise_index: Index of the noise source in the circuit.
            label: Primary QStab IR error type.
            qubit: Qubit this noise source acts on.
            tags: Additional descriptive tags.
            context: Optional metadata dict.
        """
        entry = NoiseLabel(
            noise_index=noise_index,
            label=label,
            qubit=qubit,
            tags=tags or set(),
            context=context or {},
        )

        # Remove old entry if exists
        if noise_index in self._index_to_label:
            old = self._index_to_label[noise_index]
            self._labels.remove(old)
            self._category_to_indices[old.label].remove(noise_index)
            if not self._category_to_indices[old.label]:
                del self._category_to_indices[old.label]
            for tag in old.tags:
                if tag in self._tag_to_indices:
                    self._tag_to_indices[tag].remove(noise_index)
                    if not self._tag_to_indices[tag]:
                        del self._tag_to_indices[tag]

        self._labels.append(entry)
        self._index_to_label[noise_index] = entry
        self._category_to_indices.setdefault(label, []).append(noise_index)
        for tag in entry.tags:
            self._tag_to_indices.setdefault(tag, []).append(noise_index)

    def add_tag(self, noise_index: int, tag: str) -> None:
        """Add an additional tag to an existing noise source.

        Args:
            noise_index: The noise source index.
            tag: The tag to add.

        Raises:
            KeyError: If *noise_index* is not in the map.
        """
        entry = self._index_to_label[noise_index]
        if tag not in entry.tags:
            entry.tags.add(tag)
            self._tag_to_indices.setdefault(tag, []).append(noise_index)

    # ------------------------------------------------------------------
    # Queries
    # ------------------------------------------------------------------

    def get_label(self, noise_index: int) -> str | None:
        """Return the primary label for a noise source, or ``None``."""
        entry = self._index_to_label.get(noise_index)
        return entry.label if entry else None

    def get_tags(self, noise_index: int) -> set[str]:
        """Return the tags for a noise source (empty set if not found)."""
        entry = self._index_to_label.get(noise_index)
        return entry.tags if entry else set()

    def get_entry(self, noise_index: int) -> NoiseLabel | None:
        """Return the full :class:`NoiseLabel` for a noise source."""
        return self._index_to_label.get(noise_index)

    def get_indices(self, label: str) -> list[int]:
        """Return all noise indices with the given primary label."""
        return list(self._category_to_indices.get(label, []))

    def get_indices_by_tag(self, tag: str) -> list[int]:
        """Return all noise indices that have the given tag."""
        return list(self._tag_to_indices.get(tag, []))

    @property
    def categories(self) -> list[str]:
        """All distinct primary label categories, in insertion order."""
        seen: dict[str, None] = {}
        for entry in self._labels:
            seen[entry.label] = None
        return list(seen)

    @property
    def all_tags(self) -> list[str]:
        """All distinct tags, in insertion order."""
        seen: dict[str, None] = {}
        for entry in self._labels:
            for tag in entry.tags:
                seen[tag] = None
        return list(seen)

    @property
    def total(self) -> int:
        """Number of labeled noise sources."""
        return len(self._labels)

    def category_counts(self) -> dict[str, int]:
        """Return ``{category: count}`` for all primary categories."""
        return {cat: len(idxs) for cat, idxs in self._category_to_indices.items()}

    def tag_counts(self) -> dict[str, int]:
        """Return ``{tag: count}`` for all tags."""
        return {tag: len(idxs) for tag, idxs in self._tag_to_indices.items()}

    def __len__(self) -> int:
        return len(self._labels)

    def __iter__(self):
        return iter(self._labels)

    def __contains__(self, noise_index: int) -> bool:
        return noise_index in self._index_to_label

    # ------------------------------------------------------------------
    # Serialization
    # ------------------------------------------------------------------

    def to_dict(self) -> list[dict[str, Any]]:
        """Serialize to a list of dicts."""
        return [
            {
                "noise_index": e.noise_index,
                "label": e.label,
                "qubit": e.qubit,
                "tags": sorted(e.tags),
                "context": e.context,
            }
            for e in self._labels
        ]

    @classmethod
    def from_dict(cls, data: list[dict[str, Any]]) -> NoiseLabelMap:
        """Deserialize from a list of dicts."""
        label_map = cls()
        for entry in data:
            label_map.add(
                noise_index=entry["noise_index"],
                label=entry["label"],
                qubit=entry.get("qubit", -1),
                tags=set(entry.get("tags", [])),
                context=entry.get("context", {}),
            )
        return label_map

    def to_json(self) -> str:
        """Serialize to JSON string."""
        return json.dumps(self.to_dict(), indent=2)

    @classmethod
    def from_json(cls, json_str: str) -> NoiseLabelMap:
        """Deserialize from JSON string."""
        return cls.from_dict(json.loads(json_str))

    def __repr__(self) -> str:
        counts = self.category_counts()
        parts = [f"{cat}={n}" for cat, n in counts.items()]
        return f"NoiseLabelMap({', '.join(parts)})"

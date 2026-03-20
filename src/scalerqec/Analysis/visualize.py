"""Circuit visualization with labeled noise sources.

Extends the base yquant LaTeX output from
:meth:`~scalerqec.Clifford.clifford.CliffordCircuit.get_yquant_latex`
to color-code noise boxes by their label category.

Typical usage::

    from scalerqec.Analysis.visualize import labeled_circuit_latex

    latex_str = labeled_circuit_latex(circuit, label_map)
"""

from __future__ import annotations

from typing import Any

from ..Clifford.clifford import (
    CliffordCircuit,
    Measurement,
    Reset,
    SingleQGate,
    TwoQGate,
    pauliNoise,
)
from ..Clifford.noiselabel import NoiseLabelMap


# Default color palette for the four QStab IR error types (LaTeX xcolor names)
DEFAULT_CATEGORY_COLORS: dict[str, str] = {
    "data_qubit_error": "orange!70",    # Type 0 (λ₀)
    "ghost_error": "green!60",          # Type I (λ₁)
    "hook_error": "red!80",             # Type II (λ₂)
    "measurement_error": "blue!60",     # Type III (λ₃)
}


def labeled_circuit_latex(
    circuit: CliffordCircuit,
    label_map: NoiseLabelMap,
    category_colors: dict[str, str] | None = None,
    show_legend: bool = True,
    show_noise_index: bool = False,
) -> str:
    """Generate yquant LaTeX with color-coded noise boxes by category.

    Args:
        circuit: A compiled CliffordCircuit.
        label_map: A NoiseLabelMap categorizing each noise source.
        category_colors: Optional mapping from category name to LaTeX
            color spec (e.g. ``"red!80"``).  Falls back to
            :data:`DEFAULT_CATEGORY_COLORS` for unmapped categories.
        show_legend: If ``True``, append a legend mapping colors to
            category names.
        show_noise_index: If ``True``, label noise boxes with their
            numeric index; otherwise label with the category abbreviation.

    Returns:
        A complete yquant LaTeX string.
    """
    colors = dict(DEFAULT_CATEGORY_COLORS)
    if category_colors:
        colors.update(category_colors)

    # Assign colors to any categories in the label_map not yet in the palette
    _extra_colors = [
        "teal!60", "magenta!50", "cyan!60", "lime!50",
        "pink!60", "olive!50", "violet!60",
    ]
    extra_idx = 0
    for cat in label_map.categories:
        if cat not in colors:
            colors[cat] = _extra_colors[extra_idx % len(_extra_colors)]
            extra_idx += 1

    # Build abbreviation map for category labels in boxes
    abbrevs = _build_abbreviations(label_map.categories)

    lines = []
    lines.append("\\begin{yquant}")
    lines.append("")

    lines.append("% -- Qubits and classical bits --")
    lines.append(
        "qubit {{$\\ket{{q_{{\\idx}}}}$}} q[{}];".format(circuit._qubit_num)
    )
    lines.append(
        "cbit {{$c_{{\\idx}} = 0$}} c[{}];".format(circuit._totalMeas)
    )
    lines.append("")
    lines.append("% -- Circuit Operations --")

    for gate in circuit._gatelists:
        if isinstance(gate, pauliNoise):
            entry = label_map.get_entry(gate._noiseindex)
            if entry is not None:
                color = colors.get(entry.label, "red!80")
                if show_noise_index:
                    box_label = str(gate._noiseindex)
                else:
                    box_label = abbrevs.get(entry.label, "?")
                lines.append(f"[fill={color}]")
                lines.append(
                    "box {{${}$}} q[{}];".format(box_label, gate._qubitindex)
                )
            else:
                # Unlabeled noise — use default red
                lines.append("[fill=red!80]")
                lines.append(
                    "box {{$n_{{{}}}$}} q[{}];".format(
                        gate._noiseindex, gate._qubitindex
                    )
                )

        elif isinstance(gate, TwoQGate):
            if gate._name == "CNOT":
                lines.append(
                    "cnot q[{}] | q[{}];".format(gate._target, gate._control)
                )
            elif gate._name == "CZ":
                lines.append(
                    "cz q[{}] | q[{}];".format(gate._target, gate._control)
                )

        elif isinstance(gate, SingleQGate):
            if gate._name == "H":
                lines.append("h q[{}];".format(gate._qubitindex))
            elif gate._name == "P":
                lines.append(
                    "box {{$S$}} q[{}];".format(gate._qubitindex)
                )
            elif gate._name in ("X", "Y", "Z"):
                lines.append(
                    "box {{${}$}} q[{}];".format(gate._name, gate._qubitindex)
                )

        elif isinstance(gate, Measurement):
            lines.append("measure q[{}];".format(gate._qubitindex))
            lines.append(
                "cnot c[{}] | q[{}];".format(
                    gate._measureindex, gate._qubitindex
                )
            )
            lines.append("discard q[{}];".format(gate._qubitindex))

        elif isinstance(gate, Reset):
            lines.append(
                "init {{$\\ket0$}} q[{}];".format(gate._qubitindex)
            )

    lines.append("")
    lines.append("\\end{yquant}")

    if show_legend:
        lines.append("")
        lines.append("% -- Legend --")
        lines.append("\\vspace{0.5em}")
        lines.append("\\noindent\\textbf{Noise categories:}\\\\")
        used_cats = label_map.categories
        for cat in used_cats:
            color = colors.get(cat, "red!80")
            abbrev = abbrevs.get(cat, "?")
            count = len(label_map.get_indices(cat))
            lines.append(
                f"\\colorbox{{{color}}}{{\\small ${abbrev}$}} "
                f"= {cat.replace('_', ' ')} ({count})\\quad"
            )

    return "\n".join(lines)


def noise_summary_table(label_map: NoiseLabelMap) -> str:
    """Generate a plain-text summary table of noise labels.

    Args:
        label_map: A populated NoiseLabelMap.

    Returns:
        A formatted string table.
    """
    counts = label_map.category_counts()
    total = label_map.total

    lines = []
    lines.append(f"{'Category':<25s} {'Count':>5s}  {'Fraction':>8s}")
    lines.append(f"{'-'*25} {'-'*5}  {'-'*8}")
    for cat, count in sorted(counts.items(), key=lambda x: -x[1]):
        frac = count / total if total > 0 else 0.0
        lines.append(f"{cat:<25s} {count:>5d}  {frac:>8.1%}")
    lines.append(f"{'-'*25} {'-'*5}  {'-'*8}")
    lines.append(f"{'TOTAL':<25s} {total:>5d}  {'100.0%':>8s}")

    return "\n".join(lines)


def _build_abbreviations(categories: list[str]) -> dict[str, str]:
    """Build short abbreviations for category names.

    Examples: ``hook_error`` → ``H``, ``measurement_error`` → ``M``,
    ``data_qubit_error`` → ``D``.
    """
    _known = {
        "data_qubit_error": "0",     # Type 0
        "ghost_error": "I",          # Type I
        "hook_error": "II",          # Type II
        "measurement_error": "III",  # Type III
    }

    result = {}
    used_abbrevs: set[str] = set()

    for cat in categories:
        if cat in _known and _known[cat] not in used_abbrevs:
            result[cat] = _known[cat]
            used_abbrevs.add(_known[cat])
        else:
            # Use first letter of each word
            parts = cat.replace("_", " ").split()
            abbrev = "".join(p[0].upper() for p in parts if p)
            # Ensure uniqueness
            base = abbrev
            counter = 1
            while abbrev in used_abbrevs:
                abbrev = f"{base}{counter}"
                counter += 1
            result[cat] = abbrev
            used_abbrevs.add(abbrev)

    return result

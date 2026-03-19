"""Exact symbolic LER computation via dynamic programming.

This module implements the core algorithm of ScaLERQEC: computing the exact
logical error rate (LER) as a symbolic polynomial in ``p`` (the physical
error probability) using dynamic programming over noise sources.

Algorithm overview:

1. Parse a STIM circuit and build the Quantum Error Propagation Graph (QEPG).
2. Use PyMatching to decode every possible detector outcome, building a
   lookup table of decoder predictions.
3. Identify which detector/observable outcomes correspond to logical errors
   (i.e., where the true observable disagrees with the decoder prediction).
4. Fill a 3-D DP table ``dp[i][j][v]`` where:
   - ``i`` indexes noise sources (0 .. num_noise),
   - ``j`` counts how many of those sources actually fired (weight),
   - ``v`` encodes the combined detector + observable outcome as a bitmask.
5. Sum the DP entries at outcomes that cause logical errors to obtain
   the LER polynomial.

Module-level constants:
    MAX_degree (int): Maximum polynomial degree retained in series
        expansions (default 200).
    MAX_weight (int): Maximum error weight tracked in the DP table
        (default 32).  Higher weights contribute negligibly for small ``p``.

Module-level symbols (SymPy):
    p: Symbolic physical error probability.
    q: ``1 - p``, the no-error probability.
    Px, Py, Pz: Individual Pauli error probabilities under the
        depolarizing model, each equal to ``p / 3``.
"""

from sympy import symbols, binomial, simplify, latex
from ..Clifford.clifford import *
from ..Clifford.stimparser import *
from ..Clifford.QEPGpython import *
import pymatching
from ..QEC.noisemodel import NoiseModel
from ..QEC.qeccircuit import StabCode
from tqdm.auto import tqdm

# ----------------------------------------------------------------------
# Physical-error model
p = symbols("p")
q = 1 - p  #   probability of "no error"
Px = Py = Pz = p / 3  #   probabilities of  X  Y  Z


def vec_to_idx(vec):
    """Convert a binary vector (MSB first) to its integer index.

    Interprets the vector as a big-endian binary number.

    Example::

        >>> vec_to_idx((1, 0, 1))
        5  # 1*4 + 0*2 + 1*1

    Args:
        vec (tuple[int, ...] | list[int]): Binary digits (0 or 1),
            with the most significant bit at position 0.

    Returns:
        int: Non-negative integer representation.  Returns 0 for an
        empty vector.
    """
    idx = 0
    for bit in vec:
        idx = (idx << 1) | bit
    return idx


def idx_to_vec(idx, dimension):
    """Convert an integer index to a binary vector of fixed length.

    Inverse of :func:`vec_to_idx`.  The first element of the returned
    tuple is the most significant bit (MSB).

    Example::

        >>> idx_to_vec(5, 3)
        (1, 0, 1)

    Args:
        idx (int): Non-negative integer to convert.
        dimension (int): Length of the output binary vector.  If
            ``dimension`` exceeds the number of significant bits,
            the result is zero-padded on the left.

    Returns:
        tuple[int, ...]: Binary vector of length ``dimension``.
    """
    bits = []
    temp_idx = idx
    for _ in range(dimension):
        bits.append(temp_idx & 1)
        temp_idx >>= 1
    return tuple(reversed(bits))


def idx_to_bool_list(idx, dimension):
    """Convert an integer index to a boolean list of fixed length.

    Similar to :func:`idx_to_vec` but returns booleans instead of
    integers.  Used as input to PyMatching's batch decoder.

    Example::

        >>> idx_to_bool_list(5, 3)
        [True, False, True]

    Args:
        idx (int): Non-negative integer to convert.
        dimension (int): Length of the output boolean list.

    Returns:
        list[bool]: Boolean list of length ``dimension``, MSB first.
    """
    bool_list = []
    temp_idx = idx
    for _ in range(dimension):
        bool_list.append(True if (temp_idx & 1) else False)
        temp_idx >>= 1
    return list(reversed(bool_list))


def xor_vec(vec_a, vec_b):
    """Perform element-wise XOR (GF(2) addition) on two binary vectors.

    Both vectors must have the same length.

    Args:
        vec_a (tuple[int, ...] | list[int]): First binary vector.
        vec_b (tuple[int, ...] | list[int]): Second binary vector.

    Returns:
        tuple[int, ...]: Element-wise XOR of ``vec_a`` and ``vec_b``.
    """

    return tuple(a ^ b for a, b in zip(vec_a, vec_b))


MAX_degree = 200
"""int: Maximum polynomial degree retained in series expansions."""

MAX_weight = 32
"""int: Maximum error weight tracked in the DP table.

Weights above this threshold are ignored; their contribution is
negligible for realistic physical error rates.
"""


class SymbolicLERcalc:
    """Compute the exact logical error rate as a polynomial in ``p``.

    This calculator builds a 3-D dynamic-programming table over all noise
    sources in a STIM circuit to obtain the exact LER polynomial.  The
    polynomial can then be evaluated at any physical error rate ``p`` to
    yield an exact (not sampled) logical error rate.

    The typical workflow is:

    1. Instantiate with an optional default error rate.
    2. Call :meth:`calculate_LER_from_file` or
       :meth:`calculate_LER_from_StabCode` for the full pipeline, or
       invoke the individual steps manually for finer control.

    Attributes:
        _num_detector (int): Number of detectors in the circuit.
        _num_noise (int): Number of independent noise sources.
        _error_rate (float): Physical error rate used during circuit
            compilation.
        _dp (list): 3-D DP table of shape
            ``[num_noise+1][num_noise+1][2^(num_detectors+1)]``.
            Each entry is a SymPy expression in ``p``.
        _cliffordcircuit (CliffordCircuit): The compiled Clifford circuit.
        _graph (QEPGpython): Quantum Error Propagation Graph built from
            the circuit.
        _all_predictions (numpy.ndarray): PyMatching decoder predictions
            for every possible detector outcome.
        _PROP_X (list[tuple]): Propagation vectors for X errors at each
            noise source.
        _PROP_Y (list[tuple]): Propagation vectors for Y errors.
        _PROP_Z (list[tuple]): Propagation vectors for Z errors.
        _error_row_indices (list[int]): Bitmask indices of
            detector/observable outcomes that cause logical errors.
        _subspace_LER (dict[int, sympy.Expr]): LER contribution from
            each error weight subspace.
    """

    def __init__(self, error_rate=0):
        """Initialize the symbolic LER calculator.

        Args:
            error_rate (float): Default physical error rate for
                depolarizing noise injection.  Defaults to 0.
        """
        self._num_detector = 0
        self._num_noise = 0
        self._error_rate = error_rate
        self._dp = None
        self._cliffordcircuit = CliffordCircuit(4)
        self._graph = None
        self._all_predictions = None

        self._PROP_X = None
        self._PROP_Y = None
        self._PROP_Z = None

        self._error_row_indices = []
        self._subspace_LER = {}

    def get_totalnoise(self):
        """Return the total number of independent noise sources in the circuit.

        Returns:
            int: Number of noise sources (set after parsing a circuit).
        """
        return self._num_noise

    def parse_from_file(self, filepath):
        """Load a STIM circuit from a file, compile it, and build the QEPG.

        This method reads the STIM circuit string, compiles it into a
        :class:`CliffordCircuit` with depolarizing noise injected at
        rate :attr:`_error_rate`, counts noise sources and detectors,
        and constructs the backward Quantum Error Propagation Graph.

        Args:
            filepath (str): Path to a file containing a STIM circuit
                description.
        """
        stim_str = ""
        with open(filepath, "r", encoding="utf-8") as f:
            stim_str = f.read()

        self._cliffordcircuit.compile_from_stim_circuit_str(stim_str)
        self._num_noise = self._cliffordcircuit.totalnoise
        self._num_detector = len(self._cliffordcircuit.parityMatchGroup)

        self._total_detector_outcome = 1 << (self._num_detector + 1)

        self._graph = QEPGpython(self._cliffordcircuit)
        self._graph.backword_graph_construction()

    def generate_pymatching_table(self):
        """Decode every possible detector outcome using PyMatching.

        Iterates over all :math:`2^{\\text{num\\_detectors}}` possible
        detector bit-strings, runs each through the PyMatching MWPM
        decoder, and stores the predictions in :attr:`_all_predictions`.

        This table is later used by :meth:`calc_error_row_indices` to
        determine which outcomes lead to logical errors.
        """
        # Configure a decoder using the noisy circuit.
        from ..QEC.noisemodel import SIDNoiseModel

        noisy_stim = SIDNoiseModel(self._error_rate).inject_noise(
            self._cliffordcircuit.stimcircuit
        )
        detector_error_model = noisy_stim.detector_error_model(decompose_errors=False)
        matcher = pymatching.Matching.from_detector_error_model(detector_error_model)

        all_inputs = []

        # print("Total detector outcome: ", 1<<self._num_detector)
        for i in range(0, 1 << self._num_detector):
            # print("i=",i)
            # print(1<<self._num_detector)
            # Convert the integer to a list of booleans
            bool_list = idx_to_bool_list(i, self._num_detector)
            # Print the list of booleans
            all_inputs.append(bool_list)
        # print(all_inputs)
        self._all_predictions = matcher.decode_batch(all_inputs)

    def calc_error_row_indices(self):
        """Identify detector/observable outcomes that cause logical errors.

        For each of the :math:`2^{\\text{num\\_detectors}+1}` possible
        combined (detector, observable) outcomes, compares the true
        observable bit with the decoder's predicted observable.  Outcomes
        where they disagree are logical errors; their bitmask indices
        are stored in :attr:`_error_row_indices`.

        Must be called after :meth:`generate_pymatching_table`.
        """
        self._error_row_indices = []
        for row in range(0, self._total_detector_outcome):
            full_bool_list = idx_to_bool_list(row, self._num_detector + 1)
            """
            Get the current observable outcome represented by 
            the row index.
            """
            current_obs = full_bool_list[-1]
            ind_no_obs = vec_to_idx(full_bool_list[:-1])
            """
            Return the exact prediction result from the decoder, which
            is computed by function generate_pymatching_table(self)
            """
            real_obs = self._all_predictions[ind_no_obs][0]
            if current_obs != real_obs:
                self._error_row_indices.append(row)

        # print(self._error_row_indices)

    def initialize_single_pauli_propagation(self):
        """Precompute error propagation vectors for every noise source.

        For each noise source index, computes the detector/observable
        outcome vector that results from a single X, Y, or Z error at
        that location propagating through the QEPG.  Results are stored
        in :attr:`_PROP_X`, :attr:`_PROP_Y`, and :attr:`_PROP_Z` as
        tuples of binary digits.
        """
        self._PROP_X = []
        self._PROP_Y = []
        self._PROP_Z = []
        for noiseidx in range(self._num_noise):
            self._PROP_X.append(tuple(self._graph.sample_x_error(noiseidx)))
            self._PROP_Y.append(tuple(self._graph.sample_y_error(noiseidx)))
            self._PROP_Z.append(tuple(self._graph.sample_z_error(noiseidx)))

    def initialize_dp(self):
        """Allocate and zero-initialize the 3-D DP table.

        Creates a table of shape
        ``[num_noise+1][num_noise+1][2^(num_detectors+1)]``
        with all entries set to 0.  Must be called before
        :meth:`dynamic_calculation_of_dp`.
        """
        self._dp = [
            [[0] * self._total_detector_outcome for _ in range(self._num_noise + 1)]
            for _ in range(self._num_noise + 1)
        ]

    def verify_table(self, i):
        """Verify that DP row ``i`` is a valid probability distribution.

        Asserts that the sum of ``dp[i][j][v]`` over all weights ``j``
        and all outcome vectors ``v`` equals 1 (after symbolic
        simplification).  Useful as a sanity check during development.

        Args:
            i (int): The noise-source index (row) to verify.

        Raises:
            AssertionError: If the probabilities do not sum to 1.
        """
        sum = 0
        for j in range(0, i + 1):
            for vec_index in range(self._total_detector_outcome):
                sum += self._dp[i][j][vec_index]
            # print(self._dp[i][j][vec_index])
        sum = simplify(sum)
        # print(i,sum)
        assert sum == 1

    def dynamic_calculation_of_dp(self):
        """Fill the DP table using the recurrence relation.

        Implements the core recurrence:

        .. math::

            dp[i][j][v] = q \\cdot dp[i-1][j][v]
                + P_X \\cdot dp[i-1][j-1][v \\oplus \\text{prop}_X(i)]
                + P_Y \\cdot dp[i-1][j-1][v \\oplus \\text{prop}_Y(i)]
                + P_Z \\cdot dp[i-1][j-1][v \\oplus \\text{prop}_Z(i)]

        where ``dp[i][j][v]`` is the probability that the first ``i``
        noise sources produce exactly ``j`` errors and the combined
        detector/observable outcome is ``v``.

        Uses bitmask XOR for propagation (instead of tuple-based
        ``xor_vec``) for performance.  Weights above
        :data:`MAX_weight` are skipped.

        Must be called after :meth:`initialize_single_pauli_propagation`.
        Calls :meth:`initialize_dp` internally.
        """
        MAX_I = self._num_noise

        self.initialize_dp()
        col_size = self._num_detector + 1

        # Precompute bit masks for X,Y,Z propagation if not already done
        # (do this once in __init__ ideally)
        if not hasattr(self, "_PROP_X_mask"):
            # First time: build from PROP_* vectors
            self._PROP_X_mask = [int(vec_to_idx(v)) for v in self._PROP_X]
            self._PROP_Y_mask = [int(vec_to_idx(v)) for v in self._PROP_Y]
            self._PROP_Z_mask = [int(vec_to_idx(v)) for v in self._PROP_Z]
        else:
            # If these were created earlier as numpy scalars, normalize them
            self._PROP_X_mask = [int(m) for m in self._PROP_X_mask]
            self._PROP_Y_mask = [int(m) for m in self._PROP_Y_mask]
            self._PROP_Z_mask = [int(m) for m in self._PROP_Z_mask]

        dp = self._dp
        total = self._total_detector_outcome

        PROP_X_mask = self._PROP_X_mask
        PROP_Y_mask = self._PROP_Y_mask
        PROP_Z_mask = self._PROP_Z_mask

        # (1-p)**i computed incrementally
        no_error_pow = 1  # (1-p)**0
        for i in tqdm(range(0, MAX_I + 1), desc="Dynamic DP rows", unit="row"):
            dp[i][0][0] = no_error_pow
            no_error_pow *= 1 - p

            j_max = min(i, MAX_weight)
            prev_i = i - 1

            if i == 0:
                continue

            mask_X = PROP_X_mask[prev_i]
            mask_Y = PROP_Y_mask[prev_i]
            mask_Z = PROP_Z_mask[prev_i]

            for j in range(1, j_max + 1):
                prev_j = j - 1
                dp_prev_same = dp[prev_i][j]
                dp_prev_less = dp[prev_i][prev_j]
                dp_curr = dp[i][j]

                for vec_idx in range(total):
                    # 1) no-error branch
                    acc = q * dp_prev_same[vec_idx]

                    # 2) error branches (only if j>=1, but loop starts at 1)
                    prev_idx_X = vec_idx ^ mask_X
                    prev_idx_Y = vec_idx ^ mask_Y
                    prev_idx_Z = vec_idx ^ mask_Z

                    acc += Px * dp_prev_less[prev_idx_X]
                    acc += Py * dp_prev_less[prev_idx_Y]
                    acc += Pz * dp_prev_less[prev_idx_Z]

                    # No simplify here; keep as-is
                    dp_curr[vec_idx] = acc

    # def dynamic_calculation_of_dp(self):
    #     # ----------------------------------------------------------------------
    #     # DP tables
    #     MAX_I = self._num_noise
    #     self.initialize_dp()

    #     col_size=self._num_detector+1
    #     # ----------------------------------------------------------------------
    #     # Fill   dp[i][j][·]   using the recurrence in Eq. (1)
    #     for i in tqdm(range(0, MAX_I+1),desc="Dynamic DP rows", unit = "row"):
    #         self._dp[i][0][0] = (1-p)**i

    #         for j in range(1, i+1):           # j ≤ i
    #             if j > MAX_weight:
    #                 continue

    #             #print("MAX_I=",MAX_I,"i=",i,"j=",j)
    #             for vec_idx in range(self._total_detector_outcome):

    #                 vec = idx_to_vec(vec_idx,col_size)

    #                 # 1) “no error’’ branch
    #                 acc = q * self._dp[i-1][j][vec_idx]

    #                 if j >= 1:
    #                     for (prob, prop) in ((Px, self._PROP_X[i-1]),
    #                                         (Py, self._PROP_Y[i-1]),
    #                                         (Pz, self._PROP_Z[i-1])):

    #                         prev_vec = xor_vec(prop, vec)
    #                         acc += prob * self._dp[i-1][j-1][ vec_to_idx(prev_vec) ]

    #                 self._dp[i][j][vec_idx] = simplify(acc)
    #         # if MAX_weight>=self._num_noise:
    #         #     self.verify_table(i)

    def calculate_LER_brute_force(self):
        """Compute the exact LER by brute-force enumeration.

        Enumerates all :math:`4^N` possible Pauli error patterns
        (I, X, Y, Z at each noise source) and sums the probability
        of patterns that cause a logical error.

        Note:
            Not yet implemented.  Use :meth:`dynamic_calculation_of_dp`
            followed by :meth:`calculate_LER` instead.
        """
        pass

    def calculate_LER(self):
        """Sum DP entries at error rows to obtain the LER polynomial.

        Iterates over all error weights (1 .. num_noise) and sums the
        DP entries ``dp[num_noise][weight][row]`` for each row in
        :attr:`_error_row_indices`.  The per-weight contributions are
        stored in :attr:`_subspace_LER` and the total LER polynomial
        is stored in :attr:`_ler`.

        Must be called after :meth:`dynamic_calculation_of_dp` and
        :meth:`calc_error_row_indices`.

        Returns:
            sympy.Expr: The LER as a symbolic polynomial in ``p``,
            simplified and expanded.
        """
        self._ler = 0
        for weight in range(1, self._num_noise + 1):
            subLER = 0
            for rowindex in self._error_row_indices:
                subLER += self._dp[self._num_noise][weight][rowindex]

            self._subspace_LER[weight] = simplify(subLER)
            self._ler += simplify(subLER)
        self._ler = simplify(self._ler).expand()
        # LER=LER.series(p, 0, MAX_degree).removeO()    # no .expand()

        print("LER polynomial: ", latex(self._ler))
        return self._ler

    def evaluate_LER(self, pval):
        """Evaluate the LER polynomial at a specific physical error rate.

        Args:
            pval (float): The physical error rate to substitute for ``p``.

        Returns:
            float: The numerical logical error rate at the given ``p``.
        """
        return self._ler.evalf(subs={p: pval})

    def evaluate_LER_subspace(self, pval, weight):
        """Evaluate the LER contribution from a single weight subspace.

        Args:
            pval (float): The physical error rate to substitute for ``p``.
            weight (int): The error weight subspace to evaluate.

        Returns:
            float: The numerical LER contribution from weight-``weight``
            error patterns at the given ``p``.
        """
        return self._subspace_LER[weight].evalf(subs={p: pval})

    def subspace_LER(self, weight):
        """Compute the normalized subspace LER for a given error weight.

        Divides the raw subspace LER (probability contributed by
        weight-``weight`` error patterns) by the corresponding binomial
        coefficient :math:`\\binom{N}{w} p^w (1-p)^{N-w}`, yielding
        the conditional logical error probability given that exactly
        ``weight`` noise sources fired.

        Args:
            weight (int): The error weight subspace.

        Returns:
            sympy.Expr: The normalized (conditional) LER polynomial
            for the given weight, simplified and expanded.

        Raises:
            ValueError: If :meth:`calculate_LER` has not been called
                or the given weight was not computed.
        """
        if weight in self._subspace_LER:
            bernolli_coeff = (
                binomial(self._num_noise, weight)
                * (p**weight)
                * ((1 - p) ** (self._num_noise - weight))
            )
            subspaceLER = simplify(self._subspace_LER[weight] / bernolli_coeff)
            return subspaceLER.expand()
        else:
            raise ValueError(
                "Subspace LER for weight {} is not calculated.".format(weight)
            )

    def calculate_LER_from_file(self, filepath, pvalue) -> float:
        """Full pipeline: load a STIM circuit file and compute the exact LER.

        This convenience method runs the entire symbolic LER computation
        from start to finish:

        1. Parse the circuit and inject depolarizing noise at rate ``pvalue``.
        2. Build the PyMatching decoder prediction table for all detector
           outcomes.
        3. Construct the QEPG and precompute error propagation vectors.
        4. Identify detector/observable outcomes that cause logical errors.
        5. Fill the DP table via :meth:`dynamic_calculation_of_dp`.
        6. Sum probabilities at error rows and evaluate at ``pvalue``.

        Args:
            filepath (str): Path to the STIM circuit file.
            pvalue (float): Physical error rate at which to evaluate the
                LER polynomial.

        Returns:
            float: The exact logical error rate at the given physical
            error rate.
        """
        self._error_rate = pvalue
        self.parse_from_file(filepath)
        print("---Step2: Generate the prediction table---")
        self.generate_pymatching_table()
        print("---Step2: construction QEPG--------------")
        self.initialize_single_pauli_propagation()
        print("---Step3: calculating error indices--------------")
        self.calc_error_row_indices()
        print("---Step4: dynamic algorithm--------------")
        self.dynamic_calculation_of_dp()
        self.calculate_LER()
        return self.evaluate_LER(pvalue)

    def calculate_LER_from_StabCode(
        self, qeccirc: StabCode, noise_model: NoiseModel
    ) -> float:
        """Full pipeline: compute the exact LER from a StabCode object.

        Similar to :meth:`calculate_LER_from_file` but accepts a
        :class:`~scalerqec.QEC.qeccircuit.StabCode` and
        :class:`~scalerqec.QEC.noisemodel.NoiseModel` directly,
        avoiding the need for a circuit file on disk.

        The method constructs the IR standard scheme, compiles the STIM
        circuit, injects noise via the noise model, then runs the full
        DP pipeline.

        Args:
            qeccirc (StabCode): A stabilizer code object describing the
                quantum error-correcting code.
            noise_model (NoiseModel): The noise model specifying the
                physical error rate and noise structure.

        Returns:
            float: The exact logical error rate evaluated at the noise
            model's error rate.
        """
        self._error_rate = noise_model.error_rate
        qeccirc.construct_IR_standard_scheme()
        qeccirc.compile_stim_circuit_from_IR_standard()
        self._cliffordcircuit = noise_model.reconstruct_clifford_circuit(
            qeccirc.circuit
        )

        self._num_noise = self._cliffordcircuit.totalnoise
        self._num_detector = len(self._cliffordcircuit.parityMatchGroup)
        self._total_detector_outcome = 1 << (self._num_detector + 1)
        self._graph = QEPGpython(self._cliffordcircuit)
        self._graph.backword_graph_construction()

        print("---Step1: Generate the prediction table---")
        self.generate_pymatching_table()
        print("---Step2: construction QEPG--------------")
        self.initialize_single_pauli_propagation()
        print("---Step3: calculating error indices--------------")
        self.calc_error_row_indices()
        print("---Step4: dynamic algorithm--------------")
        self.dynamic_calculation_of_dp()
        self.calculate_LER()
        self._ler = self.evaluate_LER(self._error_rate)
        print("Evaluated LER at p={} is {}".format(self._error_rate, self._ler))
        return self._ler


if __name__ == "__main__":
    tmp = SymbolicLERcalc(0.001)
    filepath = "C:/Users/username/Documents/Sampling/stimprograms/small/simple"
    print(tmp.calculate_LER_from_file(filepath, 0.001))

    num_noise = tmp._num_noise

    # for weight in range(1,11):
    #     print("LER in the subspace {} is {}".format(weight,tmp.evaluate_LER_subspace(0.001,weight)))

    for weight in range(1, 12):
        print("SubspaceLER {} is {}".format(weight, tmp.subspace_LER(weight)))

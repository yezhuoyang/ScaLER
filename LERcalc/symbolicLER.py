import sys
from sympy import symbols, binomial, Rational, simplify, latex
from LERcalc.clifford import *
from LERcalc.stimparser import *
import pymatching

# ----------------------------------------------------------------------
# Physical-error model
p = symbols('p')
q = 1 - p                     #   probability of "no error"
Px = Py = Pz = p / 3          #   probabilities of  X  Y  Z



def vec_to_idx(vec):
    """
    Converts a binary vector (A tuple or a list of 0s and 1s) to an integer index.
    The first element of the vector is treated as the most significant bit(MSB).
    Example: (1,0,1) -> 1*2^2 + 0*2^1 + 1*2^0 =5
    
    Args:
        vec: A tuple or list of binary digits(0 or 1).

    Returns:
        An integer representation of the binary vector.
        Returns 0 for an empty vector.
    """
    idx=0
    for bit in vec:
        idx = (idx<<1) | bit
    return idx


def idx_to_vec(idx, dimension):
    """
    Converts and integer index back to the binary vector of a given dimension.
    The first element of the vector is the most significant bit (MSB).
    Example: (5,3) - > (1,0,1)

    Args:
        idx: The non-negative integer index.
        dimension: The desired non-negative dimension(length) of the output binary vector

    Returns:
        A tuple of binary digits representing the index.
    """
    bits = []
    temp_idx = idx
    for _ in range(dimension):
        bits.append(temp_idx & 1)
        temp_idx >>= 1
    return tuple(reversed(bits))



def idx_to_bool_list(idx, dimension):
    """
    Converts and integer index back to the boolean vector of a given dimension.
    The first element of the vector is the most significant bit (MSB).
    Example: (5,3) - > [True,False,True]

    Args:
        idx: The non-negative integer index.
        dimension: The desired non-negative dimension(length) of the output binary vector

    Returns:
        A list of boolean representing the binary digits given the index.
    """
    bool_list = []
    temp_idx = idx
    for _ in range(dimension):
        bool_list.append(True if (temp_idx & 1) else False)
        temp_idx >>= 1
    return list(reversed(bool_list))


def xor_vec(vec_a, vec_b):
    """
    Performs element-wise XOR addition on two binary vectors.

    Args:
        vec_a: The first binary vector (tuple or list of 0s and 1s).
        vec_b: The second binary vector (tuple or list of 0s and 1s).

    Returns:
        A tuple representing the element-wise XOR of two vectors.
    """

    return tuple(a^b for a,b in zip(vec_a,vec_b))


MAX_degree=100


'''
Use symbolic algorithm to calculate the probability.
Simply enumerate all possible cases
'''
class symbolicLER:
    def __init__(self,error_rate=0):
        self._num_detector=0
        self._num_noise=0
        self._error_rate=error_rate
        self._dp=None
        self._cliffordcircuit=CliffordCircuit(4)  
        self._all_predictions=None

        self._PROP_X = None
        self._PROP_Y = None
        self._PROP_Z = None


    def parse_from_file(self,filepath):
        """
        Read the circuit, parse from the file
        """
        stim_str=""
        with open(filepath, "r", encoding="utf-8") as f:
            stim_str = f.read()
        

        self._cliffordcircuit.set_error_rate(self._error_rate)  
        self._cliffordcircuit.compile_from_stim_circuit_str(stim_str)
        self._num_noise = self._cliffordcircuit.get_totalnoise()
        self._num_detector=len(self._cliffordcircuit.get_parityMatchGroup())




    def generate_pymatching_table(self):
        """
        For all detector result, generate the prediction through pymatching.

        _all_predictions store all prediction by matching
        """
        # Configure a decoder using the circuit.
        stimcircuit=self._cliffordcircuit.get_stim_circuit()
        detector_error_model = stimcircuit.detector_error_model(decompose_errors=False)
        matcher = pymatching.Matching.from_detector_error_model(detector_error_model)

        all_inputs = []
        for i in range(0,1<<self._num_detector):
            # Convert the integer to a list of booleans
            bool_list = idx_to_bool_list(i, self._num_detector)
            # Print the list of booleans
            all_inputs.append(bool_list)
        print(all_inputs)
        self._all_predictions = matcher.decode_batch(all_inputs)




    def initialize_single_pauli_propagation(self):
        """
        Calculate and store the table of the propagation result of single pauli error
        """
        self._PROP_X = []
        self._PROP_Y = []
        self._PROP_Z = []
        for noiseidx in range(self.self._num_noise):
            pass



    def initialize_dp(self):
        """
        Given the circuit information, initialize the dp table for running the algorithm
        """
        pass





if __name__=="__main__":


    #print(idx_to_bool_list(5, 3))
    tmp=symbolicLER(0.01)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/1cnot"
    tmp.parse_from_file(filepath)

    tmp.generate_pymatching_table()

    print(tmp._all_predictions)

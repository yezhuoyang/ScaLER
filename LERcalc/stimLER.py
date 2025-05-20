from LERcalc.clifford import *
import pymatching
from LERcalc.stimparser import *


def count_logical_errors(circuit: stim.Circuit, num_shots: int) -> int:
    # Sample the circuit.
    sampler = circuit.compile_detector_sampler()
    detection_events, observable_flips = sampler.sample(num_shots, separate_observables=True)

    # Configure a decoder using the circuit.
    detector_error_model = circuit.detector_error_model(decompose_errors=True)
    matcher = pymatching.Matching.from_detector_error_model(detector_error_model)

    # Run the decoder.
    predictions = matcher.decode_batch(detection_events)

    # Count the mistakes.
    num_errors = 0
    for shot in range(num_shots):
        actual_for_shot = observable_flips[shot]
        predicted_for_shot = predictions[shot]
        if not np.array_equal(actual_for_shot, predicted_for_shot):
            num_errors += 1
    return num_errors


'''
Use stim and Monte Calo sampling method to estimate the logical error rate
'''
class stimLERcalc:
    def __init__(self):
        pass



    def calculate_LER_from_file(self,shots,filepath,pvalue):
        circuit=CliffordCircuit(2)
        circuit.set_error_rate(pvalue)

        stim_str=""
        with open(filepath, "r", encoding="utf-8") as f:
            stim_str = f.read()

        stim_circuit=rewrite_stim_code(stim_str)
        circuit.set_stim_str(stim_circuit)
        circuit.compile_from_stim_circuit_str(stim_circuit)           
        new_stim_circuit=circuit.get_stim_circuit()      

        num_ler=count_logical_errors(new_stim_circuit,shots)

        return num_ler/shots
    

if __name__ == "__main__":
    calculator=stimLERcalc()
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/surface/surface30"
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/small/1cnoth"
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/small/simpleh"
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/small/2cnot2R"
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/hexagon/hexagon3"
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/hexagon/hexagon3"
    ler=calculator.calculate_LER_from_file(1000000,filepath,0.001)

    print(ler)
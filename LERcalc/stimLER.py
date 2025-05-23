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




MIN_NUM_LE_EVENT = 100
SAMPLE_GAP = 10000000

'''
Use stim and Monte Calo sampling method to estimate the logical error rate
The sampler will finally decide how many samples to used.
Shot is the initial guess of how many samples to used.
We also need to estimate the uncertainty of the LER.
'''
class stimLERcalc:
    def __init__(self):
        self._num_LER=0
        self._sample_used=0
        self._sample_needed=0
        self._uncertainty=0
        self._estimated_LER=0
        self._samplebudget=0



    def calculate_LER_from_file(self,samplebudget,filepath,pvalue):
        circuit=CliffordCircuit(2)
        circuit.set_error_rate(pvalue)
        self._samplebudget=samplebudget

        stim_str=""
        with open(filepath, "r", encoding="utf-8") as f:
            stim_str = f.read()

        stim_circuit=rewrite_stim_code(stim_str)
        circuit.set_stim_str(stim_circuit)
        circuit.compile_from_stim_circuit_str(stim_circuit)           
        new_stim_circuit=circuit.get_stim_circuit()      

        
        sampler = new_stim_circuit.compile_detector_sampler()
        detector_error_model = new_stim_circuit.detector_error_model(decompose_errors=True)
        matcher = pymatching.Matching.from_detector_error_model(detector_error_model)        


        self._num_LER=0
        while self._num_LER<MIN_NUM_LE_EVENT:
            self._sample_used+=SAMPLE_GAP
            #self._num_LER+=count_logical_errors(new_stim_circuit,SAMPLE_GAP)


            detection_events, observable_flips = sampler.sample(SAMPLE_GAP, separate_observables=True)
            predictions = matcher.decode_batch(detection_events)
                # 3. count mismatches in vectorised form ---------------------------------
            num_errors = np.count_nonzero(observable_flips != predictions)
            self._num_LER+=num_errors

            self._estimated_LER=self._num_LER/self._sample_used
            self.calculate_standard_error()
            print("Current LER: ", self._num_LER)
            print("Current logical error rate: ", self._num_LER/self._sample_used)
            print("Current stdandard error: ", self._uncertainty)
            print("Current sample used: ", self._sample_used)

            if self._sample_used>self._samplebudget:
                print("Sample budget reached, stop sampling")
                if(self._num_LER>0):
                    self._sample_needed=int(self._sample_used*(100/self._num_LER))
                else:
                    self._sample_needed=-1
                break
            self._sample_needed=self._sample_used

        self.calculate_standard_error()
        print("Final LER: ", self._num_LER)
        print("Final logical error rate: ",  self._num_LER/self._sample_used)
        print("Final sample needed: ",  self._sample_needed)
        return self._estimated_LER
    

    def calculate_standard_error(self):
        """
        Calculate the standard error of the LER.
        """
        self._estimated_LER=self._num_LER/self._sample_used
        self._uncertainty = np.sqrt(self._estimated_LER*(1-self._estimated_LER)) / self._sample_used
        return self._uncertainty


    def get_sample_used(self):
        return self._sample_used




if __name__ == "__main__":
    calculator=stimLERcalc()
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/surface/surface30"
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/small/1cnoth"
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/small/simpleh"
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/small/2cnot2R"
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/hexagon/hexagon3"
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/hexagon/hexagon3"
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/surface/surface3"
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/surface/surface9"
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/repetition/repetition7"
    ler=calculator.calculate_LER_from_file(100000,filepath,0.001)

    print(ler)
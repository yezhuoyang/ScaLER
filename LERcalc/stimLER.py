from LERcalc.clifford import *
import pymatching
from LERcalc.stimparser import *
import time
import os
from contextlib import redirect_stdout


def count_logical_errors(circuit: stim.Circuit, num_shots: int) -> int:
    # Sample the circuit.
    sampler = circuit.compile_detector_sampler()
    detection_events, observable_flips = sampler.sample(num_shots, separate_observables=True)

    # Configure a decoder using the circuit.
    detector_error_model = circuit.detector_error_model(decompose_errors=False)
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



def format_with_uncertainty(value, std):
    """
    Format a value and its standard deviation in the form:
    1.23(±0.45)×10^k
    """
    if value == 0:
        return f"0(+{std:.2e})"
    exponent = int(np.floor(np.log10(abs(value))))
    coeff = value / (10**exponent)
    std_coeff = std / (10**exponent)
    return f"{coeff:.2f}(+{std_coeff:.2f})*10^{exponent}"




SAMPLE_GAP_INITIAL = 100
MAX_SAMPLE_GAP = 10000000


'''
Use stim and Monte Calo sampling method to estimate the logical error rate
The sampler will finally decide how many samples to used.
Shot is the initial guess of how many samples to used.
We also need to estimate the uncertainty of the LER.
'''
class stimLERcalc:
    def __init__(self,MIN_NUM_LE_EVENT=3):
        self._num_LER=0
        self._sample_used=0
        self._sample_needed=0
        self._uncertainty=0
        self._estimated_LER=0
        self._samplebudget=0
        self._MIN_NUM_LE_EVENT = MIN_NUM_LE_EVENT


    def calculate_LER_from_file(self,samplebudget,filepath,pvalue, repeat=1):
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


        Ler_list=[]
        samples_list=[]
        time_list=[]
        ler_count_list=[]
        for i in range(repeat):
            
            start = time.time()
            self._num_LER=0
            self._sample_used=0
            current_sample_gap=SAMPLE_GAP_INITIAL
            while self._num_LER<self._MIN_NUM_LE_EVENT:
                if self._num_LER==0 and self._sample_used>0:
                    current_sample_gap*=2
                    current_sample_gap=min(current_sample_gap, MAX_SAMPLE_GAP)
                elif self._num_LER>0:
                    current_sample_gap=min(int(self._MIN_NUM_LE_EVENT/self._num_LER)*self._sample_used, MAX_SAMPLE_GAP)
                self._sample_used+=current_sample_gap
                #self._num_LER+=count_logical_errors(new_stim_circuit,SAMPLE_GAP)


                detection_events, observable_flips = sampler.sample(current_sample_gap, separate_observables=True)
                predictions = matcher.decode_batch(detection_events)
                    # 3. count mismatches in vectorised form ---------------------------------
                num_errors = np.count_nonzero(observable_flips != predictions)
                self._num_LER+=num_errors

                self._estimated_LER=self._num_LER/self._sample_used
                #self.calculate_standard_error()
                # print("Current LER: ", self._num_LER)
                # print("Current logical error rate: ", self._num_LER/self._sample_used)
                # print("Current stdandard error: ", self._uncertainty)
                # print("Current sample used: ", self._sample_used)

                if self._sample_used>self._samplebudget:
                    #print("Sample budget reached, stop sampling")
                    if(self._num_LER>0):
                        self._sample_needed=int(self._sample_used*(100/self._num_LER))
                    else:
                        self._sample_needed=-1
                    break
                self._sample_needed=self._sample_used
            ler_count_list.append(self._num_LER)
            Ler_list.append(self._estimated_LER)
            samples_list.append(self._sample_used)
            elapsed = time.time() - start
            time_list.append(elapsed)

        ler_count_average=np.mean(ler_count_list)
        #print("Average number of logical errors: ", ler_count_average)
        std_ler_count=np.std(ler_count_list)

        self._estimated_LER=np.mean(Ler_list)
        self._sample_used=np.mean(samples_list)
        """
        Standard deviation
        """
        std_ler=np.std(Ler_list)
        std_sample=np.std(samples_list)
        #self.calculate_standard_error()
        time_mean=np.mean(time_list)
        time_std=np.std(time_list)
        
        print("Time(STIM): ", format_with_uncertainty(time_mean, time_std))
        print("PL(STIM): ", format_with_uncertainty(self._estimated_LER, std_ler))
        print("Nerror(STIM): ", format_with_uncertainty(ler_count_average, std_ler_count))
        print("Sample(STIM): ", format_with_uncertainty(self._sample_used, std_sample))        
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


#,0.0001,0.0000
#filepath_list=["surface/surface13","hexagon/hexagon13","square/square11","square/square13","square/square13"]
p_list=[0.00005]
#filepath_list=["repetition/repetition3","repetition/repetition5","repetition/repetition7","square/square9","hexagon/hexagon9","repetition/repetition9","surface/surface9"]


#"surface/surface9","hexagon/hexagon9","square/square9","square/square7"


# filepath_list=["repetition/repetition7","square/square9","hexagon/hexagon9","surface/surface9","surface/surface11","hexagon/hexagon11","square/square11","surface/surface13","hexagon/hexagon13","square/square13"]


if __name__ == "__main__":
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/repetition/repetition5"
    #filepath="C:/Users/yezhu/GitRepos/Sampling/stimprograms/surface/surface5"
    #filepath="C:/Users/yezhu/GitRepos/Sampling/stimprograms/surface/surface7"
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/surface/surface30"
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/small/1cnoth"
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/small/simpleh"
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/small/2cnot2R"
    base_dir = "C:/Users/yezhu/Documents/Sampling/stimprograms/"
    result_dir = "C:/Users/yezhu/Documents/Sampling/StimresultTot/"
    # for rel in filepath_list:
    #     # rel might be "surface/surface3" or "repetition/repetition5", etc.
    #     stim_path = base_dir+rel
    #     # 3) build your output filename:
    #     out_fname = stim_path + "-result.txt"     # e.g. "surface3-result.txt"
    #     # 4) redirect prints for just this file:
    #     with open(out_fname, "w") as outf, redirect_stdout(outf):
    #         print(f"---- Processing {stim_path} ----")

    #         calculator=stimLERcalc()
    #         # pass the string path into your function:
    #         ler = calculator.calculate_LER_from_file(1000000000, str(stim_path), 0.0005,5)



    # p=0.001
    # code_type="surface"
    # rel="surface/surface"
    # dlist=[3,5,7,9,11,13,15]
    # for d in dlist:
    #     stim_path = base_dir+rel+str(d)
    #     # 3) build your output filename:
    #     out_fname =  result_dir+str(p)+"-"+str(code_type)+str(d)+"-result.txt"     # e.g. "surface3-result.txt"
    #     # 4) redirect prints for just this file:
    #     with open(out_fname, "w") as outf, redirect_stdout(outf):
    #         print(f"---- Processing {stim_path} ----")

    #         calculator=stimLERcalc(100)
    #         # pass the string path into your function:
    #         ler = calculator.calculate_LER_from_file(5000000000, str(stim_path), p,5)




    p=0.001
    code_type="hexagon"
    rel="hexagon/hexagon"
    dlist=[3,5,7,9]
    for d in dlist:
        stim_path = base_dir+rel+str(d)
        # 3) build your output filename:
        out_fname =  result_dir+str(p)+"-"+str(code_type)+str(d)+"-result.txt"     # e.g. "surface3-result.txt"
        # 4) redirect prints for just this file:
        with open(out_fname, "w") as outf, redirect_stdout(outf):
            print(f"---- Processing {stim_path} ----")

            calculator=stimLERcalc(100)
            # pass the string path into your function:
            ler = calculator.calculate_LER_from_file(5000000000, str(stim_path), p,5)


    p=0.001
    code_type="square"
    rel="square/square"
    dlist=[3,5,7,9]
    for d in dlist:
        stim_path = base_dir+rel+str(d)
        # 3) build your output filename:
        out_fname =  result_dir+str(p)+"-"+str(code_type)+str(d)+"-result.txt"     # e.g. "surface3-result.txt"
        # 4) redirect prints for just this file:
        with open(out_fname, "w") as outf, redirect_stdout(outf):
            print(f"---- Processing {stim_path} ----")

            calculator=stimLERcalc(100)
            # pass the string path into your function:
            ler = calculator.calculate_LER_from_file(5000000000, str(stim_path), p,5)



    p=0.001
    code_type="repetition"
    rel="repetition/repetition"
    dlist=[3,5,7,9]
    for d in dlist:
        stim_path = base_dir+rel+str(d)
        # 3) build your output filename:
        out_fname =  result_dir+str(p)+"-"+str(code_type)+str(d)+"-result.txt"     # e.g. "surface3-result.txt"
        # 4) redirect prints for just this file:
        with open(out_fname, "w") as outf, redirect_stdout(outf):
            print(f"---- Processing {stim_path} ----")

            calculator=stimLERcalc(100)
            # pass the string path into your function:
            ler = calculator.calculate_LER_from_file(5000000000, str(stim_path), p,5)


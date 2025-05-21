from QEPG.QEPG import return_samples, return_samples_many_weights, return_detector_matrix, return_samples_many_weights_separate_obs, return_samples_numpy
from LERcalc.clifford import *
from LERcalc.stimparser import *
from time import time, perf_counter


def test_samplerate():
    distance=15
    p=0.001
    circuit=CliffordCircuit(2)
    circuit.set_error_rate(p)
    stim_circuit=stim.Circuit.generated("surface_code:rotated_memory_z",rounds=distance*3,distance=distance).flattened()
    stim_circuit=rewrite_stim_code(str(stim_circuit))
    circuit.set_stim_str(stim_circuit)
    circuit.compile_from_stim_circuit_str(stim_circuit)           
    new_stim_circuit=circuit.get_stim_circuit()      
    total_noise=circuit.get_totalnoise()



    current_time = time()
    result=return_samples_numpy(str(new_stim_circuit),5,1000000)
    print("My Time taken for return_samples: ", time()-current_time)


    current_time = time()
    sampler=new_stim_circuit.compile_sampler()
    print("Stim take {} to construct sampler".format(time()-current_time))
    print("Sampler done!")

    current_time = time()
    sampler.sample(shots=1000000)
    print("Stim take {} to sample".format(time()-current_time))

if __name__ == "__main__":
    test_samplerate()
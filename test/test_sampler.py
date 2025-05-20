#Test the correctness of the sampler in QEPG
from QEPG.QEPG import return_samples_with_noise_vector, return_samples_many_weights, return_detector_matrix
from LERcalc.clifford import *
from test.test_QEPG_by_stim import transpile_stim_with_noise_vector




error_rate=0.01
absolute_error=0.05
sample_size=10
num_subspace=3
weight=3



def convert_Floyd_sample_to_noise_vector(sample,num_noise):
    """
    Input:  Floyd samples with the form [(7, 1), (1, 2)]
    Output: The seven'th noise has type X, the first one has type Y
    """
    noise_vector=[0]*num_noise*3
    for i in range(len(sample)):
        if sample[i][1]==1:
            noise_vector[sample[i][0]]=1
        elif sample[i][1]==2:
            noise_vector[sample[i][0]+num_noise]=1
        elif sample[i][1]==3:
            noise_vector[sample[i][0]+2*num_noise]=1
    return noise_vector
    

     
def test_by_file_name(filepath):
    stim_str=""
    with open(filepath, "r", encoding="utf-8") as f:
        stim_str = f.read()
    noise_vector,samples=return_samples_with_noise_vector(stim_str,weight,sample_size)
    print("Floyd samples: ",noise_vector)

    print("samples output: ",samples)
    circuit=CliffordCircuit(3)
    circuit.compile_from_stim_circuit_str(stim_str)
    new_stim_circuit=circuit.get_stim_circuit()
    num_noise = circuit.get_totalnoise()


    detectorMatrix=np.array(return_detector_matrix(str(stim_str)))
    print("Detector matrix: ", detectorMatrix)
    detectorMatrix = detectorMatrix.T          # or: np.transpose(detector_matrix)


    noise_vector=[convert_Floyd_sample_to_noise_vector(x,num_noise) for x in noise_vector]

    print("Noise vector: ", noise_vector)


    for i in range(len(noise_vector)):
        detector_result=transpile_stim_with_noise_vector(str(new_stim_circuit),noise_vector[i],num_noise)
        print("-------------Detector result from stim: -------------")
        print(detector_result)


        #print(dectectorMatrix.shape, noise_vector.shape)
        mydetectorresult=samples[i]


        print("-------------My Detector result: -------------")
        print(mydetectorresult)

        assert((detector_result==list(mydetectorresult)))
        
        print("-------------------Pass test!-------------------------------")


if __name__ == "__main__":


    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/surface/surface3"
    test_by_file_name(filepath)
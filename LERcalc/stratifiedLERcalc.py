
from QEPG.QEPG import return_samples,return_samples_many_weights,return_detector_matrix
from LERcalc.clifford import *
import math
import pymatching

def binomial_weight(N, W, p):
    if N<5000:
        return math.comb(N, W) * ((p)**W) * ((1 - p)**(N - W))
    else:
        lam = N * p
        # PMF(X=W) = e^-lam * lam^W / W!
        # Evaluate in logs to avoid overflow for large W, then exponentiate
        log_pmf = (-lam) + W*math.log(lam) - math.lgamma(W+1)
        return math.exp(log_pmf)




'''
Use strafified sampling algorithm to calculate the logical error rate
'''
class stratifiedLERcalc:
    def __init__(self, error_rate=0, sampleBudget=10000, num_subspace=10):
        self._num_detector=0
        self._num_noise=0
        self._error_rate=error_rate
        self._cliffordcircuit=CliffordCircuit(4)  

        self._LER=0
        """
        Use a dictionary to store the estimated subspace logical error rate
        """
        self._estimated_subspaceLER={}
        self._sampleBudget=sampleBudget
        self._num_subspace=num_subspace
        self._minW=0
        self._maxW=0

        self._stim_str_after_rewrite=""


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
        self._stim_str_after_rewrite=stim_str

        # Configure a decoder using the circuit.
        self._detector_error_model = self._cliffordcircuit.get_stim_circuit().detector_error_model(decompose_errors=False)
        self._matcher = pymatching.Matching.from_detector_error_model(self._detector_error_model)


    def subspace_sampling(self):
        """
        Sample around the subspaces.
        """
        ave_error_weight=self._error_rate*self._num_noise

        """
        Evenly distribute the sample budget across all subspace

        [ave_error_weight-self._num_subspace//2, ave_error_weight+self._num_subspace//2 ]
        """
        
        
        if(ave_error_weight-self._num_subspace//2<0):
            self._minW=0
        else:
            self._minW=int(ave_error_weight-self._num_subspace//2)

        if(ave_error_weight+self._num_subspace//2>self._num_noise):
            self._maxW=self._num_noise
        else:
            self._maxW=int(ave_error_weight+self._num_subspace//2)    

        wlist = list(range(self._minW, self._maxW + 1))
        slist=[self._sampleBudget//self._num_subspace]*len(wlist)


        result=return_samples_many_weights(self._stim_str_after_rewrite,wlist,slist)

        # print("wlist: ",wlist)
        # print("slist: ",slist)
        # print("Result shape is ",len(result)," ",len(result[0])," ",len(result[0][0]))
        for i in range(len(wlist)):
            states, observables = [], []

            for j in range(0,slist[i]):
                states.append(result[i][j][:-1])
                observables.append([result[i][j][-1]])


            shots=len(states)
            predictions =self._matcher.decode_batch(states)
            num_errors = 0
            for shot in range(shots):
                actual_for_shot = observables[shot]
                predicted_for_shot = predictions[shot]
                if not np.array_equal(actual_for_shot, predicted_for_shot):
                    num_errors += 1

            self._estimated_subspaceLER[wlist[i]]=num_errors/shots
            print("Logical error rate when w={} ".format(wlist[i])+str(self._estimated_subspaceLER[wlist[i]]))
        




    # ----------------------------------------------------------------------
    # Calculate logical error rate
    # The input is a list of rows with logical errors
    def calculate_LER(self):
        self._LER=0
        for weight in range(1,self._num_noise+1):
            if weight in self._estimated_subspaceLER.keys():
                self._LER+=self._estimated_subspaceLER[weight]*binomial_weight(self._num_noise, weight,self._error_rate)
        return self._LER    


    def get_LER_subspace(self,weight):
        return self._estimated_subspaceLER[weight]*binomial_weight(self._num_noise, weight,self._error_rate)


    def calculate_LER_from_file(self,filepath,pvalue):
        pass



if __name__ == "__main__":
    tmp=stratifiedLERcalc(0.001,sampleBudget=500000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/surface/surface9"
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/small/1cnot"
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/small/surface3r1"
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/small/cnot01h01"
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/small/1cnoth"
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/small/simpleh"
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/small/1cnot1R"
    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/small/2cnot2R"
    tmp.parse_from_file(filepath)
    tmp.subspace_sampling()


    LER=tmp.calculate_LER()

    print(LER)

    num_noise=tmp._num_noise

    for weight in range(4,14):
        print("LER in the subspace {} is {}".format(weight,tmp.get_LER_subspace(weight)))    



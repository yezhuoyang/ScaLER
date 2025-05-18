
from QEPG.QEPG import return_samples,return_samples_many_weights,return_detector_matrix
from LERcalc.clifford import *
import math
import pymatching
from scipy.optimize import curve_fit
from scipy.stats import norm

def binomial_weight(N, W, p):
    if N<200:
        return math.comb(N, W) * (p**W) * ((1 - p)**(N - W))
    else:
        lam = N * p
        # PMF(X=W) = e^-lam * lam^W / W!
        # Evaluate in logs to avoid overflow for large W, then exponentiate
        log_pmf = (-lam) + W*math.log(lam) - math.lgamma(W+1)
        return math.exp(log_pmf)



def scurve_function(x, mu, sigma):
    cdf_values = 0.5*norm.cdf(x, loc=mu, scale=sigma)
    return cdf_values


'''
Use strafified sampling + Scurve fitting  algorithm to calculate the logical error rate
'''
class stratified_Scurve_LERcalc:
    def __init__(self, error_rate=0, sampleBudget=10000, num_subspace=5):
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
        self._minw=0
        self._maxw=0
        self._estimated_wlist=[]

        self._stim_str_after_rewrite=""




    def calc_logical_error_rate_with_fixed_w(self, shots, w):
        """
        Calculate the logical error rate with fixed w
        """
        result=return_samples(self._stim_str_after_rewrite,w,shots)
        states, observables = [], []

        for j in range(0,shots):
            states.append(result[j][:-1])
            observables.append([result[j][-1]])

        shots=len(states)
        predictions =self._matcher.decode_batch(states)
        num_errors = 0
        for shot in range(shots):
            actual_for_shot = observables[shot]
            predicted_for_shot = predictions[shot]
            if not np.array_equal(actual_for_shot, predicted_for_shot):
                num_errors += 1

        return num_errors/shots


    '''
    Use binary search to determine the 
    exact number of errors that give non-zero logical error rate
    We just try 10 samples
    '''
    def binary_search_zero(self,low,high,shots):
        left=low
        right=high
        while left<right:
            print(left,right)
            mid=(left+right)//2
            er=self.calc_logical_error_rate_with_fixed_w(shots,mid)
            print("er: ",er)
            if er>0:
                right=mid
            else:
                left=mid+1
        return left


    '''
    Use binary search to determine the exact number of errors 
    that give saturate logical error rate
    We just try 10 samples
    '''
    def binary_search_half(self,low,high, shots, epsilon=0.05):
        left=low
        right=high
        while left<right:
            print(left,right)
            mid=(left+right)//2
            er=self.calc_logical_error_rate_with_fixed_w(shots,mid)
            if er>(0.5-epsilon):
                right=mid
            else:
                left=mid+1
        return left



    def determine_w(self,shots=100):
        """
        Use binary search to determine the minw and maxw
        """
        self._minw=0
        self._maxw=self.binary_search_half(0,self._num_noise,shots)
        print("minw: ",self._minw)
        print("maxw: ",self._maxw)



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
        After we determine the minw and maxw, we generate an even distribution of points 
        between minw and maxw
        """

        gap=int((self._maxw-self._minw)/self._num_subspace)


        wlist =list(np.arange(self._minw, self._maxw, gap, dtype=int))
        self._estimated_wlist=wlist
        print("wlist: ",wlist)
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


    '''
    Fit the distribution by 1/2-e^{alpha/W}
    '''
    def fit_Scurve(self):
        # Initial guess for alpha
        initial_guess = [10,1.0]


        sigmas=[x*2 for x in self._estimated_wlist]

        # Perform the curve fit with the bounds
        popt, pcov = curve_fit(
            scurve_function, 
            self._estimated_wlist, 
            [self._estimated_subspaceLER[x] for x in self._estimated_wlist], 
            p0=initial_guess, 
            sigma=sigmas,
        )

        # Extract the best-fit parameter (alpha)
        mu,sigma = popt[0],popt[1]
        return mu,sigma


    def calc_logical_error_rate_by_curve_fitting(self):
        mu,sigma=self.fit_Scurve()
        self._LER=0
        for weight in range(1,self._num_noise+1):
            """
            If the weight is in the estimated list, we use the estimated value
            Else, we use the curve fitting value
            """
            if weight in self._estimated_wlist:
                self._LER+=self._estimated_subspaceLER[weight]*binomial_weight(self._num_noise,weight,self._error_rate)  
            else:
                self._LER+=scurve_function(weight,mu,sigma)*binomial_weight(self._num_noise,weight,self._error_rate)
        return self._LER


if __name__ == "__main__":
    tmp=stratified_Scurve_LERcalc(0.0001,sampleBudget=10000,num_subspace=8)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/surface/surface3"
    tmp.parse_from_file(filepath)


    #result=tmp.calc_logical_error_rate_with_fixed_w(3,2)
    #print(result)
    #tmp.subspace_sampling()


    tmp.determine_w(1000)

    tmp.subspace_sampling()

    LER=tmp.calculate_LER()
    print("LER: ",LER)






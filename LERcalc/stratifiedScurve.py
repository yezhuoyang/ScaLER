
from QEPG.QEPG import return_samples,return_samples_many_weights,return_detector_matrix,return_samples_many_weights_numpy,return_samples_many_weights_separate_obs
from LERcalc.clifford import *
import math
import pymatching
from scipy.optimize import curve_fit
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


def binomial_weight(N, W, p):
    if N<200:
        return math.comb(N, W) * (p**W) * ((1 - p)**(N - W))
    else:
        lam = N * p
        # PMF(X=W) = e^-lam * lam^W / W!
        # Evaluate in logs to avoid overflow for large W, then exponentiate
        log_pmf = (-lam) + W*math.log(lam) - math.lgamma(W+1)
        return math.exp(log_pmf)



# def scurve_function(x, mu, sigma):
#     cdf_values = 0.5*norm.cdf(x, loc=mu, scale=sigma)
#     return cdf_values


def scurve_function(x, center, sigma):
    return 0.5/(1+np.exp(-(x - center) / sigma))
    #return 0*x



def scurve_function_with_distance(x, cd, mu, sigma):
    """
    Piece-wise S-curve:
        0                          for x < cd
        0.5 * Φ((x - μ) / σ)       for x ≥ cd
    where Φ is the standard normal CDF.
    """
    x = np.asarray(x)                      # ensure array
    y = 0.5 * norm.cdf(x, loc=mu, scale=sigma)
    return np.where(x < cd, 0.0, y)        # vectorised “if”




def r_squared(y_true, y_pred, clip=False):
    """
    Compute the coefficient of determination (R²).

    Parameters
    ----------
    y_true : array-like
        Observed data.
    y_pred : array-like
        Model-predicted data (same length as y_true).
    clip : bool, default False
        If True, negative R² values are clipped to 0 so the
        score lies strictly in the interval [0, 1].

    Returns
    -------
    float
        The R² statistic.
    """
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)

    if y_true.shape != y_pred.shape:
        raise ValueError("y_true and y_pred must have the same shape")

    ss_res = np.sum((y_true - y_pred) ** 2)        # residual sum of squares
    ss_tot = np.sum((y_true - y_true.mean()) ** 2) # total sum of squares

    # Handle the degenerate case where variance is zero
    if ss_tot == 0.0:
        return 1.0 if ss_res == 0.0 else 0.0

    r2 = 1.0 - ss_res / ss_tot
    return max(0.0, r2) if clip else r2


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
        Use a dictionary to store the estimated subspace logical error rate,
        how many samples have been used in each subspace
        """
        self._estimated_subspaceLER={}
        self._subspace_LE_count={}
        self._estimated_subspaceLER_second={}
        self._subspace_sample_used={}

        self._sampleBudget=sampleBudget
        self._sample_used=0
        self._circuit_level_code_distance=1
        self._num_subspace=num_subspace
        self._minw=1
        self._maxw=10000000000000
        self._saturatew=10000000000000       
        self._estimated_wlist=[]

        self._stim_str_after_rewrite=""

        self._mu=0
        self._sigma=0

        #In the area we are interested in, the maximum value of the logical error rate
        self._rough_value_for_subspace_LER=0




    def calc_logical_error_rate_with_fixed_w(self, shots, w):
        """
        Calculate the logical error rate with fixed w
        TODO: Optimize this function
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
    Use binary search to determine the exact number of errors 
    that give saturate logical error rate
    We just try 10 samples
    '''
    def binary_search_half(self,low,high, shots):
        left=low
        right=high
        epsilon=0.5-10*self._rough_value_for_subspace_LER
        while left<right:
            print(left,right)
            mid=(left+right)//2
            er=self.calc_logical_error_rate_with_fixed_w(shots,mid)
            if er>(0.5-epsilon):
                right=mid
            else:
                left=mid+1
        return left



    def determine_w(self,shots=10000):
        """
        Use binary search to determine the minw and maxw
        """
        self._minw=1
        self._maxw=self.binary_search_half(1,self._num_noise,shots)
        self._saturatew=self._maxw
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



    def get_rough_subspace_LER(self,maximum_shots):
        """
        Input: the maximum number of shots
        Output: the rough value of the subspace logical error rate
        """
        max_estimated_subspaceLER=0
        shots=1000
        while max_estimated_subspaceLER==0:
            self._sampleBudget=shots
            wlist,slist=self.distribute_samples(equal=True)
            for w,s in zip(wlist,slist):
                if not w in self._subspace_LE_count.keys():
                    self._subspace_LE_count[w]=0
                    self._subspace_sample_used[w]=s
                else:
                    self._subspace_sample_used[w]+=s
            detector_result,obsresult=return_samples_many_weights_separate_obs(self._stim_str_after_rewrite,wlist,slist)
            predictions_result = self._matcher.decode_batch(detector_result)

            begin_index=0
            for w_idx, (w, quota) in enumerate(zip(wlist, slist)):

                observables =  np.asarray(obsresult[begin_index:begin_index+quota])                    # (shots,)
                predictions = np.asarray(predictions_result[begin_index:begin_index+quota]).ravel()

                # 3. count mismatches in vectorised form ---------------------------------
                num_errors = np.count_nonzero(observables != predictions)

                self._subspace_LE_count[w]+=num_errors
                self._estimated_subspaceLER_second[w] = self._subspace_LE_count[w] / self._subspace_sample_used[w]
                self._estimated_subspaceLER[w]=self._subspace_LE_count[w]/self._subspace_sample_used[w]
                
                if self._estimated_subspaceLER[w]>max_estimated_subspaceLER:
                    max_estimated_subspaceLER=self._estimated_subspaceLER[w]

                begin_index+=quota
            shots*=10
            if(shots>maximum_shots):
                break

        self._rough_value_for_subspace_LER=max_estimated_subspaceLER



    def subspace_sampling_first_round(self,sampleBudget):

        """
        After we determine the minw and maxw, we generate an even distribution of points 
        between minw and maxw.

        The goal is for the curve fitting in the next step to get more accurate.
        """

        gap=int((self._maxw-self._minw)/self._num_subspace)
        if gap==0:
            gap=1

        wlist =list(np.arange(self._minw, self._maxw, gap, dtype=int))
        self._estimated_wlist=wlist
        slist=[sampleBudget//self._num_subspace]*len(wlist)


        detector_result,obsresult=return_samples_many_weights_separate_obs(self._stim_str_after_rewrite,wlist,slist)
        predictions_result = self._matcher.decode_batch(detector_result)

        for w,s in zip(wlist,slist):
            if not w in self._subspace_LE_count.keys():
                self._subspace_LE_count[w]=0
                self._subspace_sample_used[w]=s
                self._estimated_subspaceLER[w]=0
            else:
                self._subspace_sample_used[w]+=s        

        begin_index=0
        for w_idx, (w, quota) in enumerate(zip(wlist, slist)):

            observables =  np.asarray(obsresult[begin_index:begin_index+quota])                    # (shots,)
            predictions = np.asarray(predictions_result[begin_index:begin_index+quota]).ravel()

            # 3. count mismatches in vectorised form ---------------------------------
            num_errors = np.count_nonzero(observables != predictions)

            self._subspace_LE_count[w]+=num_errors
            self._estimated_subspaceLER[w]=self._subspace_LE_count[w]/self._subspace_sample_used[w]
            #print("Logical error rate when w={} ".format(w)+str(self._estimated_subspaceLER[w]))
            begin_index+=quota
        print("---------------After first round sampling-------------------")
        print("Subspace LE count: ",self._subspace_LE_count)
        print("self._subspace_sample_used: ",self._subspace_sample_used)
        print("self._estimated_subspaceLER: ",self._estimated_subspaceLER)


    def calculate_R_square_score(self):
        y_observed = [self._estimated_subspaceLER[x] for x in self._estimated_wlist]
        y_predicted = [scurve_function(x,self._mu,self._sigma) for x in self._estimated_wlist]
        r2 = r_squared(y_observed, y_predicted)
        print("R^2 score: ", r2)
        return r2


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


    '''
    Fit the distribution by 1/2-e^{alpha/W}
    '''
    def fit_Scurve(self):
        # Initial guess for alpha
        self._saturatew
        # sigma_guess = self._saturatew / 3.2898
        # mu_guess    = self._saturatew/2          # centre in the middle of that span
        # initial_guess  = (0,mu_guess, sigma_guess)

        center = self._saturatew /2 
        sigma    = self._saturatew/7          # centre in the middle of that span
        initial_guess  = (center, sigma )


        # initial_guess  = (mu_guess, sigma_guess)
        # Perform the curve fit with the bounds
        # popt, pcov = curve_fit(
        #     scurve_function, 
        #     self._estimated_wlist, 
        #     [self._estimated_subspaceLER[x] for x in self._estimated_wlist], 
        #     p0=initial_guess
        # )

        popt, pcov = curve_fit(
            scurve_function, 
            self._estimated_wlist, 
            [self._estimated_subspaceLER[x] for x in self._estimated_wlist], 
            p0=initial_guess
        )

        self._codedistance = 0
        # Extract the best-fit parameter (alpha)
        self._mu,self._sigma = popt[0] , popt[1]
        return self._codedistance,self._mu,self._sigma


    def calc_logical_error_rate_after_curve_fitting(self):
        #self.fit_Scurve()
        self._LER=0
        for weight in range(self._minw,self._num_noise+1):
            """
            If the weight is in the estimated list, we use the estimated value
            Else, we use the curve fitting value

            If the weight is less than the minw, we just declare it as 0
            """
            if weight in self._estimated_wlist:
                self._LER+=self._estimated_subspaceLER[weight]*binomial_weight(self._num_noise,weight,self._error_rate)  
            else:
                self._LER+=scurve_function(weight,self._mu,self._sigma)*binomial_weight(self._num_noise,weight,self._error_rate)
                #self._LER+=scurve_function(weight,self._mu,self._sigma)*binomial_weight(self._num_noise,weight,self._error_rate)
        return self._LER


    """
    After the first round of sampling, we use the fitted curve value to determine the 
    number of samples needed for each subspace
    """
    def distribute_samples(self, equal=False):
        proportion_list=[]
        """
        Sample around the subspaces.
        """
        ave_error_weight=self._error_rate*self._num_noise

        """
        Evenly distribute the sample budget across all subspace
        [ave_error_weight-self._num_subspace//2, ave_error_weight+self._num_subspace//2 ]
        """
        if(ave_error_weight-self._num_subspace//2<0):
            self._minw=1
        else:
            self._minw=int(ave_error_weight-self._num_subspace//2)

        if(ave_error_weight+self._num_subspace//2>self._num_noise):
            self._maxw=self._num_noise
        else:
            self._maxw=int(ave_error_weight+self._num_subspace//2)    

        if self._maxw<self._num_subspace:
            self._maxw=self._num_subspace

        wlist = list(range(self._minw, self._maxw + 1))       

        if equal:
            proportion_list = [1 for i in range(len(wlist))]
        else:
            for i in range(len(wlist)):
                if wlist[i] in self._estimated_wlist:
                    proportion_list.append(self._estimated_subspaceLER[wlist[i]])
                else:
                    proportion_list.append(scurve_function(wlist[i],self._mu,self._sigma)) 

        maxelement=max(proportion_list)


        """
        If all the elements are 0, we set all of them to 1
        """
        if(maxelement==0):
            proportion_list=[1 for i in range(len(proportion_list))]  
        else:
            """
            If some of the elements are 0, we set them to the minimum/10 of the other elements
            """
            proportion_list = [1 / (x if x else min(y for y in proportion_list if y)/10) for x in proportion_list]

        #We Normalize the proportions
        sum_proportion = sum(proportion_list)
        proportion_list = [x / sum_proportion for x in proportion_list]

        slist = [int(x * self._sampleBudget) or 10  for x in proportion_list]


        return wlist,slist


    def subspace_sampling_second_round(self):
        """
        The second round of subspace stratified sampling
        We already know the rough value of the subspace logical error rate
        We also have get the fitted curve.
        The goal is to get more accurate value of the logical error rate, specially for these subspaces
        """
        wlist,slist=self.distribute_samples()     
        detector_result,obsresult=return_samples_many_weights_separate_obs(self._stim_str_after_rewrite,wlist,slist)
        predictions_result = self._matcher.decode_batch(detector_result)


        for w,s in zip(wlist,slist):
            if not w in self._subspace_LE_count.keys():
                self._subspace_LE_count[w]=0
                self._subspace_sample_used[w]=s
            else:
                self._subspace_sample_used[w]+=s

        begin_index=0
        for w_idx, (w, quota) in enumerate(zip(wlist, slist)):

            observables =  np.asarray(obsresult[begin_index:begin_index+quota])                    # (shots,)
            predictions = np.asarray(predictions_result[begin_index:begin_index+quota]).ravel()

            # 3. count mismatches in vectorised form ---------------------------------
            num_errors = np.count_nonzero(observables != predictions)

            self._subspace_LE_count[w]+=num_errors
            self._estimated_subspaceLER_second[w] = self._subspace_LE_count[w] / self._subspace_sample_used[w]
            self._estimated_subspaceLER[w]=self._subspace_LE_count[w]/self._subspace_sample_used[w]
            begin_index+=quota


        print("---------------After second round sampling-------------------")
        print("Subspace LE count: ",self._subspace_LE_count)
        print("self._subspace_sample_used: ",self._subspace_sample_used)
        print("Samples used:{}".format(self._sample_used))
        print("circuit level code distance:{}".format(self._circuit_level_code_distance))



    def plot_scurve(self, filename,title="S-curve"):
        """Plot the S-curve and its discrete estimate."""


        print("LER count: ",self._subspace_LE_count)
        keys   = list(self._estimated_subspaceLER.keys())
        values = [self._estimated_subspaceLER[k] for k in keys]

        print("keys: ",keys)
        print("values: ",values)
        fig=plt.figure()
        # bars ── discrete estimate
        plt.bar(keys, values,
                color='tab:blue',         # pick any color you like
                alpha=0.8,
                label='Estimated subspace LER by sampling')


        keys_second = list(self._estimated_subspaceLER_second.keys())
        values_second = [self._estimated_subspaceLER_second[k] for k in keys_second]


        print("keys_second: ",keys_second)
        print("values_second: ",values_second)
        plt.bar(keys_second, values_second,
                color='tab:red',         # pick any color you like
                alpha=0.8,
                label='Second round sampling LER by sampling')        


        plt.axvline(x=self._error_rate*self._num_noise, color="red", linestyle="--", linewidth=1.2, label="Average Error number") # vertical line at x=0.5

        print("Saturate w: ",self._saturatew)

        # smooth curve ── fitted S-curve
        x = np.linspace(0, self._saturatew, 1000)
        y = scurve_function(x, self._mu, self._sigma)
        plt.plot(x, y,
                color='tab:orange',
                linewidth=2.0,
                label='Fitted S-curve')

        plt.xlabel('Weight')
        plt.ylabel('Logical Error Rate in subspace')
        plt.title(title+" (PL={})".format(self._LER))
        plt.legend()                     # <- shows the two labels

        # ── Force integer ticks on the X axis ──────────────────────────
        plt.gca().xaxis.set_major_locator(mticker.MaxNLocator(integer=True))


        plt.tight_layout()               # optional: nicely fit everything
        plt.savefig(filename, dpi=300)
        #plt.show()
        plt.close(fig)



    def calculate_LER_from_file(self,filepath,pvalue,figname,titlename):
        self._error_rate=pvalue
        self.parse_from_file(filepath)
        self.determine_w()
        self.subspace_sampling_first_round(int(self._sampleBudget//4))
        self.fit_Scurve()
        self.subspace_sampling_second_round()
        self.fit_Scurve()        
        self.calc_logical_error_rate_after_curve_fitting()
        #self.calculate_LER()
        print("Final LER: ",self._LER)
        self.plot_scurve(figname,titlename)





def generate_all_surface_code_figure():
    p=0.001

    tmp=stratified_Scurve_LERcalc(p,sampleBudget=100000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/surface/surface3"
    ler=tmp.calculate_LER_from_file(filepath,p,"S3.png","Surface3")



    tmp=stratified_Scurve_LERcalc(p,sampleBudget=500000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/surface/surface5"
    ler=tmp.calculate_LER_from_file(filepath,p,"S5.png","Surface5")


    tmp=stratified_Scurve_LERcalc(p,sampleBudget=800000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/surface/surface7"
    ler=tmp.calculate_LER_from_file(filepath,p,"S7.png","Surface7")



    tmp=stratified_Scurve_LERcalc(p,sampleBudget=100000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/surface/surface9"
    ler=tmp.calculate_LER_from_file(filepath,p,"S9.png","Surface9")



    tmp=stratified_Scurve_LERcalc(p,sampleBudget=120000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/surface/surface11"
    ler=tmp.calculate_LER_from_file(filepath,p,"S11.png","Surface11")



    tmp=stratified_Scurve_LERcalc(p,sampleBudget=140000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/surface/surface13"
    ler=tmp.calculate_LER_from_file(filepath,p,"S13.png","Surface13")

    tmp=stratified_Scurve_LERcalc(p,sampleBudget=160000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/surface/surface15"
    ler=tmp.calculate_LER_from_file(filepath,p,"S15.png","Surface15")



    tmp=stratified_Scurve_LERcalc(p,sampleBudget=180000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/surface/surface17"
    ler=tmp.calculate_LER_from_file(filepath,p,"S17.png","Surface17")



    tmp=stratified_Scurve_LERcalc(p,sampleBudget=200000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/surface/surface19"
    ler=tmp.calculate_LER_from_file(filepath,p,"S19.png","Surface19")



    tmp=stratified_Scurve_LERcalc(p,sampleBudget=220000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/surface/surface21"
    ler=tmp.calculate_LER_from_file(filepath,p,"S21.png","Surface21")




def generate_all_repetition_code_figure():
    p=0.001

    tmp=stratified_Scurve_LERcalc(p,sampleBudget=100000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/repetition/repetition3"
    ler=tmp.calculate_LER_from_file(filepath,p,"R3.png","Repetition3")



    tmp=stratified_Scurve_LERcalc(p,sampleBudget=500000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/repetition/repetition5"
    ler=tmp.calculate_LER_from_file(filepath,p,"R5.png","Repetition5")


    tmp=stratified_Scurve_LERcalc(p,sampleBudget=800000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/repetition/repetition7"
    ler=tmp.calculate_LER_from_file(filepath,p,"R7.png","Repetition7")



    tmp=stratified_Scurve_LERcalc(p,sampleBudget=100000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/repetition/repetition9"
    ler=tmp.calculate_LER_from_file(filepath,p,"R9.png","Repetition9")



    tmp=stratified_Scurve_LERcalc(p,sampleBudget=120000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/repetition/repetition11"
    ler=tmp.calculate_LER_from_file(filepath,p,"R11.png","Repetition11")



    tmp=stratified_Scurve_LERcalc(p,sampleBudget=140000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/repetition/repetition13"
    ler=tmp.calculate_LER_from_file(filepath,p,"R13.png","Repetition13")

    tmp=stratified_Scurve_LERcalc(p,sampleBudget=160000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/repetition/repetition15"
    ler=tmp.calculate_LER_from_file(filepath,p,"R15.png","Repetition15")



    tmp=stratified_Scurve_LERcalc(p,sampleBudget=180000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/repetition/repetition17"
    ler=tmp.calculate_LER_from_file(filepath,p,"R17.png","Repetition17")




def generate_all_hexagon_code_figure():
    p=0.001

    tmp=stratified_Scurve_LERcalc(p,sampleBudget=100000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/hexagon/hexagon3"
    ler=tmp.calculate_LER_from_file(filepath,p,"H3.png","hexagon3")



    tmp=stratified_Scurve_LERcalc(p,sampleBudget=500000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/hexagon/hexagon5"
    ler=tmp.calculate_LER_from_file(filepath,p,"H5.png","hexagon5")


    tmp=stratified_Scurve_LERcalc(p,sampleBudget=800000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/hexagon/hexagon7"
    ler=tmp.calculate_LER_from_file(filepath,p,"H7.png","hexagon7")



    tmp=stratified_Scurve_LERcalc(p,sampleBudget=100000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/hexagon/hexagon9"
    ler=tmp.calculate_LER_from_file(filepath,p,"H9.png","hexagon9")



    tmp=stratified_Scurve_LERcalc(p,sampleBudget=120000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/hexagon/hexagon11"
    ler=tmp.calculate_LER_from_file(filepath,p,"H11.png","hexagon11")



    tmp=stratified_Scurve_LERcalc(p,sampleBudget=140000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/hexagon/hexagon13"
    ler=tmp.calculate_LER_from_file(filepath,p,"H13.png","hexagon13")

    tmp=stratified_Scurve_LERcalc(p,sampleBudget=160000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/hexagon/hexagon15"
    ler=tmp.calculate_LER_from_file(filepath,p,"H15.png","hexagon15")




def generate_all_square_code_figure():
    p=0.001

    tmp=stratified_Scurve_LERcalc(p,sampleBudget=100000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/square/square3"
    ler=tmp.calculate_LER_from_file(filepath,p,"Sq3.png","square3")



    tmp=stratified_Scurve_LERcalc(p,sampleBudget=500000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/square/square5"
    ler=tmp.calculate_LER_from_file(filepath,p,"Sq5.png","square5")


    tmp=stratified_Scurve_LERcalc(p,sampleBudget=800000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/square/square7"
    ler=tmp.calculate_LER_from_file(filepath,p,"Sq7.png","square7")



    tmp=stratified_Scurve_LERcalc(p,sampleBudget=100000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/square/square9"
    ler=tmp.calculate_LER_from_file(filepath,p,"Sq9.png","square9")



    tmp=stratified_Scurve_LERcalc(p,sampleBudget=120000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/square/square11"
    ler=tmp.calculate_LER_from_file(filepath,p,"Sq11.png","square11")



    tmp=stratified_Scurve_LERcalc(p,sampleBudget=140000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/square/square13"
    ler=tmp.calculate_LER_from_file(filepath,p,"Sq13.png","square13")

    tmp=stratified_Scurve_LERcalc(p,sampleBudget=160000,num_subspace=10)
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/square/square15"
    ler=tmp.calculate_LER_from_file(filepath,p,"Sq15.png","square15")




if __name__ == "__main__":
    #generate_all_repetition_code_figure()
    #generate_all_square_code_figure()
    #generate_all_hexagon_code_figure()
    p=0.001

    tmp=stratified_Scurve_LERcalc(p,sampleBudget=10000000,num_subspace=10)

    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/surface/surface7"
    figname="S7.png"
    titlename="Surface7"

    tmp.parse_from_file(filepath)
    tmp.get_rough_subspace_LER(10000)
    tmp.determine_w()

    tmp.subspace_sampling_first_round(int(10000000//4))
    tmp.fit_Scurve()
    tmp._sampleBudget=10000000
    tmp.subspace_sampling_second_round()
    tmp.fit_Scurve()        
    tmp.calc_logical_error_rate_after_curve_fitting()
    #self.calculate_LER()
    print("Final LER: ",tmp._LER)
    tmp.plot_scurve(figname,titlename)
    #ler=tmp.calculate_LER_from_file(filepath,p,"S5.png","Surface5")
    #tmp.calculate_R_square_score()

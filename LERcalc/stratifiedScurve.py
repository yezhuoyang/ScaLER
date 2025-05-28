
from QEPG.QEPG import return_samples,return_samples_many_weights,return_detector_matrix,return_samples_many_weights_numpy,return_samples_many_weights_separate_obs
from LERcalc.clifford import *
import math
import pymatching
from scipy.optimize import curve_fit
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.stats import binom

# def binomial_weight(N, W, p):
#     if N<200:
#         return math.comb(N, W) * (p**W) * ((1 - p)**(N - W))
#     else:
#         lam = N * p
#         # PMF(X=W) = e^-lam * lam^W / W!
#         # Evaluate in logs to avoid overflow for large W, then exponentiate
#         log_pmf = (-lam) + W*math.log(lam) - math.lgamma(W+1)
#         return math.exp(log_pmf)



def binomial_weight(N, W, p):
    return binom.pmf(W, N, p)
    #return math.comb(N, W) * (p**W) * ((1 - p)**(N - W))

def linear_function(x, a, b):
    """
    Linear function for curve fitting.
    """
    return a * x + b


def modified_linear_function(x, a, b,c,d):
    """
    Linear function for curve fitting.
    """
    return a * x + b+c/(x-d)


def quadratic_function(x, a, b,c):
    """
    Linear function for curve fitting.
    """
    return a * x**2+b*x + c


def poly_function(x, a, b,c,d):
    """
    Linear function for curve fitting.
    """
    return a * x**3+b*x**2 + c*x+d


MIN_NUM_LE_EVENT = 200
SAMPLE_GAP=100

def subspace_size(num_noise, weight):
    """
    Calculate the size of the subspace
    """
    return math.comb(num_noise, weight)

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




def evenly_spaced_ints(minw, maxw, N):
    if N == 1:
        return [minw]
    if N > (maxw - minw + 1):
        return list(range(minw, maxw + 1))
     
    # Use high-resolution linspace, round, then deduplicate
    raw = np.linspace(minw, maxw, num=10 * N)
    rounded = sorted(set(map(int, raw)))
    
    # Pick N evenly spaced indices from the unique set
    indices = np.linspace(0, len(rounded) - 1, num=N, dtype=int)
    return [rounded[i] for i in indices]


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
        """
        minw and maxw store the range of subspace we need to fit.
        This is determined by the uncertainty value
        """
        self._minw=1
        self._maxw=10000000000000
        """
        self._saturatew is the weight of the subspace where the 
        logical error get satureated
        """
        self._saturatew=10000000000000       
        self._estimated_wlist=[]

        self._stim_str_after_rewrite=""

        self._mu=0
        self._sigma=0

        #In the area we are interested in, the maximum value of the logical error rate
        self._rough_value_for_subspace_LER=0

        self._stratified_succeed=False




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

    TODO: Restructure the function.
    Add the threshold as an input parameter.
    '''
    def binary_search_half(self,low,high, shots):
        left=low
        right=high
        epsion=0.45
        while left<right:
            print(left,right)
            mid=(left+right)//2
            er=self.calc_logical_error_rate_with_fixed_w(shots,mid)
            if er>epsion:
                right=mid
            else:
                left=mid+1
        return left



    def determine_saturated_w(self,shots=1000):
        """
        Use binary search to determine the minw and maxw
        """
        self._saturatew=self.binary_search_half(self._minw,self._num_noise,shots)
        print("Self._saturatew: ",self._saturatew)



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
        self._detector_error_model = self._cliffordcircuit.get_stim_circuit().detector_error_model(decompose_errors=True)
        self._matcher = pymatching.Matching.from_detector_error_model(self._detector_error_model)



    def determine_range_to_sample(self,k=3):
        """
        We need to be exact about the range of w we want to sample. 
        We don't want to sample too many subspaces, especially those subspaces with tiny binomial weights.
        This should comes from the analysis of the weight of each subspace.

        We use the standard deviation to approimxate the range
        """
        sigma=int(np.sqrt(self._error_rate*(1-self._error_rate)*self._num_noise))
        if sigma==0:
            sigma=1
        ep=int(self._error_rate*self._num_noise)
        self._minw=max(1,ep-3*sigma)
        self._maxw=max(self._num_subspace,ep+k*sigma)

        print("Needed Minw:    ")
        print(self._minw)
        print("Needed Maxw:    ")
        print(self._maxw)




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
        print("Rough value for subspace LER:")
        print(max_estimated_subspaceLER)



    def subspace_sampling(self):
        """
        wlist store the subset of weights we need to sample and get
        correct logical error rate.
        
        In each subspace, we stop sampling until 100 logical error events are detected, or we hit the total budget.
        """
        #wlist_need_to_sample = list(range(self._minw, self._maxw + 1))
        wlist_need_to_sample=evenly_spaced_ints(self._minw, self._maxw, self._num_subspace)

        for weight in wlist_need_to_sample:
            if not weight in self._estimated_wlist:
                self._estimated_wlist.append(weight)
                self._subspace_LE_count[weight]=0
                self._subspace_sample_used[weight]=0

        print(wlist_need_to_sample)

        min_non_zero_weight=1e9
        self._sample_used=0
        while True:
            if(self._sample_used>self._sampleBudget):
                break
            slist=[]
            wlist=[]
            """
            Case 1 to end the while loop: We have consumed all of our sample budgets
            """
            if(self._sample_used>self._sampleBudget):
                break

            for weight in wlist_need_to_sample:
                if weight<=self._circuit_level_code_distance:
                    self._estimated_subspaceLER[weight]=0
                    continue
                if(self._subspace_LE_count[weight]==0):
                    if(subspace_size(self._num_noise, weight)<2*self._subspace_sample_used[weight]):
                        self._circuit_level_code_distance=weight
                        self._estimated_subspaceLER[weight]=0
                        print("Circuit level code distance: ",self._circuit_level_code_distance)
                        continue         
                if(self._subspace_LE_count[weight]<MIN_NUM_LE_EVENT):
                    if(self._subspace_LE_count[weight]>=1):
                        sample_num_required=int(MIN_NUM_LE_EVENT/self._subspace_LE_count[weight])* self._subspace_sample_used[weight]
                        slist.append(sample_num_required)
                        self._subspace_sample_used[weight]+=sample_num_required  
                        self._sample_used+=sample_num_required
                    else:                   
                        slist.append(max(SAMPLE_GAP,self._subspace_sample_used[weight]*2))
                        self._subspace_sample_used[weight]+=max(SAMPLE_GAP,self._subspace_sample_used[weight]*2)
                        self._sample_used+=max(SAMPLE_GAP,self._subspace_sample_used[weight]*2)
                    wlist.append(weight)
            """
            Case 2 to end the while loop: We have get 100 logical error events for all these subspaces
            """
            if(len(wlist)==0):
                break

            print("wlist: ",wlist)
            print("slist: ",slist)
            detector_result,obsresult=return_samples_many_weights_separate_obs(self._stim_str_after_rewrite,wlist,slist)
            predictions_result = self._matcher.decode_batch(detector_result)
            print("Result get!")

            
            begin_index=0
            for w_idx, (w, quota) in enumerate(zip(wlist, slist)):

                observables =  np.asarray(obsresult[begin_index:begin_index+quota])                    # (shots,)
                predictions = np.asarray(predictions_result[begin_index:begin_index+quota]).ravel()

                # 3. count mismatches in vectorised form ---------------------------------
                num_errors = np.count_nonzero(observables != predictions)

                self._subspace_LE_count[w]+=num_errors
                self._estimated_subspaceLER[w] = self._subspace_LE_count[w] / self._subspace_sample_used[w]


                print(f"Logical error rate when w={w}: {self._estimated_subspaceLER[w]*binomial_weight(self._num_noise, w,self._error_rate):.6g}")

                begin_index+=quota
            print(self._subspace_LE_count)
            print(self._subspace_sample_used)


        total_ler_count=np.sum(list(self._subspace_LE_count.values()))
        if(total_ler_count>0.8*MIN_NUM_LE_EVENT*self._num_subspace):
            self._stratified_succeed=True
        else:
            self._stratified_succeed=False

        
        print("Samples used:{}".format(self._sample_used))
        print("circuit level code distance:{}".format(self._circuit_level_code_distance))
        print(self._estimated_wlist)
        print(self._estimated_subspaceLER)



    def subspace_sampling_to_fit_curve(self,sampleBudget):

        """
        After we determine the minw and maxw, we generate an even distribution of points 
        between minw and maxw.

        The goal is for the curve fitting in the next step to get more accurate.
        """

        wlist=evenly_spaced_ints(self._minw,self._saturatew,self._num_subspace)
        for weight in wlist:
            if not (weight in self._estimated_wlist):
                self._estimated_wlist.append(weight)
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


    def fit_linear(self):
        
        x_list = [x for x in self._estimated_subspaceLER.keys() if 0 < self._estimated_subspaceLER[x] < 0.5]
        y_list = [np.log(0.5/self._estimated_subspaceLER[x]-1) for x in x_list]
        #x_list=np.log(x_list)  # Take the log of the weights
        popt, pcov = curve_fit(
            modified_linear_function, 
            x_list, 
            y_list
        )

        self._codedistance = 0
        # Extract the best-fit parameter (alpha)
        a,b,c,d= popt[0] , popt[1], popt[2], popt[3]


        #Plot the fitted line
        x_fit = np.linspace(min(x_list), max(x_list), 1000)

        y_fit = modified_linear_function(x_fit, a, b,c,d)
        plt.figure()
        plt.scatter(x_list, y_list, label='Data points', color='blue')
        plt.plot(x_fit, y_fit, label='Fitted line', color='orange')
        plt.xlabel('Weight')
        plt.ylabel('Log(0.5 / LER - 1)')
        plt.title('Linear Fit of S-curve')
        plt.legend()
        plt.savefig("linear_fit.png", dpi=300)
        plt.close()



    '''
    Fit the distribution by 1/2-e^{alpha/W}
    '''
    def fit_Scurve(self):
        # Initial guess for alpha
        # sigma_guess = self._saturatew / 3.2898
        # mu_guess    = self._saturatew/2          # centre in the middle of that span
        # initial_guess  = (0,mu_guess, sigma_guess)
        if self._stratified_succeed:
            self._saturatew=self._maxw

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
        if not self._stratified_succeed:
            sigma=int(np.sqrt(self._error_rate*(1-self._error_rate)*self._num_noise))
            if sigma==0:
                sigma=1
            ep=int(self._error_rate*self._num_noise)
            self._minw=max(1,ep-5*sigma)
            self._maxw=max(2,ep+5*sigma)

        for weight in range(self._minw,self._maxw+1):
            """
            If the weight is in the estimated list, we use the estimated value
            Else, we use the curve fitting value

            If the weight is less than the minw, we just declare it as 0
            """
            if weight in self._estimated_subspaceLER.keys():
                self._LER+=self._estimated_subspaceLER[weight]*binomial_weight(self._num_noise, weight,self._error_rate)
            else:
                self._LER+=scurve_function(weight,self._mu,self._sigma)*binomial_weight(self._num_noise,weight,self._error_rate)
            #self._LER+=scurve_function(weight,self._mu,self._sigma)*binomial_weight(self._num_noise,weight,self._error_rate)
        return self._LER




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
    p=0.0005
    sample_buget=1000000
    tmp=stratified_Scurve_LERcalc(p,sampleBudget=sample_buget,num_subspace=20)

    #filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/surface/surface3"
    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/surface/surface9"
    figname="S5.png"
    titlename="Surface5"
    tmp.parse_from_file(filepath)
    tmp.determine_saturated_w()
    tmp.subspace_sampling_to_fit_curve(sample_buget)
    # k=3
    # while not tmp._stratified_succeed:
    #     print("------------------------k=",k,"------------------------")
    #     tmp.determine_range_to_sample(k)
    #     tmp.subspace_sampling()
    #     k=k*2
    #     # 
    #     # tmp.subspace_sampling_to_fit_curve(sample_buget)

    tmp.fit_linear()
    #tmp.fit_Scurve()
    #tmp.calc_logical_error_rate_after_curve_fitting()
    #print("Final LER: ",tmp._LER)
    #tmp.plot_scurve(figname,titlename)


    #tmp.calculate_R_square_score()

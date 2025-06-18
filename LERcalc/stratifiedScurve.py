
from QEPG.QEPG import return_samples,return_samples_many_weights,return_detector_matrix,return_samples_many_weights_numpy,return_samples_many_weights_separate_obs, compile_QEPG, return_samples_many_weights_separate_obs_with_QEPG, return_samples_with_fixed_QEPG
from LERcalc.clifford import *
import math
import pymatching
from scipy.optimize import curve_fit
from scipy.stats import norm
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from scipy.stats import binom
from contextlib import redirect_stdout

import warnings
from scipy.optimize import OptimizeWarning

# Suppress only OptimizeWarning
#warnings.filterwarnings("ignore", category=OptimizeWarning)
# def binomial_weight(N, W, p):
#     if N<200:
#         return math.comb(N, W) * (p**W) * ((1 - p)**(N - W))
#     else:
#         lam = N * p
#         # PMF(X=W) = e^-lam * lam^W / W!
#         # Evaluate in logs to avoid overflow for large W, then exponentiate
#         log_pmf = (-lam) + W*math.log(lam) - math.lgamma(W+1)
#         return math.exp(log_pmf)


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


def binomial_weight(N, W, p):
    return binom.pmf(W, N, p)
    #return math.comb(N, W) * (p**W) * ((1 - p)**(N - W))

def linear_function(x, a, b):
    """
    Linear function for curve fitting.
    """
    return a * x + b


# def modified_linear_function(x, a, b,c,d):
#     """
#     Linear function for curve fitting.
#     """
#     return a * x + b+c/(x-d)


def modified_linear_function_with_d(x, a, b, c, d):
    eps   = 1e-12
    delta = (x - d)**0.5
    delta = np.where(np.abs(delta) < eps, np.sign(delta)*eps, delta)
    return a * x + b + c / delta



# Strategy A: keep the model safe near the pole
def modified_linear_function(d):
    def tempfunc(x,a,b,c,d=d):
        return modified_linear_function_with_d(x, a, b, c, d)
    return tempfunc



def modified_sigmoid_function(x, a, b,c,d):
    """
    Modified sigmoid function for curve fitting.
    This function is used to fit the S-curve.
    """
    z = a*x + b + c/((x - d)**0.5)
    # ignore overflows in exp → exp(z) becomes np.inf, so 0.5/(1+inf) = 0.0
    with np.errstate(over='ignore'):
        y = 0.5 / (1 + np.exp(z))
    return y

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



# Redefine turning point where the 2nd term is still significant in dy/dw
def refined_sweat_spot(alpha, beta, t, ratio=2):
    # We define turning point by solving: 1/alpha = ratio * (1/2) * beta / (w - t)^{3/2}
    # => (w - t)^{3/2} = (ratio * beta * alpha) / 2
    # => w = t + [(ratio * beta * alpha / 2)]^{2/3}
    return t + ((ratio * beta * alpha / 2) ** (2 / 3))






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

    def __init__(self, error_rate=0, sampleBudget=10000, k_range=3, num_subspace=5,beta=4):
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
        self._t=1
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
        self._has_logical_errorw=0
        self._estimated_wlist=[]

        self._stim_str_after_rewrite=""

        self._mu=0
        self._sigma=0

        #In the area we are interested in, the maximum value of the logical error rate
        self._rough_value_for_subspace_LER=0

        self._stratified_succeed=False

        self._k_range=k_range
        self._QEPG_graph=None

        self._R_square_score=0
        self._beta=beta

        self._sweat_spot=None


        self._MIN_NUM_LE_EVENT = 100
        self._SAMPLE_GAP=100
        self._MAX_SAMPLE_GAP=1000000
        self._MAX_SUBSPACE_SAMPLE=5000000


    def set_sample_bound(self, MIN_NUM_LE_EVENT,SAMPLE_GAP, MAX_SAMPLE_GAP, MAX_SUBSPACE_SAMPLE):
        """
        Set the sample bound for the subspace sampling
        """
        self._MIN_NUM_LE_EVENT=MIN_NUM_LE_EVENT
        self._SAMPLE_GAP=SAMPLE_GAP
        self._MAX_SAMPLE_GAP=MAX_SAMPLE_GAP
        self._MAX_SUBSPACE_SAMPLE=MAX_SUBSPACE_SAMPLE


    def clear_all(self):
        self._estimated_subspaceLER={}
        self._subspace_LE_count={}
        self._estimated_subspaceLER_second={}
        self._subspace_sample_used={}
        self._sample_used=0
        self._LER=0        
        self._estimated_wlist=[]
        self._saturatew=10000000000000               
        self._minw=self._t+1
        self._maxw=10000000000000
        self._cliffordcircuit=CliffordCircuit(4)  
        self._R_square_score=0


    def calc_logical_error_rate_with_fixed_w(self, shots, w):
        """
        Calculate the logical error rate with fixed w
        TODO: Optimize this function
        """
        result= return_samples_with_fixed_QEPG(self._QEPG_graph,w,shots)
        self._sample_used+=shots
        if w not in self._subspace_LE_count.keys():
            self._subspace_LE_count[w]=0
            self._subspace_sample_used[w]=shots
        else:
            self._subspace_sample_used[w]+=shots
        arr=np.asarray(result)
        states=arr[:,:-1]
        observables=arr[:,-1]
        predictions =np.squeeze(self._matcher.decode_batch(states))
        num_errors = np.count_nonzero(observables != predictions)
        self._subspace_LE_count[w]+=num_errors
        return num_errors/shots


    '''
    Use binary search to determine the exact number of errors 
    that give saturate logical error rate
    We just try 10 samples

    TODO: Restructure the function.
    Add the threshold as an input parameter.
    '''
    def binary_search_upper(self,low,high, shots):
        left=low
        right=high
        epsion=0.3
        while left<right:
            mid=(left+right)//2
            er=self.calc_logical_error_rate_with_fixed_w(shots,mid)
            if er>epsion:
                right=mid
            else:
                left=mid+1
        return left


    def binary_search_lower(self,low,high, shots=10000):
        left=low
        right=high
        epsion=3e-3
        while left<right:
            print("left: ",left)
            print("right: ",right)
            mid=(left+right)//2
            er=self.calc_logical_error_rate_with_fixed_w(shots,mid)
            if er>epsion:
                right=mid
            else:
                left=mid+1
        return left


    def determine_lower_w(self):
        self._has_logical_errorw=self.binary_search_lower(self._t+1,self._num_noise//10)


    def determine_saturated_w(self,shots=100):
        """
        Use binary search to determine the minw and maxw
        """
        #self._saturatew=self._num_detector//30
        self._saturatew=self.binary_search_upper(self._minw,self._num_noise//10,shots)
        #print("Self._saturatew: ",self._saturatew)



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

        self._QEPG_graph=compile_QEPG(stim_str)

    def determine_range_to_sample(self):
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
        self._minw=max(self._t+1,ep-self._k_range*sigma)
        self._maxw=max(self._num_subspace,ep+self._k_range*sigma)
        self._maxw=min(self._maxw,self._num_noise)
        # print("Needed Minw:    ")
        # print(self._minw)
        # print("Needed Maxw:    ")
        # print(self._maxw)
        # print("Average number of noise:    ")
        # print(ep)



    def subspace_sampling(self):
        """
        wlist store the subset of weights we need to sample and get
        correct logical error rate.
        
        In each subspace, we stop sampling until 100 logical error events are detected, or we hit the total budget.
        """
        #wlist_need_to_sample = list(range(self._minw, self._maxw + 1))
        wlist_need_to_sample=evenly_spaced_ints(self._sweat_spot,self._saturatew,self._num_subspace)
        #print("wlist_need_to_sample: ",wlist_need_to_sample)
        for weight in wlist_need_to_sample:
            if not weight in self._estimated_wlist:
                self._estimated_wlist.append(weight)
                self._subspace_LE_count[weight]=0
                self._subspace_sample_used[weight]=0

        #print(wlist_need_to_sample)
        self._sample_used=0
        total_LE_count=0
        while True:
            x_list = [x for x in self._estimated_subspaceLER.keys() if (self._estimated_subspaceLER[x] < 0.5 and self._estimated_subspaceLER[x]>0)]

            slist=[]
            wlist=[]
            """
            Case 1 to end the while loop: We have consumed all of our sample budgets
            """
            if(self._sample_used>self._sampleBudget):
                break

            for weight in wlist_need_to_sample:
                """
                When we declare the circuit level code distance, we don't need to sample these subspaces
                """
                if(self._subspace_sample_used[weight]>self._MAX_SUBSPACE_SAMPLE):
                    continue


                if(self._subspace_LE_count[weight]<self._MIN_NUM_LE_EVENT):
                    if(self._subspace_LE_count[weight]>=1):
                        """
                        For larger subspaces, when we have already get some logical error, 
                        we can estimate how many we still need to sample
                        """
                        sample_num_required=int(self._MIN_NUM_LE_EVENT/self._subspace_LE_count[weight])* self._subspace_sample_used[weight]
                        if sample_num_required>self._MAX_SAMPLE_GAP:
                            sample_num_required=self._MAX_SAMPLE_GAP
                        slist.append(sample_num_required)
                        self._subspace_sample_used[weight]+=sample_num_required  
                        self._sample_used+=sample_num_required
                    else:        
                        """
                        For larger subspaces, if we have not get any logical error, then we double the sample size
                        """
                        sample_num_required=max(self._SAMPLE_GAP,self._subspace_sample_used[weight]*10)
                        if sample_num_required>self._MAX_SAMPLE_GAP:
                            sample_num_required=self._MAX_SAMPLE_GAP
                        slist.append(sample_num_required)
                        self._subspace_sample_used[weight]+=sample_num_required
                        self._sample_used+=sample_num_required
                    wlist.append(weight)
            """
            Case 2 to end the while loop: We have get 100 logical error events for all these subspaces
            """
            if(len(wlist)==0):
                break

            print("wlist: ",wlist)
            print("slist: ",slist)
            #detector_result,obsresult=return_samples_many_weights_separate_obs(self._stim_str_after_rewrite,wlist,slist)
            detector_result,obsresult=return_samples_many_weights_separate_obs_with_QEPG(self._QEPG_graph,wlist,slist)
            predictions_result = self._matcher.decode_batch(detector_result)
            # print("Result get!")

            
            begin_index=0
            for w_idx, (w, quota) in enumerate(zip(wlist, slist)):

                observables =  np.asarray(obsresult[begin_index:begin_index+quota])                    # (shots,)
                predictions = np.asarray(predictions_result[begin_index:begin_index+quota]).ravel()

                # 3. count mismatches in vectorised form ---------------------------------
                num_errors = np.count_nonzero(observables != predictions)
                total_LE_count+=num_errors
                self._subspace_LE_count[w]+=num_errors
                self._estimated_subspaceLER[w] = self._subspace_LE_count[w] / self._subspace_sample_used[w]


                # print(f"Logical error rate when w={w}: {self._estimated_subspaceLER[w]*binomial_weight(self._num_noise, w,self._error_rate):.6g}")

                begin_index+=quota
            print(self._subspace_LE_count)
            print("Subspace LE count: ",self._subspace_LE_count)
            print("self._subspace_sample_used: ",self._subspace_sample_used)

        
        # print("Samples used:{}".format(self._sample_used))
        #print("circuit level code distance:{}".format(self._circuit_level_code_distance))
        #print(self._subspace_LE_count)



    def subspace_sampling_to_fit_curve(self,sampleBudget):

        """
        After we determine the minw and maxw, we generate an even distribution of points 
        between minw and maxw.

        The goal is for the curve fitting in the next step to get more accurate.
        """

        wlist=evenly_spaced_ints(self._has_logical_errorw,self._saturatew,self._num_subspace)
        for weight in wlist:
            if not (weight in self._estimated_wlist):
                self._estimated_wlist.append(weight)
        slist=[sampleBudget//self._num_subspace]*len(wlist)


        #detector_result,obsresult=return_samples_many_weights_separate_obs(self._stim_str_after_rewrite,wlist,slist)
        detector_result,obsresult=return_samples_many_weights_separate_obs_with_QEPG(self._QEPG_graph,wlist,slist)
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
        # print("---------------After first round sampling-------------------")
        # print("Subspace LE count: ",self._subspace_LE_count)
        # print("self._subspace_sample_used: ",self._subspace_sample_used)
        # print("self._estimated_subspaceLER: ",self._estimated_subspaceLER)


    def calculate_R_square_score(self):
        y_observed = [self._estimated_subspaceLER[x] for x in self._estimated_wlist]
        y_predicted = [scurve_function(x,self._mu,self._sigma) for x in self._estimated_wlist]
        #y_predicted = [scurve_function_with_distance(x,self._mu,self._sigma) for x in self._estimated_wlist]
        r2 = r_squared(y_observed, y_predicted)
        #print("R^2 score: ", r2)
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




    def fit_linear_area(self):
        x_list = [x for x in self._estimated_subspaceLER.keys() if (self._estimated_subspaceLER[x] < 0.5 and self._estimated_subspaceLER[x]>0)]
        y_list = [np.log(0.5/self._estimated_subspaceLER[x]-1) for x in x_list]

        popt, pcov = curve_fit(
            linear_function,
            x_list,
            y_list
        )


        sigma=int(np.sqrt(self._error_rate*(1-self._error_rate)*self._num_noise))
        if sigma==0:
            sigma=1
        ep=int(self._error_rate*self._num_noise)
        self._minw=max(self._t+1,ep-self._k_range*sigma)
        self._maxw=max(2,ep+self._k_range*sigma)



        self._a,self._b= popt[0] , popt[1]

        #Plot the fitted line
        x_fit = np.linspace(min(x_list), max(x_list), 1000)

        y_fit = linear_function(x_fit, self._a, self._b)

        #print("Fitted parameters: a={}, b={}".format(self._a,self._b))
        plt.figure()
        plt.plot(x_fit, y_fit, label='Fitted line', color='orange')
        plt.scatter(x_list, y_list, label='Data points', color='blue')


        # ── NEW: highlight linear-fit window ───────────────────────────────────
        plt.axvline(self._minw, color='red',  linestyle='--', linewidth=1.2, label=r'$w_{\min}$')
        plt.axvline(self._maxw, color='green', linestyle='--', linewidth=1.2, label=r'$w_{\max}$')
        plt.axvspan(self._minw, self._maxw, color='gray', alpha=0.15)  # translucent fill
        # ───────────────────────────────────────────────────────────────────────


        plt.xlabel('Weight')
        plt.ylabel('Linear')
        plt.title('Linear Fit of S-curve')
        plt.legend()
        plt.savefig("total_linear_fit.png", dpi=300)
        plt.close()


    def fit_log_S_model(self,filename):
        
        #print("circuit d:",self._circuit_level_code_distance)
        x_list = [x for x in self._estimated_subspaceLER.keys() if (self._estimated_subspaceLER[x] < 0.5 and self._estimated_subspaceLER[x]>0)]

        y_list = [np.log(0.5/self._estimated_subspaceLER[x]-1) for x in x_list]
        #x_list=np.log(x_list)  # Take the log of the weights

        print("x_list: ",x_list)
        print("y_list: ",y_list)


        non_zero_indices=[x for x in x_list if self._estimated_subspaceLER[x]>0]
        upper_bound_code_distance=min(non_zero_indices) if len(non_zero_indices)>0 else self._circuit_level_code_distance*10


        #print("Upper bound code distance: ",upper_bound_code_distance)


        center = self._saturatew /2 
        sigma    = self._saturatew/7          # centre in the middle of that span
        b=self._b
        a=self._a
        c=self._beta
        initial_guess  = (a, b ,c)

        #print("Initial guess d: ",int((self._circuit_level_code_distance+upper_bound_code_distance)/2))
        sigma=int(np.sqrt(self._error_rate*(1-self._error_rate)*self._num_noise))
        if sigma==0:
            sigma=1
        ep=int(self._error_rate*self._num_noise)
        self._minw=max(self._t+1,ep-self._k_range*sigma)
        self._maxw=max(2,ep+self._k_range*sigma)


        self._num_detector
        self._num_noise

        # ── lower bounds for [param1, param2, param3, param4]
        lower = [ min(self._a*0.5,self._a*1.5), min(self._b*0.5,self._b*1.5),  0]

        # ── upper bounds for [param1, param2, param3, param4]
        upper = [ max(self._a*0.5,self._a*1.5), max(self._b*0.5,self._b*1.5)  , np.inf]



        popt, pcov = curve_fit(
            modified_linear_function(self._t),
            x_list,
            y_list,
            p0=initial_guess,          # len(initial_guess) must be 4 and within the bounds above
            bounds=(lower, upper),     # <-- tuple with two arrays
            maxfev=50_000              # or max_nfev in newer SciPy
        )

        self._codedistance = 0
        # Extract the best-fit parameter (alpha)
        self._a,self._b,self._c= popt[0] , popt[1], popt[2]


        #print("circuit d:",self._circuit_level_code_distance)
        y_list = [np.log(0.5/self._estimated_subspaceLER[x]-1) for x in x_list]
        y_predicted = [modified_linear_function_with_d(x,self._a,self._b,self._c,self._t) for x in x_list]
        #y_predicted = [scurve_function_with_distance(x,self._mu,self._sigma) for x in self._estimated_wlist]
        self._R_square_score = r_squared(y_list, y_predicted)
        #print("R^2 score: ", self._R_square_score)

        #Plot the fitted line
        x_fit = np.linspace(self._t+1, max(x_list), 1000)

        y_fit = modified_linear_function_with_d(x_fit, self._a, self._b,self._c,self._t)

        self.calc_logical_error_rate_after_curve_fitting()
        
        alpha= -1/self._a
        self._sweat_spot = refined_sweat_spot(alpha, self._c, self._t, ratio=0.2)
        sweat_spot_y = modified_linear_function_with_d(self._sweat_spot, self._a, self._b, self._c, self._t)

        #print("Fitted parameters: a={}, b={}, c={}, d={}".format(self._a, self._b, self._c, self._d))
        plt.figure()
        plt.scatter(x_list, y_list, label='Data points', color='blue')
        plt.plot(x_fit, y_fit, label=f'Fitted line, R2={self._R_square_score}', color='orange')

        # ── NEW: highlight sweat spot ────────────────────────────────────────
       
        plt.scatter(self._sweat_spot, sweat_spot_y, color='purple', marker='o', s=50, label='Sweat Spot')

        # ground_x=self._ground_estimated_subspaceLER.keys()
        # ground_y=self._ground_estimated_subspaceLER.values()
        # ground_x=[x for x in ground_x if self._ground_estimated_subspaceLER[x]>0]
        # ground_y=[np.log(0.5/x-1) for x in ground_y if x > 0]
        # plt.bar(ground_x, ground_y, width=0.4, color='red', alpha=0.6, label='Ground truth')

        # ── NEW: highlight linear-fit window ───────────────────────────────────
        plt.axvline(self._minw, color='red',  linestyle='--', linewidth=1.2, label=r'$w_{\min}$')
        plt.axvline(self._maxw, color='green', linestyle='--', linewidth=1.2, label=r'$w_{\max}$')
        plt.axvspan(self._minw, self._maxw, color='gray', alpha=0.15)  # translucent fill
        # ───────────────────────────────────────────────────────────────────────


        exp_str = "{0:.2e}".format(self._LER)
        base, exp = exp_str.split('e')
        exp = int(exp)  # will be negative or positive

        pl_formatted = r'$P_L={0}\times 10^{{{1}}}$'.format(base, exp)

        textstr = '\n'.join((
            r'$\alpha=%.4f$' % alpha,
            r'$\mu =%.4f$' % (alpha*self._b),
            r'$\beta=%.4f$' % self._c,
            r'$w_{\min}=%d$' % self._minw,
            r'$w_{\max}=%d$' % self._maxw,
            r'$w_{sweet}=%d$' % self._sweat_spot,         
            r'$\#\mathrm{detector}=%d$' % self._num_detector,
            r'$\#\mathrm{noise}=%d$' % self._num_noise,
            pl_formatted
        ))

        fig = plt.gcf()
        fig.subplots_adjust(right=0.75)  # Make room on the right
        fig.text(0.78, 0.5, textstr,
                fontsize=10,
                va='center', ha='left',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.95))
        # ──────────────────────────────────────────────────────────────────

        plt.xlabel('Weight')
        plt.ylabel(r'$\log\left(\frac{0.5}{\mathrm{LER}} - 1\right)$')
        plt.title('Linear Fit of S-curve')
        plt.legend()
        plt.tight_layout()
        plt.savefig(filename, dpi=300)
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


    def ground_truth_subspace_sampling(self):
        """
        Sample around the subspaces.
        This is the ground truth value to test the accuracy of the curve fitting.
        """
        sigma=int(np.sqrt(self._error_rate*(1-self._error_rate)*self._num_noise))
        if sigma==0:
            sigma=1
        ep=int(self._error_rate*self._num_noise)
        minw=max(self._t+1,ep-self._k_range*sigma)
        maxw=max(self._num_subspace,ep+self._k_range*sigma)
        maxw=min(maxw,self._num_noise)
        wlist_need_to_sample = list(range(minw, maxw + 1))
        self._ground_sample_used=0
        self._ground_estimated_subspaceLER={}
        self._ground_subspace_LE_count={}
        self._ground_subspace_sample_used={}
        for weight in wlist_need_to_sample:
            self._ground_subspace_LE_count[weight]=0
            self._ground_subspace_sample_used[weight]=0        

        while True:
            slist=[]
            wlist=[]        

            if(self._ground_sample_used>self._sampleBudget):
                break
            
            for weight in wlist_need_to_sample:
                if(self._ground_subspace_sample_used[weight]>self._MAX_SUBSPACE_SAMPLE):
                    continue
                if(self._ground_subspace_LE_count[weight]==0):
                    wlist.append(weight)
                    sample_num_required=max(self._MAX_SAMPLE_GAP,self._ground_subspace_sample_used[weight]*10)
                    self._ground_subspace_sample_used[weight]+=sample_num_required  
                    self._ground_sample_used+=sample_num_required
                    slist.append(sample_num_required)
                    continue
                if(self._ground_subspace_LE_count[weight]<self._MIN_NUM_LE_EVENT):
                    sample_num_required=int(self._MIN_NUM_LE_EVENT/self._ground_subspace_LE_count[weight])* self._ground_subspace_sample_used[weight]
                    if sample_num_required>self._MAX_SAMPLE_GAP:
                        sample_num_required=self._MAX_SAMPLE_GAP
                    self._ground_subspace_sample_used[weight]+=sample_num_required  
                    self._ground_sample_used+=sample_num_required
                    wlist.append(weight)
                    slist.append(sample_num_required)

            if(len(wlist)==0):
                break
            #detector_result,obsresult=return_samples_many_weights_separate_obs(self._stim_str_after_rewrite,wlist,slist)
            
            print("Ground truth wlist: ",wlist)
            print("Ground truth slist: ",slist)
            
            detector_result,obsresult=return_samples_many_weights_separate_obs_with_QEPG(self._QEPG_graph,wlist,slist)
            predictions_result = self._matcher.decode_batch(detector_result)

            
            begin_index=0
            for w_idx, (w, quota) in enumerate(zip(wlist, slist)):

                observables =  np.asarray(obsresult[begin_index:begin_index+quota])                    # (shots,)
                predictions = np.asarray(predictions_result[begin_index:begin_index+quota]).ravel()

                # 3. count mismatches in vectorised form ---------------------------------
                num_errors = np.count_nonzero(observables != predictions)

                self._ground_subspace_LE_count[w]+=num_errors
                self._ground_estimated_subspaceLER[w] = self._ground_subspace_LE_count[w] / self._ground_subspace_sample_used[w]


                print(f"Logical error rate when w={w}: {self._ground_estimated_subspaceLER[w]*binomial_weight(self._num_noise, w,self._error_rate):.6g}")

                begin_index+=quota
            print(self._ground_subspace_LE_count)
            print(self._ground_subspace_sample_used)
        print("Samples used:{}".format(self._ground_sample_used))




    def calc_logical_error_rate_after_curve_fitting(self):
        #self.fit_Scurve()
        self._LER=0

        sigma=int(np.sqrt(self._error_rate*(1-self._error_rate)*self._num_noise))
        if sigma==0:
            sigma=1
        ep=int(self._error_rate*self._num_noise)
        self._minw=max(self._t+1,ep-self._k_range*sigma)
        self._maxw=max(2,ep+self._k_range*sigma)

        for weight in range(self._minw,self._maxw+1):
            """
            If the weight is in the estimated list, we use the estimated value
            Else, we use the curve fitting value

            If the weight is less than the minw, we just declare it as 0
            """
            if weight in self._estimated_subspaceLER.keys():
                self._LER+=self._estimated_subspaceLER[weight]*binomial_weight(self._num_noise, weight,self._error_rate)
                #print("Weight: ",weight," LER: ",self._estimated_subspaceLER[weight]*binomial_weight(self._num_noise, weight,self._error_rate))
            else:
                fitted_subspace_LER=modified_sigmoid_function(weight,self._a,self._b,self._c,self._t)
                self._LER+=fitted_subspace_LER*binomial_weight(self._num_noise,weight,self._error_rate)
                #print("Weight: ",weight," LER: ",fitted_subspace_LER*binomial_weight(self._num_noise, weight,self._error_rate))
            #self._LER+=scurve_function(weight,self._mu,self._sigma)*binomial_weight(self._num_noise,weight,self._error_rate)
        return self._LER




    def plot_scurve(self, filename,title="S-curve"):
        """Plot the S-curve and its discrete estimate."""
        keys   = list(self._estimated_subspaceLER.keys())
        values = [self._estimated_subspaceLER[k] for k in keys]

        fig=plt.figure()
        # bars ── discrete estimate
        plt.bar(keys, values,
                color='tab:blue',         # pick any color you like
                alpha=0.8,
                label='Estimated subspace LER by sampling')
    
        plt.axvline(x=self._error_rate*self._num_noise, color="red", linestyle="--", linewidth=1.2, label="Average Error number") # vertical line at x=0.5

        #print("Saturate w: ",self._saturatew)

        # smooth curve ── fitted S-curve
        x = np.linspace(1.001, self._saturatew, 1000)
        #x = np.linspace(self._d*1.01, self._saturatew, 1000)
        y = modified_sigmoid_function(x, self._a, self._b,self._c,self._t)

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


    def set_t(self, t):
        """
        Set the t value for the S-curve fitting.
        This is used to determine the range of subspace we need to sample.
        """
        self._t = t


    def calculate_LER_from_file(self,filepath,pvalue,codedistance,figname,titlename, repeat=1):
        self._error_rate=pvalue
        self._circuit_level_code_distance=codedistance
        ler_list=[]
        sample_used_list=[]
        r_squared_list=[]
        Nerror_list=[]
        time_list=[]
        for i in range(repeat):
            self.clear_all()
            self.parse_from_file(filepath)
            start = time.time()
            #self.determine_range_to_sample()
            #self.subspace_sampling()
            self.determine_lower_w()

            #self.ground_truth_subspace_sampling()
            #self._has_logical_errorw=self._t+50
            self.determine_saturated_w()


            self.subspace_sampling_to_fit_curve(10000*self._num_subspace)
            '''
            Fit the curve first time just to get the estimated sweat spot
            '''
            self.fit_linear_area()
            self.fit_log_S_model(figname+"First.png")
            '''
            Second round of samples
            '''
            self.subspace_sampling()
            self.fit_linear_area()
            self.fit_log_S_model(figname+"Final.png")

            self.calc_logical_error_rate_after_curve_fitting()
            end = time.time()
            time_list.append(end - start)
            self.plot_scurve(figname,titlename)
            r_squared_list.append(self._R_square_score)
            self._sample_used=np.sum(list(self._subspace_sample_used.values()))
            # print("Final LER: ",self._LER)
            # print("Total samples used: ",self._sample_used)
            ler_list.append(self._LER)
            sample_used_list.append(self._sample_used)
            Nerror_list.append(sum(self._subspace_LE_count.values()))

        # Compute means
        self._LER = np.mean(ler_list)
        self._sample_used = np.mean(sample_used_list)

        # Compute standard deviations
        ler_std = np.std(ler_list)
        sample_used_std = np.std(sample_used_list)
        r2_mean = np.mean(r_squared_list)
        r2_std = np.std(r_squared_list)
        Nerror_mean = np.mean(Nerror_list)
        Nerror_std = np.std(Nerror_list)

        time_mean = np.mean(time_list)
        time_std = np.std(time_list)

        # Print with scientific ± formatting
        print("k: ", self._k_range)
        print("beta: ",self._beta)
        print("Subspaces: ", self._num_subspace)
        print("R2: ", format_with_uncertainty(r2_mean, r2_std))
        print("Samples(ours): ", format_with_uncertainty(self._sample_used, sample_used_std))
        print("Time(our): ", format_with_uncertainty(time_mean, time_std))
        print("PL(ours): ", format_with_uncertainty(self._LER, ler_std))
        print("Nerror(ours): ", format_with_uncertainty(Nerror_mean, Nerror_std))




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


import time
p_list=[0.005,0.001,0.0005,0.0001,0.00005]





surface_file_list=["surface/surface"+str(d) for d in [15,17,19,21]]
surface_code_name_list=["surface"+str(d) for d in [15,17,19,21]]
surface_t_list=[7,8,9,10]
surface_k_list=[16,20,24,28]
surface_subspace_list=[8,8,8,8]  # All surface codes have num_subspace=10


hexagon_file_list=["hexagon/hexagon"+str(d) for d in [9,11,13,15,17,19,21]]
hexagon_code_name_list=["hexagon"+str(d) for d in [9,11,13,15,17,19,21]]
hexagon_t_list=[4,5,6,7,8,9,10]  # All hexagon codes have beta=3
hexagon_k_list=[10,12,16,20,24,28,32]
hexagon_subspace_list=[8,8,8,8,8,8,8]  # All surface codes have num_subspace=10

square_file_list=["square/square"+str(d) for d in [9,11,13,15,17,19,21]]
square_code_name_list=["square"+str(d) for d in [9,11,13,15,17,19,21]]
square_t_list=[4,5,6,7,8,9,10]  # All square codes have beta=3
square_k_list=[10,12,16,20,24,28,32]
square_subspace_list=[8,8,8,8,8,8,8]  # All surface codes have num_subspace=10

# repetition_file_list=["repetition/repetition"+str(d) for d in [9,11,13,15,17,19,21]]
# repetition_code_name_list=["repetition9","repetition11","repetition13","repetition15"]
# repetition_t_list=[4,5,6,7]
# repetition_k_list=[8,12,12,12]
# repetition_subspace_list=[12,15,15,15]  # All surface codes have num_subspace=10

all_file_list=surface_file_list+hexagon_file_list+square_file_list
all_k_list=surface_k_list+hexagon_k_list+square_k_list
all_t_list=surface_t_list+hexagon_t_list+square_t_list
all_subspace_list=surface_subspace_list+hexagon_subspace_list+square_subspace_list
all_code_name_list=surface_code_name_list+hexagon_code_name_list+square_code_name_list


if __name__ == "__main__":
    #generate_all_repetition_code_figure()
    #generate_all_square_code_figure()
    #generate_all_hexagon_code_figure()


    # base_dir = "C:/Users/yezhu/Documents/Sampling/stimprograms/"
    # result_dir = "C:/Users/yezhu/Documents/Sampling/myresultTot/"
    # p=0.001
    # index=0
    # for file in all_file_list:
    #     # 1) build your input filename:
    #     stim_path = base_dir + file
    #     # 2) build your output filename:
       
    #     out_fname =  result_dir + all_code_name_list[index]     # e.g. "surface3-result.txt"
    #     # 4) redirect prints for just this file:
    #     with open(out_fname, "w") as outf, redirect_stdout(outf):
    #         print("---------------Processing code: ", all_code_name_list[index], " with p: ", p, "-----------------")
    #         sample_buget=5*10**8
    #         tmp=stratified_Scurve_LERcalc(p,sampleBudget=sample_buget,k_range=all_k_list[index],num_subspace=all_subspace_list[index])
    #         figname="tmp.png"
    #         titlename="tmp"
    #         tmp.set_t(all_t_list[index])
    #         tmp.calculate_LER_from_file(stim_path,p,0,figname,titlename,5)
    #         index+=1


      

    # p_list=[0.005,0.001,0.0005,0.0001,0.00005]
    # code_type_list=["square"]
    # rel_list=["square/square7"]
    # for code_type, rel in zip(code_type_list, rel_list):
    #     for p in p_list:
    #         # rel might be "surface/surface3" or "repetition/repetition5", etc.
    #         stim_path = base_dir+rel
    #         # 3) build your output filename:
    #         out_fname =  result_dir+str(p)+"-"+str(code_type)+"-result.txt"     # e.g. "surface3-result.txt"
    #         # 4) redirect prints for just this file:
    #         with open(out_fname, "w") as outf, redirect_stdout(outf):
    #             print("---------------Processing code type: ", code_type, " with p: ", p, "-----------------")
    #             sample_buget=1000000000
    #             tmp=stratified_Scurve_LERcalc(p,sampleBudget=sample_buget,k_range=5,num_subspace=12,beta=21)
    #             tmp.set_t(3)
    #             figname="tmp.png"
    #             titlename="tmp"
    #             tmp.calculate_LER_from_file(stim_path,p,0,figname,titlename,5)


    # code_type_list=["hexagon"]
    # rel_list=["hexagon/hexagon7"]
    # for code_type, rel in zip(code_type_list, rel_list):
    #     for p in p_list:
    #         # rel might be "surface/surface3" or "repetition/repetition5", etc.
    #         stim_path = base_dir+rel
    #         # 3) build your output filename:
    #         out_fname =  result_dir+str(p)+"-"+str(code_type)+"-result.txt"     # e.g. "surface3-result.txt"
    #         # 4) redirect prints for just this file:
    #         with open(out_fname, "w") as outf, redirect_stdout(outf):
    #             print("---------------Processing code type: ", code_type, " with p: ", p, "-----------------")
    #             sample_buget=1000000000
    #             tmp=stratified_Scurve_LERcalc(p,sampleBudget=sample_buget,k_range=5,num_subspace=12,beta=21)
    #             tmp.set_t(3)
    #             figname="tmp.png"
    #             titlename="tmp"
    #             tmp.calculate_LER_from_file(stim_path,p,0,figname,titlename,5)



    # code_type_list=["square"]
    # rel_list=["square/square7"]
    # for code_type, rel in zip(code_type_list, rel_list):
    #     for p in p_list:
    #         # rel might be "surface/surface3" or "repetition/repetition5", etc.
    #         stim_path = base_dir+rel
    #         # 3) build your output filename:
    #         out_fname =  result_dir+str(p)+"-"+str(code_type)+"-result.txt"     # e.g. "surface3-result.txt"
    #         # 4) redirect prints for just this file:
    #         with open(out_fname, "w") as outf, redirect_stdout(outf):
    #             print("---------------Processing code type: ", code_type, " with p: ", p, "-----------------")
    #             sample_buget=1000000000
    #             tmp=stratified_Scurve_LERcalc(p,sampleBudget=sample_buget,k_range=5,num_subspace=12,beta=12)
    #             tmp.set_t(3)
    #             figname="tmp.png"
    #             titlename="tmp"
    #             tmp.calculate_LER_from_file(stim_path,p,0,figname,titlename,5)



    # code_type_list=["repetition"]
    # rel_list=["repetition/repetition7"]
    # for code_type, rel in zip(code_type_list, rel_list):
    #     for p in p_list:
    #         # rel might be "surface/surface3" or "repetition/repetition5", etc.
    #         stim_path = base_dir+rel
    #         # 3) build your output filename:
    #         out_fname =  result_dir+str(p)+"-"+str(code_type)+"-result.txt"     # e.g. "surface3-result.txt"
    #         # 4) redirect prints for just this file:
    #         with open(out_fname, "w") as outf, redirect_stdout(outf):
    #             print("---------------Processing code type: ", code_type, " with p: ", p, "-----------------")
    #             sample_buget=1000000000
    #             tmp=stratified_Scurve_LERcalc(p,sampleBudget=sample_buget,k_range=5,num_subspace=12,beta=4)
    #             tmp.set_t(3)
    #             figname="tmp.png"
    #             titlename="tmp"
    #             tmp.calculate_LER_from_file(stim_path,p,0,figname,titlename,5)



    p = 0.001
    sample_budget = 100_000_000

    for d in range(3, 28, 2):
        t = (d - 1) // 2

        stim_path = f"C:/Users/yezhu/Documents/Sampling/stimprograms/surface/surface{d}"
        figname = f"Surface{d}"
        titlename = f"Surface{d}"
        output_filename = f"surface{d}.txt"

        tmp = stratified_Scurve_LERcalc(p, sampleBudget=sample_budget, k_range=5, num_subspace=8, beta=4)
        tmp.set_t(t)
        tmp.set_sample_bound(
            MIN_NUM_LE_EVENT=100,
            SAMPLE_GAP=100_000,
            MAX_SAMPLE_GAP=1_000_000,
            MAX_SUBSPACE_SAMPLE=1_000_000
        )

        with open(output_filename, "w") as f:
            with redirect_stdout(f):
                tmp.calculate_LER_from_file(stim_path, p, 0, figname, titlename, 5)


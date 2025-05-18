#Test our stratified algorithm by comparing it with dp algorithm
#This test only work for small scale

from LERcalc.stratifiedLERcalc import stratifiedLERcalc
from LERcalc.stimLER import stimLERcalc
from LERcalc.symbolicLER import symbolicLER



error_rate=0.01
absolute_error=0.05
sample_size=100000
num_subspace=3


def test_by_file_name(filepath):
    symbolic_calculator=symbolicLER(error_rate)
    ground_truth=symbolic_calculator.calculate_LER_from_file(filepath,error_rate)
    print("Exact ground truth: ",ground_truth)


    tmp=stratifiedLERcalc(error_rate,sampleBudget=sample_size,num_subspace=num_subspace)
    tmp.parse_from_file(filepath)
    tmp.subspace_sampling()

    LER=tmp.calculate_LER()
    print(LER)
    num_noise=tmp._num_noise

    for weight in range(1,3):
        print("LER in the subspace {} is {}".format(weight,tmp.get_LER_subspace(weight)))    



if __name__ == "__main__":


    filepath="C:/Users/yezhu/Documents/Sampling/stimprograms/small/simpleMultiObs"
    test_by_file_name(filepath)
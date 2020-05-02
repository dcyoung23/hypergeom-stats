__author__ = 'dcyoung23'

import math
from scipy import stats

'''
ConfidenceInterval = ci
ConfidenceLevel = cl
SampleErrors = X
SampleSize = n
Population = N
SampleErrorRate = phat
'''

def wilson_ci(ci,cl,X,n,N):

    phat = float(X)/n

    alpha = 1-cl
    z = stats.norm.ppf(1-alpha/ci)

    # Break down calculation in components
    a = phat+((z**2)/(2*n))
    b = math.sqrt((((phat*(1-phat))/n)*((N-n)/(float(N)-1)))+(z**2)/(4*n**2))
    c = 1+((z**2)/n)

    lower = (a-(z*b))/c
    upper = (a+(z*b))/c

    return lower, upper

if __name__ == "__main__":

    wilson_lower, wilson_upper = wilson_ci(2,.95,0,351,9999999)
    print(wilson_lower, wilson_upper)
    wilson_mid = (wilson_lower + wilson_upper)/2
    print(wilson_mid)
    print(int(wilson_mid*9999999))
    est_errors = int(wilson_lower*9999999)
    #print est_errors
    inc_factor = int(min(est_errors/2.0,9999999*.005))
    #print inc_factor

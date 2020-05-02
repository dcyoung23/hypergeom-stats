__author__ = 'dcyoung23'

import math
from time import time
from scipy import stats
from wilsonscore_stats import wilson_ci

# Normal approximation to the binomial
def normapprox_ss(ci,cl,N,p,B):

    alpha = 1-cl
    z = stats.norm.ppf(1-alpha/ci)

    n = (N*(z**2)*(p*(1-p)))/(((N-1)*(B**2))+((z**2)*(p*(1-p))))
    n = int(round(n,0))

    return n

# Increment factor update
def inc_factor_update(loop_change,inc_factor):
    if loop_change == 1:
        if abs(inc_factor) == 1:
            inc_factor = inc_factor*-1
        else:
            # Make sure to apply ceiling when increment factor is positive
            if inc_factor < 0:
                inc_factor = int(math.ceil((inc_factor/2.0)*-1))
            else:
                inc_factor = int((math.ceil(inc_factor/2.0))*-1)
    return inc_factor

# Confidence interval loop variables
def loop_ci(ci_bounds,fract,cdf,inc_ci_factor):
    # Condition for loop exit
    ci_bounds = ci_bounds[:1].upper()
    if ci_bounds == "L" and fract <= cdf and abs(inc_ci_factor) == 1:
        loop_ci_cond = 1
    elif ci_bounds == "U" and fract >= cdf and abs(inc_ci_factor) == 1:
        loop_ci_cond = 1
    else:
        loop_ci_cond = 0

    # Condition for loop direction change
    if fract <= cdf and inc_ci_factor < 0:
        loop_ci_change = 1
    elif fract >= cdf and inc_ci_factor > 0:
        loop_ci_change = 1
    else:
        loop_ci_change = 0

    return loop_ci_cond, loop_ci_change

# Hypergeometric confidence limits
def hypergeom_ci(ci_bounds,ci,cl,X,n,N):
    # Increment adjustment factor
    inc_adj = .01
    # Set Wilson Score Interval Lower/Upper bounds
    wilson_lower, wilson_upper = wilson_ci(ci,cl,X,n,N)
    wilson_mid = (wilson_lower + wilson_upper)/2

    ci_bounds = ci_bounds[:1].upper()
    alpha = 1-cl
    
    if ci_bounds == "L":
        wilson_cl = wilson_lower
        fract = 1-round(alpha/ci,3)
        inc_ci_factorinit = -1
    elif ci_bounds == "U":
        wilson_cl = wilson_upper
        fract = round(alpha/ci,3)
        inc_ci_factorinit = 1

    # Use Wilson Score Interval for starting point for estimated errors
    est_errors = int(wilson_cl*N)

    # Set initial increment factor
    inc_ci_factor = int(N*inc_adj)
    
    # Set initial direction
    inc_ci_factor = inc_ci_factor * inc_ci_factorinit

    # Set initial Hypergeometric CDF
    cdf = stats.hypergeom.cdf(X,N,est_errors,n)
    loop_ci_cond, loop_ci_change = loop_ci(ci_bounds,fract,cdf,inc_ci_factor)
        
    conf_limit = float(est_errors)/N
    previnc_ci_factor = inc_ci_factor
    loop_ci_cnt = 0

    all_ci_results = []

    while loop_ci_cond == 0:
        ci_results = {}
        
        # Infinite loop break
        if loop_ci_cnt > 100:
            break

        # Update increment factor
        inc_ci_factor = inc_factor_update(loop_ci_change,inc_ci_factor)

        # Update new estimated errors
        est_errors = est_errors + inc_ci_factor

        # Ensure estimated errors does not increment beyond 0 and population
        if est_errors < 0:
            est_errors = 0
        elif est_errors > N:
            est_errors = N
            
        cdf = stats.hypergeom.cdf(X,N,est_errors,n)
        conf_limit = float(est_errors)/N
        previnc_ci_factor = inc_ci_factor
        loop_ci_cond, loop_ci_change = loop_ci(ci_bounds,fract,cdf,inc_ci_factor)
        loop_ci_cnt += 1

        # Capture results
        ci_results["bounds"] = ci_bounds
        ci_results["est_errors"] = est_errors
        ci_results["inc_ci_factor"] = inc_ci_factor
        ci_results["cdf"] = cdf
        ci_results["conf_limit"] = conf_limit
        ci_results["loop_cond"] = loop_ci_cond
        ci_results["loop_change"] = loop_ci_change
        ci_results["loop_cnt"] = loop_ci_cnt

        # Append loop to results
        all_ci_results.append(ci_results)

    return ci_results, all_ci_results

# Sample Size loop variables
def loop_ss(inc_ss_factor,prec_desired,prec_exact,true_ss_cnt):
    # Condition for sample size condition
    if prec_exact <= prec_desired:
        loop_ss_cond = 1
    else:
        loop_ss_cond = 0

    # Condition for loop direction change
    if inc_ss_factor == 1:
        loop_ss_change = 0
    elif prec_exact > prec_desired and inc_ss_factor > 1:
        loop_ss_change = 0
    elif prec_exact <= prec_desired and inc_ss_factor < 0:
        loop_ss_change = 0
    else:
        loop_ss_change = 1

    # Condition for loop exit - Loop for at least 10 iterations to avoid stopping early
    if loop_ss_cond == 1 and true_ss_cnt >= 10 and inc_ss_factor == 1:
        loop_ss_exit = 1
    else:
        loop_ss_exit = 0

    return loop_ss_cond, loop_ss_change, loop_ss_exit

# Hypergeometric sample size
def hypergeom_ss(ss_bounds,cl,N,p,prec_desired):
    # Set rounding digits
    r = 6

    # Set confidence interval (2 or 1 tail)
    ss_bounds = ss_bounds[:1].upper()
    if ss_bounds == "B":
        ci = 2
    else:
        ci = 1

    # Use normal approximation for binomial as starting point
    n = normapprox_ss(ci,cl,N,p,prec_desired)

    # Set final sample size from starting n
    n_final = n

    # Set starting loop count and increment factor
    loop_ss_cnt = 0
    inc_ss_factor = 50

    # True sample size count will count the number of consecutive loops that condition is true
    true_ss_cnt = 0

    # Initiate variables
    loop_ss_cond = 0
    loop_ss_change = 0
    loop_ss_exit = 0

    # All results list
    all_ss_results = []

    while loop_ss_exit == 0:
        # Current results dict
        ss_results = {}

        # Update increment factor
        inc_ss_factor = inc_factor_update(loop_ss_change,inc_ss_factor)
 
        # Increment sample size
        n = n + inc_ss_factor

        # Set estimated sample size errors based on expected error rate
        est_ss_errors = round(n*p,0)
        est_ss_p = round(est_ss_errors/n,r)

        # Get confidence interval results
        if ss_bounds != "U":
            lower_results, lower_all_results = hypergeom_ci("lower",ci,cl,est_ss_errors,n,N)
            # Set lower confidence limit
            lower_conf_limit = round(lower_results["conf_limit"],r)
            # Set lower precision exact for lower only and two-tailed
            lower_prec_exact = est_ss_p-lower_conf_limit
        if ss_bounds != "L":
            upper_results, upper_all_results = hypergeom_ci("upper",ci,cl,est_ss_errors,n,N)
            # Set upper confidence limit
            upper_conf_limit = round(upper_results["conf_limit"],r)
            # Set upper precision exact for lower only and two-tailed
            upper_prec_exact = upper_conf_limit-est_ss_p

        # Set precision exact
        if ss_bounds == "L":
            prec_exact = lower_prec_exact       
        elif ss_bounds == "U":
            prec_exact = upper_prec_exact
        else:
            prec_exact = max(est_ss_p-lower_conf_limit,upper_conf_limit-est_ss_p)

        # Increment loop
        loop_ss_cnt += 1

        loop_ss_cond, loop_ss_change, loop_ss_exit = loop_ss(inc_ss_factor,prec_desired,prec_exact,true_ss_cnt)

        if loop_ss_cond == 1 and true_ss_cnt == 0:
            # Set final sample size for the first occurrence of the loop condition true
            n_final = n
            prec_exact_final = prec_exact
            true_ss_cnt += 1
        elif loop_ss_cond == 1 and true_ss_cnt > 0:
            true_ss_cnt += 1
        elif loop_ss_cond == 0:
            # Reset true condition sample size count
            true_ss_cnt = 0

        # Capture results
        ss_results["bounds"] = ss_bounds
        ss_results["n"] = n
        ss_results["inc_factor"] = inc_ss_factor
        ss_results["prec_exact"] = prec_exact
        ss_results["prec_desired"] = prec_desired
        ss_results["loop_cond"] = loop_ss_cond
        ss_results["loop_change"] = loop_ss_change
        ss_results["loop_exit"] = loop_ss_exit
        ss_results["loop_cnt"] = loop_ss_cnt
        ss_results["true_cnt"] = true_ss_cnt

        # Append loop to results
        all_ss_results.append(ss_results)
        
    return n_final, prec_exact_final, all_ss_results


if __name__ == "__main__":
    start_time = time()

    # Lets run a test
    
    # Upper/Lower/Both
    tail = "both"
    # Confidence Level
    cl = .95
    # Population
    N = 9999999
    # Expected Error Rate
    p = .01
    # Precision
    B = .02

    # 95/0/1 2 tail - 368
    n_final, all_ss_results = hypergeom_ss(tail,cl,N,p,B)
    
    #print "Sample Size calculation took %f seconds" % (time()-start_time)
    print(f"Sample Size equals {n_final}")

    # Confidence Interval
    ci = 2
    # Number of errors
    X = 3
    # Sample Size
    n = 284

    # Lower confidence limit
    lower_results, lower_all_results = hypergeom_ci("lower",ci,cl,X,n,N)

    lower_conf_limit = lower_results["conf_limit"]

    # Upper confidence limit
    upper_results, upper_all_results = hypergeom_ci("upper",ci,cl,X,n,N)

    upper_conf_limit = upper_results["conf_limit"]

    ss_p = round(float(X)/n,6)
    
    # Precision exact
    prec_exact = round(max(ss_p-lower_conf_limit,upper_conf_limit-ss_p),6)
    print(f"Exact precision equals {prec_exact:.2f}")


    



    

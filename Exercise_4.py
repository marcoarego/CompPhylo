# -*- coding: utf-8 -*-
"""
Created on Sun Feb 01 18:03:17 2015

@author: Marcola
"""

# **** CODE BELOW TO BE POSTED TO GITHUB BY THURSDAY, FEB. 5TH ****


"""
Sometimes it will not be feasible or efficient to calculate the likelihoods for every
value of a parameter in which we're interested. Also, that approach can lead to large
gaps between relevant values of the parameter. Instead, we'd like to have a 'hill
climbing' function that starts with some arbitrary value of the parameter and finds
values with progressively better likelihood scores. This is an ML optimization
function. There has been a lot of work on the best way to do this. We're going to try
a fairly simple approach that should still work pretty well, as long as our likelihood 
surface is unimodal (has just one peak). Our algorithm will be:
(1) Calculate the likelihood for our starting parameter value (we'll call this pCurr)
(2) Calculate likelihoods for the two parameter values above (pUp) and below (pDown)
our current value by some amount (diff). So, pUp=pCurr+diff and pDown=pCurr-diff. To
start, set diff=0.1, although it would be nice to allow this initial value to be set
as an argument of our optimization function.
(3) If either pUp or pDown has a better likelihood than pCurr, change pCurr to this
value. Then repeat (1)-(3) until pCurr has a higher likelihood than both pUp and
pDown.
(4) Once L(pCurr) > L(pUp) and L(pCurr) > L(pDown), reduce diff by 1/2. Then repeat
(1)-(3).
(5) Repeat (1)-(4) until diff is less than some threshold (say, 0.001).
(6) Return the final optimized parameter value.
Write a function that takes some starting p value and observed data (k,n) for a
binomial as its arguments and returns the ML value for p.
To write this function, you will probably want to use while loops. The structure of
these loops is
while (someCondition):
    code line 1 inside loop
    code line 2 inside loop
    
As long as the condition remains True, the loop will continue executing. If the
condition isn't met (someCondition=False) when the loop is first encountered, the 
code inside will never execute.
If you understand recursion, you can use it to save some lines in this code, but it's
not necessary to create a working function.
"""

from exer2funct import *
from scipy import stats
from scipy.stats import rv_discrete
from scipy.stats import binom
import matplotlib.pyplot as plt

def binomPMF(n,k,p): #formely called bernou
    """this function has 3 parameters: n, K and p. The objects were created in 
    order to simplify the funtion, making it more readable"""
    a=binCoef2(n,k)
    b=pow(1-p,n-k)
    return a*(p**k)*b

data = 12
numTrials = 20
 
# Write a function that finds the ML value of p for a binomial, given k and n.

def P_ml(numTrials,data,pCurr,diff=0.1):  #pCurr is the first parameter and also the last, once it has to be replaced by a best value
    """
    this function will return the p value associated with the max likelihood value 
    """
    pUp=pCurr+diff
    pDown=pCurr-diff    
    L_pCurr=binomPMF(numTrials,data,pCurr)    
    L_pUp=binomPMF(numTrials,data,pUp)    
    L_pDown=binomPMF(numTrials,data,pDown)
    if (L_pCurr < L_pUp):
        while (L_pCurr < L_pUp):
            pCurr=pUp
            L_pCurr=binomPMF(numTrials,data,pCurr)
            pUp=pCurr+diff
            L_pUp=binomPMF(numTrials,data,pUp)
        return pCurr
    else:
        while (L_pCurr < L_pDown):
            pCurr=pDown   
            L_pCurr=binomPMF(numTrials,data,pCurr)
            pDown=pCurr-diff
            L_pDown=binomPMF(numTrials,data,pDown)
        return pCurr
    
    
test=P_ml(20,12,0.9,0.0001)
print test


def ml_value (numTrials,data,pCurr,diff=0.1):
    """This function will return the max likelihood value itself"""
    pUp=pCurr+diff
    pDown=pCurr-diff    
    L_pCurr=binomPMF(numTrials,data,pCurr)    
    L_pUp=binomPMF(numTrials,data,pUp)    
    L_pDown=binomPMF(numTrials,data,pDown)
    if (L_pCurr < L_pUp):
        while (L_pCurr < L_pUp):
            pCurr=pUp
            L_pCurr=binomPMF(numTrials,data,pCurr)
            pUp=pCurr+diff
            L_pUp=binomPMF(numTrials,data,pUp)
        return L_pCurr
    else:
        while (L_pCurr < L_pDown):
            pCurr=pDown   
            L_pCurr=binomPMF(numTrials,data,pCurr)
            pDown=pCurr-diff
            L_pDown=binomPMF(numTrials,data,pDown)
        return L_pCurr


"""
In the exercise above, you tried to find an intuitive cutoff for likelihood ratio
scores that would give you a reasonable interval in which to find the true value of 
p. Now, we will empirically determine one way to construct such an interval. To do 
so, we will ask how far away from the true value of a parameter the ML estimate 
might stray. Use this procedure: (1) start with a known value for p, (2) simulate
a bunch of datasets, (3) find ML parameter estimates for each simulation, and then 
(4) calculate the likelihood ratios comparing the true parameter values and the ML
estimates. When you do this, you will be constructing a null distribution of
likelihood ratios that might be expected if the value of p you picked in (1)
was true. Note that the ML values for these replicates are very often greater than
L(true value of P), because the ML value can only ever be >= L(true value). Once 
you have this distribution, find the likelihood ratio cutoff you need to ensure 
that the probability of seeing an LR score that big or greater is <= 5%. 
"""

# Set a starting, true value for p

trueP = 0.3

# Simulate 1,000 datasets of 200 trials from a binomial with this p

"""
def DiscrSample(x,p, sam_siz):
    test = stats.rv_discrete(name='test', values=(x, p))
    result = list(test.rvs(size=sam_siz))
    return result
"""

# If you haven't already done so, you'll want to import the binom class from scipy:
# from scipy.stats import binom
# binom.rvs(n,p) will then produce a draw from the corresponding binomial.

x=[0,1]
p=[0.7,0.3]
sam_siz=200
k_list1=[]

for y in range (1000):
    test = DiscrSample(x=[0,1],p=[0.7,0.3],sam_siz=200)
    count1 = test.count(1)
    k_list1.append(count1)
print k_list1 #list containing results of 1000 tests where the number of trials (sam_siz) in each test is 200 with p~0.3


# Now find ML parameter estimates for each of these trials(tests)
P_ml_list1=[]
for value in k_list1:
    q= P_ml(200,value,pCurr=0.5,diff=0.001)    
    P_ml_list1.append(q)
print P_ml_list1 #list containing p values related to each number in 'k_list1'

ml_list1=[]
for value in k_list1:
    j=ml_value(200,value,pCurr=0.5,diff=0.001)
    ml_list1.append(j)
print ml_list1 #this list contains the ML for each test result


# Calculate likelihood ratios comparing L(trueP) in the numerator to the maximum
# likelihood (ML) in the denominator. Sort the results and find the value
# corresponding to the 95th percentile.

TrueP_Like_list=[]
for u in k_list1:
    r=binomPMF(200,u,0.3)
    TrueP_Like_list.append(r)
print TrueP_Like_list #list containing Likelihood values related to each value in k_list1, using the True P

L_Rts = [float(b) / float(m) for b,m in zip(TrueP_Like_list, ml_list1)]# in this line I am dividing the TrueP ML list by the ML list. This will provide a list of ML ratios
print L_Rts

import numpy as np
print np.percentile(L_Rts,95) # this is the 95th percentile
print np.percentile(L_Rts,50) # this is the 50th percentile
print np.percentile(L_Rts,25) # this is the 25th percentile

# Now, convert the likelihood ratios (LRs) to -2ln(LRs) values.
log_L_Rts=(-2)*(np.log(L_Rts))
print log_L_Rts

# Find the 95th percentile of these values. Compare these values to this table:
# https://people.richland.edu/james/lecture/m170/tbl-chi.html. In particular, look at the 0.05 column. Do any of these values seem similar to the one you calculated?

print np.percentile(log_L_Rts,95) ##Result --> 3.593. This number is close to the first value in the 0.05 column which is 3.841

print np.percentile(log_L_Rts,75)
print np.percentile(log_L_Rts,50)


# Any idea why that particular cell would be meaningful?
#####Because in this exercise we are working with a Binomial scenario, which presents only one degree of freedom


# Based on your results (and the values in the table), what LR statistic value 
# [-2ln(LR)] indicates that a null value of p is far enough away from the ML value
# that an LR of that size is <=5% probable if that value of p was true?

####Values of [-2ln(LR)] larger than 3.593 in my case, or larger than 3.841 according to the Chi Square table can help me to exclude the respective values of p from my confidence interval.


rel_p = []
for i in range(0,105,5):
    x=i/100.0
    rel_p.append(x)
print rel_p

L_scores=[]
for p in rel_p:
    like=binomPMF(5,4,p)
    L_scores.append(like)
print L_scores

maximun=max(L_scores)
print maximun

like_ratios=[]
for y in L_scores:
    t=y/maximun
    like_ratios.append(t)
print like_ratios

logRts=(-2)*(np.log(like_ratios))
print logRts


data = 12
numTrials = 20

rel_p2 = []
for i in range(0,105,5):
    x=i/100.0
    rel_p2.append(x)
print rel_p2

L_scores2=[]
for p in rel_p2:
    like1=binomPMF(numTrials,data,p)
    L_scores2.append(like1)
print L_scores2 

maximun2=max(L_scores2)
print maximun2

like_ratios2=[]
for zz in L_scores2:
    aa=zz/maximun2
    like_ratios2.append(aa)
print like_ratios2

logRts2=(-2)*(np.log(like_ratios2))
print logRts2

# We've talked in previous classes about two ways to interpret probabilities. Which
# interpretation are we using here to define these intervals?
###This is a frequentist interpretation in my opinion
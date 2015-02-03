# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 11:21:05 2015

@author: Marcola
"""

# -*- coding: utf-8 -*-
"""
An Introduction to Likelihood
@author: jembrown
"""

"""
There are two primary ways to use probability models. Given what we know to be 
true about an experimental setup, we can make predictions about what we expect 
to see in an upcoming trial. For this purpose, probability functions are what 
we need. If R is the outcome of an experiment (i.e., an event) and p is a 
parameter of the probability function defined for these experimental outcomes, 
we can write these expectations or predictions as: P(R|p).

This conditional probability tells us what to expect about the outcomes of our 
experiment, given knowledge of the underlying probability model. Thus, it is a 
probability function of R given a value for p (and the model itself).

However, we might also wish to ask what we can learn about p itself, given 
outcomes of trials that have already been observed. This is the purview of the 
likelihood. Likelihoods are functions of parameters (or hypotheses, more 
generally) given some observations. The likelihood function of a parameter 
value is defined as: L(p;R) = P(R|p)

Note that this is the same probability statement we saw above. However, in this
context we are considering the outcome (R) to be fixed and we're interested in 
learning about p. Note that the likelihood is sometimes written in several 
different ways: L(p;R) or L(p) or L(p|R). P(R|p) gives a probability when R is 
discrete or a probability density when R is continuous. Since likelihoods are 
only compared for some particular R, we do not need to worry about this 
distinction. Technically speaking, likelihoods are just said to be proportional
to P(R|p), with the constant of proportionality being arbitrary.

There are some very important distinctions between likelihoods and probabilities.
First, likelihoods do NOT sum (or integrate) to 1 over all possible values of p.
Therefore, the area under a likelihood curve is not meaningful, as it is for 
probability.
It does not make sense to compare likelihoods across different R. For instance,
smaller numbers of observations generally produce higher values of P(R|p), 
because there are fewer total outcomes.
Likelihood curves provide useful information about different possible values of
p. When we are interested in comparing discrete hypotheses instead of 
continuous parameters, the likelihood ratio is often used:

L(H1;R)     P(R|H1)
-------  =  -------
L(H2;R)     P(R|H2)

Now, let's try using likelihoods to learn about unknown aspects of the process 
that's producing some data.

---> Inferring p for a binomial distribution <---
First, we'll start by trying to figure out the unknown probability of success 
associated with a Binom(5,p) random variable. If you want to try this on your 
own later, the following code will perform draws from a binomial with 5 trials.
You can simply change the associated value of p to whatever you'd like. To make
the inference blind, have a friend set this value and perform the draws from 
the Binomial for you, without revealing the value of p that they used.
"""
from exer2funct import *
from scipy.stats import binom
import matplotlib.pyplot as plt

n = 5
p = 0.5 # Change this and repeat

data = binom.rvs(n,p)

print data

"""
For the in-class version of this exercise, I'm going to perform a manual draw 
from a binomial using colored marbles in a cup. We'll arbitrarily define dark 
marbles as successes and light marbles as failures.
Record the outcomes here:
Draw 1: D
Draw 2: D
Draw 3: D
Draw 4: D
Draw 5: W
Number of 'successes': 4
Now record the observed number of succeses as in the data variable below.
"""

data =   4 # Supply observed number of successes here.
numTrials = 5


"""
Since we are trying to learn about p, we define the likelihood function as;
L(p;data) = P(data|p)
If data is a binomially distributed random variable [data ~ Binom(5,p)]
P(data=k|p) = (5 choose k) * p^k * (1-p)^(n-k)
So, we need a function to calculate the binomial PMF. Luckily, you should have
just written one and posted it to GitHub for your last exercise. Copy and paste 
your binomial PMF code below. For now, I will refer to this function as binomPMF(). 
"""

def binomPMF(n,k,p): #formely called bernou
    """this function has 3 parameters: n, K and p. The objects were created in 
    order to simplify the funtion, making it more readable. . We used """
    a=binCoef2(n,k)
    b=pow(1-p,n-k)
    return a*(p**k)*b

from JMB_binomial import * # I've stored my function in the file JMB_binomial.py
import scipy

"""
Now we need to calculate likelihoods for a series of different values for p to compare likelihoods. There are an infinite number of possible values for p, so let's confine ourselves to steps of 0.05 between 0 and 1.
""""
# Set up a list with all relevant values of p
rel_p = []
for i in range(0,105,5):
    x=i/100.0
    rel_p.append(x)
print rel_p

# Calculate the likelihood scores for these values of p, in light of the data you've collected
L_scores=[]
for p in rel_p:
    like=binomPMF(5,4,p)
    L_scores.append(like)
print L_scores

# Find the maximum likelihood value of p (at least, the max in this set)
maximun=max(L_scores)
print maximun

# What is the strength of evidence against the most extreme values of p (0 and 1)?
SE_val_0 = maximun/L_scores[0]
SE_val_1 = maximun/L_scores[-1]

# Calculate the likelihood ratios comparing each value (in the numerator) to the max value (in the denominator)
like_ratios=[]
for y in L_scores:
    t=y/maximun
    like_ratios.append(t)
print like_ratios

plt.scatter(rel_p,L_scores) # the maximum value of returned by the graph is 0.8. This reflects the 'p' of the maximum likelihood

"""
Now let's try this all again, but with more data. This time, we'll use 20 draws from our cup of marbles.
"""

data = 12  # Supply observed number of successes here.
numTrials = 20


# Calculate the likelihood scores for these values of p, in light of the data you've collected
rel_p2 = []
for i in range(0,1001,1):
    x=i/1000.0
    rel_p2.append(x)
print rel_p2

L_scores2=[]
for p in rel_p2:
    like1=binomPMF(numTrials,data,p)
    L_scores2.append(like1)
print L_scores2 #List with likelihood scores


# Find the maximum likelihood value of p (at least, the max in this set)
maximun2=max(L_scores2)
print maximun2 # this is the maximum likelihood value inside the list 'L_scores2'


# What is the strength of evidence against the most extreme values of p (0 and 1)?
SE_val_2 = maximun/L_scores2[0]
SE_val_3 = maximun/L_scores2[-1]


# Calculate the likelihood ratios comparing each value (in the numerator) to the max value (in the denominator)
like_ratios2=[]
for zz in L_scores2:
    aa=zz/maximun2
    like_ratios2.append(aa)
print like_ratios2 #this is a list of likelihood ratios

hist2 = dict(zip(rel_p2,L_scores2))
print hist2

plt.scatter(rel_p2,L_scores2) # the maximum value of returned by the graph is 0.6. This reflects the 'p' of the maximum likelihood
# When is the ratio small enough to reject some values of p?
#0.05
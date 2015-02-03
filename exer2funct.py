# -*- coding: utf-8 -*-
"""
Created on Thu Jan 29 11:46:40 2015

@author: Marcola
"""

from scipy import stats
from scipy.stats import rv_discrete
from scipy.stats import binom

def fact(h1,h0):
    """
    The main purpose of this function is to execute factorial calculations.
    It presents two arguments. The first one, 'h1' will be the integer that we
    intend to calculate its factorial. In the second argument, you need to place
    the number which you want the multiplication to stop. For instance, if you want
    the full factorial of the number 6, you will put the number 6 as the first
    argument of this function and the number 1 as the second argument: fact (6,1).
    If you want to multiply just 6*5*4, you can only change h0's value to 4 like
    this: fact(6,4).
    """
    x = 1
    while (h1 > h0-1 or 0): #The 'while' loop executes a statement while it is true
        x = x * h1 # in these two idented formulas the factorial is executed by 
        h1 = h1-1 # automatically replacing the values of 'h1' and 'x'
    return x # the final value of 'x' is the result of the factorial

def binCoef (n,k):
    """
    This function represents the function 2a where all factorial are fully
    calculated. It has two arguments that represent the same arguments as the 
    ones used in the Binomial Coefficient formula, named 'n'and 'k'. See chapter
    01, item 1.4.13, for the formula.
    """
    n1 = fact(n,1) #we called the previous created function 'fact' in order to
    k1 = fact(k,1) #reduce the number of lines we needed to write.
    c1 = fact (n-k,1)
    return n1/(c1*k1)

def binCoef2 (n2,k2):
    """
    This function does exactly the same thing as 'bioCoef(n,k)', although in 
    this new function the factorials aren't fully calculated. This saves a lot of
    time!!!
    """
    f = 1
    minus = n2-k2
    while n2 > minus:
        f = f*n2
        n2 = n2 - 1
    k3 = fact(k2,1)
    return f/(k3)
    
def binomPMF(n,k,p):
    """this function has 3 parameters: n, K and p. The objects were created in 
    order to simplify the funtion, making it more readable. . We used """
    a=binCoef2(n,k)
    b=pow(1-p,n-k)
    return a*(p**k)*b
    
def DiscrSample(x,p, sam_siz):
    test = stats.rv_discrete(name='test', values=(x, p))
    result = list(test.rvs(size=sam_siz))
    return result
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 24 23:04:02 2015

@author: Marcola
"""
import random
from scipy import stats
from scipy.stats import rv_discrete
import matplotlib.pyplot as plt

# Item 1
"""
(1) Write a function that multiplies all consecutively decreasing numbers 
between a maximum and a minimum supplied as arguments. (Like a factorial, but 
not necessarily going all the way to 1). This calculation would look like: 
max * max-1 * max-2 * ... * min
"""

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
fact(6,1)
fact(6,4)

#----------------------------------------------------------------------------#

# Item 2
"""
(2) Using the function you wrote in (1), write a function that calculates the 
binomial coefficient (see Definition 1.4.12 in the probability reading). 
Actually, do this twice. The first time (2a) calculate all factorials fully. 
Now re-write the function and cancel as many terms as possible so you can avoid
unnecessary multiplication (see the middle expression in Theorem 1.4.13).
"""
# This is the first function where all the multiplication are done (2a).
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
int(binCoef(7,4))

#Following next is the second function (2b) where we avoid unnecessary multiplication 
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
binCoef2 (7,4)

#----------------------------------------------------------------------------#

# Item 3

"""
(3) Try calculating different binomial coefficients using both the functions 
from (2a) and (2b) for different values of n and k. Try some really big values 
there is a noticeable difference in speed between the (2a) and (2b) function. 
Which one is faster? By roughly how much?


Now we will test both functions created to execute the Binomial Coefficient 
formula. The first formula with the fully calculated factorials and the second
formula where most of the multiplication process were cutted off.
"""
#This is the first function where the factorials are fully calculated
test1=binCoef(50000,3000)
#Up next is the second function, where most of the multiplication where cutted off
test2=binCoef2(50000,3000)
#Bellow is a logical test to see if the results of both functions are the same
test1==test2
"""
As a result, the function with less multiplication process was 5 times faster than
the one with all the multiplication.
"""

#----------------------------------------------------------------------------#

# Item 4
"""
(4) Use either function (2a) or (2b) to write a function that calculates the 
probability of k successes in n Bernoulli trials with probability p. This is 
called the Binomial(n,p) distribution. See Theorem 3.3.5 for the necessary equation. 
[Hint: pow(x,y) returns x^y (x raised to the power of y).]
"""
"""
Bernouli test - I used the function binCoef2 inside the function 'bernou'. The
function bernou calculates the probability of k successes in n Bernoulli trials
with probability p 
"""

def bernou(n,k,p):
    """this function has 3 parameters: n, K and p. The objects were created in 
    order to simplify the funtion, making it more readable. . We used """
    a=binCoef2(n,k) 
    b=pow(1-p,n-k)  
    return a*(p**k)*b #This is the final equation (see theorem 3.3.5)
bernou (5,2,0.4)      #the double ** as a replacement to the function pow(x,y)

#----------------------------------------------------------------------------#

#Item 5
"""
(5) Now write a function to sample from an arbitrary discrete distribution. 
This function should take two arguments. The first is a list of arbitrarily 
labeled events and the second is a list of probabilities associated with these 
events. Obviously, these two lists should be the same length.
"""
"""For this task, I used an existing function, rv_discrete, to simplify my
function. The rv_discrete function was obtained through the scipy module.
"""
def DiscrSample(x,p, sam_siz):
    test = stats.rv_discrete(name='test', values=(x, p))
    result = list(test.rvs(size=sam_siz))
    return result
    
#Testing the function with x being a list of numbers from 0 to 5
x = range(5)
#p represents a list with the probability of each number to occur
p = [0.1, 0.3,0.2,0.2,0.1,0.1]
DiscrSample(x,p,10)#10 is the number of times we'll run the test
    
#----------------------------------------------------------------------------#    
    
# Item 6
    
"""
(6) For an alignment of 400 sites, with 200 sites of type 1 and 200 of
type 2, sample a new alignment (a new set of site pattern counts) with 
replacement from the original using your function from (5). Print out the 
counts of the two types."""


types = range(2)# the types 1 and 2 will be 0 and 1
prob = [0.5,0.5]#0 and 1 have the same probability (50%)
Item6 = list(DiscrSample(types,prob,400))#generate 400 outcomes
print Item6
countA = Item6.count(0)
print countA
countB = Item6.count(1)
print countB

#----------------------------------------------------------------------------#

#Item 7 and 8

"""
(7) Repeat (6) 100 times and store the results in a list.

(8) Of those 100 trials, summarize how often you saw particular proportions 
    of type 1 vs. type 2. 
"""
## Setting the objects
lista0=[]
lista1=[]
prop0=[]

types = range(2)
prob = [0.5,0.5]

## for loops
for y in range (100):
    test = DiscrSample(types,prob,sam_siz=400)
    count0 = test.count(0)
    lista0.append(count0)
print lista0

prop=[float(y) for y in lista0]

prop1=[]
for d in prop:
    val = d/400
    prop1.append(val)
print prop1

#histograms displaying the frequencies of each value in the two proportion lists
plt.hist(lista0)
plt.hist(prop1)


#----------------------------------------------------------------------------#

#Item 9

""" (9) Calculate the probabilities of the proportions you saw in (8) 
using the binomial probability mass function (PMF) from (4).
"""

pmf_list=[]
for w in lista0:
    pmf= bernou(n=400,k=w,p=0.5)
    pmf_list.append(pmf)
print pmf_list

plt.hist(pmf_list)
    
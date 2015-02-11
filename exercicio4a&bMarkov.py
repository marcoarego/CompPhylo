# -*- coding: utf-8 -*-
"""
Created on Sun Feb 08 15:28:39 2015

@author: Marcola
"""

"""
Exercise 4
Discrete-time Markov chains
@author: jembrown
"""

"""
In this exercise, we will explore Markov chains that have discrete state spaces
and occur in discrete time steps. To set up a Markov chain, we first need to 
define the states that the chain can take over time, known as its state space.
To start, let's restrict ourselves to the case where our chain takes only two
states. We'll call them A and B.
"""

# Create a tuple that contains the names of the chain's states

state_names = ("A","B")
print state_names
type(state_names)


"""
The behavior of the chain with respect to these states will be determined by 
the probabilities of taking state A or B, given that the chain is currently in 
A and B. Remember that these are called conditional probabilities (e.g., the 
probability of going to B, given that the chain is currently in state A is 
P(B|A).)
We record all of these probabilities in a transition matrix. Each row
of the matrix records the conditional probabilities of moving to the other
states, given that we're in the state associated with that row. In our example
row 1 will be A and row 2 will be B. So, row 1, column 1 is P(A|A); row 1, 
column 2 is P(B|A); row 2, column 1 is P(A|B); and row 2, column 2 is P(B|B). 
All of the probabilities in a ROW need to sum to 1 (i.e., the total probability
associated with all possibilities for the next step must sum to 1, conditional
on the chain's current state).
In Python, we often store matrices as "lists of lists". So, one list will be 
the container for the whole matrix and each element of that list will be 
another list corresponding to a row, like this: mat = [[r1c1,r1c2],[r2c1,r2c2]]. 
We can then access individual elements use two indices in a row. For instance,
mat[0][0] would return r1c1. Using just one index returns the whole row, like
this: mat[0] would return [r1c1,r1c2].
Define a transition matrix for your chain below. For now, keep the probabilties
moderate (between 0.2 and 0.8).
"""

# Define a transition probability matrix for the chain with states A and B

rowA=[0.2,0.8]
rowB=[0.2,0.8]
matrix=[rowA,rowB]

# Try accessing a individual element or an individual row 
# Element
print matrix [0] [0] # print 0.2
print matrix [0] [1] # print 0.8
print matrix [1] [0] # print 0.2
print matrix [1] [1] # print 0.8

# Row
print matrix [0]
print matrix [1]

"""
Now, write a function that simulates the behavior of this chain over n time
steps. To do this, you'll need to return to our earlier exercise on drawing 
values from a discrete distribution. You'll need to be able to draw a random
number between 0 and 1 (built in to scipy), then use your discrete sampling 
function to draw one of your states based on this random number.
"""

# Import scipy U(0,1) random number generator
import scipy
from scipy import stats
from scipy import random
from numpy.random import uniform
from scipy.stats import rv_discrete


# Paste or import your discrete sampling function

def DiscrSample(x,p, sam_siz):
    test = stats.rv_discrete(name='test', values=(x, p))
    result = int(test.rvs(size=sam_siz)) # I changed the original function so it woul return me an integer and not a list
    return result

# Write your Markov chain simulator below. Record the states of your chain in 
# a list. Draw a random state to initiate the chain.
rowA=[0.4,0.6]
rowB=[0.7,0.3]
mat=[rowA,rowB]
print mat
len(mat)


def Markov2events (matrix,steps):
    iv = uniform(0.0,1.0)
    currValue = DiscrSample(x=[0,1],p=[iv,1-iv],sam_siz=1)
    list1 = [currValue]
    for i in range(steps-1):
        currValue = DiscrSample(x=[currValue,1-currValue],p=matrix[currValue],sam_siz=1)
        list1.append(currValue)      
    return list1

   
# Run a simulation of 10 steps and print the output.

test=Markov2events(matrix=mat,steps=10)
print test

# ----> Try to finish the above lines before Tues, Feb. 10th <----

# Now try running 100 simulations of 100 steps each. How often does the chain
# end in each state? How does this change as you change the transition matrix?
rowA=[0.4,0.6]
rowB=[0.7,0.3]
mat=[rowA,rowB]
print mat

lastNumber = []
for i in range (100):
    a = Markov2events(matrix=mat,steps=100)
    b = a[-1] # appending the last value to the list  
    lastNumber.append(b)

print lastNumber
lastNumber.count(0) # returned a proportion near 40/60
lastNumber.count(1)

rowAa=[0.5,0.5]
rowBb=[0.5,0.5]
mat2=[rowAa,rowBb]

lastNumber2 = []
for h in range (100):
    r = Markov2events(matrix=mat2,steps=100)
    s = r[-1] # appending the last value to the list  
    lastNumber2.append(s)

lastNumber2.count(0) # returned a proportion near 50/50
lastNumber2.count(1)


# Try defining a state space for nucleotides: A, C, G, and T. Now define a 
# transition matrix with equal probabilities of change between states.
row1 = [0.25,0.25,0.25,0.25]
row2 = [0.25,0.25,0.25,0.25]
row3 = [0.25,0.25,0.25,0.25]
row4 = [0.25,0.25,0.25,0.25]
mainMat = [row1,row2,row3,row4]
print mainMat

    
def Markov4events (matrix,steps):
    """This function works only with integers as states. Since it has 4 states, it will present 0, 1, 2 and 3 as states"""
    currValue = int(random.randint(0,4)) 
    list1 = [currValue]
    for i in range(steps-1):
        if currValue == 0:
            currValue = DiscrSample(x=[currValue,currValue+1,currValue+2,currValue+3],p=matrix[currValue],sam_siz=1)
        elif currValue == 1:
            currValue = DiscrSample(x=[currValue-1,currValue,currValue+1,currValue+2],p=matrix[currValue],sam_siz=1)
        elif currValue == 2:
            currValue = DiscrSample(x=[currValue-2,currValue-1,currValue,currValue+1],p=matrix[currValue],sam_siz=1)
        elif currValue == 3:
            currValue = DiscrSample(x=[currValue-3,currValue-2,currValue-1,currValue],p=matrix[currValue],sam_siz=1)
        list1.append(currValue)
    list2 = [str(x) for x in list1]
    list2 = [w.replace('0', 'A') for w in list2] # Converting the numbers to letters
    list2 = [p.replace('1', 'C') for p in list2]
    list2 = [u.replace('2', 'G') for u in list2]
    list2 = [n.replace('3', 'T') for n in list2]
    return list2

test2=Markov4events(mainMat,40)
print test2
   
# Again, run 100 simulations of 100 steps and look at the ending states. Then
# try changing the transition matrix.
   
test3=Markov4events(mainMat,100)
print test3

mainMat2 = [[0.2,0.3,0.4,0.1],[0.25,0.15,0.3,0.3],[0.2,0.2,0.3,0.3],[0.7,0.1,0.1,0.1]]

lastNumber3=[]
for m in range(100):
    r = Markov4events(matrix=mainMat2,steps=100)
    s = r[-1] # appending the last value to the list
    lastNumber3.append(s)
print lastNumber3

lastNumber3.count('A')
lastNumber3.count('T')
lastNumber3.count('C')
lastNumber3.count('G')
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 10:36:48 2015

@author: Marcola
"""

import scipy 
import random 
from numpy.random import uniform
from math import log 


### Creating the object that will contain a 'Q' and a function that will generate a Continuous Markov Chain simulation.

class ContMarkSim(object): # to define an object we need to use the 'class' statement and the word 'object' inside the parentheses
    
    '''Defining the object variables'''    
    T_end=50 #Ending Time, in this case 50. I will use a while statement that will execute the function untill that time peiod is reached
    Q_matrix=[[-1.916,0.541,0.787,0.588],[0.148,-1.069,0.417,0.506],[0.286,0.170,-0.591,0.135],[0.525,0.236,0.594,-1.355]] # this is the Q matrix, or Rates matrix. The negative numbers are in the positions [0][0],[1][1],[2][2]and [3][3]. By doing this we are assigning them to the main diagonal of the matrix.
    
    '''Defining the function that will make my simulation'''
    def MarkSim(self):
        tup=(0,1,2,3) # State spaces for the simulation
        state=random.choice(tup) # using the random.choice function to sort the initial state
        Time=0 # Defining the starting point in the branch
        list1=[] # This will be the list with the waiting times
        list2=[] # this list will hold the probabilities of going to any particular event based on the present one
        list3=[] # list containing the possible events to come in the next change in the chain
        list4=[state] # This will be the list with the states. Note that thi list already has the first state choosed from the random.choice function
       
        '''definning a function inside a function'''
        def discSamp(events,probs): # This is the Discrete sampling function used in other exercises
            ranNum = scipy.random.random()
            cumulProbs = []
            cumulProbs.extend([probs[0]])
            for t in range(1,len(probs)):
                cumulProbs.extend([probs[t]+cumulProbs[-1]])
            for t in range(0,len(probs)):
                if ranNum < cumulProbs[t]:
                    return events[t]
            return None
        
        '''Now back to the main function body'''
        while Time <= self.T_end: # executing the function untill it reaches the ending time
            u1=uniform(low=0.0, high=1.0, size=None)
            Wtime=-(1/-(self.Q_matrix[state][state]))*log(u1) # Calculating the waiting time
            list1.append(Wtime)
            Time=sum(z for z in list1)
            for g in self.Q_matrix[state]:
                if g >= 0:
                    p=-(g/(self.Q_matrix[state][state]))##calculating the probabilities of the next events according to the current state. We divide the positive rates of the row, by the diagonal value
                    list2.append(p)
            for next_state in [0,1,2,3]: # list with possible next states
                if next_state != state: # for next state being different from current
                    list3.append(next_state) # appending following states to list3
            i=discSamp(events=list3,probs=list2)
            list4.append(i)
            
            '''Converting the original numbers to letters corresponding the nucleotide bases'''
            list5 = [str(x) for x in list4]
            list5 = [A.replace('0', 'A') for A in list5] # Converting the numbers to letters
            list5 = [C.replace('1', 'C') for C in list5] # Converting the numbers to letters
            list5 = [G.replace('2', 'G') for G in list5] # Converting the numbers to letters
            list5 = [T.replace('3', 'T') for T in list5] # Converting the numbers to letters           
            
        return list5,list1, sum(list1) # this function is returning the list of events (list5), the waiting time for each event to happen, and the final sum of the waiting time so we can check last value
        
d= ContMarkSim()
print d.T_end

d.T_end = (500) # changing the value of the variable inside the objec

print d.T_end
print d.Q_matrix
print d.MarkSim()


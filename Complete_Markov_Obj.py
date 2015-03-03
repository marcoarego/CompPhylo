# -*- coding: utf-8 -*-
"""
Created on Mon Mar 02 12:15:54 2015

@author: Marcola
"""

# Continuous time Markov Simulation
''' Importing modules'''
import scipy  
from numpy.random import uniform
from math import log 
from math import exp
from itertools import tee, islice, chain, izip
import operator
import functools
from scipy import linalg
import numpy as np


'''Defining a class'''
class ContMarkSim(object):
    
    def __init__ (self,v=10,Q=[[-1.916,0.541,0.787,0.588],[0.148,-1.069,0.417,0.506],[0.286,0.170,-0.591,0.135],[0.525,0.236,0.594,-1.355]],waittimes=[],sta_numb=[], states=[],TotalProbs=[],MargProb=[]): # variable defaults were set in the constructor funcion
        '''the self variable represents the instance of the object itself'''
        self.v=v # Ending Time, in this case I chose 10
        self.Q=Q # Q matrix, or Rates matrix
        ##Tuple with my state space
        self.waittimes=waittimes##This will be the list with the waiting times
        self.sta_numb=sta_numb##This will be the list with the states (in numbers where 0=A,1=C,2=G and 3=T)
        self.states=states##This will give you the states as nucleotides.
        self.TotalProbs=TotalProbs # probability of any given chain (this is usually very low)
        self.MargProb=MargProb # probability of starting with 'x' and ending with 'y' (for example)
    # Markov Simulation Function
    def MarkSim(self):
        statespace=[0,1,2,3]
        # defining the function of discrete sample, to sort the events to change
        def discSamp(events,probs):
            ranNum = scipy.random.random()
            cumulProbs = []
            cumulProbs.extend([probs[0]])
            for t in range(1,len(probs)):
                cumulProbs.extend([probs[t]+cumulProbs[-1]])
            for t in range(0,len(probs)):
                if ranNum < cumulProbs[t]:
                    return events[t]
            return None
        A = np.squeeze(np.asarray(self.Q)) # this transforms the Q matrix in an array A
        StationaryProbs=linalg.expm(A*1000) # we need an array to input in the linalg function
        i=discSamp(events=statespace,probs=StationaryProbs[1]) # Sorting the first state
        E=0 # Branch starting point
        self.sta_numb=[i] # list with the numbers representing the states.
        list6=[] # this list will keep the exponPDF values
        list7=[]
        list10=[] 
        M4=StationaryProbs[1][i]
        # this 'for loop' calculates the transition matrix values for each state based on the Q matrix
        for d in range(len(self.Q)):
            lista=[]    
            for g in self.Q[d]:
                if g > 0: 
                    p=-(g/(self.Q[d][d]))
                else:
                    p = 0
                lista.append(p)
        list10.append(lista)
        T= list10 # this is the transition matrix
        
        def previous_and_next(some_iterable):
            '''this function will be used to set two lists, one with the current
            sampled states, and another with the respective next states. This will 
            be used to calculate the probabilities of each state change.'''
            prevs, items, nexts = tee(some_iterable, 3)
            prevs = chain([None], prevs)
            nexts = chain(islice(nexts, 1, 0), [0])
            return izip(prevs, items, nexts)
        while E <= self.v:##this loop will occur while the sum of the waiting times does not reach the ending time
            p_list=[] # this list will hold the probabilities of next events, given the previous events
            u1=uniform(low=0.0, high=1.0, size=None) # sorting an uniform value to give to the next equation
            Wtime=-(1/-(self.Q[i][i]))*log(u1) # calculating the waiting time, by assuming that lambda is equal to the -diagonal (-Qii) value
            self.waittimes.append(Wtime) # putting the waiting times in list 1
            E=sum(z for z in self.waittimes) # adding the waiting times to make the simulation stop when reaching the ending time
            for g in self.Q[i]:
                if g > 0: 
                    p=-(g/(self.Q[i][i])) # calculating the probabilities of the next events according to the current state. We divide the positive rates of the row, by the diagonal value
                else:
                    p = 0
                p_list.append(p) # puting the probabilities in a list
            for s in self.waittimes:
                exponPDF = -(self.Q[i][i]) * exp(- (-(self.Q[i][i]))*s) # this function calculates the exponential probability density function for each waiting time according to the current state
                list6.append(exponPDF) # this list keeps exponPDF values
                M1=functools.reduce(operator.mul, list6, 1) # this list has the multiplication of all exponPDF values
            i=discSamp(events=[0,1,2,3],probs=p_list)##sampling the next state according to the current state
            self.sta_numb.append(i) # appending the states in a list
            lista1=[]
            lista2=[]
            for previous, item, nxt in previous_and_next(self.sta_numb):
                lista1.append(item) # this list will keep the previous sampled nucleotide
                lista2.append(nxt) # this list will keep the subsequent nucleotide
            for c,m in zip (lista1,lista2): # this loop will give the right values to transition matrix values
                if m != lista2[-1]:
                    statesprobs=T[c][m] # this will give the probabilities of going from the previous state to the next
                    list7.append(statesprobs) # this will append the probabilities of each state in list7
                else:
                    pass
                    M2=functools.reduce(operator.mul, list7, 1) # this will multiply the probabilities of sampling each state
            MargProbs=linalg.expm(A*self.v) # this calculates the marginal probabilities
            finaltime=self.waittimes[-1]
            finalstate=self.sta_numb[-1]
            cdffinaltime=1-(2.71826**(-(self.Q[finalstate][finalstate])*finaltime))
            M3=1-cdffinaltime
            self.MargProb=MargProbs[self.sta_numb[0]][self.sta_numb[-1]]
            self.TotalProbs=M4*M1*M2*M3 # this multiplies the probabilities for each waiting time and for the each nucleotide, in other words, the probability of all the events that occurred in my simulation
            self.states = [str(h) for h in self.sta_numb] # transforming the numbers into strings so we can replace them for their respective letters
            self.states = [z.replace('0', 'A') for z in self.states] # Converting the numbers to letters representing nucleotides
            self.states = [a.replace('1', 'C') for a in self.states] # Converting the numbers to letters representing nucleotides
            self.states = [u.replace('2', 'G') for u in self.states] # Converting the numbers to letters representing nucleotides
            self.states = [n.replace('3', 'T') for n in self.states] # Converting the numbers to letters representing nucleotides
        return self.states,self.waittimes,self.TotalProbs,self.MargProb # returning the states(substitutions) and waiting times until reaching an ending time




# Giving a name to the object        
d= ContMarkSim()
# Getting the results of the Continuous time Markov Simultaion
print d.MarkSim()
print d.states
print d.waittimes
print d.TotalProbs
print d.MargProb
# Testing changes in the waiting time
banana=ContMarkSim(v=200)
print banana.MarkSim()
print banana.states
print banana.waittimes
print banana.TotalProbs
print banana.MargProb

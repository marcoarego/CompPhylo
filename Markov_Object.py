# -*- coding: utf-8 -*-
"""
Created on Wed Feb 18 21:51:15 2015

@author: Marcola
"""

import scipy as sp
def discSamp(events,probs):
    ranNum = sp.random.random()
    cumulProbs = []
    cumulProbs.extend([probs[0]])
    for i in range(1,len(probs)):
        cumulProbs.extend([probs[i]+cumulProbs[-1]])
    for i in range(0,len(probs)):
        if ranNum < cumulProbs[i]:
            return events[i]
    return None


class MarkovObj(object):
    characters=("A","T","C","G")
    probabilities=[[0.25,0.25,0.25,0.25],[0.25,0.25,0.25,0.25],[0.25,0.25,0.25,0.25],[0.25,0.25,0.25,0.25]]
    num=100

    def dmcSim(self):

    # Define list to hold chain's states
        final_chain = []    

    # Draw a state to initiate the chain
        currState = discSamp(events=self.characters,probs=[1.0/len(self.characters) for x in self.characters])
        final_chain.extend(currState)

    # Simulate the chain over n-1 steps following the initial state
        for step in range(1,self.num):
            pro = self.probabilities[self.characters.index(currState)] # Grabbing row associated with currState
            currState = discSamp(self.characters,pro) # Sample new state
            final_chain.extend(currState)        
        
        return final_chain
banana=MarkovObj()
print banana.characters
print banana.probabilities
print banana.dmcSim()
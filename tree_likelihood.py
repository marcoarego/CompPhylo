# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 20:54:03 2015

@author: Marcola
"""
# import modules

from scipy import linalg
import numpy as np
import re

# creating the class Node, an object capable of hold the relationships of different nodes of a tree
class Node:
    
    def __init__(self,name="",parent=None,children=None):
        self.name = name
        self.parent = None
        if children is None:
            self.children = []
        else:
            self.children = children

extVar = None
    
def ReadNewick (data,base="root",space = None, extVar = None):
    """ This function gets the string of a tree, like (A:1,(B:0.1,(C:0.2,D:1):1):0.6) for example,
    and construct the relationship between the nodes. The PC memory stores their relationship thanks
    to object class Node
    """
    if base == "root" :# the content of this if statement will deal with the sons of the root
        root = Node("root")
        root.brl = "0" # the branch length of the root is zero
        extVar = root # this sentence grab the memory space that has the root of the tree. With this, we can access the tree and all the relationship of its nodes outside the function
        son1 = data.partition('(')[-1].rpartition(')')[0] # this will take off the external parentheses
        regEx1 = re.compile(r'(.*?)\(.*\)') # we used the regex module to read in and outside parentheses of the input string tree
        result1 = re.findall(regEx1, son1)
        result1=result1.pop(0)
        result2 = son1.replace(result1,"")
        result1 = result1[:-1]
        result1=result1.split(",") 
        result1.append(result2) # the list "result1" contains the sons of the root
        for obj in result1:
            branchname=obj.rpartition(":")[0] # separates the name of the son from his branch length
            child=Node(name=branchname,parent="root")
            print "root son: ", child.name
            root.children.append(child)
            brl=obj.rpartition(":")[-1]
            child.brl=brl
            print "son's branch length: ",child.brl
            child.parent = root
            print "son's parent: ", child.parent.name,"\n"
        for item in root.children:
            if item.name[0]=="(":
                ReadNewick(data=item.name,base="x",space=item) # this recursion will take us to other nodes other than the root
        return extVar
        
    else: # if the node is not the "root"
        """ There are several if and else statement in this part of the function because we
        tried to deal with all tree structure that we could        
        """
        son2 = data.partition('(')[-1].rpartition(')')[0] # this will take off the external parentheses
        regEx1 = re.compile(r'(.*?)\(.*\)')
        result3 = re.findall(regEx1, son2)
        if result3 != [] and result3 != ['']:
            result3=result3.pop(0)
            result4 = son2.replace(result3,"")
            result3=result3.split(",")
            result3 = result3[:-1]
            result3.append(result4) # the list "result3" contains the sons of a given node, depending on the structure of original tree string
            for obj in result3:
                branchname=obj.rpartition(":")[0]
                child=Node(name=branchname,parent=data)
                print "internal node child: ",child.name # this will give the child of each internal node
                space.children.append(child)
                brl=obj.rpartition(":")[-1]
                child.brl=brl
                print "internal node branch length: ",child.brl # this will give the branch length for each internal node
                child.parent=space
                print "son's parent: ", child.parent.name,"\n" # parent
                if child.name[0]=="(": # if the next node starts with "(", repeat this 'if' statement; the function will go to the next "else" otherwise
                    ReadNewick(data=child.name,base="x",space=child)
        else:
            if son2.count("(") > 0:
                son3=son2.split(",(")
                print son3
                for value in son3:
                    if value[0]!="(":
                        value1 = son3.pop(son3.index(value))
                        value2 = "(" + value1
                        son3.append(value2)
                print son3
                for obj in son3:
                    branchname=obj.rpartition(":")[0]
                    child=Node(name=branchname,parent=data)
                    space.children.append(child)
                    brl=obj.rpartition(":")[-1]
                    child.brl=brl
                    print "node child: ",child.name
                    print "node branch length: ",child.brl
                    child.parent=space 
                    print "son's parent: ",child.parent.name,"\n" # parent
                    ReadNewick(data=child.name,base="x",space=child)
            else:
                son3 = son2.split(",") # son3 is now a list containing terminal nodes, hopefully... :p
                for obj in son3:
                    branchname=obj.rpartition(":")[0]
                    child=Node(name=branchname,parent=data)
                    space.children.append(child)
                    brl=obj.rpartition(":")[-1]
                    child.brl=brl
                    print "node child: ",child.name
                    print "node branch length: ",child.brl
                    child.parent=space
                    print "son's parent: ",child.parent.name,"\n" #parent
                    
                    
SomeTree="(A:1,(B:0.1,(C:0.2,D:1):1):0.6)" 

x=ReadNewick(data=SomeTree)


for y in x.children:
    print y.name
  
def newick(node): # this newick was made by Subir. 
        """
        A method of a Tree object that will print out the Tree as a 
        parenthetical string (Newick format).
        """
        parenthetical = "(" 
        if len(node.children) == 0:
            return node.name + ":" + str(node.brl)
        else:
            for child in node.children:
                if node.children[-1] == child: 
                    parenthetical += newick(child)
                else:
                    parenthetical += newick(child) + ","
            if node.brl is not None:
                parenthetical += "):" + str(node.brl)
            else:
                parenthetical += ")"
            return parenthetical   


def printNames(node):
        """
        A method of a Tree object that will print out the names of its
        terminal nodes.
        """
        if len(node.children) > 0:
            node.nucl = []
            for child in node.children:
                printNames(child)
        else:
            if node.name == "A":
                node.nucl=[1,0,0,0]
            elif node.name == "B":
                node.nucl=[1,0,0,0]
            elif node.name == "C":
                node.nucl=[1,0,0,0]
            elif node.name == "D":
                node.nucl=[1,0,0,0]
            print node.nucl
        return node
                
            
banana = printNames(x)

def TreeLikelihood (node,Q): # listA,listC,listG,listT):
    """this function calculates the tree likelihood for one nucleotide site"""
    array = np.squeeze(np.asarray(Q)) # this transforms my Q matrix in an array A    
    if len(node.children) > 0 and node.nucl == []:
       for child in node.children:
           TreeLikelihood(child,Q) # listA,listC,listG,listT  
    else:
        # lists that will hold the probabilities of having each base on each internal nodes
        listA=[]
        listC=[]
        listG=[]
        listT=[]
        if node.name != "root":
            for x in node.parent.children:
                if x.nucl != []:
                    print "node name: ", x.name,"(","branch length: ",x.brl,"/","nucleotide probs(A,C,T,G): ",x.nucl,")","\n"
                
                if x.nucl != []: # do the following calculations if the node already has information on probabilities for each nucleotide
                    
                    MP=linalg.expm(array*float(x.brl)) # getting the marginal probs according to each branch length size
                    
                    #calculating the probabilities of going from one base to others
                    x.resA=MP[0][0]*x.nucl[0]+MP[0][1]*x.nucl[1]+MP[0][2]*x.nucl[2]+MP[0][3]*x.nucl[3] # moving from A
                    x.resC=MP[1][0]*x.nucl[0]+MP[1][1]*x.nucl[1]+MP[1][2]*x.nucl[2]+MP[1][3]*x.nucl[3] # .....  ...  C
                    x.resG=MP[2][0]*x.nucl[0]+MP[2][1]*x.nucl[1]+MP[2][2]*x.nucl[2]+MP[2][3]*x.nucl[3] # .....  ...  T
                    x.resT=MP[3][0]*x.nucl[0]+MP[3][1]*x.nucl[1]+MP[3][2]*x.nucl[2]+MP[3][3]*x.nucl[3] # .....  ...  G
                    
                    # lists for the results of each nucleotide
                    listA.append(x.resA)
                    listC.append(x.resC)
                    listG.append(x.resG)
                    listT.append(x.resT)
                    
                if len(listA)>1:
                    squareA=reduce(lambda z, y: z*y,listA)
                if len(listC)>1:
                    squareC=reduce(lambda z, y: z*y,listC)
                if len(listG)>1:                
                    squareG=reduce(lambda z, y: z*y,listG)
                if len(listT)>1:
                    squareT=reduce(lambda z, y: z*y,listT)
                    lista = [squareA,squareC,squareG, squareT]
                    x.parent.nucl=lista
                    print "Parent node with their probabilities for each nucleotide: ",x.parent.name,"(nucleotide probs (A,C,T,G): ",x.parent.nucl,")\n"
                    if x.parent.nucl != []:
                            TreeLikelihood (node.parent,Q)
                    else:
                        pass                        
                        
        else:            
            StatProb=linalg.expm(array*100) # calculating the stationary probability matrix to use it on the likelihood calculation
            Likelihood = node.nucl[0]*StatProb[0][0]+node.nucl[1]*StatProb[0][1]+node.nucl[2]*StatProb[0][2]+node.nucl[3]*StatProb[0][3]
            print "Likelihood: ", Likelihood, "\n"        
            return Likelihood # returning the likelihood value so it can be used afterwards    
                
    

TreeLikelihood(node=banana,Q=[[-1.916,0.541,0.787,0.588],[0.148,-1.069,0.417,0.506],[0.286,0.170,-0.591,0.135],[0.525,0.236,0.594,-1.355]])

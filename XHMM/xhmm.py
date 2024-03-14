"""
NAME: Caroline Cutter

CSCI 0321A - Bioinformatics

ASSIGNMENT: XHMM Implementation 
            

DUE: Friday, April 22

"""

# IMPORTS
import numpy as np
import math
from scipy.stats import norm


# GLOBALS

Q = 1/6            # numbers for the transition matrix 
P = 0.00000001

transition = np.array([[1-Q, Q, 0],[P, 1-2*P, P],[0, Q, 1-Q]], dtype=np.longdouble)   # table 1 from paper 





# FUNCTIONS

# Read File
def read_xhmm(filename):
    """
    Reads in the XHMM File and returns arrays of samples, targets and data
    
    ** This is modified from Prof Linderman's implementation suggestions**
    """
    
    with open(filename) as file:
        targets = file.readline().strip().split()[1:]
        
        samples = np.genfromtxt(filename, skip_header=1, usecols=0, dtype=str)
        
        data = np.genfromtxt(filename, skip_header=1, usecols=range(1,1+len(targets)), dtype=np.longdouble)
    
    return samples, targets, data

    
    

def xhmm_viterbi(samples, targets, data):
    '''
    The viterbi algorithm using xhmm implementation
    
    Args:
        samples: matrix of samples 
        targets: matrix of targets
        data: data from file 
    
    Returns:
        NONE but prints samples in SAMPLE CNV START END Q_EXACT form 
    '''


    # For each example 
    for i in range(len(samples)):
        
        # Make Emission Matrix - A 3 x |targets| matrix computed from the observations in the input file
            # Probability density function of a normal distribution with variance 
            # of 1 and mean of -M, 0, and M for DEL, DIP, and DUP respectively (where M=3)
        emission = np.zeros((3, len(data[i,:])), dtype=float)
        emission[0,:] = norm.pdf((data[i,:]), -3, 1)                         # -M, 1
        emission[1,:] = norm.pdf((data[i,:]), 0, 1)                          #  0, 1
        emission[2,:] = norm.pdf((data[i,:]), 3, 1)                          # +M, 1
        
        

        # VITIRBI - finding most likely path 
    
        n = 3
        m = 265
    
    
        # Initialize Matrices
        score = np.zeros((n, m), dtype=float)           
        back = np.zeros((n, m), dtype=int)
    

        # BASE CASES:
        score[0,0] = emission[0,0] * P         # using fixed probabilities from transition table        
        score[1,0] = emission[1,0] * 1-(2*P)
        score[2,0] = emission[2,0] * P 
    
        # RECURSIVE CASE: 
        for d in range(1, m):
            for k in range(n):
            
                # Get incoming edges with matrix multiplication:
                incoming_edges = (score[:,d-1] * transition[:,k] * emission[k,d])
                score[k,d] = max(incoming_edges)
                back[k,d] = np.argmax(incoming_edges)
     

     
        # BACKTRACK - get sequence
 
        end_state = np.argmax(score[:,-1])
        path = []

        for v in range(m - 1, -1, -1):
            path = [end_state] + path 
            end_state = back[end_state, v]
        
        # Keep track of start and end indicies 
        start = -1
        end = -1
     
        for x in range(len(path)):
            if path[x] != 1:
                count = 0
                start = x
                while path[start + count] != 1:
                    count += 1
                end = x + count -1
                break
         
        # Skip the rest if start = -1 
        if start == -1:
            continue

        start_target = targets[start]           # Find the tart and end targets 
        end_target = targets[end]
        
    
        # FORWARD ALGORITHM 
        forward = np.zeros((n, m), dtype= float)

        # BASE CASES 
        forward[0,0] = emission[0,0] * P
        forward[1,0] = (1-(2 * P)) * emission[1,0]
        forward[2,0] = emission[2,0] * P 

        # RECURSIVE CASES 
        for s in range(1, m):
            for j in range(n):
                edges = (forward[:,s - 1] * transition[:,j] * emission[j,s])
                forward[j,s] = sum(edges)
        
        
        
        # BACKWARD ALGORITHM 
        backward = np.zeros((n, m), dtype= float)
    
        # BASE CASE:
        for b in range(n):
            backward[b,-1] = 1
        
      
        # RECURSIVE CASE: 
        for c in range(m-2, -1, -1):
            for k in range(n):
                edges = (backward[:,c+1] * transition[k,:] * emission[:,(c+1)])
                backward[k,c] = sum(edges)
    
    
        
        # Now get answer by summing forward * backward / forward(sink)
        weight = 1
        for l in range(start + 1, end + 1):
            weight = weight * transition[path[l-1], path[l]] * emission[path[l], l]
        
        forward_sum = forward[path[start], start]
        backward_sum = backward[path[end], end]
        
        foward_sink = sum(forward[:,-1])
        
        answer = (weight * forward_sum * backward_sum)/ foward_sink 
    
    

        # XHMM Specifics
        
        q = -10*(math.log10(1 - (answer)))        # Phred Score using forward-backward probabilities 

        CNV = ""                                # Find whether it is DEL or DUP (excluding DIP) 
        
        if path[start] == 0:
            CNV = "DEL"
        elif path[start] == 2:
            CNV = "DUP"

        


        
      
        print(samples[i], CNV, start_target, end_target, q)
      








# MAIN
if __name__ == "__main__":
    
    samples, targets, data = read_xhmm("XHMM.in.txt")
    
    print("SAMPLE CNV START END Q_EXACT")
   
    xhmm_viterbi(samples, targets, data)
    

    
    
    
    
    
    
    


"""
NAME: Caroline Cutter

CSCI 0321A - Bioinformatics

ASSIGNMENT: Rosalind Randomized Algorithms #7
            Compute the Probability of a String Emitted by an HMM
            

DUE: Friday, April 15

"""


# IMPORTS
import numpy as np
import viterbi
import outcome_likelihood



# FUNCTIONS

# NOTE:
    # Forward algorithm is the "outcome_likelihood" function
    # Backward algorithm is in the "soft decoding" problem 


def soft_decoding(sequence, alphabet, states, transition, emission):
    '''
    Solves the soft decoding problem 
     
    Args:
        sequence: a string of characters
        alphabet: list of characters in which sequence was constucted from
        states: list of states
        transition: 2D transition matrix
        emission: 2D emission matrix 
    
    Returns:
        Integer probability Pr(πi = k|x) that the HMM was in state k at step i (for each state k and each step i).
    '''
    
    # FORWARD MATRIX
    forward, outcome = outcome_likelihood.outcome_likelihood(sequence, alphabet, states, transition, emission)
    
    
    
    # BACKWARD MATRIX 
    
    
    n = len(states)
    m = len(sequence)
    
    
    # Initialize Matrices
    backward = np.zeros((n, m), dtype=float)           
    
    
    # BASE CASE: (Backward(k,n) = 1)
    for k in range(n):
         backward[k][-1] = 1.0
         
    
    
    # RECURSIVE CASE: backward(k,i) = sum all states l ( backward(l,i+1 * transition(k, l) * emission-l(x(i+1)) ) 
    for i in range(m - 1, 0, -1):
        for k in range(n):
            for l in range(n):
                
                symbol_at_x1 = alphabet.index(sequence[i])
                
                backward[k][i - 1] = backward[k][i-1] + (backward[l][i] * transition[k][l] * emission[l][symbol_at_x1])
                
               

    # SOFT DECODE
    decode_matrix = np.zeros((m, n), dtype = float)      
    
    
    # Calculate Forward*Backward / Forward-sink for each state k and position i
    
    for i in range(m):
        for k in range(n):
            forward_sink = 0
            for l in range(n):
                forward_sink =  forward_sink + forward[l][-1]

                
            decode_matrix[i][k] = forward[k][i] * backward[k][i] / forward_sink 
    
    return decode_matrix









# MAIN

if __name__ == "__main__":
    
    sequence, alphabet, states, transition, emission = viterbi.read10C("rosalind_ba10j.txt")
    
    answer = soft_decoding(sequence, alphabet, states, transition, emission)
    
    
    # turn into string
        
    n = len(states)
    m = len(sequence)
    
    
    str_array = [["" for i in range(n)] for j in range(m)]
                         
    for i in range(m):
        for k in range(n):
            str_array[i][k] = str(answer[i][k])



    # Print out answer
    print('\t'.join(states))
    
    for line in str_array:
        print('\t'.join(line))
        
    
        

    


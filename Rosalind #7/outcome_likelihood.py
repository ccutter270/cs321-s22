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



# FUNCTIONS


def outcome_likelihood(sequence, alphabet, states, transition, emission):
    '''
    Given HMM values (Σ, States, Transition, Emission), and gives
    Pr(x) that a string x was emitted by the HMM
     
    Args:
        sequence: a string of characters
        alphabet: list of characters in which sequence was constucted from
        states: list of states
        transition: 2D transition matrix
        emission: 2D emission matrix 
    
    Returns:
        Integer probability Pr(x) that the HMM emits x
    '''
    
    n = len(states)
    m = len(sequence)
    
    
    # Initialize Matrices
    score = np.zeros((n, m), dtype=float)           
    
    
    # BASE CASE: (Score(source) = 1) * Emission(k, i) * transition(l, k) <- transition is 1/num_states
    for k in range(n):
         first_symbol = alphabet.index(sequence[0])
         score[k][0] = 1.0 * emission[k][first_symbol] * (1.0/ n)

    
    # RECURSIVE CASE: score(k,i) = sum all states l ( scorel,i-1 * weight(l,k, i-1) ) 
    for i in range(1, m):
        for k in range(n):
            for l in range(n):
                
                symbol_at_i = alphabet.index(sequence[i])
                
                score[k][i] = score[k][i] + (score[l][i-1] * transition[l][k] * emission[k][symbol_at_i])   # sum all previous states
                
 
    # Get final sum by adding all values in the last column 
    outcome = 0
    for i in range(n):
        outcome += score[i][-1]
    
    return score, outcome








# MAIN


if __name__ == "__main__":
    
    sequence, alphabet, states, transition, emission = viterbi.read10C("rosalind_ba10d.txt")
    
    score, answer = outcome_likelihood(sequence, alphabet, states, transition, emission)
    
    print(answer)
    

    
    
    
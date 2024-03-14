"""
NAME: Caroline Cutter

CSCI 0321A - Bioinformatics

ASSIGNMENT: Rosalind Randomized Algorithms #7
            Implement the Viterbi Algorithm
            

DUE: Friday, April 15

"""


# IMPORTS
import numpy as np



# FUNCTIONS

# Read File
def read10C(filename):
    '''
    Takes in a file containing A string x, followed by
    the alphabet Σ from which x was constructed,
    followed by the states States, transition matrix Transition,
    and emission matrix Emission of an HMM (Σ, States, Transition, Emission).
     
    Args:
        filename: a string of file name to be read
    
    Returns:
        sequence, alphabet, states, transition[], emission[]
    '''
    with open(filename, "r") as file:
        
        # Sequence
        sequence = file.readline().strip()
        file.readline()
        
        # Alphabet
        alphabet = file.readline().split()
        file.readline()
         
        # States
        states = file.readline().split()
        file.readline()
        file.readline()
     
         # Transition
        transition = [] 
        for i in range(len(states)):
            transition.append(list(map(float, file.readline().split()[1:])))  
        file.readline()
        file.readline()

        # Emission
        emission = []
        for i in range(len(states)):
            emission.append(list(map(float, file.readline().split()[1:])))

        return sequence, alphabet, states, transition, emission






def viterbi(sequence, alphabet, states, transition, emission):
    '''
    Given HMM values (Σ, States, Transition, Emission), and gives
    a string that maximized the unconditional probability
    Pr(x, π) over all possible paths π
     
    Args:
        sequence: a string of characters
        alphabet: list of characters in which sequence was constucted from
        states: list of states
        transition: 2D transition matrix
        emission: 2D emission matrix 
    
    Returns:
        A string of characters that maximizes the (unconditional) probability Pr(x, π) over all possible paths π.
    '''
    
    n = len(states)
    m = len(sequence)
    
    
    # Initialize Matrices
    score = np.zeros((n, m), dtype=float)           
    back = np.zeros((n, m), dtype=int)
    
    
    # BASE CASE: (Score(source) = 1) * Emission(k, i) * transition(l, k) <- transition is 1/num_states
    for k in range(n):
         first_symbol = alphabet.index(sequence[0])
         score[k][0] = 1.0 * emission[k][first_symbol] * (1.0/ n)

    
    # RECURSIVE CASE: score(k,i) = max all states l ( scorel,i-1 * weight(l,k, i-1) ) 
    for i in range(1, m):
        for k in range(n):
            for l in range(n):
                
                symbol_at_i = alphabet.index(sequence[i])
                
                scr = score[k][i]
                
                if scr < score[l][i-1] * transition[l][k] * emission[k][symbol_at_i]:              
                    score[k][i] = score[l][i-1] * transition[l][k] * emission[k][symbol_at_i]    # Set score to max of previous states
                    back[k][i] = l                                                               # Set backtrack to max of states
                     
 
    # BACKTRACK - get sequence
 
    # Find largest value in last column - this is the ending sequence
    end_state = 0            
    max_val = 0
    
    for i in range(n):
        if score[i][m-1] > max_val:
            end_state = i
            max_val = score[i][m-1]
    
    
    # Start from largest value in last column, follow backtrace grid for rest of sequence 
    cur_state = end_state

    path = states[end_state]
    
    # Iterate backwards through the backtrace to get path 
    for i in range(m - 1, 0, -1):
        prev_state = back[cur_state][i]
        path = states[prev_state] + path
        cur_state = prev_state
        
        
    return path








# MAIN
if __name__ == "__main__":
    
    sequence, alphabet, states, transition, emission = read10C("rosalind_ba10c.txt")
    
    answer = viterbi(sequence, alphabet, states, transition, emission)
    
    print(answer)
    

    
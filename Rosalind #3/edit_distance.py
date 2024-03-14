"""
NAME: Caroline Cutter

CSCI 0321A - Bioinformatics

ASSIGNMENT: Rosalind Randomized Algorithms #3
            Compute the Edit Distance Between Two Strings

DUE: Friday, March 11, 2022

"""

# IMPORTS
import itertools
import random
import numpy as np
import global_alignment
from Bio.Align import substitution_matrices
from enum import IntEnum



# FUNCTIONS

def edit_distance(sequence1, sequence2):
    '''
    ** This code is modified from info in 'Bioinformatics' and
       Professor Linderman's Lectures **
    
    Finds the edit distance between two strings
    
    Args:
        sequence1: a string of amino acids
        sequence2: a string of amino acids
    
    Returns:
        an integer of the edit distance between sequence1 and sequence2
        
    '''  
    n = (len(sequence1))
    m = (len(sequence2))

    grid = np.zeros(((n+1), (m+1)), dtype=int)             # initialize grid array
    
    # BASE CASES - one of the strings is empty
    for i in range(n):     # second string is empty - requires i insertions       
        grid[i][0] = i
     
     
    for j in range(m):     # first string is empty - requires j insertions 
        grid[0][j] = j 

    
    # RECURRENCE RELATION
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            
            if(sequence1[i-1] == sequence2[j-1]):      # Matches - don't add one, use diagonal
                grid[i][j] = grid[i-1][j-1]
            
            else:                                      # No Match - add one, plus the minimum score from...
                grid[i][j] = 1 + min(grid[i][j-1],     # vertical score             
                                     grid[i-1][j],     # horizontal score
                                     grid[i-1][j-1])   # diagonal score

  
    # the score is the bottom right (end of both strings)
    return grid[n][m]











# MAIN
if __name__ == "__main__":
    
    sequence1, sequence2 = global_alignment.read5E("rosalind_ba5g.txt")

    answer = edit_distance(sequence1, sequence2)

    print(answer)





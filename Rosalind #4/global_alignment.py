"""
NAME: Caroline Cutter

CSCI 0321A - Bioinformatics

ASSIGNMENT: Rosalind Randomized Algorithms #3
            Find a Highest-Scoring Alignment of Two Strings
            (Global Alignment)

DUE: Friday, March 11, 2022

"""

# IMPORTS
import itertools
import random
import numpy as np
from Bio.Align import substitution_matrices
from enum import IntEnum


# GLOBALS 
sub_mat = substitution_matrices.load("BLOSUM62")
penalty = 5


# CLASSES 
class Back(IntEnum):
    MAT = 0
    VRT = 1
    HRZ = 2
    
    
    

# READ FILE 
def read5E(filename):
    '''
    Takes in a file containing two strings and returns them as seprate strings
     
    Args:
        filename: a string of file name to be read
    
    Returns:
        Two strings - sequence1, sequence2
    '''  
    with open(filename, "r") as file:

        sequence1 = str(file.readline().strip())   # reads first sequence
        sequence2 = str(file.readline().strip())   # reads second string

    return sequence1, sequence2
        
        
  
  
  
  
def global_alignment(sequence1, sequence2):
    '''
    ** This code is modified from info in 'Bioinformatics' **
    
    Finds the highest-scoring alignment between two strings using a scoring matrix
    
    Args:
        sequence1: a string of amino acids
        sequence2: a string of amino acids
    
    Returns:
        a list of highest score of alignments between sequence 1 & 2 and the alignment
        [score, alignment of sequence1, alignment of sequence2]
        
    '''
    n = (len(sequence1))
    m = (len(sequence2))
    
    
    # MAKE THE GRID FOR THE ALIGNMENTS 
    
    grid = np.zeros(((n+1), (m+1)), dtype=int)             # initialize grid and backtrack array
    back = np.zeros((n + 1, m + 1), dtype=int)
    
    # BASE CASES
    for i in range(n):            
        grid[i+1][0] = grid[i][0] - penalty                # account for gaps moving down
     
     
    for j in range(m):
        grid[0][j+1] = grid[0][j] - penalty                # account for gaps moving right
     

   # RECURRENCE RELATION
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            
            # Find scores of diagonal, vertical and horizontal
            diag = grid[i-1][j-1] + sub_mat[(sequence1[i-1], sequence2[j-1])]
            vert = grid[i-1][j] - penalty
            horz = grid[i][j-1] - penalty
            
            
            incoming = [diag, vert, horz]
            
            back[i][j] = max_edge = np.argmax(incoming)      # set backtrack as either diag, vert, horz
            grid[i][j] = incoming[max_edge]       # set grid as the max of the scores

            


    # FIND SCORE - at bottom right of grid
    score = grid[n][m] 
    
    # Initialize variables for backtracking
    i = n 
    j = m    
    align1 = ""
    align2 = ""
    
    # Iterate through backtracking array 
    while(i != 0 and j != 0):

       # DIAGONAL BACKTRACK
        if(back[i][j] == 0):
            i -= 1
            j -= 1
            align1 = sequence1[i] + align1
            align2 = sequence2[j] + align2
        
        # VERTICAL BACKTRACK
        elif(back[i][j] == 1):
            i -= 1
            align1 = sequence1[i] + align1
            align2 = "-" + align2
          
        # HORIZONTAL BACKTRACK
        else:
            j -= 1
            align1 = "-" + align1
            align2 = sequence2[j] + align2
           
    
    # LEFTOVER DASHES 
     
    # Vertical leftover dashes
    while(i != 0):
        i -= 1
        align1 = sequence1[i] + align1
        align2 = "-" + align2
    
    # Horizontal leftover dashes 
    while(j != 0):   
        j -= 1
        align1 = "-" + align1
        align2 = sequence2[j] + align2
        

    return score, align1, align2





# MAIN

if __name__ == "__main__":
    
    sequence1, sequence2 = read5E("rosalind_ba5e.txt")

    answer = global_alignment(sequence1, sequence2)

    print(answer[0])
    print(answer[1])
    print(answer[2])

        
        
        
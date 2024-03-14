"""
NAME: Caroline Cutter

CSCI 0321A - Bioinformatics

ASSIGNMENT: Rosalind Randomized Algorithms #4
            Align Two Strings Using Affine Gap Penalties
            (Alignment with Affine Gap Penalties Problem)

DUE: Friday, March 18, 2022

"""

# IMPORTS
import numpy as np
import global_alignment
from Bio.Align import substitution_matrices
import itertools


# GLOBALS 
sub_mat = substitution_matrices.load("BLOSUM62")
opening_penalty = 11
extension_penalty = 1
infinity = float('inf')



def affine_align(sequence1, sequence2):
    '''
    ** This code is modified from info in 'Bioinformatics' and
       Professor Linderman's Lectures **
    
    Construct a highest-scoring global alignment of two strings (with affine gap penalties).
    
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
    
    lower = np.zeros(((n+1), (m+1)), dtype=float)          # initialize upper, middle, lower grids
    middle = np.zeros(((n+1), (m+1)), dtype=float)           
    upper = np.zeros(((n+1), (m+1)), dtype=float)            
    
    # MAKE BACKTRACK ALIGNMENTS
    back_lower = []
    back_middle = []
    back_upper = []
    

    # BASE CASES - for all three grids 
    for i in range(n + 1):
        back_lower.append([0] * (m + 1))
        back_middle.append([0] * (m + 1))
        back_upper.append([0] * (m + 1))
    
    
    
    # RECURSIVE CASE
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            
            # SET LOWER & UPPER  
            
            # Lower 
            lower[i][j] = max(lower[i-1][j] - extension_penalty, middle[i-1][j] - opening_penalty)
            
            # Lower backtrack
            if lower[i][j] == lower[i-1][j] - extension_penalty:
                back_lower[i][j] = (0, i-1, j)
            
            else:
                back_lower[i][j] = (1, i-1, j)
            
            
            
            # Upper 
            upper[i][j] = max(upper[i][j-1] - extension_penalty, middle[i][j-1] - opening_penalty)
            
            # Upper backtrack 
            if upper[i][j] == upper[i][j-1] - extension_penalty:
                back_upper[i][j] = (2, i, j-1)
            else:
                back_upper[i][j] = (1, i, j-1)

            
            
            # MIDDLE
            score = sub_mat[(sequence1[i-1], sequence2[j-1])]      # match or mismatch 
            
            
            middle[i][j] = max(lower[i][j],
                                middle[i-1][j-1] + score,
                                upper[i][j])
            
            
            # Middle Backtrack
            if middle[i][j] == lower[i][j]:
                back_middle[i][j] = (0, i, j)
            elif middle[i][j] == upper[i][j]:
                back_middle[i][j] = (2, i, j)
            else:
                back_middle[i][j] = (1, i-1, j-1)
                    

    # FIND SOLUTION
    
    # NOTE:
        # Upper = 0
        # Middle = 1
        # Upper = 2
    
    align1 = ""
    align2 = ""
    
    trans = back_middle[n][m]  # Start at middle matrix
    current = 1                     
    
    score = int(middle[-1][-1])     # score is bottom right of middle matrix 
    
    while trans != 0:
            
        # Lower Matrix 
        if current == 0:
            if trans[0] == 0:
                align1 = sequence1[trans[1]-1] + align1
                align2 = "-" + align2
                trans = back_lower[trans[1]][trans[2]]
            else:
                trans = back_middle[trans[1]][trans[2]]
                current = 1
                
        # Middle Matrix 
        elif current == 1:
            if trans[0] == 1:
                align1 = sequence1[trans[1]] + align1
                align2 = sequence2[trans[2]] + align2
                trans = back_middle[trans[1]][trans[2]]
            elif trans[0] == 0:
                align1 = sequence1[trans[1]-1] + align1
                align2 = "-" + align2
                trans = back_lower[trans[1]][trans[2]]
                current = 0
            else:
                align1 = "-" + align1
                align2 = sequence2[trans[2]-1] + align2
                trans = back_upper[trans[1]][trans[2]]
                current = 2
        
        # Upper Matrix
        else:
            if trans[0] == 2:
                align1 = "-" + align1
                align2 = sequence2[trans[2]-1] + align2
                trans = back_upper[trans[1]][trans[2]]
            else:
                trans = back_middle[trans[1]][trans[2]]
                current = 1
                
                
    return score, align1, align2

    
  
  
  
  
  
  
  
  
# MAIN
if __name__ == "__main__":
    
    sequence1, sequence2 = global_alignment.read5E("rosalind_ba5j.txt")

    answer = affine_align(sequence1, sequence2)
    
    print(answer[0])
    print(answer[1])
    print(answer[2])



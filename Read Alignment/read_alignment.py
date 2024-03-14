"""
NAME: Caroline Cutter

CSCI 0321A - Bioinformatics

ASSIGNMENT: Read Alignment 

DUE: Friday, March 18, 2022

"""



# IMPORTS
import numpy as np
import pysam



# FUNCTIONS


# Note - this is the same as the rosalind assignment, but added an output for position of
#        where the alignment starts 

def fitting_alignment(sequence1, sequence2):
    '''
    ** This code is modified from info in 'Bioinformatics' and
       Professor Linderman's Lectures **
    
    Construct a highest-scoring fitting alignment between two strings.
    
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
    grid = np.zeros(((n+1), (m+1)), dtype=int)             # initialize grid and backtrack arrays
    back = np.zeros((n + 1, m + 1), dtype=int)
 
    # BASE CASES - n = 0 and m = 0 --> grid = 0 (taken care of in code above)

    # RECURRENCE RELATION
    for i in range(1, n + 1):
        for j in range(1, m + 1):
 
            # Find scores of diagonal, vertical and horizontal
            
            
            if (sequence1[i-1] == sequence2[j-1]):     # match (+1)
                diag = grid[i-1][j-1] + 1
            else:                                      # mismatch (-1) 
                diag = grid[i-1][j-1] - 1

            vert = grid[i-1][j] - 1                    # indel (-1)
            horz = grid[i][j-1] - 1
                  
            
            incoming = [diag, vert, horz]
            
            back[i][j] = max_edge = np.argmax(incoming)      # set backtrack as either diag, vert, horz
            grid[i][j] = incoming[max_edge]                  # set grid as the max of the scores
            
            
    # FIND SCORE - the highest score from last column 
    
    # Search for max in last column
    row = 0
    max_val = 0
    
    for i in range(n):
        if grid[i][m] >= max_val:
            row = i
            max_val = grid[i][m]
    
    
    score = grid[row][m]               # score = max of last column 
    
    i = row                            # Initialize location to start backtracking 
    j = m
    

    # CREATE SEQUENCES (using backtracking  
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
           
    
    return align2, i     # i returns the starting position for the alignment of the shorter sequence relative to the longer sequence.






def read_align(sequence1, sequence2):
    answer, start = fitting_alignment(sequence1, sequence2)
    
    n = len(sequence1)
    m = len(answer)
    
    filler = n - m - start

    answer = (start * "-") + answer + filler * "-"
    
    return answer






# MAIN
if __name__ == "__main__":
    
    # TESTER PART 1
#     answer = read_align('TTAGTAGGCTTAAGGTTA', 'TAGATA')
#     print(answer)
    
    sequence1 = "AATCTGCCTTGTGGTCGGCTCTTACCTTCAGGCTGCTCTGAGCCCAGAGCAGAATGGTCATCACAGCTCTCCTCAACTTGGCATTGCCTGAGATCAGGATGGCTGCATGCCCAGAGGGACAAGCTGCCATTATCCCAACACAAACCATCACCCCTATTTTGTCGCGCCACAGAATCAGTAGGGGCACAGAGATGAAGGCAGC"
    
    # PART 2: WORKING WITH REAL DATA
    with pysam.AlignmentFile("NA12878-ready.sample.bam", mode="rb") as bam_reader:
        for read in bam_reader:
            # Operations on each read
            answer = read_align(sequence1, read.query_sequence)
            print(answer)





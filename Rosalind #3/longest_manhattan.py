"""
NAME: Caroline Cutter

CSCI 0321A - Bioinformatics

ASSIGNMENT: Rosalind Randomized Algorithms #3
            Find the Length of a Longest Path in a Manhattan-like Grid

DUE: Friday, March 11, 2022

"""

# IMPORTS
import itertools
import random
import numpy as np



# READ FILE 
def read5B(filename):
    '''
    ** This code was modified from Professor Linderman's implemetation suggestions **
    Takes in a file containing n, m, n × (m+1) matrix Down and an (n+1) × m matrix Right
     
    Args:
        filename: a string of file name to be read
    
    Returns:
        A list [n, m, [Down 2D matrix], [Right 2D matrix]] 
    '''   
    with open(filename, "r") as file:

        numbers = str(file.readline().strip())   # this block splits up the first line
        split_numbers = numbers.split(" ")       # and assigns integers to n and m accordingly 
        n = int(split_numbers[0])
        m = int(split_numbers[1])

        # Create "down" array 
        down = []
        for i in range(n):                                                 # down array is n lines
            line = [int(num) for num in file.readline().split()]           # create an list of numbers 
            down.append(line)                                              # adds array of numbers to down array

        file.readline() # reads the "-" out
                    
        # Create the "right" array
        right = []
        for i in range(n+1):                                               # right array has n+1 lines
            line = [int(num) for num in file.readline().split()]
            right.append(line)

        return n, m, down, right


# Functions

def ManhattanTourist(n, m, down, right):
    '''
    ** This code is modified from pseudocode from pg 245 of 'Bioinformatics' and
        Professor Linderman's pesudocode from lecture (which helped with indexing **
    
    Finds the length of a longest path in a rectangular city
    
    Args:
        n: integer of rows of grid
        m: integer of columns of rid
        down: 2D array of lengths of down paths 
        right: 2D array of lengths of right paths 
    
    Returns:
        the length of the longest path in a grid 
    '''
    
    # 1. S(0,0) <- 0       S is 2D array
    S = np.zeros((n+1, m+1), dtype=int)
    
    # 2. for i <- 1 to n
    for i in range(1, n + 1):
    
        # 3. S(i, 0) <- S(i - 1, 0) + down(i - 1 , 0)
        S[i][0] = S[i - 1][0] + down[i - 1][0]
        
    # 4. for j <- 1 to m
    for j in range(1, m + 1):
    
        # 5.  S(0, j) = S(0, j - 1) + right(0, j -1 )
        S[0][j] = S[0][j - 1] + right[0][j - 1]
        
    # 6. for i <- 1 to n
    for i in range(1, n + 1):
    
        # 7. for j <- 1 to m
        for j in range(1, m + 1):
        
            # 8. S(i, j) = max(S(i-1, j) + down(i, j), S(i, j - 1) + right(i, j)
            S[i][j] = max((S[i-1][j] + down[i-1][j]), (S[i][j-1] + right[i][j-1]))
            
    # 9. Return S(n, m)
    return S[n][m]








# MAIN

if __name__ == "__main__":
    n, m, down, right = read5B("rosalind_ba5b.txt")

    answer = ManhattanTourist(n, m, down, right)
    
    print(answer)

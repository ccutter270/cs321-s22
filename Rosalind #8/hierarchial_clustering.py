"""
NAME: Caroline Cutter

CSCI 0321A - Bioinformatics

ASSIGNMENT: Rosalind Randomized Algorithms #8
            Implement Hierarchical Clustering
            

DUE: Friday, April 29

"""


# IMPORTS
import numpy as np
import lloyd
from scipy.spatial import distance


# READ FILE
def read8E(filename):
    '''
    An integer n, followed by an nxn distance matrix
     
    Args:
        filename: a string of file name to be read
    
    Returns:
        an integer n and a numpy matrix of size n x n
    '''
    with open(filename, "r") as file:
       
        n = int(file.readline())         # get n
        matrix = np.loadtxt(file)        # get n x n matrix
 
        return n, matrix
 
 
 
 
 
 
 def hierarchial_clustering(n, matrix):
     
     # 1.  Clusters ← n single-element clusters labeled 1, ... , n
     
     # 2. construct a graph T with n isolated nodes labeled by single elements 1, ... , n
 
 






# MAIN
if __name__ == "__main__":
    
    n, matrix = read8E("rosalind_ba8e.txt")

    print(n)
    print(matrix)
    
    


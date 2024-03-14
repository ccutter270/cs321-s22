"""
NAME: Caroline Cutter

CSCI 0321A - Bioinformatics

ASSIGNMENT: Rosalind Randomized Algorithms #8
            Implement the Soft k-Means Clustering Algorithm
            

DUE: Friday, April 29

"""


# IMPORTS
import numpy as np
import lloyd
from scipy.spatial import distance



# READ FILE
def read8D(filename):
    '''
    Integers k and m, followed by a stiffness parameter β, followed by a set of points Data in m-dimensional space.
     
    Args:
        filename: a string of file name to be read
    
    Returns:
        integer k, integer m, float B, list of points
    '''
    with open(filename, "r") as file:
       
       
        integers = file.readline().strip().split()
        
        k = int(integers[0])
        m = int(integers[1])
        
        b = float(file.readline())
        
        points = [] 
        
        for line in file:
            points.append([float(coord) for coord in line.strip().split()]) 
 
    
        return k, m, b, points
    




def soft_clustering(k, m, b, points):
    '''
    Implement the Soft k-Means Clustering Algorithm
    
    Args:
        k: number of points
        m: dimensions
        b: the stiffness parameter
        points: list of points 
    
    Returns:
        a list of centers after runnin soft k-means clustering algorithm 
    
    


    '''
    
    # Generate k random centers
    centers = []
    for i in range(k):
        centers.append(points[i])
        

    # Create hidden Matrix 
    hidden_matrix = np.zeros((k, len(points)))
    
    
    # Repeat 100 times (specified by rosalind) 
    for iterator in range(100): 
    
        # initialize new centers
        new_centers = np.zeros(m)
    
        # CENTERS TO SOFT CLUSTERS 
        for i in range(len(points)):
            
            distances = []
            
            for n in range(k):
                
                distances.append(distance.euclidean(points[i], centers[n]))
            
    
            data = np.exp(-b*np.array(distances))
            data /= np.sum(data)
        
            hidden_matrix[:,i] = data
            
    
        # SOFT CLUSTERS TO CENTERS      
        new_centers = np.matmul(hidden_matrix, points)
        denom =  np.sum(hidden_matrix, axis = 1)
        new_centers /= denom[:, np.newaxis]
    
        centers = new_centers

    # Return best centers after 100 iterations 
    return centers

    
    
        




# MAIN
if __name__ == "__main__":
    
    k, m, b, points = read8D("rosalind_ba8d.txt")

    centers = soft_clustering(k, m, b, points)
    
    for center in centers:
        print(' '.join(str(point) for point in center))
   
    
    
    


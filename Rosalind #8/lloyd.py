"""
NAME: Caroline Cutter

CSCI 0321A - Bioinformatics

ASSIGNMENT: Rosalind Randomized Algorithms #8
            Implement the Lloyd Algorithm for k-Means Clustering
            

DUE: Friday, April 29

"""


# IMPORTS
import numpy as np



# READ FILE
def read8C(filename):
    '''
    Integers k and m followed by a set of points Data in m-dimensional space.
     
    Args:
        filename: a string of file name to be read
    
    Returns:
        integer k, integer m, list of points
    '''
    with open(filename, "r") as file:
       
       
        integers = file.readline().strip().split()
        
        k = int(integers[0])
        m = int(integers[1])
        
        points = [] 
        
        for line in file:
            points.append([float(coord) for coord in line.strip().split()]) 
 

    
        return k, m, points


 
def dist(point1, point2, m):
    '''
    Finds distance between points in m dimensional space 
     
    Args:
        point1: list of coordinates of point 1
        point2: list of coordinates of point 2
        m: integer dimensional space
    
    Returns:
        float of euclidean distance between point1 and point2
    '''
    total = 0
    
    for i in range(m):
        total += (point1[i] - point2[i])**2
    
    return total



def closest_point(point, centers, k, m):
    '''
    Finds closest point from a point to centers
     
    Args:
        point: list of coordinates of a point
        centers: a list of points
        k: number of centers
        m: integer dimensional space
    
    Returns:
        list of the coordinates for the best point
    
    '''
    
    best_coord = -1
    best_dist = float("inf")
    
    
    for i in range(k):
        distance = dist(point, centers[i], m)
        
        if distance < best_dist:
            best_coord = i
            best_dist = distance
            
    print(best_coord)
            
    return best_coord


def cluster_to_center(index, indices, m, points):
    
    
    '''
    Transforms clusters to centers
     
    Args:
        index: index of a point
        indices: indicies of centers
        m: integer dimensional space
        points: list of points 
    
    Returns:
        returns the max point of the clusters to make new centers

    '''
    
    total = 0
    coords = np.zeros(m)

    
    for i in range(len(points)):
        
        if indices[i] == index:
            total += 1
            for j in range(m):
                coords[j] += points[i][j]
    
    max_point = [c/max(total,1) for c in coords]
       
    return max_point
        
    



def lloyds(k, m, points):
    '''
    Implements lloyds clustering algorithm
     
    Args:
        k: number of points 
        m: integer dimensional space
        points: list of points 
    
    Returns:
        the centers of the clusters 

    '''

    # Select k arbitrary data points as Centers
    centers = []
    for i in range(k):
        centers.append(points[i])

    
    difference = float("inf")
    
    while difference > 0:
        
        indices = [closest_point(point, centers, k, m) for point in points]
    
        # Centers to Clusters
        new_centers = [cluster_to_center(index, indices, m, points) for index in range(k)]
        
        # Clusters to Centers
        difference = sum([dist(new_centers[i], centers[i], m) for i in range(k)] )
        
        centers = new_centers[:]

    
    return centers
    
    




# MAIN
if __name__ == "__main__":
    
    k, m, points = read8C("rosalind_ba8c.txt")
    
    centers = lloyds(k, m, points)
    
    for center in centers:
         print(' '.join(str(point) for point in center))
    
    
    
    
    

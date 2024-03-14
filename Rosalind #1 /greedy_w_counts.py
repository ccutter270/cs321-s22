"""
NAME: Caroline Cutter
COLLABORATORS: None

CSCI 0321A - Bioinformatics

ASSIGNMENT: Rosalind Randomized Algorithms #1
            Greedy Algorithm with Pseudocounts

DUE: Friday, Feb 25, 2022

"""

# IMPORTS
import itertools
import numpy as np
import median_string
import profile_most
import greedy


# GLOBALS
infinity = float('inf')
DNA = "ACGT"



def profile_matrix_pseudocounts(motifs):
    '''
    Creates a profile matrix (with pseudocounts) from a given list of motifs
     
    Args:
        motifs: list of dna sequences
    
    Returns:
        A 2D list profile of probabilities from given motifs 
    '''
    
    k = len(motifs[0])          # length of sequences
    t = len(motifs)             # number of motifs 

    
    profile = np.full((4,k), 1, dtype=np.float_)     # initialize profile array with float 1's  (for pseudo)
    

    # CREATE COUNT MATRIX - with the number of nucleotides in sequences
    
    for motif in motifs:                             # runs through each sequence in motif                        
        
        i = 0
        for nucleotide in motif:                     # runs through each nucleotide in sequence
            
            profile[DNA.find(nucleotide)][i] += 1    # adds count to [nucleotide #][positition in sequence] 
            
            i += 1 

     # CHANGE COUNTS TO PROBABILITIES  
    for i in range(4):
        for j in range(k):
            profile[i][j] = profile[i][j] / (t + t)  # divide count accounting for pseudocounts
            
    return profile 
    


def greedy_pesudocounts(dna, k, t):
    '''
    ** NOTE: this is from a combination of pesuedocodes from Rosalind & 'Bioinformatic Algorithms'**
    Uses a greedy approach to find the collection of best motifs from a given set of DNA with pseudocounts 
     
    Args:
        dna: list of dna sequences
        k: integer length of desired kmers
        t: integer of sequences in DNA 
    
    Returns:
        A list of best motifs from a given set of DNA using the greedy approach 
    '''
    
    # 1. BestMotifs ← motif matrix formed by first k-mers in each string from Dna
    best_motifs = greedy.create_best_motifs(dna, k)
    
    
    # 2. BestScore <- INF  & Initialize motifs matrix
    best_score = infinity
    motifs = []
    
   
    # 3. For each k-mer in DNA-0 
    kmers = median_string.unique_kmer_composition(dna[0], k)      # create kmer list from DNA-0
    
    for kmer in kmers:
        

        # 4. Motifs-0 <- k-mer
        motifs = [kmer]
        
        # 5. For i in 1 to t
        for i in range(1, t):
            
            # 6. form Profile from Motifs0 .... Motifs i-1
            profile = profile_matrix_pseudocounts(motifs)

            # 7. Motifs-i <- Prfile-most probable k-mer in DNA-i  
            motifs.append(profile_most.profile_most_probable(dna[i], k, profile))
       
       
        # 8. if Score(Motifs < BestScore
        if greedy.score(motifs, profile) < best_score:
           
            # 9. BestScore <- Score(Motifs) &  BestMotifs <- Motifs
            best_score = greedy.score(motifs, profile)
            best_motifs = motifs
    
    return best_motifs
    


# MAIN 

if __name__ == "__main__":
    k, t, dna = greedy.read2D("rosalind_ba2e.txt")
    answer = greedy_pesudocounts(dna, k, t)
    for string in answer:
        print(string)

"""
NAME: Caroline Cutter
COLLABORATORS: None

CSCI 0321A - Bioinformatics

ASSIGNMENT: Rosalind Randomized Algorithms #1
            Greedy Algorithm

DUE: Friday, Feb 25, 2022

"""

# IMPORTS
import itertools
import numpy as np
import median_string
import profile_most


# GLOBALS
infinity = float('inf')
DNA = "ACGT"




# READ FILES 
def read2D(filename):
    '''
    ** This code was taken from Professor Linderman's implemetation suggestions **
    Takes in a file containing k (kmer length), t (number of sequences), DNA (dna sequences)
     
    Args:
        filename: a string of file name to be read
    
    Returns:
        A list [k, t, [sequneces]] 
    '''
    with open(filename, "r") as file:

        numbers = str(file.readline().strip())   # this block splits up the first line
        split_numbers = numbers.split(" ")       # and assigns integers to k and t accordingly
        k = int(split_numbers[0])
        t = int(split_numbers[1])

        kmers = [line.strip() for line in file.readlines()]  # this line takes the sequences and
                                                             # adds them to a list of kmers
        return k, t, kmers




def create_best_motifs(dna, k):
    '''
    Takes a list of dna sequences and makes a list of the kmers of length
    k from the beginning of each sequence 
     
    Args:
        dna: list of dna seuqences
        k: integer of desired length of kmers
    
    Returns:
        A list of kmers of length k from the beginning of each sequence 
    '''
    
    motif_matrix = []
    
    for sequence in dna:                    # creates kmer from each sequence 
        motif_matrix.append(sequence[0:k])  # adds kmer to matrix to be returned 
        
    return motif_matrix




def profile_matrix(motifs):
    '''
    Creates a profile matrix (without pseudocounts) from a given list of motifs
     
    Args:
        motifs: list of dna sequences
    
    Returns:
        A 2D list profile of probabilities from given motifs 
    '''
    
    k = len(motifs[0])          # length of sequences
    t = len(motifs)             # number of motifs 

    
    profile = np.full((4,k), 0, dtype=np.float_)     # initialize profile array with float 0's 
    

    # CREATE COUNT MATRIX - with the number of nucleotides in sequences
    
    for motif in motifs:                             # runs through each sequence in motif                        
        
        i = 0
        for nucleotide in motif:                     # runs through each nucleotide in sequence
            
            profile[DNA.find(nucleotide)][i] += 1    # adds count to [nucleotide #][positition in sequence] 
            
            i += 1 

     # CHANGE COUNTS TO PROBABILITIES  
    for i in range(4):
        for j in range(k):
            profile[i][j] = profile[i][j] / t       # divide count by how many sequences there are (t) 
            
    return profile 
    

    
def score(motifs, profile):
    '''
    Scores motifs from a given profile compared to the consensus string
     
    Args:
        motifs: list of dna sequences
        profile: 2D list of probabilities of each nucleotide at a position 
    
    Returns:
        An integer score of motifs compared to the consensus 
    '''
    
    # FIND CONSENSUS STRING - this was from Prof. Lindermans implementation suggestions
    consensus = ""
    for i in np.argmax(profile, axis=0):
        consensus += DNA[i]
        
        
    # SCORE - compare each sequence to consensus and get the sum total
    score = median_string.distance_to_all_sequences(consensus, motifs)
    
    return score
  

def greedy(dna, k, t):
    '''
    ** NOTE: this is from a combination of pesuedocodes from Rosalind & 'Bioinformatic Algorithms'**
    Uses a greedy approach to find the collection of best motifs from a given set of DNA 
     
    Args:
        dna: list of dna sequences
        k: integer length of desired kmers
        t: integer of sequences in DNA 
    
    Returns:
        A list of best motifs from a given set of DNA using the greedy approach 
    '''
    
    # 1. BestMotifs ← motif matrix formed by first k-mers in each string from Dna
    best_motifs = create_best_motifs(dna, k)
    
    
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
            profile = profile_matrix(motifs)

            # 7. Motifs-i <- Prfile-most probable k-mer in DNA-i  
            motifs.append(profile_most.profile_most_probable(dna[i], k, profile))
       
       
        # 8. if Score(Motifs < BestScore
        if score(motifs, profile) < best_score:
           
            # 9. BestScore <- Score(Motifs) &  BestMotifs <- Motifs
            best_score = score(motifs, profile)
            best_motifs = motifs
    
    return best_motifs
    

# MAIN 

if __name__ == "__main__":
    k, t, dna = read2D("rosalind_ba2d.txt")
    answer = greedy(dna, k, t)
    for string in answer:
        print(string)

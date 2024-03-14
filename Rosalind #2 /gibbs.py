"""
NAME: Caroline Cutter

CSCI 0321A - Bioinformatics

ASSIGNMENT: Rosalind Randomized Algorithms #2
            Gibbs Sampler

DUE: Friday, March 4, 2022

"""

# IMPORTS
import itertools
import random
import numpy as np
import median_string
import profile_most
import greedy
import greedy_w_counts as gr
import randomized



# GLOBALS
infinity = float('inf')
DNA = "ACGT"




# READ FILES
def read2G(filename):
    '''
    ** This code was modified from Professor Linderman's implemetation suggestions **
    Takes in a file containing k (kmer length), t (number of sequences), N (
    DNA (dna sequences)
     
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
        N = int(split_numbers[2])

        kmers = [line.strip() for line in file.readlines()]  # this line takes the sequences and
                                                             # adds them to a list of kmers
        return k, t, N, kmers





# FUNCTIONS

def profile_randomly_generated(sequence, k, profile):
    '''
    ** This code uses modified snippets from Professor Linderman's implementation suggestions  **
    Generates a profile-random kmer from a dna sequece
     
    Args:
        sequence: a string of dna sequence
        k: integer length of desired kmer
        profile: a 2D array of probabilities for nucleotides in a sequence of length k
    
    Returns:
        a string of the profile-random kmer from the sequence 
    '''
    
    possible_kmers = median_string.unique_kmer_composition(sequence, k)   # create all possible kmers from sequence

    kmer_probs = []
         
    for kmer in possible_kmers:                                           # find the probability of each kmer and add to a list
        kmer_probs.append(profile_most.kmer_probability(kmer, profile))
        
    kmer_probs /= np.sum(kmer_probs)                                      # normalize the probabilities so the sum = 1
    
    rand = np.digitize(np.random.random(), np.cumsum(kmer_probs))         # pick a random value based off the probability distribution
    
    return possible_kmers[rand]                                           # return the random sequnece 



def gibbs_sampler(dna, k, t, N):
    '''
    ** NOTE: this is from the pesuedocode from Rosalind & 'Bioinformatic Algorithms'**
    Uses the gibbs sampler approach to find the collection of best motifs from a given set of DNA 
     
    Args:
        dna: list of dna sequences
        k: integer length of desired kmers
        t: integer of sequences in DNA
        N: number of iterations 
    
    Returns:
        A list of best motifs from a given set of DNA using the gibbs sampler approach 
    '''

    # 1. Randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
    motifs = randomized.randomly_select_kmers(dna, k)

    # 2. BestMotifs ← Motifs
    best_motifs = motifs[:] 
    
    # 3. for j ← 1 to N
    for j in range(1, N):
    
        # 4. i ← Random(t)
        i = random.randint(0, t-1)
        
        # 5. Profile ← profile matrix constructed from all strings in Motifs except for Motif-i
        motifs_without_i = motifs[:i] + motifs[i+1:] 
        profile = gr.profile_matrix_pseudocounts(motifs_without_i)
        
        # 6. Motif-i ← Profile-randomly generated k-mer in the i-th sequence
        motifs[i] = profile_randomly_generated(dna[i], k, profile)

        # 7. if Score(Motifs) < Score(BestMotifs)
        if randomized.score(motifs, k) < randomized.score(best_motifs, k): 

            # 8. BestMotifs ← Motifs
            best_motifs = motifs[:]
              
    #9. return BestMotifs
    return best_motifs
      
 
 


# MAIN

if __name__ == "__main__":
    k, t, N, dna = read2G("rosalind_ba2g.txt")
    
    
    best_m = gibbs_sampler(dna, k, t, N)
    best_score = randomized.score(best_m, k)
    
    i = 0
    
    while i < 19:
        
        new_motifs = gibbs_sampler(dna, k, t, N)
        scr = randomized.score(new_motifs, k)
        
        if scr < best_score:
            best_m = new_motifs
            best_score = scr 
        
        i += 1
      
    for m in best_m:
        print(m) 

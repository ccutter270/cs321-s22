"""
NAME: Caroline Cutter

CSCI 0321A - Bioinformatics

ASSIGNMENT: Rosalind Randomized Algorithms #2
            Randomized Motif Search

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




# GLOBALS
infinity = float('inf')
DNA = "ACGT"




# IMPORTS
import itertools
import random
import numpy as np
import median_string
import profile_most
import greedy
import greedy_w_counts as gr




# GLOBALS
infinity = float('inf')
DNA = "ACGT"



# FUNCTIONS



def randomly_select_kmers(dna, k): 
    '''
    Takes in a list of sequences from dna and randomly selects k-mers
    Motifs = (Motif1, …, Motift) in each string from Dna
     
    Args:
        dna: list of dna seuqences
        k: integer of desired length of kmers
    
    Returns:
        A list of random kmers of length k from each of the DNA sequences
    '''
    motif_matrix = []
    
    for sequence in dna:                    # creates kmer from each sequence
        
        rand = random.randrange(0, len(sequence) - k) # ** This could be a different
        
        motif_matrix.append(sequence[rand:rand + k])  # adds kmer to matrix to be returned 
        
    return motif_matrix




"""
def profile_random_kmers(dna, k):
    '''
    Takes in a list of sequences from dna and randomly selects k-mers
    Motifs = (Motif1, …, Motift) in each string from Dna from a probability
    distribution of the kmers 
     
    Args:
        dna: list of dna seuqences
        k: integer of desired length of kmers
    
    Returns:
        A list of random kmers of length k from each of the DNA sequences
    '''
    
    random_kmers = []
    
    for sequence in dna:
          
        # get possible kmers from this sequnce 
        possible_kmers = median_string.unique_kmer_composition(sequence, k)
        
        profile = gr.profile_matrix_pseudocounts(possible_kmers)

        kmer_probs = []
        

        for kmer in possible_kmers:
            kmer_probs.append(profile_most.kmer_probability(kmer, profile))
        
        kmer_probs /= np.sum(kmer_probs)
          
        rand = np.random.choice(range(len(kmer_probs)), p=kmer_probs)

        random_kmers.append(possible_kmers[rand])
    
    return random_kmers
    
"""




def profile_to_motifs(k, dna, profile):
    '''
    Takes in a profile and creates a list of profile most probable sequences
    of length k from a list of dna sequences 
     
    Args:
        k: integer of desired length of kmers
        dna: list of dna seuqences
        profile: a 2D array of provavilities of nucleotide in a dna sequences

    Returns:
        A list of profile most probable motifs 
    '''
    
    motifs = []                   # inititalize motifs array 
    
    for sequence in dna:          # for each dna sequence
                 
                                  # adds the profile most probable string to motifs  
        motifs.append(profile_most.profile_most_probable(sequence, k, profile))  
        
    return motifs





def consensus(motifs, k):
    '''
    ** This function was modified from Professor Linderman's implementation suggestions**
    
    Finds the consensus string of length k from a list of motifs 
     
    Args:
        motifs: list of dna motifs
        k: integer of desired length of consensus string

    Returns:
        A consensus string of length k 
    '''
    counts = np.zeros((len(DNA), k), dtype=np.int_)
    for motif in motifs:
        for i in range(len(motif)):
            counts[DNA.find(motif[i]), i] += 1

    consensus = ""
    for i in np.argmax(counts, axis=0):
        consensus += DNA[i]

    return consensus



def score(motifs, k):
    '''
    Finds the score of a motif matrix with dna strings of length k 
     
    Args:
        motifs: list of dna motifs
        k: integer of desired length of consensus string

    Returns:
        The score of motifs 
    '''
    con = consensus(motifs, k)                 # find consensus string 
    
                                               # calculate score compared to consensus 
    score = median_string.distance_to_all_sequences(con, motifs)  
    
    return score




def randomized_motif_search(dna, k, t):
    '''
    ** NOTE: this is from the pesuedocode from Rosalind & 'Bioinformatic Algorithms'**
    Uses a randomized approach to find the collection of best motifs from a given set of DNA 
     
    Args:
        dna: list of dna sequences
        k: integer length of desired kmers
        t: integer of sequences in DNA 
    
    Returns:
        A list of best motifs from a given set of DNA using a randomized approach 
    '''
    
    # 1. randomly select k-mers Motifs = (Motif1, …, Motift) in each string from Dna
    motifs = randomly_select_kmers(dna, k)
    
    # 2. BestMotifs ← Motifs
    best_motifs = motifs
    
    # 3. while forever
    while True:
        
        # 4. Profile ← Profile(Motifs)
        profile = gr.profile_matrix_pseudocounts(motifs)
        
        # 5. Motifs ← Motifs(Profile, Dna)
        motifs = profile_to_motifs(k, dna, profile)
        
        
        #best_profile = gr.profile_matrix_pseudocounts(best_motifs)
        
        # 6. if Score(Motifs) < Score(BestMotifs)
        if score(motifs, k) < score(best_motifs, k): 
              
              # 7. BestMotifs ← Motifs
              best_motifs = motifs
        
        #  8. else return BestMotifs
        else:
            
            return best_motifs
        
            
              
        
    
    





# MAIN

if __name__ == "__main__":
    k, t, dna = greedy.read2D("rosalind_ba2f.txt")
    
    best_m = randomized_motif_search(dna, k, t)
    best_score = score(best_m, k)
    
    i = 0
    
    while i < 1000:
        
        new_motifs = randomized_motif_search(dna, k, t)
        scr = score(new_motifs, k)
        
        if scr < best_score:
            best_m = new_motifs
            best_score = scr 
        
        i += 1
      
    for m in best_m:
        print(m) 






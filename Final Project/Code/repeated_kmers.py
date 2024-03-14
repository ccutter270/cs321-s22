"""
NAMES: Julia Fairbank, Mia Tarantola, Caroline Cutter

CSCI 0321A - Bioinformatics

ASSIGNMENT: Final Project
            Analyzing Frequencies of repeated k-mers 
            

DUE: Monday, May 16th

________________________________________________________

Explanation:
     Running this file will print out the number of repeated kmers for k
     values 1 - 12. The data will print in the console and write to a
     file called "repeated_kmer_data.csv". The DNA sequence being read
     is in the file named "constructed_genome.txt".

"""


# IMPORTS
import csv



# READ FILE
def read_file(filename):
    '''
    A file containing genome sequences 
     
    Args:
        filename: a string of file name to be read
    
    Returns:
        a string sequence of the reconstructed genome 
    '''
    with open(filename, "r") as file:
       
        sequence = file.readline()
        
    
        return sequence




def unique_kmer_composition(sequence, k):
    '''
    Return set of unique k-mers in sequence
    
    Args:
        sequence: dna sequence string to find k-mers from 
        k: integer of desired length of k-mers
    
    Returns:
        A list of all the unique kmers of lenght k in the given sequence 
    '''
    i = 0
    kmers = []                           # initialize empty kmer list 
    
    while i <= len(sequence) - k:        # run through length of sequence
        kmers.append(sequence[i:i+k])    # add kmer at that location to kmer list 
        i += 1
        
    return kmers 






def repeated_kmers(k, sequence):
    '''
    Return number of repeated sequences of length k 
    
    Args:
        sequence: dna sequence string to find k-mers from 
        k: integer of desired length of k-mers
    
    Returns:
        An integer of number of repeated kmers in sequence 
    ''' 

    kmers = unique_kmer_composition(sequence, k)          # Create array of all possible kmers
    
    kmer_set = list(set(kmers))                           # Remove repeats with set function
    
    total_repeats = len(kmers) - len(kmer_set)            # Subtract lenght to get repeated numbers
  
    return total_repeats
    
    
    
    

# MAIN
if __name__ == "__main__":
    
    sequence = read_file("constructed_genome.txt")
    
    
    # WRITE DATA TO CSV FILE 
    with open('repeated_kmer_data.csv', 'w') as f:
    
        writer = csv.writer(f)

        writer.writerow(['K', 'Repeats'])                                   # write header

        
        for i in range(1, 13):
            writer.writerow([i, repeated_kmers(i, sequence)])    # write data
        
    
    
    
    # PRINT DATA IN CONSOLE 
    print("Number of repeated 1-mers: ")
    print(repeated_kmers(1, sequence))
    
    print("Number of repeated 2-mers: ")
    print(repeated_kmers(2, sequence))
    
    print("Number of repeated 3-mers: ")
    print(repeated_kmers(3, sequence))
    
    print("Number of repeated 4-mers: ")
    print(repeated_kmers(4, sequence))
    
    print("Number of repeated 5-mers: ")
    print(repeated_kmers(5, sequence))
    
    print("Number of repeated 6-mers: ")
    print(repeated_kmers(6, sequence))
    
    print("Number of repeated 7-mers: ")
    print(repeated_kmers(7, sequence))
    
    print("Number of repeated 8-mers: ")
    print(repeated_kmers(8, sequence))
    
    print("Number of repeated 9-mers: ")
    print(repeated_kmers(9, sequence))
    
    print("Number of repeated 10-mers: ")
    print(repeated_kmers(10, sequence))
    
    print("Number of repeated 11-mers: ")
    print(repeated_kmers(11, sequence))
    
    print("Number of repeated 12-mers: ")
    print(repeated_kmers(12, sequence))
    

    

    
                   
                   
                   
                   
                   




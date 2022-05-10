#!/usr/bin/env python3

''' usage: kmerCounter.py inFilename
For each line in a specified input file of DNA strings (inFilename), a dataframe of k values and kmer counts is written to a csv output files and the DNA string, output filename, and linguistic complexity are printed to the command line'''

import sys # import sys library in order to take command line input arguments

# For use in later functions: A function to return a list of all possible k values for a specified string
def k_vals(DNAstr):
    '''(str) -> list
    Return a list containing all possible kmer sizes (k values) for DNAstr)'''
    return list(range(1, len(DNAstr)+1))

# Define a function to count possible kmers of size k from a string, where k and the string are specified as arguments:
def possible_kmers(k,DNAstr):
    '''(int, str) -> int
    Return the number of possible substrings (kmers) of length k (k > 0) in DNAstr'''
    return(min(len(DNAstr)-k+1,4**k))

# Define a function to count observed kmers of size k from a string, where k and the string are specified as arguments.
def observed_kmers(k,DNAstr):
    '''(int, str) -> int
    Return the number of observed substrings (unique kmers) of size k (k > 0) in DNAstr'''
    l=[] #initialize list of unique kmers
    # for each observed kmer, if it is not in the unique kmer list, add it to the list
    i=0 #initiate an index
    while i <= len(DNAstr)-k:
        kmer = DNAstr[i:i+k]
        i+=1 #increment the index
        if kmer not in l:
            l.append(kmer)
    return(len(l)) # return the length of the unique kmer list

# Define a function to calculate linguistic complexity. The linguistic complexity of the string is defined as the number of kmers observed for all possible k-mer lengths, divided by the total number that are theoretically possible.
def ling_complexity(DNAstr):
    '''(str) -> float
    Return the linguistic complexity of DNAstr'''
    sum_obs = 0 # initiate sum of all observed kmers
    sum_pos = 0 # initiate sum of all possible kmers
    for k in k_vals(DNAstr): # iterate over all possible kmer lengths for DNAstr
        sum_obs += observed_kmers(k,DNAstr) # calc observed kmers and add to sum_obs
        sum_pos += possible_kmers(k,DNAstr) # calc possible kmers and add to sum_pos
    return sum_obs/sum_pos # return linguistic complexity (sum of obs kmers / sum of possible kmers)

# Define a function to create a pandas data frame containing all possible k and the associated number of observed and expected kmers (see above table).
def kmer_df(DNAstr):
    '''(str) -> pandas.core.frame.DataFrame
    Return a pandas dataframe containing all possible valuse of k
    and the associated number of observed and expected kmers for DNAstr'''
    import pandas as pd # import pandas library
    obs_kmers = [] # initiate an observed kmer list
    pos_kmers = [] # initiate a possible kmer list
    all_k = k_vals(DNAstr) # generate list of all possible kmer lengths (k) for DNAstr
    for k in all_k: # iterate over all possible k
        obs_kmers.append(observed_kmers(k,DNAstr)) # generate number of observed kmers and add to observed kmer list
        pos_kmers.append(possible_kmers(k,DNAstr)) # generate number of possible kmers and add to possible kmer list
    # use pandas to generate three-column dataframe containing the all_k, obs_kmers, and pos_kmers lists:
    kmer_df = pd.DataFrame({
        'k' : all_k,
        'Observed kmers' : obs_kmers,
        'Possible kmers' : pos_kmers
        }
    )
    return(kmer_df) # return kmer dataframe

# The main function will read the specified input file and, for each string in the input file, write a kmer dataframe to a separate output file and print the linguistic complexity to the command line
def main(inFilename):
    ''' (str) -> commandline output, .csv file(s)
    For each line in a specified input file of DNA strings (inFilename)
    main() will write a kmer dataframe to a separate csv output file
    and print the string, output filename, and linguistic complexity to the command line'''
    inFile = open(inFilename, 'r') # open the input file for reading
    lines = inFile.readlines() # read all lines of input file and return each line as a string element in a list (lines)
    i = 1 # initiate an iterator for naming outFiles
    for line in lines:
        line = line.strip(";\n") # remove end line and ; characters from line
        outFilename = "kmer_df_" + str(i) + ".csv" # generate outFile name
        kmer_df(line).to_csv(outFilename, index=None) # generate kmer dataframe and write to outFile as csv
        print("DNA sequence " + str(i) + ": " + line)
        print(" - linguistic complexity: " + str(ling_complexity(line)))
        print(" - kmer dataframe saved to: " + outFilename + "\n")
        i += 1 # increment i for naming next outFile


if __name__ == "__main__":
    ## script execution commands:
    inFilename = sys.argv[1] # save command line argument as inFile
    main(inFilename)
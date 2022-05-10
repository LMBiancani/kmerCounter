# kmerCounter
## ReadMe

For each line in a specified input file of DNA strings `kmerCounter.py` produces a dataframe of all possible k values, as well as the potential and observed kmer counts for each k value. Each dataframe is written to a separate .csv file and the DNA string, output filename, and linguistic complexity are printed to the command line.

### Usage:
* `$ python3 kmerCounter.py inFilename`
* `inFilename` is the name (or path and name) of an input file passed as a command line argument

### Other Files:
* `sampledna.txt` sample input file
* `kmer_df_1.csv`, `kmer_df_2.csv`, `kmer_df_3.csv` output files generated using `sampledna.txt`
* `test_kmerCounter.py` test file for unit testing with `py.test`

### Functions
* `k_vals(DNAstr)` Returns a list containing all possible kmer sizes (k values) for a specified string (DNAstr)
* `possible_kmers(k,DNAstr)` Returns the number of possible substrings (kmers) of length k (k > 0) in a specified string (DNAstr)
* `observed_kmers(k,DNAstr)` Returns the number of observed substrings (unique kmers) of length k (k > 0) in a specified string (DNAstr)
* `ling_complexity(DNAstr)` Returns the linguistic complexity of a specified string (DNAstr)
* `kmer_df(DNAstr)` Returns a pandas dataframe containing all possible values of k and the associated number of observed and expected kmers a specified string (DNAstr)
* `main(inFilename)` For each line in a specified input file of DNA strings (inFilename) this function will write a kmer dataframe to a separate csv output file and print the string, output filename, and linguistic complexity to the command line

### Background
A string is a sequence of characters that we store as an object. This string can be divided into substrings, where the possible number of substrings of length k (termed kmers) is the number of possible characters to the k. The linguistic complexity of the string is defined as the number of kmers that are observed for all possible kmer lengths, divided by the total number that are theoretically possible.

As an example, consider the string ATTTGGATT. From the following table you can see that the linguistic complexity is 35 / 40 = 0.875. Note that the possible number of kmers (usually 4^k) is limited by the length of the sequence. Thus Possible Kmers is calculated as the minimum of (1) the length of the string minus k plus 1, and (2) 4^k (i.e. the number of possible k-mers of length 9 in the sequence is 1, not 4^9).

| k |Observed kmers|Possible kmers|
|:-:|:------------:|:------------:|
| 1 |      3       |      4       |
| 2 |      5       |      8       |
| 3 |      6       |      7       |
| 4 |      6       |      6       |
| 5 |      5       |      5       |
| 6 |      4       |      4       |
| 7 |      3       |      3       |
| 8 |      2       |      2       |
| 9 |      1       |      1       |


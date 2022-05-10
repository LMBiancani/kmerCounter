#!/usr/bin/env python3

from kmerCounter import *

def test_example_str1():
    DNAstr = "ATTTGGATT"
    assert k_vals(DNAstr) == [1, 2, 3, 4, 5, 6, 7, 8, 9]
    assert ling_complexity(DNAstr) == 0.875
    k = 1
    assert possible_kmers(k,DNAstr) == 4 # only 4 unique kmers possible with 4-letter alphabet
    assert observed_kmers(k,DNAstr) == 3
    k = 2
    assert possible_kmers(k,DNAstr) == 8
    assert observed_kmers(k,DNAstr) == 5
    k = 9 # full length of input string
    assert possible_kmers(k,DNAstr) == 1
    assert observed_kmers(k,DNAstr) == 1
    df = kmer_df(DNAstr)
    assert list(df['k']) == [1, 2, 3, 4, 5, 6, 7, 8, 9]
    assert list(df['Observed kmers']) == [3, 5, 6, 6, 5, 4, 3, 2, 1]
    assert list(df['Possible kmers']) == [4, 8, 7, 6, 5, 4, 3, 2, 1]
    
def test_example_str2():
    DNAstr = "ACTGCAGCGCGATGATGAGAGAGATTTCAGGACACACATTGCCAAATTGAGGCAT"
    assert k_vals(DNAstr) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55]
    assert ling_complexity(DNAstr) == 0.9738111647139903
    k = 1
    assert possible_kmers(k,DNAstr) == 4 # only 4 unique kmers possible with 4-letter alphabet
    assert observed_kmers(k,DNAstr) == 4
    k = 5
    assert possible_kmers(k,DNAstr) == 51
    assert observed_kmers(k,DNAstr) == 47
    k = 55 # full length of input string
    assert possible_kmers(k,DNAstr) == 1
    assert observed_kmers(k,DNAstr) == 1
    df = kmer_df(DNAstr)
    assert list(df['k']) == k_vals(DNAstr)
    assert list(df['Observed kmers']) == [4, 14, 31, 43, 47, 49, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
    assert list(df['Possible kmers']) == [4, 16, 53, 52, 51, 50, 49, 48, 47, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33, 32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18, 17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
    

def test_example_str3():
    DNAstr = "ATATATATATATATATA"
    assert k_vals(DNAstr) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
    assert ling_complexity(DNAstr) == 0.2357142857142857
    k = 1
    assert possible_kmers(k,DNAstr) == 4 # only 4 unique kmers possible with 4-letter alphabet
    assert observed_kmers(k,DNAstr) == 2
    k = 2
    assert possible_kmers(k,DNAstr) == 16
    assert observed_kmers(k,DNAstr) == 2
    k = 17 # full length of input string
    assert possible_kmers(k,DNAstr) == 1
    assert observed_kmers(k,DNAstr) == 1
    df = kmer_df(DNAstr)
    assert list(df['k']) == [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17]
    assert list(df['Observed kmers']) == [2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1]
    assert list(df['Possible kmers']) == [4, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1]
    
def test_mono_str():
    DNAstr = "AAA"
    assert k_vals(DNAstr) == [1, 2, 3]
    assert ling_complexity(DNAstr) == 0.5
    k = 1
    assert possible_kmers(k,DNAstr) == 3
    assert observed_kmers(k,DNAstr) == 1
    k = 2
    assert possible_kmers(k,DNAstr) == 2
    assert observed_kmers(k,DNAstr) == 1
    k = 3 # full length of input string
    assert possible_kmers(k,DNAstr) == 1
    assert observed_kmers(k,DNAstr) == 1
    df = kmer_df(DNAstr)
    assert list(df['k']) == [1, 2, 3]
    assert list(df['Observed kmers']) == [1, 1, 1]
    assert list(df['Possible kmers']) == [3, 2, 1]
    

def test_single_char():
    DNAstr = "G"
    assert k_vals(DNAstr) == [1]
    assert ling_complexity(DNAstr) == 1
    k = 1 # full length of input string
    assert possible_kmers(k,DNAstr) == 1
    assert observed_kmers(k,DNAstr) == 1
    df = kmer_df(DNAstr)
    assert list(df['k']) == [1]
    assert list(df['Observed kmers']) == [1]
    assert list(df['Possible kmers']) == [1]
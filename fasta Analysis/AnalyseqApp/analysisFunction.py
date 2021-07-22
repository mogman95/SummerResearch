"""Import Necessary Libraries"""
from Bio import SeqIO
import time
# from fp_vari import file_path

def main(file_path):
    """Define Functions"""
    '''Complement Function:'''
    ### Takes a DNA sequences as a string and returns its complementary strand sequence
    def comp(seq):
        ### Initialize base pair dictionary and 'comp' string
        bp = {'A': 'T',
               'T': 'A',
               'C': 'G',
               'G': 'C'}
        comp = ''
        ### Sequence creation
        for b in seq:     # For each base in the sequence...
            comp += bp[b] # use the dictionary to add the next complementary base in the sequence...
        return comp       # output the new complementary sequence

    '''Reverse Function'''
    ### Takes any string and returns its reverse string
    def rev(seq):
        ### Initialize 'rev' string
        rev = ''
        ### Reversing mechanism
        for i in range(len(seq)): # For each index in the length of the given string...
            rev += seq[-i-1]      # index the given string from the last charater to the first,
                                  # adding each character to the new string
        return rev                # output the now reversed string

    '''Counting Function'''
    ### Scrolls across a string one charater at a time counting
    ### how many times a certain target phrase shows up
    def nuc_count(nuc, seq):
        ### Initialize nuc_locs list to later be counted, my_bool useful boolian variable, and start variable for indexing
        nuc_locs = []
        my_bool = True
        start = 0
        ### Counting sequence
        while my_bool:                         # While my_bool is true
            a = seq.find(nuc, start, len(seq)) # Find the first place the target
            if a != -1:                        # If search doesn't not find a value (finds a value)...
                nuc_locs.append(a)             # Append location of match to list
                start = a+1                    # Iterate starting location for search
            else:                              # If search does find a value...
                my_bool = False                # Break from the loop
        return len(nuc_locs)                   # Return the length of the string that holds each occurance
                                               # of the target phrase (i.e. number of times the phrase appears)

    '''Nucleotide Frequency Function'''
    ### Counts the frequencies of all 4 nucleotides in a sequence
    def nuc_freq(seq):
        ### Initialize the length of the sequence, l, and the dictionary that will hold values
        l = len(seq)
        nuc_f_dict = {}

        ### Perform count
        for char in "ATCG":                          # For A, T, C, & G...
            nuc_f_dict[char] = nuc_count(char,seq)/l # Count how many times that base shows up in the sequence,
                                                     # divide by the length of the sequence,
                                                     # & store the value with its corresponding base in the dictionary
        return nuc_f_dict                            # Return the dictionary with the nuc freq values

    '''Nucleotide Frequency Function'''
    ### Counts the frequencies of all 16 dinucleotides in a sequence
    def dinuc_freq(seq):
        ### Initialize the length of the sequence, l, and the dictionary that will hold values
        l = len(seq)
        dinuc_f_dict = {}
        ### Perform count
        for char1 in "ATCG":                                      # >
            for char2 in "ATCG":                                  # >
                dinuc = char1 + char2                             # > For A, T, C, & G add A, T, C, & G to make each dinucleotide (AA, AT, AC, AG, TA, TT, etc.)
                dinuc_f_dict[dinuc] = nuc_count(dinuc, seq)/(l-1) # Count how many times each dinuc shows up in the sequence,
                                                                  # divide by 1 minus the length of the sequence (e.g. 'AAA' -> length: 3| total dinuc: 2 = length-1),
                                                                  # & store the value with its corresponding dinuc in the dictionary
        return dinuc_f_dict                                       # Return the dictionary with the dinuc freq values

    """Define file inputs"""
    ### Commenting out allows rapid analysis of different files
    # file = "SARS-CoV-2.fasta"
    # file = "baccoa.fna"

    """Running Main program"""
    # tic = time.perf_counter() # Start timer
    """Parse File"""
    ### Initialize string to create the sequence in case there are multiple chromosomes/scaffolds/contigs in the FASTA file
    ref_seq = ""
    ### File parsing
    records = SeqIO.parse(file_path, "fasta")    # Create iterator for FASTA file
    for record in records:                  # For each chromosome/scaffold/contig in the FASTA file...
        seq = record.seq                    # Initializes a sequence
        seq = seq.__str__().replace("N","") # Converts sequence into a string & removes any N's from genome
        seq = seq.upper()                   # Makes sure sequence is uppercase so each charater is counted correctly
        ref_seq += seq                      # Stitches seperate sequences together into one

    # GC_cont = (nuc_count("G",ref_seq) + nuc_count("C",ref_seq))/len(ref_seq) # Add the number of Gs and Cs in the sequence and divide by the length of the sequence
    # nuc_dict = nuc_freq(ref_seq)                                             # Get nucleotide freqs
    # dinuc_dict = dinuc_freq(ref_seq)                                         # Get dinucleotide freqs
    # comp_seq = comp(ref_seq)                                                 # Get the complement of the sequence (3' to 5' if the sequence is 5' to 3')
    # rev_comp = rev(comp(ref_seq))    # Reverse the complementary sequence (to 5' to 3' from 3' to 5' or vice versa)
    def sequence():
        return ref_seq

    def comp_seq():
        return comp(ref_seq)

    def rev_comp():
        return rev(comp(ref_seq))

    def GC_cont():
        GC_cont = str((nuc_count("G", ref_seq) + nuc_count("C", ref_seq)) / len(ref_seq)*100)
        return GC_cont

    def length():
        return str(len(ref_seq))

    def nuc_results():
        nuc_dict = nuc_freq(ref_seq)
        results = ""
        for nuc in nuc_dict:  # Print each nucleotide and its corresponding frequency value
            results += nuc + ": " + str(nuc_dict[nuc]) + "\n"
        return results

    def dinuc_results():
        dinuc_dict = dinuc_freq(ref_seq)
        results = ""
        for nuc in dinuc_dict:  # Print each nucleotide and its corresponding frequency value
            results += nuc + ": " + str(dinuc_dict[nuc]) + "\n"
        return results

    # toc = time.perf_counter() # Stop timer

    return sequence(), comp_seq(), rev_comp(), GC_cont(), length(), nuc_results(), dinuc_results()

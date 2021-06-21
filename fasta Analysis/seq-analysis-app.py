"""Import Necessary Libraries"""
from Bio import SeqIO
from Bio.Seq import Seq
import time
import sys
from PyQt5.QtCore import QSize, Qt
from PyQt5.QtWidgets import (
    QApplication,
    QLabel,
    QMainWindow,
    QPushButton,
    QVBoxLayout,
    QWidget,
    QAction,
    QCheckBox,
    QMainWindow,
    QStatusBar,
    QToolBar,
    QTextEdit,
)
from PyQt5.QtGui import QIcon

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
file = "SARS-CoV-2.fasta"
# file = "baccoa.fna"


"""Running Main program"""
tic = time.perf_counter() # Start timer
"""Parse File"""
### Initialize string to create the sequence in case there are multiple chromosomes/scaffolds/contigs in the FASTA file
ref_seq = ""
### File parsing
records = SeqIO.parse(file, "fasta")    # Create iterator for FASTA file
for record in records:                  # For each chromosome/scaffold/contig in the FASTA file...
    seq = record.seq                    # Initializes a sequence
    seq = seq.__str__().replace("N","") # Converts sequence into a string & removes any N's from genome
    seq = seq.upper()                   # Makes sure sequence is uppercase so each charater is counted correctly
    ref_seq += seq                      # Stitches seperate sequences together into one

# GC_cont = (nuc_count("G",ref_seq) + nuc_count("C",ref_seq))/len(ref_seq) # Add the number of Gs and Cs in the sequence and divide by the length of the sequence
# nuc_dict = nuc_freq(ref_seq)                                             # Get nucleotide freqs
# dinuc_dict = dinuc_freq(ref_seq)                                         # Get dinucleotide freqs
comp_seq = comp(ref_seq)                                                 # Get the complement of the sequence (3' to 5' if the sequence is 5' to 3')
rev_comp = rev(comp(ref_seq))                                            # Reverse the complementary sequence (to 5' to 3' from 3' to 5' or vice versa)
def GC_cont():
    GC_cont = str((nuc_count("G", ref_seq) + nuc_count("C", ref_seq)) / len(ref_seq)*100)
    return GC_cont + "\n"

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

toc = time.perf_counter() # Stop timer
""" Output results to a TXT file """
# ### Initialize output file name
# output_file = "analysis.txt"
# with open(output_file,'w+') as f: # Open a new file, or override and old file, that can be written to as 'f'
#     f.write("GC Content:\n" + str(GC_cont*100) + " %\n")
#     f.write("\nNucleotide Frequencies:\n")
#     for nuc in nuc_dict: # Print each nucleotide and its corresponding frequency value
#         f.write(nuc + ": " + str(nuc_dict[nuc]) + "\n")
#     f.write("\nDinucleotide Frequencies:\n")
#     for dinuc in dinuc_dict: # Print each dinucleotide and its corresponding frequency value
#         f.write(dinuc + ": " + str(dinuc_dict[dinuc]) + "\n")
#     f.write("\nSequence:\n" + ref_seq)
#     f.write("\n\nComplement:\n" + comp_seq)
#     f.write("\n\nReverse Complement:\n" + rev_comp)
#     f.write("\nCode executed in %0.4f seconds" % (toc-tic))
# f.close()


class AnotherWindow(QWidget):
    """
    This "window" is a QWidget. If it has no parent,
    it will appear as a free-floating window.
    """

    def __init__(self):
        super().__init__()
        layout = QVBoxLayout()
        self.setLayout(layout)
        self.button = QPushButton('close')
        layout.addWidget(self.button)
        self.setFixedSize(QSize(500,300))

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("Sequence Analysis")
        self.resize(1000,1000)
        self.window1 = AnotherWindow()
        self.window1.button.clicked.connect(self.toggle_window1)

        menu = self.menuBar()

        save_comp_analysis = QAction("Save Complete Analysis", self)
        save_comp_analysis.triggered.connect(self.save_comp_analysis)
        save = QAction("Save", self)
        save.triggered.connect(self.save_comp_analysis)
        save_as = QAction("Save as...", self)
        save_as.triggered.connect(self.save_comp_analysis)
        # self.setStatusBar(QStatusBar(self))
        file_menu = menu.addMenu("&File")
        file_menu.addAction(save_comp_analysis)
        file_menu.addAction(save)
        file_menu.addAction(save_as)

        newana_action = QAction('New Analysis', self)
        newana_action.triggered.connect(self.save_comp_analysis)
        editana_action = QAction('Edit Analysis', self)
        editana_action.triggered.connect(self.save_comp_analysis)
        analysis_menu = menu.addMenu('&Analysis')
        analysis_menu.addAction(newana_action)
        analysis_menu.addAction(editana_action)

        tutorial_action = QAction('Tutorial', self)
        tutorial_action.triggered.connect(self.save_comp_analysis)
        help_menu = menu.addMenu('&Help')
        help_menu.addAction(tutorial_action)

        l = QVBoxLayout()
        self.text = QTextEdit()
        self.text.setText("GC Content:\n" + GC_cont()
                        + "\nNucleotide Frequencies:\n" + nuc_results()
                        + "\nNucleotide Frequencies:\n" + dinuc_results()
                        + "\nSequence:\n" + ref_seq
                        + "\n\nComplement:\n" + comp_seq
                        + "\n\nReverse Complement:\n" + rev_comp
                        + "\nCode executed in %0.4f seconds" % (toc-tic)
                                               )
        self.text.setDisabled(False)
        l.addWidget(self.text)

        w = QWidget()
        w.setLayout(l)
        self.setCentralWidget(w)

    def toggle_window1(self, checked):
        if self.window1.isVisible():
            self.window1.hide()

        else:
            self.window1.show()

    def save_comp_analysis(self):
        output_file = "analysis.txt"
        with open(output_file,'w+') as f: # Open a new file, or override and old file, that can be written to as 'f'
            f.write("GC Content:\n" + GC_cont())
            f.write("\nNucleotide Frequencies:\n" + nuc_results())
            f.write("\nDinucleotide Frequencies:\n" + dinuc_results())
            f.write("\nSequence:\n" + ref_seq)
            f.write("\n\nComplement:\n" + comp_seq)
            f.write("\n\nReverse Complement:\n" + rev_comp)
            f.write("\nCode executed in %0.4f seconds" % (toc-tic))
        f.close()

app = QApplication(sys.argv)

window = MainWindow()
window.show()

app.exec()
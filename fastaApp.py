from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import GC
import tkinter as tk
from tkinter import filedialog

root = tk.Tk()
root.withdraw()
file_path = filedialog.askopenfilename()

record_iterator = SeqIO.parse(file_path, "fasta")
first_record = next(record_iterator)
print(first_record,"\n")

for seq_rec in SeqIO.parse(file_path, "fasta"):
    seq = seq_rec.seq
    l = len(seq)
    if l < 50000:
        print("\nSequence:",seq)
        print("\nComplement:",seq.complement())
        print("\nReverse Complement:",seq.reverse_complement())
    else:
        print("\nFull sequence too large to display")
    print(f"\nLength: {l} bp")
    print(f"\nGC Content: {GC(seq)}%")

A = seq.count("A")
T = seq.count("T")
C = seq.count("C")
G = seq.count("G")
Ap = A/l
Tp = T/l
Cp = C/l
Gp = G/l
print(f"\nNucleotide Frequency:\nA: {Ap}\nT: {Tp}\nC: {Cp}\nG: {Gp}")

AA = 0
AT = 0
AC = 0
AG = 0
TA = 0
TT = 0
TC = 0
TG = 0
CA = 0
CT = 0
CC = 0
CG = 0
GA = 0
GT = 0
GC = 0
GG = 0
bp_num = 0

for bp in seq:
    if seq[bp_num]+seq[bp_num+1] == "AA":
        AA += 1
    elif seq[bp_num]+seq[bp_num+1] == "AT":
        AT += 1
    elif seq[bp_num]+seq[bp_num+1] == "AC":
        AC+= 1
    elif seq[bp_num]+seq[bp_num+1] == "AG":
        AG += 1
    elif seq[bp_num]+seq[bp_num+1] == "TA":
        TA += 1
    elif seq[bp_num]+seq[bp_num+1] == "TT":
        TT += 1
    elif seq[bp_num]+seq[bp_num+1] == "TC":
        TC+= 1
    elif seq[bp_num]+seq[bp_num+1] == "TG":
        TG += 1
    elif seq[bp_num]+seq[bp_num + 1] == "CA":
        CA += 1
    elif seq[bp_num]+seq[bp_num + 1] == "CT":
        CT += 1
    elif seq[bp_num]+seq[bp_num + 1] == "CC":
        CC += 1
    elif seq[bp_num]+seq[bp_num + 1] == "CG":
        CG += 1
    elif seq[bp_num]+seq[bp_num+1] == "GA":
        GA += 1
    elif seq[bp_num]+seq[bp_num+1] == "GT":
        GT += 1
    elif seq[bp_num]+seq[bp_num+1] == "GC":
        GC+= 1
    elif seq[bp_num]+seq[bp_num+1] == "GG":
        GG += 1
    bp_num += 1
    if bp_num == l - 1:
        break

dinuc_count = AA+AT+AC+AG+TA+TT+TC+TG+CA+CT+CC+CG+GA+GT+GC+GG
AAp = AA/dinuc_count
ATp = AT/dinuc_count
ACp = AC/dinuc_count
AGp = AG/dinuc_count
TAp = TA/dinuc_count
TTp = TT/dinuc_count
TCp = TC/dinuc_count
TGp = TG/dinuc_count
CAp = CA/dinuc_count
CTp = CT/dinuc_count
CCp = CC/dinuc_count
CGp = CG/dinuc_count
GAp = GA/dinuc_count
GTp = GT/dinuc_count
GCp = GC/dinuc_count
GGp = GG/dinuc_count
print(f"\nDinucleotide Frequency:\nAA: {AAp}\nAT: {ATp}\nAC: {ACp}\nAG: {AGp}\nTA: {TAp}\nTT: {TTp}\nTC: {TCp}\nTG: {TGp}\nCA: {CAp}\nCT: {CTp}\nCC: {CCp}\nCG: {CGp}\nGA: {GAp}\nGT: {GTp}\nGC: {GCp}\nGG: {GGp}")
input("\n\nPress 'enter' to exit\n->")
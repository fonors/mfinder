#!/usr/bin/env python3
'''
    Motif Finder - A FASTA parser for motifs
    Copyright (C) 2023  fonors & Wil-s0n

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''

from sys import argv, stderr, exit
from re import match

# Opens a file, returns an error and prints it in case it doesn't find a valid file
try:
    fasta_file = open(argv[1], "r")
except IndexError as err:
    exit("ERROR: No file specified. Please specify a valid FASTA file.")
else:
    try:
        if argv[1].endswith(".fasta") == False:
            raise RuntimeError("File provided is not a valid FASTA file.")
    except RuntimeError as err:
        exit(err)

# Defines the veriables necessary for finding the motifs in the sequences and stores them
motif = input("Type the motif you wish to find: ")
motif = motif.upper()
seq_dict = {}
revcomp_seq_dict = {}

nucleotide_check = match('^[ATGC]+$', motif)
if len(motif) >= 5 and nucleotide_check:
    for line in fasta_file:
        line = line.strip()
        if line.startswith(">"):
            seq_name = line
            seq_name = seq_name[1:len(seq_name)]
            seq_dict[seq_name] = ""
        else:
            seq_dict[seq_name] += line
elif len(motif) < 5 and nucleotide_check:
    print("The motif you have typed is too short!", file=stderr)
    exit()
elif len(motif) >= 5 and not nucleotide_check:
    print("The motif you have typed doesn't have the permitted nucleotides!", file=stderr)
    exit()
else:
    print("The motif you have typed is too short and doesn't have the permitted nucleotides!", file=stderr)
    exit()
fasta_file.close()

# Creates the reverse-complement motif for additional analysis
revcomp_motif_temp = motif.replace("A", "%temp%").replace("T", "A").replace("%temp%", "T")
revcomp_motif = revcomp_motif_temp.replace("C", "%temp%").replace("G", "C").replace("%temp%", "G")
revcomp_motif = revcomp_motif[::-1]

# Checks for multiple motifs in the same sequence, and its reverse-complement, and prints out en error if it finds multiple ones in the same sequence. It then checks the position and the frame of the motifs in the sequences it's found in
for seq in seq_dict:
    seq_len = len(seq_dict[seq])
    if seq_dict[seq].count(motif) > 0 or seq_dict[seq].count(revcomp_motif) > 0:
        motif_count = seq_dict[seq].count(motif) + seq_dict[seq].count(revcomp_motif)
        for each_motif in range(motif_count):
            if motif_count > 1:
                if motif in seq_dict[seq] and revcomp_motif in seq_dict[seq]:
                    temp = 0
                    motif_index = temp + seq_dict[seq].index(motif)
                    if motif_index % 3 == 0:
                        frame = 3
                    elif motif_index % 3 == 1:
                        frame = 2
                    elif motif_index % 3 == 2:
                        frame = 1
                    del temp
                print(f"{seq}\t{seq_len}\t{motif_index}\t{frame}")
            elif motif_count == 1:
                if motif in seq_dict[seq]:
                    motif_index = seq_dict[seq].index(motif) + 1
                    if motif_index % 3 == 0:
                        frame = 3
                    elif motif_index % 3 == 1:
                        frame = 2
                    elif motif_index % 3 == 2:
                        frame = 1
                elif revcomp_motif in seq_dict[seq]:
                    motif_index = seq_dict[seq].index(revcomp_motif) + 1
                    if motif_index % 3 == 0:
                        frame = -1
                    elif motif_index % 3 == 1:
                        frame = -2
                    elif motif_index % 3 == 2:
                        frame = -3
                print(f"{seq}\t{seq_len}\t{motif_index}\t{frame}")
        if motif_count > 1:
            print(f"Multiple matches found for the inputed motif in {seq}", file=stderr)

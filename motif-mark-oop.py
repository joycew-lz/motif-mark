#!/usr/bin/env python

#--------------
# import modules
#--------------

import argparse
import re
import cairo
import math
import IPython.display

#--------------
# argparse
#--------------

def get_args():
    parser = argparse.ArgumentParser(description="to __")
    parser.add_argument("-f", help="designates file path to the fasta file", type=str, required=True)
    parser.add_argument("-m", help="designates file path to the motifs text file", type=str, required=True)
    # parser.add_argument("-r", help="designates file path to __", type=str, required=False)
    return parser.parse_args()

args = get_args()

# Set global argparse variables:
fasta_file_path = args.f
motif_file_path = args.m
# report = args.r

#--------------
# classes with instance attributes and methods
#--------------

class FastaSequence:
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence 
        self.length = len(sequence)
        self.exons, self.introns = self._find_regions()  # list of (start, end) tuples for exons and introns in the sequence

    def _find_regions(self):
        '''
        Takes a FastaSequence object and returns a list of tuples for the start and end positions of
        the exons and of the introns in the sequence.
        Exons are always "sandwiched" between introns, so if the first region is an exon, it starts at position 0 and ends at the start of the first intron. 
        If the first region is an intron, it starts at position 0 and ends at the start of the first exon. 
        '''
        exons = []
        introns = []

        if self.length == 0:
            return exons, introns
        
        # look at the first character of the sequence to determine if the first region is an exon or intron:
        if self.sequence[0].isupper():
            region_type = "exon"
        else:
            region_type = "intron"
        

        return exons, introns
    
def read_fasta(fasta_file: str):
    '''
    Reads a FASTA file and return a list of FastaSequence objects, where each object contains the header beginning with ">", 
    the sequence in the following lines, the sequence length, and exons and introns ranges for a single gene sequence.
    Handles multiple-line sequences and preserves capitalization to allow for exon and intron detection.
    '''
    sequences = []
    header = None
    sequence_lines = ""

    with open(fasta_file, 'r') as fh:
        for line in fh:
            line = line.strip('\n')
        
            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    sequences.append(FastaSequence(header, sequence_lines))
                header = line[1:]
                sequence_lines = "" 
            else:
                sequence_lines += line
        
        # final sequence after the loop ends
        if header is not None:
            sequences.append(FastaSequence(header, sequence_lines))

    return sequences

fasta_sequences = read_fasta(fasta_file_path) # List of FastaSequence objects for all sequences in the fasta file


class Motif:
    def __init__(self, name, raw_sequence, color):
        self.name = name  # Instance attribute: store motif name
        self.raw_sequence = raw_sequence  # Instance attribute: store motif raw sequence
        self.length = len(raw_sequence)  # Instance attribute: store motif length
        self.regex_pattern = self._convert_to_regex()  # Instance attribute: store motif regex pattern
        self.color = color # Instance attribute: store visualization color

    # Methods:
    def _convert_to_regex(self):
        '''
        Converts a motif raw sequence to a regex pattern using an IUPAC code dictinary.
        Converts RNA to DNA if necessary by converting U's to T's. 
        For example, if the raw sequence is "AUGC", it will be converted to "ATGC".
        '''
        # Convert motif raw sequence to regex pattern
        # Example: "ATCG" -> "ATCG"
        return self.raw_sequence

class MotifLocation:
    def __init__(self, motif, start, end, sequence_header):
        self.motif = motif  # Instance attribute: Motif object
        self.sequence_header = sequence_header  # Instance attribute: header string of the sequence where the motif is found
        self.start = start  # Instance attribute
        self.end = end  # Instance attribute

class MotifFinder:
    def __init__(self, motif_file):
        self.motif_file = motif_file  # Instance attribute: store motif file path

    # Methods:
    # def find_motifs_in_sequence(self, fasta_sequence) -> returns list of MotifLocation objects for the given motif in the given fasta sequence, using the motif's regex pattern Motif.regex_pattern.
        # deal with overlapping mathes
    
class FigureBuilder:
    def __init__(self, width=1000, height=900, gene_sequences, motif_locations):
        self.width = width  # Instance attribute: store figure width
        self.height = height  # Instance attribute: store figure height
        self.surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, self.width, self.height)  # Instance attribute: create cairo surface
        self.context = cairo.Context(self.surface)  # Instance attribute: create cairo context

        self.gene_sequences = gene_sequences   # Instance attribute: list of FastaSequence objects
        self.motifs = motif_locations # Instance attribute: dict mapping sequence header to list of MotifLocation objects for that sequence

    # Methods:
    # def draw_figure(output_png_path)


# List of Classes for OOP Motif Mark

# 1. **FastaSequence**: Represents a single gene sequence, including exons, introns, and motifs.
# 2. **Motif**: Represents a single motif location on a gene as well as the motif sequence, length, and color in the figure.
# 3. **MotifLocation**: Represents the location of a motif on a gene sequence, including the motif, start and end positions, and the header of the sequence it is found on.
# 4. **MotifFinder**: Contains methods to find motifs in the gene sequences.
# 5. **FigureBuilder**: Uses pycairo to generate the figure.
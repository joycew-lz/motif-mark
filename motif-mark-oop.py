#!/usr/bin/env python

#--------------
# import modules
#--------------

import argparse
import os
import re
import cairo
import math

#--------------
# argparse
#--------------

def get_args():
    parser = argparse.ArgumentParser(description="to visualize sequence motifs using pycairo")
    parser.add_argument("-f", help="designates file path to the fasta file containing sequences", type=str, required=True)
    parser.add_argument("-m", help="designates file path to the motifs text file, with one motif per line", type=str, required=True)
    return parser.parse_args()

#--------------
# classes with instance attributes and methods
#--------------

# Class to represent a FASTA record.
class FastaSequence:
    '''
    Represents one FASTA record's sequence, including the header, raw sequence, length, and exon/intron regions.
    '''
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence 
        self.length = len(sequence)
        self.exons, self.introns = self._find_regions()  # list of (start, end) tuples for exons and introns in the sequence

    # Internal method to find exon and intron regions based on capitalization of the sequence
    def _find_regions(self):
        '''
        Returns two lists of (start, end) tuples for the start and end positions of
        the exons and of the introns in the sequence.
        Exons are always "sandwiched" between introns, so if the first region is an exon, it starts at position 0 and ends at the start of the first intron. 
        If the first region is an intron, it starts at position 0 and ends at the start of the first exon. 
        '''
        exons = []
        introns = []

        if self.length == 0:
            return exons, introns
        
        # Determine the first region type based on the capitalization of the first base in the sequence
        if self.sequence[0].isupper():
            current_type = "exon"
        else:
            current_type = "intron"

        # Start at position 0 for the first base
        start = 0

        # Loop through the base positions starting at position 1,
        # and check for changes in capitalization to determine exon and intron boundaries
        for i in range(1, len(self.sequence)):

            # Determine the i'th base position's region type
            if self.sequence[i].isupper():
                new_type = "exon"
            else:
                new_type = "intron"
            
            # If we have reached the end of the current region type and we're changing from exon --> intron, or intron --> exon:
            if new_type != current_type:
                if current_type == "exon":
                    exons.append((start, i))  # add the exon region's position to the list of exons
                if current_type == "intron":
                    introns.append((start, i))  # add the intron region's position to the list of introns
                
                start = i # update the start position to the current position for the new region type
                current_type = new_type # update the current region type to the new region type

        # After the loop ends, add the final region to the appropriate region type list      
        if current_type == "exon":
            exons.append((start, len(self.sequence)))
        else:
            introns.append((start, len(self.sequence)))

        return exons, introns

# Function to read in a FASTA file.   
def read_fasta(fasta_file: str):
    '''
    Reads a FASTA file and return a list of FastaSequence objects, where each object contains the header beginning with ">", 
    the sequence in the following lines, the sequence length, and exons and introns ranges for a single gene sequence.
    Handles multiple-line sequences and preserves capitalization to allow for exon and intron detection.
    '''
    records = []
    header = None
    sequence_lines = []

    with open(fasta_file, 'r') as fh:
        for line in fh:
            line = line.strip('\n')
        
            if not line:
                continue

            # If this is a header line:
            if line.startswith(">"):

                # If there already is a header, finalize the previous record and add it to the list of records
                if header is not None:
                    records.append(FastaSequence(header, ''.join(sequence_lines)))
                
                # Update header to this new header line and reset sequence lines for the new record
                header = line[1:]
                sequence_lines = []
            
            # If this is not a header line, store it as a sequence line
            else:
                sequence_lines.append(line)
        
        # After the loop ends, append the final record
        if header is not None:
            records.append(FastaSequence(header, ''.join(sequence_lines)))

    return records

# Map IUPAC codes to regex character classes.
IUPAC_REGEX = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "U": "T", # U in RNA should be converted to a T
    "N": "[ACGT]", # Any base
    "Y": "[CT]", # Pyrimidine
    "R": "[AG]", # Purine
    "K": "[GT]", # Ketone
    "M": "[AC]", # Amino
    "S": "[CG]", # Strong interaction
    "W": "[AT]", # Weak interaction
    "B": "[CGT]", # Not A
    "D": "[AGT]", # Not C
    "H": "[ACT]", # Not G
    "V": "[ACG]" # Not T
}

# Class to represent a motif.
class Motif:
    '''
    Represents __.
    '''
    def __init__(self, raw_sequence, color):
        self.raw_sequence = raw_sequence.strip() # Store motif raw sequence
        self.color = color # Store motif color
        self.length = len(self.raw_sequence)  # Store motif length
        self.regex_pattern = self._convert_to_regex()  # Store motif regex pattern as a string

    # Internal method to convert the motif raw sequence to a regex pattern using the IUPAC code dictionary.
    def _convert_to_regex(self):
        '''
        Converts a motif raw sequence to a regex pattern using an IUPAC code dictionary:
        For example, convert a motif like "YGCY" into a regex like "[CT]GC[CT]".
        Also, U's in RNA gets converted to a T.
        '''
        # Convert raw sequence to all uppercase.
        sequence = self.raw_sequence.upper()

        pattern = ""

        for base in sequence:
            if base in IUPAC_REGEX:
                pattern += IUPAC_REGEX[base]
            else:
                pattern += base

        return pattern

# Function to read in a motif file.
def read_motifs(motif_file: str):
    '''
    Reads a motif file and returns a list of Motif objects, where each object contains the motif raw sequence, color, and regex pattern.
    The motif file should have one motif per line.
    '''
    motifs = []
    # Create a color palette with RGB values ranging from 0.00 to 1.00
    color_palette = [
        (0.15, 0.40, 1.00),
        (0.60, 0.80, 1.00),
        (0.70, 0.40, 0.25),
        (0.85, 0.60, 0.45),
        (0.95, 0.90, 0.70),
    ]
    color_index = 0

    with open(motif_file, 'r') as fh:
        for line in fh:
            line = line.strip('\n')
        
            if not line:
                continue
            
            # Cycle through the color palette and pick one color for each motif
            color = color_palette[color_index % len(color_palette)]
            # Increment the color index for the next motif
            color_index += 1

            motifs.append(Motif(line, color))
    
    return motifs

# Function to detect motif matches in a FASTA record's sequence.
def find_motif_locations(record, motifs):
    '''
    Detects motif matches in a FASTA record's sequence, and returns a list of MotifLocation objects, where each object contains the motif match's motif, start position, and end position.
    Also uses lookahead regex to find overlapping matches.
    '''
    # Convert raw sequence to all uppercase
    sequence = record.sequence.upper()

    # Store hit
    motif_locations = []

    for motif in motifs:
        # Use lookahead to find overlapping matches
        pattern = f"(?=({motif.regex_pattern}))"

        for match in re.finditer(pattern, sequence):
            # Start is the position of the first base in the motif match
            start = match.start()
            # End is start + motif length
            end = start + motif.length
            # MotifLocation object stores the motif match's motif, start position, and end position
            motif_locations.append(MotifLocation(motif, start, end))

    return motif_locations


# Class to represent a motif match location on a gene sequence
class MotifLocation:
    '''
    Represents one motif match on a gene, with a start position and end position.
    '''
    def __init__(self, motif, start, end):
        self.motif = motif # Store motif object for this match
        self.start = start # Store match start (0-based)
        self.end = end # Store match end (end-exclusive)
    



class FigureBuilder:
    def __init__(self, gene_sequences, motif_locations, width=1000, height=900):
        # Store figure width and height
        self.width = width 
        self.height = height

        # Create cairo surface and context
        self.surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, self.width, self.height)  # Instance attribute: create cairo surface
        self.context = cairo.Context(self.surface)  # Instance attribute: create cairo context

        # Store gene sequences and motif locations
        self.gene_sequences = gene_sequences 
        self.motifs = motif_locations

    # Methods:
    # def draw_figure(output_png_path)

# Main function:
def main():
    # Parse arguments and set global argparse variables:
    args = get_args()
    fasta_file_path = args.f
    motif_file_path = args.m

    # Read FASTA records
    fasta_sequences = read_fasta(args.f) # List of FastaSequence objects for all sequences in the fasta file

    # Read motifs
    motifs = read_motifs(args.m) # List of Motif objects for all motifs in the motif file

    motif_locations = {}
    for seq in fasta_sequences:
        motif_locations[seq.header] = find_motif_locations(seq, motifs)


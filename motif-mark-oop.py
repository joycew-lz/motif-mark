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
                    exons.append((start, i)) # add the exon region's position to the list of exons
                if current_type == "intron":
                    introns.append((start, i)) # add the intron region's position to the list of introns
                
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

# Convert FASTA filename to a figure PNG name
def figure_name(fasta_path):
    '''
    Given a FASTA file path, return a PNG filename for the figure
    with the same file name as the FASTA file (but with a .png extension instead of .fasta or .fa).
    '''
    # Get the base file name (without the directory path)
    output_name = os.path.basename(fasta_path)
    prefix, ext = os.path.splitext(output_name)
    return prefix + ".png"

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
        (0.733, 0.843, 0.949), # light blue
        (0.949, 0.882, 0.729), # yellow
        (0.784, 0.612, 0.475), # brown
        (0.949, 0.729, 0.788), # pink
        (0.796, 0.969, 0.643) # light green

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

    # Return all hits
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

# Class to build the figure for the gene sequences and motif marks.
class FigureBuilder:
    def __init__(self, gene_sequences, motif_locations, width=1500, height=900):
        # Store figure width and height
        self.width = width 
        self.height = height

        # Create cairo surface and context
        self.surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, self.width, self.height)
        self.context = cairo.Context(self.surface)

        # Store gene sequences
        self.gene_sequences = gene_sequences
        self.motif_locations = motif_locations

        # Fixed layout constants
        self.left_margin = 200 # Left margin is larger to allow for gene labels
        self.right_margin = 60 # Right margin is smaller to allow for more space for the gene sequences
        self.vertical_margin = 65 # Margin on top and bottom of the figure
        self.gene_spacing = 150
        self.exon_height = 50

    # Internal method to draw out the backbone gene line, exon rectangles, and gene labels for the figure,
    # and then draws motif marks on the gene structure.
    # Note, this is a public interface method (no underscore) of FigureBuilder that can be called from outside the class.
    def draw_figure(self, output_png_path):
        '''
        Draws gene lines to scale, exon rectangles to scale with the length of an exon, and gene labels.
        Then, draws motif marks on the gene structure.
        '''
        ctx = self.context

        # Set a white background color
        ctx.set_source_rgb(1, 1, 1)
        ctx.paint()

        # Set a font
        ctx.select_font_face("Sans", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        ctx.set_font_size(20)

        # Horizontal scaling factor for drawing gene lines
        max_gene_length = max(sequence.length for sequence in self.gene_sequences) # Longest gene sequence length
        available_width = self.width - self.left_margin - self.right_margin # Available width for drawing (total width minus margins on both sides)
        px_per_base = available_width / max_gene_length # Scaling factor to convert from gene sequence length to pixels

        # Draw the gene structure:
        current_y = self.vertical_margin

        for sequence in self.gene_sequences:

            # Draw intron line, aka backbone gene line, to scale
            x1 = self.left_margin
            x2 = self.left_margin + sequence.length * px_per_base
            ctx.set_source_rgb(0, 0, 0) # black
            ctx.set_line_width(5)
            ctx.move_to(x1, current_y) # sets line starting point at the left margin and current y position
            ctx.line_to(x2, current_y) # sets line ending point at the left margin plus gene length and current y position
            ctx.stroke()

            # Draw exons as rectangles, to scale
            for exon_start, exon_end in sequence.exons:
                exon_x = self.left_margin + exon_start * px_per_base # exon x position is based on exon start position and scaling factor, starting from the left margin
                exon_y = current_y - self.exon_height/2 # exon y position centered vertically on backbone line
                exon_width = (exon_end - exon_start) * px_per_base # exon width is based on exon length and scaling factor
                ctx.set_source_rgb(0, 0, 0) # black
                ctx.rectangle(exon_x, exon_y, exon_width, self.exon_height)
                ctx.fill()
        
            # Draw the gene label, which is the FASTA sequence header, above the gene structure
            ctx.set_source_rgb(0, 0, 0) # black
            ctx.move_to(10, current_y - 32) # position the label to the left of the gene backbone and slightly above
            ctx.show_text(sequence.header)

            # Draw motifs on the gene structure using the internal method _draw_motif_marks()
            self._draw_motif_marks(sequence, px_per_base, current_y)

            # Increment y position for the next gene sequence
            current_y += self.gene_spacing
        
        # Draw the figure legend for the motifs on the left, bottom side of the figure using the internal method _draw_legend()
        # self._draw_legend()

        # Save PNG with a specific name
        self.surface.write_to_png(output_png_path)
    
    # Internal method to draw motif marks on the figure.
    def _draw_motif_marks(self, sequence, px_per_base, y):
        '''
        Draws motif rectangles for a single gene and assign hits to lanes so motif overlaps do not overlap visually on the figure.
        seq = FastaSequence object
        px_per_base = scaling factor
        y = vertical position of the gene line
        '''
        ctx = self.context

        # Sort motif matches (hits) by start and then end
        hits_sorted = sorted(self.motif_locations[sequence.header], key = lambda x: x.start)

        # Store the end position of the last motif mark in each vertical lane, or track
        tracks = []

        # Constants to later determine the y position for this motif mark, based on stacking it below previous motif marks if there are overlaps!
        motif_base_height =  self.exon_height # a motif height matches exon rectangle height
        motif_stacked_height = self.exon_height * 0.5
        motif_vertical_spacing = 6

        for hit in hits_sorted:
            # Track if this motif has been placed in a lane yet
            placed = False

            # Go through each lane:
            for track_index in range(len(tracks)):
                if hit.start >= tracks[track_index]:
                    tracks[track_index] = hit.end
                    level = track_index
                    placed = True
                    break
            
            if placed == False:
                tracks.append(hit.end)
                level = len(tracks) - 1

            motif_x = self.left_margin + hit.start * px_per_base # motif x position is based on motif start position and scaling factor, starting from the left margin
            motif_width = (hit.end - hit.start) * px_per_base # motif width is based on motif length and scaling factor
            
            if level == 0: # motif centered on the gene bar
                motif_height = motif_base_height
                motif_y = y - motif_height/2
            else:
                motif_height = motif_stacked_height
                motif_y = (y + motif_height/2 + (level - 1))

            ctx.set_source_rgb(*hit.motif.color) # set color for this motif, using the color that was assigned to the motif earlier
            ctx.rectangle(motif_x, motif_y, motif_width, motif_height)
            ctx.fill()

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

    # Create figure builder instance and draw figure
    png_name = figure_name(fasta_file_path) 
    builder = FigureBuilder(fasta_sequences, motif_locations)
    builder.draw_figure(png_name) # Note, this is a public interface method of FigureBuilder 

# Run main if script is executed
if __name__ == "__main__":
    main()

# Run the script as:
# ./motif-mark-oop.py -f Figure_1.fasta -m Fig_1_motifs.txt
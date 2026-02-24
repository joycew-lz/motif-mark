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
    parser = argparse.ArgumentParser(description="to visualize motif marks along a gene sequence using pycairo")
    parser.add_argument("-f", help="designates file path to the fasta file containing sequences", type=str, required=True)
    parser.add_argument("-m", help="designates file path to the motifs text file, with one motif per line", type=str, required=True)
    return parser.parse_args()

#--------------
# classes with instance attributes and internal methods
#--------------

# Class to represent a FASTA record.
class FastaSequence:
    '''
    Represents a FASTA record, including its header, raw sequence, length, and exon/intron regions.
    '''
    def __init__(self, header, sequence):
        self.header = header
        self.sequence = sequence 
        self.length = len(sequence)
        self.exons, self.introns = self._find_regions()  # list of (start, end) tuples for exons and introns in the sequence

    # Internal method to find exon and intron regions based on capitalization of the sequence
    def _find_regions(self):
        '''
        Returns one list each (for exons as well as introns) of (start, end) tuples for the positions of
        exons and introns in the sequence.
        Exons are always "sandwiched" between introns, so if the first region is an exon, it starts at position 0 and ends at the start of the first intron. 
        If the first region is an intron, it starts at position 0 and ends at the start of the first exon.
        Exons are distinguished by uppercase bases. 
        Introns are lowercase bases.
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

        # Loop through all other base positions in the sequence starting at position 1,
        # and check for changes in capitalization to determine exon and intron boundaries
        for i in range(1, len(self.sequence)):

            # Determine the i'th base position's region type
            if self.sequence[i].isupper():
                new_type = "exon"
            else:
                new_type = "intron"
            
            # If we reach the end of the current region type and we're changing from exon --> intron, or intron --> exon:
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

# Function to read in a FASTA file and return a list of FastaSequence class objects.   
def read_fasta(fasta_file: str):
    '''
    Reads a FASTA file and return a list of FastaSequence objects, where each object contains the header beginning with ">", 
    the sequence, the sequence length, and exons and introns ranges.
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
                    # Create a FastaSequence object using the header and sequence;
                    # The FastaSequence constructor also computes the sequence length and exon/intron ranges internally
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

# Function to convert a FASTA filename to a figure PNG name.
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
    "Y": "[CTU]", # Pyrimidine
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
    Represents a motif and its raw sequence, assigned color, length, and regex pattern.
    '''
    def __init__(self, raw_sequence, color):
        self.raw_sequence = raw_sequence.strip() # Store motif raw sequence
        self.color = color # Store motif color
        self.length = len(self.raw_sequence.strip()) # Store motif length
        self.regex_pattern = self._convert_to_regex()  # Store motif's regex pattern as a string

        # print("Raw motif:", self.raw_sequence)
        # print("Regex pattern:", self.regex_pattern)
        # print()

    # Internal method to convert the motif raw sequence to a regex pattern using the IUPAC code dictionary.
    def _convert_to_regex(self):
        '''
        Converts a motif raw sequence to a regex pattern using an IUPAC code dictionary:
        For example, convert a motif like "YGCY" into a regex like "[CT]GC[CT]".
        Also, U's in RNA motifs are mapped to T via the IUPAC dictionary.
        This is used for finding motif marks along a sequence.
        '''
        # Convert raw sequence to all uppercase.
        sequence = self.raw_sequence.upper().strip()

        pattern = ""

        for base in sequence:
            if base in IUPAC_REGEX:
                pattern += IUPAC_REGEX[base]
            else:
                continue # skip unknown characters not in the IUPAC dictionary

        return pattern

# Function to read in a motif file and return a list of Motif class objects.
def read_motifs(motif_file: str):
    '''
    Reads a motif file and returns a list of Motif objects, where each object contains the motif raw sequence, assigned color,
    length, and regex pattern.
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

            # Create a Motif object using the raw sequence and assigned color;
            # The Motif constructor also computes the motif length and regex pattern internally
            motifs.append(Motif(line, color))
    
    return motifs

# Function to detect motif matches in a FASTA record's sequence.
def find_motif_locations(record, motifs):
    '''
    Detects motif matches in a FASTA record's sequence, and returns a list of MotifLocation objects, where each object contains 
    the motif match's Motif object, start position, and end position.
    Also uses lookahead regex to find overlapping matches.
    '''
    # Convert raw sequence to all uppercase and replace all U's with T's, since U in RNA should be treated as T for motif matching
    sequence = record.sequence.upper().replace("U", "T")

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

            # Create a MotifLocation object that stores the motif match's motif, start position, and end position
            motif_locations.append(MotifLocation(motif, start, end))

    # Return all hits
    return motif_locations

# Class to represent a motif match's Motif object and location on a gene sequence.
class MotifLocation:
    '''
    Represents one motif match on a gene, storing the Motif object and the match's start and end positions.
    '''
    def __init__(self, motif, start, end):
        self.motif = motif # Store motif object for this match
        self.start = start # Store match start (0-based)
        self.end = end # Store match end (end-exclusive)

# Class to build the figure for the gene sequences and motif marks.
class FigureBuilder:
    '''
    Represents a figure builder for gene sequences and motif marks along the gene.
    Introns are represented as a gene line, and exons are represented as black rectangles.
    Motifs are marked as colored rectangles on the gene structure, with overlapping motifs stacked below the gene line.
    '''

    def __init__(self, gene_sequences, motifs, motif_locations, width=1500, height=900):
        # Store figure width and height
        self.width = width 
        self.height = height

        # Create cairo surface and context
        self.surface = cairo.ImageSurface(cairo.FORMAT_ARGB32, self.width, self.height)
        self.context = cairo.Context(self.surface)

        # Store gene sequences, motifs, motif hits locations
        self.gene_sequences = gene_sequences
        self.motifs = motifs
        self.motif_locations = motif_locations

        # Fixed entire png file layout constants
        self.left_margin = 200 # Left margin is larger to allow for gene labels
        self.right_margin = 60 # Right margin is smaller to allow for more space for the gene sequences
        self.vertical_margin = 65 # Margin on top and bottom of the figure
        self.gene_spacing = 180 # Space between each gene line
        self.exon_height = 50

    # Internal method to draw out the backbone gene line, exon rectangles, and gene labels for the figure,
    # and then draws motif marks with overlapping motifs as well.
    # Note, this is a public interface method (no underscore) of FigureBuilder that can be called from outside the class.
    def draw_figure(self, output_png_path):
        '''
        Draws gene lines to scale, exon rectangles to scale with the length of an exon, and gene labels.
        Then, draws motif marks on the gene structure with overlapping motifs stacked below the gene line.
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
        current_y = self.vertical_margin # Initialize the current y position to start at the vertical margins

        # Loop to draw one gene structure for each record in the FASTA file!:
        for sequence in self.gene_sequences:

            # -----Draw intron line, aka backbone gene line, to scale-----
            x1 = self.left_margin
            x2 = self.left_margin + sequence.length * px_per_base
            ctx.set_source_rgb(0, 0, 0) # black
            ctx.set_line_width(5)
            ctx.move_to(x1, current_y) # sets line starting point at the left margin and current y position
            ctx.line_to(x2, current_y) # sets line ending point at the left margin plus gene length, and current y position
            ctx.stroke()

            # -----Draw exons as rectangles, to scale-----
            for exon_start, exon_end in sequence.exons:
                exon_x = self.left_margin + exon_start * px_per_base # exon x position is based on exon start position and scaling factor, starting from the left margin
                exon_y = current_y - self.exon_height/2 # exon y position centered vertically on backbone line
                exon_width = (exon_end - exon_start) * px_per_base # exon width is based on exon length and scaling factor
                ctx.set_source_rgb(0, 0, 0) # black
                ctx.rectangle(exon_x, exon_y, exon_width, self.exon_height)
                ctx.fill()
        
            # -----Draw the gene label, which is the FASTA sequence header, above the gene structure-----
            ctx.set_source_rgb(0, 0, 0) # black
            ctx.move_to(10, current_y - 32) # position the label to the left of the gene backbone and slightly above
            ctx.show_text(sequence.header)

            # -----Draw motifs on the gene structure using the internal method _draw_motif_marks()-----
            self._draw_motif_marks(sequence, px_per_base, current_y)

            # Increment y position by the space between each gene line for the next gene sequence!
            current_y += self.gene_spacing
        
        # -----Outside the loop, draw the figure legend for the motifs on the left, bottom side of the figure using the internal method _draw_legend()-----
        self._draw_legend(self.motifs)

        # -----Save PNG with a specific name-----
        self.surface.write_to_png(output_png_path)
    
    # Internal method to draw motif marks on the figure.
    def _draw_motif_marks(self, sequence, px_per_base, y):
        '''
        Draws motif marks for a given gene sequence and shows overlapping motifs stacked below the gene line.
        seq = FastaSequence object
        px_per_base = scaling factor
        y = vertical position of the gene line
        '''
        ctx = self.context

        # list of MotifLocation objects for this sequence, which contains the motif match's motif, start position, and end position
        hits = self.motif_locations[sequence.header]

        # -----Draw all motif hits on the gene line directly-----
        full_height = self.exon_height

        for hit in hits:
            motif_x = self.left_margin + hit.start * px_per_base # motif x position is based on motif start position and scaling factor
            motif_y = y - full_height/2 # motif y position is centered vertically on backbone line
            motif_width = (hit.end - hit.start) * px_per_base # motif width is based on motif length and scaling factor

            ctx.set_source_rgb(*hit.motif.color)
            ctx.rectangle(motif_x, motif_y, motif_width, full_height)
            ctx.fill()
        
        # -----Lane assignment so that overlapping motifs are stacked properly-----

        # Sort motif matches (hits) by start position
        hits_sorted = sorted(self.motif_locations[sequence.header], key = lambda x: x.start)
       
        lanes = [] # list of lanes, where each lane is a list of motif hits that DO NOT overlap with each other
        lane_ends = [] # list of the last end position for each lane

        for hit in hits_sorted:
            placed = False

            # Try to place this hit into an existing lane.
            # A hit can go into lane i only if its start position is at or after the end of the last hit already in that lane.
            for i, last_end in enumerate(lane_ends):
                # If the hit does not overlap with the last hit in this lane, place it in this lane
                if hit.start >= last_end:
                    lanes[i].append(hit)
                    lane_ends[i] = hit.end
                    placed = True
                    break
            
            # If the hit was not placed in any existing lane, it overlaps with all of them.
            # Create a new lane for it
            if not placed:
                lanes.append([hit])
                lane_ends.append(hit.end)

        # -----Draw smaller bars of overlapping motifs under the gene line-----

        # Parameters for drawing smaller bars below the gene line
        smallmotif_height = 8 # height for lower lanes below the gene line
        lane_gap = 8 # vertical spacing between lanes

        for lane_index, lane in enumerate(lanes):
            # y position for this lane, starting 35 pixels below the gene line (y)
            # stack each lane evenly downward with the spacing of the height of one WHOLE lane block
            smallmotif_y_uncentered = y + 35 + lane_index * (smallmotif_height + lane_gap)
            smallmotif_y = smallmotif_y_uncentered - smallmotif_height/2 # motif y position is centered vertically on the lane's y position

            for hit in lane: 
                smallmotif_x = self.left_margin + hit.start * px_per_base # motif x position is based on motif start position and scaling factor
                smallmotif_width = (hit.end - hit.start) * px_per_base # motif width is based on motif length and scaling factor

                ctx.set_source_rgb(*hit.motif.color) 
                ctx.rectangle(smallmotif_x, smallmotif_y, smallmotif_width, smallmotif_height) 
                ctx.fill()

    # Internal method to draw the figure legend for the motifs on the left, bottom side of the figure.
    def _draw_legend(self, motifs):
        '''
        Draws the figure legend for the motifs on the left, bottom side of the figure.
        '''
        ctx = self.context

        # Parameters for legend layout
        legend_x = 20 # left margin for legend
        legend_y = self.height - 150 # legend y position (top)
        box_size = 15 # size of each motif color square
        legend_spacing = 15 # spacing between motifs

        # Draw the title 
        ctx.set_source_rgb(0, 0, 0) # black
        ctx.move_to(legend_x, legend_y - box_size) # leave vertical space between title and first motif text
        ctx.show_text("Motifs (legend)")

        # Draw out motifs
        for motif in motifs:
            color = motif.color

            # Draw the colored square 
            ctx.set_source_rgb(*color) 
            ctx.rectangle(legend_x, legend_y, box_size, box_size) 
            ctx.fill()

            # Draw the motif name
            ctx.move_to(legend_x + box_size + 10, legend_y + box_size)
            ctx.show_text(motif.raw_sequence)

            # Move down for next motif square and name entry
            legend_y += box_size + legend_spacing

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

    # Map each sequence header to its list of MotifLocation objects, or motif matches
    motif_locations = {}
    for record in fasta_sequences:
        motif_locations[record.header] = find_motif_locations(record, motifs)

    # Create figure builder instance and draw the figure
    png_name = figure_name(fasta_file_path)
    builder = FigureBuilder(fasta_sequences, motifs, motif_locations)
    builder.draw_figure(png_name) # Note, draw_figure() is a public interface method of FigureBuilder

# Run main if script is executed
if __name__ == "__main__":
    main()

# conda activate my_pycairo (An environment with pycairo and necessary packages installed)
# Run the script as:
# ./motif-mark-oop.py -f Figure_1.fasta -m Fig_1_motifs.txt
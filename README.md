# Overview
`motif-mark` is a visualization tool for DNA sequences with exon/intron structure and highlights where protein-binding motifs occur and overlap.

Given a FASTA file of DNA sequences and a file containing the sequences of known motifs, this object-oriented code identifies where each motif occurs within a DNA sequence's exons and introns, assigns each motif a different color, and draws out a diagram of motif marks along the sequence using `pycairo`.

# Repository structure
- `Fig_1_motifs.txt`: A test motif text file
- `Figure_1.fasta`: A test FASTA file
- `Figure_1.png`: The expected output from the test motif and FASTA files after running `motif-mark-oop.py`
- `motif-mark-oop.py`: Object-oriented program that visualizes motifs along the gene sequences in a FASTA file

# Requirements
`argparse`
`os`
`re`
`pycairo`

# Inputs
This program requires the following inputs:

- FASTA file
    - The FASTA file (.fasta or .fa) is expected to contain one or more DNA sequences.
    - Exons are indicated as uppercase nucleotide bases. Introns are indicated as lowercase.

- Motif file
     - The motif file (.txt) has one motif per line.
     - The IUPAC dictionary indicates how ambiguous nucleotide codes (e.g., Y) are expanded into other character classes for motif matching.

# Running motif-mark-oop.py
Use the help option to view all arguments and usage details: `./motif-mark-oop.py -h`

The following command is an example of how to run this script: 
`./motif-mark-oop.py -f Figure_1.fasta -m Fig_1_motifs.txt`
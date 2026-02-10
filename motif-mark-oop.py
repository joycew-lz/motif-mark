class Sequence:
    def __init__(self, sequence):
        self.sequence = sequence  # Instance attribute

class Exon:
    def __init__(self, location, length):
        self.location = location  # Instance attribute
        self.length = length  # Instance attribute

class Intron:
    def __init__(self, location, length):
        self.location = location  # Instance attribute
        self.length = length  # Instance attribute

class Motif:
    type = "Moti_on_gene" # Class attribute
    
    def __init__(self, location, length):
        self.location = location  # Instance attribute
        self.length = length  # Instance attribute

class FASTAParser:
    def __init__(self, input_file):
        self.input_file = input_file  # Instance attribute


# read in a fasta file's gene sequence, where the sequence can be either DNA, RNA, or protein
def oneline_fasta(input_file: str, output_file: str):
    '''
    Takes a fasta file and writes a new file where each sequence is on a single line
    following its  ">" header.
    '''
    with open(input_file, 'r') as i_fh:
        with open(output_file, 'w') as o_fh:
            sequence_lines = ''
            for line in i_fh:
                line = line.strip('\n')
                if ">" in line:
                    if sequence_lines == '':
                        o_fh.write(line + '\n')
                    else:
                        o_fh.write(sequence_lines + '\n')
                        sequence_lines = ''
                        o_fh.write(line + '\n')
                else:
                    sequence_lines += line
            o_fh.write(sequence_lines + '\n')



# Class (like cat) has attributes (like color) and can do things (like functions like meow)


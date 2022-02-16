from random import choice
import random

def create_genome(length):
    DNA=""
    for count in range(length):
        DNA+=choice("CGTA")
    return DNA

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    #for k,v in complement.items():
    #    seq = seq.replace(k,v)
    return reverse_complement

contigs_sizes = [
	4000,4000,4000,4000,
	50000,50000,50000,50000
]
gap_sizes = [
	1000,2000,3000,4000,
	-100,0,5000,0 ##last gap is neglectible but just wanted this array to be same size
]
## True = forward, False = reverse
contig_orientations = [
	True, True, False, False,
	True, True, False, False
]
# C1 --1000-- C2 --2000-- rC3 --3000-- rC4 --4000-- C5 --(-100)-- C6 --0-- rC7 --5000-- rC8

genome_size = 1000000 
assert genome_size > sum(contigs_sizes) + sum(gap_sizes)
genome = create_genome(genome_size) ## genome size

last_genome_index_written = 0
contigs_file = open("simulated_contigs.fa", "w")
for index in range(len(contigs_sizes)):
    ## write contig id
    contigs_file.write(">" + str(index + 1) + "\n")

    if contig_orientations[index]:
        contigs_file.write( genome[last_genome_index_written
                            :last_genome_index_written+contigs_sizes[index]] + "\n")
    else:
        contigs_file.write( reverse_complement(genome[last_genome_index_written
                            :last_genome_index_written+contigs_sizes[index]]) + "\n")

    last_genome_index_written += contigs_sizes[index] + gap_sizes[index] 
contigs_file.close()

## no error simulation in genome(used as read for tests)
draft_genome_file = open("simulated_genome.fa", "w")
draft_genome_file.write(">1\n" + genome)
draft_genome_file.close()

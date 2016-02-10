'''This program will identify internal stop codons (TAA, TAG, TGA)
	in a nucleotide alignment and print a list of genes with
	internal stop codons. 

Note: This assumes prior filtering has been done so that:
   -the reference sequence is a conserved open reading frame.
   -each nucleotide sequence is divisible by three evenly.'''

from sys import argv
from glob import glob
from Bio import AlignIO

def openFiles(n, path, outfile):
    '''Open all input files in the directory'''
    inpath = path + "*_anoCar2"
    files = glob(inpath)
    for file in files:
        count = False
        filename = file.split("/")[-1]
        geneid = filename.split("_")[0]
        with open(file, "r") as infile:
            count, seqs = seqDict(infile, n)
            if count == True:
                countStops(geneid, seqs, outfile)
            else:
                pass

def seqDict(infile, n):
    '''Converts fasta into separate sequence objects, determine sequence names
    and create dictionary entries for each set of codons'''
    seqs = {}
    try:
        alignment = AlignIO.parse(infile, "fasta", seq_count=n)
        for item in alignment:
            for i in range(0, n):
                codons = []
                seq = str(item[i].seq)
                species = str(item[i].id.split("_")[-1])
                for j in range(0, len(seq), 3):
                    codons.append(seq[j:j +3])
                    j += 3
                # Replace terminal stop codons so the program can identify
                # remaining internal stops
                codons[-1] = "---"
                seqs[species] = codons
        return True, seqs
    except ValueError:
        return False, seqs
    
def countStops(geneid, seqs, outfile):
    '''This will remove stop codons from the sequences'''
    with open(outfile, "a") as output:
        for species in seqs:
            for idx,codon in enumerate(seqs[species]):
                codon = codon.upper()
                if codon == "TAA" or codon == "TAG" or codon == "TGA":
                    output.write(geneid + "\t" + species + "\t" + str(idx) +
                                 "\t" + str(len(seqs[species])) + "\n")
                    break
                 
def main():
    if argv[1] == "-h" or argv[1] == "--help":
        print("Usage: python countStopCodons.py <number of species> \
<path to inut directory> <output file name>")
        quit()
    else:
        n = int(argv[1])
        path = argv[2]
        outfile = argv[3]
        openFiles(n, path, outfile)

if __name__ == "__main__":
    main()

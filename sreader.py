#! /usr/bin/env python3
# -*- coding: utf-8 -*-
""" Read-in DNA sequence - one chromosome at a time. Soft repeat-masked fasta file: repeat sequences are in lower case. """

__author__  = "Ray Stefancsik"
__version__ = "0.2"

#################################################################
# Import command-line parsing module from the Python standard library.
from argparse import ArgumentParser, RawTextHelpFormatter
#import argparse
# Import Biopython modules
#from Bio import Seq, SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
# Module to open compressed files
import gzip
#################################################################


parser = ArgumentParser( description='Read-in DNA sequence - one chromosome at a time.', formatter_class=RawTextHelpFormatter )
#parser = argparse.ArgumentParser( description='Read-in DNA sequence - one chromosome at a time.', formatter_class=argparse.RawTextHelpFormatter )

# positional argument
parser.add_argument( 'input_file', help='Path to your input file.' )

# Optional argument
parser.add_argument('-z', '--gzip', help='use gzip compressed input file', action='store_true')


args = parser.parse_args()

# Is is a gzip compressed fasta file? 
compressedFasta = args.gzip

######################################################################

# input filename and path
fname = args.input_file

######################################################################
########                    Functions                        #########
######################################################################

def readFasta( pathToFile, isCompressed):
    """Read-in DNA sequence from fasta format file and return a record as a tuple of two strings, the title line (everything after the > character) and the sequence (as a plain string. Note that only the first sequence record in the input file is processed."""

    # initialise a tuple for chromosome sequence
    seqRecord = tuple()
#   print('UNCOMPRESSED?', isCompressed)

    try:
        if isCompressed:
            with gzip.open(pathToFile, 'rt') as handle:
#for testing#   print('Compressed')
                count = 0
                total_len = 0
                for title, seq in SimpleFastaParser(handle):
                    if count < 1:
                        seqRecord = (title, seq.upper()) 
                        total_len += len(seq)
                    count += 1
                    if count > 1:
                        print('WARNING: more than one sequence in fasta file!')
                return(seqRecord)
        else:
            with open(fname, 'rt') as handle:
#for testing#   print('Uncompressed')
                count = 0
                total_len = 0
                for title, seq in SimpleFastaParser(handle):
                    if count < 1:
                        seqRecord = (title, seq.upper()) 
                        total_len += len(seq)
                    count += 1
                    if count > 1:
                        print('WARNING: more than one sequence in fasta file!')
                return(seqRecord)
    except Exception as e:
        print(e)
        raise SystemExit


def fetchCompareSeq( refSequence, givePosition, yourSequence):
    """Fetch a short subset from your reference DNA sequence starting from the start coordinate you provide and compare the fetched reference sequence to the short DNA sequence that you provide."""

    seq = refSequence
    pos1 = givePosition
    nt = yourSequence
    yourSeqLength = len(yourSequence)
    if nt == seq[ (pos1-1): (pos1-1+yourSeqLength) ]:
        print('OK, match at position:', pos1, '\t', nt, seq[ (pos1-1): (pos1-1+yourSeqLength) ] )
    else:
        print('mismatch at position:', pos1, '\t', nt, seq[ (pos1-1): (pos1-1+yourSeqLength) ] )

#################################################################
# Read in a data file.
#################################################################
testSeq = readFasta(fname, compressedFasta)
print ('Sequence name:', testSeq[0], sep = '\t')
if len(testSeq[1]) > 1000000:
    print ('length:', len(testSeq[1])/1000000, 'Mbp')
else:
    print ('\tlength:', len(testSeq[1]), 'bp')

#################################################################
### validating sequence position
#################################################################
fetchCompareSeq( testSeq[1], 5, 'ACG')
fetchCompareSeq( testSeq[1], 6, 'CG')
fetchCompareSeq( testSeq[1], 7, 'G')
fetchCompareSeq( testSeq[1], 8, 'G')

#################################################################
#################################################################
#################################################################

if len(testSeq[1]) > 70:
    print(''.join( [testSeq[1][:70], '...']) )
else:
    print(testSeq[1])

### validating sequence position


#pos1 = 5
#nt = 'ACG'
#insLen = len(nt)
#seq = testSeq[1]
#if nt == seq[ (pos1-1): (pos1-1+insLen) ]:
#    print('OK, match at position:', pos1, '\t', nt, seq[ (pos1-1): (pos1-1+insLen) ] )
#else:
#    print('mismatch at position:', pos1, '\t', nt, seq[ (pos1-1): (pos1-1+insLen) ] )
#
#
#################################################################
#################################################################
#################################################################

# Store your data in a list of lists
data = list()
# define variables
chrm = ''
pos1 = ''
pos2 = ''
ref = ''
alt = ''
mtype = '' # mutation type
hgvs = '' # hgvs syntax for mutation

### Dictionary of RefSeq chromosome accession numbers
#   from https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh37.p13
# Genome version #assembly=GCF_000001405.25
GRCh37 = {
'1' : 'NC_000001.10',
'2' : 'NC_000002.11',
'3' : 'NC_000003.11',
'4' : 'NC_000004.11',
'5' : 'NC_000005.9',
'6' : 'NC_000006.11',
'7' : 'NC_000007.13',
'8' : 'NC_000008.10',
'9' : 'NC_000009.11',
'10' : 'NC_000010.10',
'11' : 'NC_000011.9',
'12' : 'NC_000012.11',
'13' : 'NC_000013.10',
'14' : 'NC_000014.8',
'15' : 'NC_000015.9',
'16' : 'NC_000016.9',
'17' : 'NC_000017.10',
'18' : 'NC_000018.9',
'19' : 'NC_000019.9',
'20' : 'NC_000020.10',
'21' : 'NC_000021.8',
'22' : 'NC_000022.10',
'X' : 'NC_000023.10',
'Y' : 'NC_000024.9'
}

### Dictionary of RefSeq chromosome accession numbers
#   from https://www.ncbi.nlm.nih.gov/grc/human/data?asm=GRCh38.p10
# Genome version: GRCh38
GRCh38 = {
'1' : 'NC_000001.11',
'2' : 'NC_000002.12',
'3' : 'NC_000003.12',
'4' : 'NC_000004.12',
'5' : 'NC_000005.10',
'6' : 'NC_000006.12',
'7' : 'NC_000007.14',
'8' : 'NC_000008.11',
'9' : 'NC_000009.12',
'10' : 'NC_000010.11',
'11' : 'NC_000011.10',
'12' : 'NC_000012.12',
'13' : 'NC_000013.11',
'14' : 'NC_000014.9',
'15' : 'NC_000015.10',
'16' : 'NC_000016.10',
'17' : 'NC_000017.11',
'18' : 'NC_000018.10',
'19' : 'NC_000019.10',
'20' : 'NC_000020.11',
'21' : 'NC_000021.9',
'22' : 'NC_000022.11',
'X' : 'NC_000023.11',
'Y' : 'NC_000024.10'
}




#################################################################
#################################################################
#################################################################
## Close file if not using 'with' statement.
#handle.close()

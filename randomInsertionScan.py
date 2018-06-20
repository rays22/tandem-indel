#! /usr/bin/env python3
# -*- coding: utf-8 -*-
""" Take a tabulated file that contains columns chromosome identifier, genomic position, and DNA sequence of corresponding insertion alleles.  Check if the inserted sequences are related to the genomic sequences surrounding the insertion sites by simulating insertions at randomised insertion sites. In other words, use the insertion alleles from the input file, but shuffle the corresponding insertion sites observed for the set of insertions also from the file. An insertion allele is defined as related to the sequence surrounding the insertion site if it creates or expands a tandem repeat. Report the proportion of repeat-altering insertions for the set in the input file along with the results of a number of randomised iterations."""

__author__  = "Ray Stefancsik"
__version__ = "0.1"

#################################################################
# Module to open compressed files
import gzip
# Import Biopython modules
from Bio.SeqIO.FastaIO import SimpleFastaParser
# Import command-line parsing module from the Python standard library.
from argparse import ArgumentParser, RawTextHelpFormatter
#import argparse
## Core tools for working with streams
#from io import StringIO
#import io
from io import StringIO
#################################################################

# parse user input

parser = ArgumentParser( description='Read-in CSV FILE that contains three columns in the following order:\n 1. chromosome identifier, 2. genomic position, and 3. DNA sequence of corresponding insertion alleles.\n\n %(prog)s allows you to check if the inserted sequences are related to the genomic sequences surrounding the insertion sites. You can do this with the original insertions in the input file or by simulating insertions at randomised insertion sites based on the ones in your input. Here an insertion allele is defined as being related to the sequence surrounding the insertion site if it creates or expands a tandem repeat. Reference genome file must be either a fasta or a gzip compressed fasta.gz file.', formatter_class=RawTextHelpFormatter )

# Limit the maximum number of allowed iterations:
imax = 10

# positional argument
parser.add_argument( 'reference_genome', help='Path to reference genome fasta file.' )
parser.add_argument( 'csv', help='Path to input CSV FILE.' )
parser.add_argument( 'iterations', type=int, default=imax, help='Enter number of iterations to run (maximum %(default)s).' )

# Optional argument
#parser.add_argument('-u', '--uncompressed', help='use uncompressed fasta file', action='store_true')
#parser.add_argument('-z', '--gzip', help='use gzip compressed fasta file', action='store_true')
parser.add_argument('-i', '--inputScore', help='compute score for original input in addition to simulation scores', action='store_false')

#IMPLEMENT OPTION: compute score for original input or not? (default = no)

args = parser.parse_args()

# Is fasta file gzip compressed ? 
#isFastaCompressed = args.gzip

######################################################################

genomeFilePath = args.reference_genome # get reference genome fasta file path from user input

######################################################################
### Read-in input file
# The input file is expected to contain the full path to the CSV file.
# The comma separated values should correspond to these comlumns:
#1. chromosome identifier, 2. genomic position, and 3. DNA sequence of corresponding insertion alleles
# PATH/INPUT_FILE_CSV, FILENAME AND PATH
######################################################################

iterations = args.iterations

if iterations > imax:
    print('The number of requested iterations is too large. The maximum allowed value is:', imax ,  '.')
    raise SystemExit

if iterations < 1:
    print('Iteration must be a natural number.')
    raise SystemExit

fname = args.csv

delimiter = ',' # comma-separated values (csv)
inputdata = []

try:
    with open(fname, 'rt') as fhandle:
        for line in fhandle:
            if line.startswith('#'): # ignore comment lines
                continue
            else:
                words = line.split(delimiter)
                chrm = words[0].strip() # chromosome
                posn = words[1].strip() # genomic position
                alll = words[2].strip().upper() # upper case allele sequence
                inputdata.append( [chrm, posn, alll] )
except Exception as eCSV:
    print('ERROR:', eCSV)

### just for testing ###
for i in inputdata:
    print(i)

######################################################################
########                    Functions                        #########
######################################################################

def readRefGenome( pathToFile, isCompressed):
    """Read-in reference genome DNA sequence from fasta format file and return a dictionary of chromosome or contig names as keys and the corresponding sequences as values (uppercase strings)."""

    def readFasta(fileHandle):

        """This function is doing the actual work for readFasta."""

        chrDict1 = dict()

        try:
            for title, seq in SimpleFastaParser(fileHandle):
                kromosome = title.lstrip().split()[0] # this works only if first word in title is the chromosome name
                chrDict1[kromosome] = seq.upper() # NOTE: convert sequence to uniform upper case
            if len(chrDict1) >0:
                return(chrDict1)
            else:
                print('ERROR: Check your fasta file!')
#       except Exception as e1:
        except Exception as e1:
            print('e1:', e1)


    if isCompressed:
        try:
            with gzip.open(pathToFile, mode='rt') as handle:
                return(readFasta(handle))
        except Exception as e2:
            print('e2:', e2)
    else:
            try:
                with open(pathToFile, mode='rt') as handle:
                    return(readFasta(handle))
            except Exception as e3:
                print('e3:', e3)

#####################################################
# Try to determine if the fasta file gzip compressed.
#####################################################
try:
    fastaFile = gzip.open(genomeFilePath, mode='rb')
    fastaFile.peek(1) # check quickly if file is indeed a gzip compressed file
    isFastaCompressed = True
#   print('found a gzip compressed file')
except OSError as e2b:
    isFastaCompressed = False # found an uncompressed file
#   print('e2b:', e2b)



### just for testing ###
print('isFastaCompressed:', isFastaCompressed)

# Load reference genome into memory - subsequent functions use it repeatedly
refGenomeRecord = readRefGenome( genomeFilePath, isFastaCompressed ) # get genomeFilePath from user input

def vcfPosValidator ( chromosomeID, genomicPosition, refAllele ):

    """Validate reference allele sequence and position. The refGenomeRecord argument should refer to a dictionary of chromosome or contig names as keys and the corresponding sequences as values (uppercase strings). Nucleotide numbering: 1st base having position 1."""

    validated = False # set default value

    refAllele = refAllele.upper() # convert string to uppercase
    yourSeqLength = len(refAllele)
    fetchedGenomic = refGenomeRecord[chromosomeID][(genomicPosition-1) : (genomicPosition - 1 + yourSeqLength) ] # this is the fetched genomic sequence
    if fetchedGenomic == refAllele:
        validated = True

    return(validated) 

######################################################################

def repeatScan(chromosomeID, genomicPosition, repeatUnit):

    """Scan insertion or deletion site for tandem repeats in both directions. The genomic position should be such that the repeatUnit completely overlaps with any segment of the tandem repeat region otherwise the function will not find it. It requires a reference genome chromosome/contig identifier and the genomic position based on numbering system where the first nucleotide is 1. The repeatUnit is a short nt sequence (without context nucleotides). Because it is intended to be used with VCF files (that provide 5-prime context nucleotides), it expects that the insertion and deletion site positions are NOT shifted to the 3-prime-most position."""

    repeatUnit = repeatUnit.strip().upper()

    def scan5prime(chromosomeID, genomicPosition, repeatUnit):

        """Scan insertion or deletion site for tandem repeats in the 5-prime direction. It requires a reference genome chromosome/contig identifier and the genomic position based on numbering system where the first nucleotide is 1. The repeatUnit is a short nt sequence."""

        if vcfPosValidator(chromosomeID, genomicPosition, repeatUnit):
            genomicPosition -= len(repeatUnit)
            return(scan5prime( chromosomeID, genomicPosition, repeatUnit ))
        else:
            return( genomicPosition + len(repeatUnit) ) # returns the position of the first nucleotide of a repeat using a 1-start based nucleotide numbering system
            

    def scan3prime(chromosomeID, genomicPosition, repeatUnit):

        """Scan insertion or deletion site for tandem repeats in the 3-prime direction. It requires a reference genome chromosome/contig identifier and the genomic position based on numbering system where the first nucleotide is 1. The repeatUnit is a short nt sequence."""


        if vcfPosValidator(chromosomeID, genomicPosition, repeatUnit):
            genomicPosition += len(repeatUnit)
            return(scan3prime( chromosomeID, genomicPosition, repeatUnit ))
        else:
            return( genomicPosition -1 ) # returns the position of the last nucleotide of a repeat using a 1-start based nucleotide numbering system.

    genomicPos1 = scan5prime(chromosomeID, genomicPosition, repeatUnit) - 1 # position of the nucleotide preceding the repeat-region
    if scan3prime(chromosomeID, genomicPosition, repeatUnit) > 0:
        genomicPos2 = scan3prime(chromosomeID, genomicPosition, repeatUnit) + len(repeatUnit) -1
    else: # do not go beyond 5-prime end
        genomicPos2 = scan3prime(chromosomeID, genomicPosition, repeatUnit) + len(repeatUnit)
    repeatSequence = refGenomeRecord[chromosomeID][genomicPos1 : genomicPos2]
    repeatSequence2 = refGenomeRecord[chromosomeID][genomicPos1 - len(repeatUnit) : genomicPos2] # look further upstream
    repeatCount = repeatSequence.count(repeatUnit)
    repeatCount2 = repeatSequence2.count(repeatUnit)
    length = len(repeatUnit)
    for i in range(length):
        if repeatCount2 > repeatCount:
            genomicPos1 -= 1
            repeatSequence = refGenomeRecord[chromosomeID][genomicPos1 : genomicPos2]
            repeatCount = repeatSequence.count(repeatUnit)
        else:
            break
        
# genomicPos1 + repeatUnit, '*', repeatCount should reconstruct the reference allele:
    mySyntax = ''.join( ['g.', chromosomeID, ':', str(genomicPos1), repeatUnit, '[', str(repeatCount), ']'] )
    results = ( mySyntax, chromosomeID, genomicPos1, repeatUnit, repeatCount  )
    return( results ) # genomicPos1 + repeatUnit, '*', repeatCount should reconstruct the reference allele

######################################################################

######################################################################
##############   works up to this point   ############################
######################################################################

def parseVCF( vcfFile, sampleID='testSample' ):

    """Read-in VCF file and parse it. It can handle uncompressed and gzip compressed VCF files as well."""

    try:
        vcf_reader = vcf.Reader(filename=vcfFile)
    except Exception as errVCF1:
        print('VCF error 1:', errVCF1)

    output= StringIO() # buffered stream as output
    delimiter = '\t' # for output formatting
    header = ( '#Sample identifier'
             , 'CHR'
             , 'START'
             , 'STOP'
             , 'REF'
             , 'ALT'
             , 'MUT_TYPE'
             , 'Genomic Position Validation'
             , 'Repeat-syntax'
             , 'Repeat on chromosome'
             , 'Repeat position'
             , 'Repeat unit'
             , 'Repeat count'
             , 'Repeat change')
    print( delimiter.join(header), file=output )
    
    try:
       for record in vcf_reader:
            for n in record.ALT:
                ref = str(record.REF).upper()
                if ref.find(',') > -1 :
#                   raise ValueError ('ERROR: can not handle multiple reference alleles in VCFs', str(record) ) # The VCF specification allows multiple reference alleles, but these were ambiguous with regards to which allele is changed to what.
                    continue # The VCF specification allows multiple reference alleles, but these were ambiguous with regards to which allele is changed to what. Skip record, because it is not informative for the scope of the study this programme is part of.
                alt = str(n).upper()
#               chrom = record.CHROM
                chrom = str(record.CHROM).upper().lstrip('CHR') # remove 'chr' before all chromosome values
                pos = record.POS
                mutType = mutTyperVCF(ref, alt)
                if vcfPosValidator(chrom, pos, ref): # genomic position validation
                    validation = 'PASS' # True # 'PASS'
                else:
                    validation = 'FAIL' # False # 'FAIL'
                repeatChange = False # The effect of the variant on repeat-unit count ( True for gain or loss, False otherwise).
                ####################
                if mutType == 'insertion':
                ####################
                    if alt.startswith(ref): # Trim common suffix.
                        alt = alt[len(ref):]
                        insertionStart = pos+len(ref) - 1
                        insertionEnd = insertionStart + 1
                        ref = '-'
                        repeatScanResults = repeatScan(chrom, insertionStart + 1, alt)
                        if repeatScanResults[4] == 0:
                            repeatScanResults = repeatScan(chrom, insertionStart - len(ref), alt) # look upstream
                        pos = repeatScanResults[2] # NOTE: do a 5-prime shift as for deletions.
                        end = pos + 1
                        repeats = repeatScanResults[4]
                        if validation: # ignore mutation records that have failed genomic position validation
                            if repeats > 0: # 'repeats == 1' means duplication
                                repeatChange = True
                        cosmicFormatted = ( chrom, pos, end, ref, alt, mutType )
                        cosmicFormatted = map( str, cosmicFormatted ) # convert to strings for printing
                        repeatScanResults = map( str, repeatScanResults ) # convert to strings for printing
                        print( sampleID, delimiter.join(cosmicFormatted), validation, delimiter.join(repeatScanResults), repeatChange, sep=delimiter, file=output )
                    else:
                        raise ValueError ('ERROR: incorrect mutation type.')
                ####################
                elif mutType == 'deletion': 
                ####################
                    if ref.startswith(alt): # Trim common suffix.
                        ref0 = ref
                        pos0 = pos # genomic position of REF allele with context nucleotides
                        ref = ref[len(alt):]
                        pos1 = pos + len(alt) # genomic position of REF allele without context nucleotides
                        repeatScanResults = repeatScan(chrom, pos1, ref)
                        repeats = repeatScanResults[4]
                        if repeats > 1: # if deletion affects repeat region
                            repeatChange = True
                            pos = repeatScanResults[2] + 1
                            end = repeatScanResults[2] + len(ref) # NOTE: do a 5-prime shift in repeat region so that equivalent alleles resulting from deletions at different positions can be identified. The 5-prime shift works better for VCF files than 3-prime shift, because VCF alleles are already 5-prime shifted. For example, compare the two shift methods when even number of nucleotides are deleted in regions that contain odd number of homo-nucleotides.
                        else:
                            pos = pos1
                            end = pos + len(ref) - 1
                        alt = '-'
                        cosmicFormatted = ( chrom, pos, end, ref, alt, mutType )
                        cosmicFormatted = map( str, cosmicFormatted ) # convert to strings for printing
                        repeatScanResults = map( str, repeatScanResults ) # convert to strings for printing
                        print( sampleID, delimiter.join(cosmicFormatted), validation, delimiter.join(repeatScanResults), repeatChange, sep=delimiter, file=output )
                    else:
                        raise ValueError ('ERROR: incorrect mutation type.')
                ####################
                else:
                ####################
                    end = pos + len(ref) - 1
                    cosmicFormatted = ( chrom, pos, end, ref, alt, mutType )
                    cosmicFormatted = map( str, cosmicFormatted ) # convert to strings for printing
                    rslts = ( 'NA','NA','NA','NA','NA' )
                    print( sampleID, delimiter.join(cosmicFormatted), validation, delimiter.join(rslts), 'NA', sep=delimiter, file=output )

    except Exception as errVCF3:
        print('VCF error 3:', errVCF3)

    contents = output.getvalue()
    output.close()
    outfileName = sampleID + '-out.tsv'
    with open(outfileName, 'w') as f_out:
        f_out.write(contents)

#################################################################

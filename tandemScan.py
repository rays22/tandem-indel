#! /usr/bin/env python3
# -*- coding: utf-8 -*-
""" Read-in VCF files and check genomic sequences upstream and downstream of insertion or deletion sites for tandem repeat expansion or reduction, respectively. Reference genome files must be fasta files. You can use compressed fasta.gz files. The results are written out into a tsv file in the current working directory."""

__author__  = "Ray Stefancsik"
__version__ = "0.5"

#################################################################
# Module to open compressed files
import gzip
# Import Biopython modules
from Bio.SeqIO.FastaIO import SimpleFastaParser
# Module to parse VCF files
import vcf
# Import command-line parsing module from the Python standard library.
from argparse import ArgumentParser, RawTextHelpFormatter
#import argparse
## Core tools for working with streams
#from io import StringIO
#import io
from io import StringIO
#################################################################

# parse user input

parser = ArgumentParser( description='Read-in VCF files (one for each tumour sample) from a csv file, convert mutation records to COSMIC format. In addition to these, check genomic sequences upstream and downstream of insertion or deletion sites for tandem repeat expansion or reduction, respectively. Reference genome files must be fasta files. You can use compressed fasta.gz files.', formatter_class=RawTextHelpFormatter )

# positional argument
parser.add_argument( 'reference_genome', help='Path to reference genome fasta file.' )
parser.add_argument( 'vcf_list', help='Path to VCF input file.' )

# Optional argument
parser.add_argument('-z', '--gzip', help='use gzip compressed fasta file', action='store_true')

args = parser.parse_args()

# Is fasta file gzip compressed ? 
isFastaCompressed = args.gzip

######################################################################

genomeFilePath = args.reference_genome # get reference genome fasta file path from user input

######################################################################
### Read-in list of VCF files.
# The input file is expected to contain the full path to the VCF file and the
# corresponding sample name, one comma-separated filename-sample pair per line:
# PATH/INPUT_FILE_VCF, FILENAME AND PATH
######################################################################

fname = args.vcf_list

vcflist = list()
delimiter = ',' # comma-separated values (csv)


try:
    with open(fname, 'rt') as fhandle:
        for line in fhandle:
            if line.startswith('#'): # ignore comment lines
                continue
            else:
                words = line.split(delimiter)
                vcfpath = words[0].strip()
                sampleID = words[1].strip()
                if len(sampleID) < 1:
                    print('ERROR: Missing sample names!')
                    raise SystemExit
                vcflist.append( (vcfpath, sampleID) )
except Exception as eVCF:
    print('vcf_list ERROR:', eVCF)

if len(vcflist) < 1: # check if VCF list is empty
    print('ERROR: Check your list of VCF files!')
    raise SystemExit

######################################################################

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
                chrDict1[kromosome] = seq.upper() # NOTE: convert sequence to uniform upper case:
            if len(chrDict1) >0:
                return(chrDict1)
            else:
                print('ERROR: Check your fasta file!')
        except Exception as e1:
            print('e1:', e1)


    if isCompressed:
        try:
            with gzip.open(pathToFile, 'rt') as handle:
                return(readFasta(handle))
        except Exception as e2:
            print('e2:', e2)
    else:
        try:
            with open(pathToFile, 'rt') as handle:
                return(readFasta(handle))
        except Exception as e3:
            print('e3:', e3)

def mutTyperVCF(refA, altA):

    """Return mutation type for allele pair from VCF file."""

    ### Sequence Ontology mutation type categories ###
    ######### and their COSMIC translations  #########
    snv = 'SUB' #		SO:0001483
    mnv = 'COMPLEX' #	SO:0002007
    insertion = 'INS' #	SO:0000667
    deletion = 'DEL' #	SO:0000159
    indel = 'COMPLEX' #	SO:1000032

    mut_type = 'ERROR: unexpected data format'

    if refA == altA:
        mut_type = 'ERROR: unexpected data format' # this is not a mutation
    elif altA == 'NONE':
        mut_type = 'NA: monomorphic reference' # this is not a mutation
    elif len(refA) == len(altA):
        if len(refA) > 1:
            mut_type = mnv # multiple-nucleotide substitution
        elif len(refA) == 1:
            mut_type = snv # or single-nucleotide substitution
        else:
            mut_type = 'ERROR: unexpected data format' # unexpected case
    elif altA.startswith(refA): # VCF-style insertion
        mut_type = insertion
    elif refA.startswith(altA): # VCF-style deletion
            mut_type = deletion
    else:
        mut_type = indel
    return(mut_type)

# Load reference genome into memory - subsequent functions use it repeatedly
refGenomeRecord = readRefGenome( genomeFilePath, isFastaCompressed ) # get genomeFilePath and isFastaCompressed from user input

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
                chrom = str(record.CHROM)
                pos = record.POS
                mutType = mutTyperVCF(ref, alt)
                if vcfPosValidator(chrom, pos, ref): # genomic position validation
                    validation = 'PASS' # True # 'PASS'
                else:
                    validation = 'FAIL' # False # 'FAIL'
                repeatChange = False # The effect of the variant on repeat-unit count ( True for gain or loss, False otherwise).
                ####################
                if mutType == 'INS':
                ####################
                    if alt.startswith(ref): # Trim common suffix.
                        alt = alt[len(ref):]
                        insertionStart = pos+len(ref) - 1
                        insertionEnd = insertionStart + 1
                        ref = '-'
                        repeatScanResults = repeatScan(chrom, insertionStart + 1, alt)
                        if repeatScanResults[4] == 0:
                            repeatScanResults = repeatScan(chrom, insertionStart - len(ref), alt) # look upstream
                        pos = repeatScanResults[2] # NOTE: do a 5-prime shift as for DEL.
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
                elif mutType == 'DEL': 
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
######################## Read in a VCF file(s).

# process VCF files
for item in vcflist:
    parseVCF( item[0], item[1] )

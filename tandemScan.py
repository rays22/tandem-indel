#! /usr/bin/env python3
# -*- coding: utf-8 -*-
""" Read-in VCF file and check genomic sequences upstream and downstream of insertion or deletion sites for tandem repeat expansion or reduction, respectively. Reference genome files must be fasta files. You can use compressed fasta.gz files."""

__author__  = "Ray Stefancsik"
__version__ = "0.1"

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
#import io
#################################################################

# parse user input

parser = ArgumentParser( description='Read-in VCF file and check genomic sequences upstream and downstream of insertion or deletion sites for tandem repeat expansion or reduction, respectively. Reference genome files must be fasta files. You can use compressed fasta.gz files..', formatter_class=RawTextHelpFormatter )

# positional argument
parser.add_argument( 'reference_genome', help='Path to reference genome fasta file.' )
parser.add_argument( 'vcf_file', help='Path to VCF input file.' )

# Optional argument
parser.add_argument('-z', '--gzip', help='use gzip compressed fasta file', action='store_true')

args = parser.parse_args()

# Is fasta file gzip compressed ? 
isFastaCompressed = args.gzip

######################################################################

# input filename and path
fname = args.vcf_file
genomeFilePath = args.reference_genome # get reference genome fasta file path from user input

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
                raise SystemExit
        except Exception as e1:
            print('e1:', e1)
            raise SystemExit


    if isCompressed:
        try:
            with gzip.open(pathToFile, 'rt') as handle:
                return(readFasta(handle))
        except Exception as e2:
            print('e2:', e2)
            raise SystemExit
    else:
        try:
            with open(pathToFile, 'rt') as handle:
                return(readFasta(handle))
        except Exception as e3:
            print('e3:', e3)
            raise SystemExit

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

#   refGenomeRecord = readRefGenome( genomeFilePath, isFastaCompressed ) # use global variable

    refAllele = refAllele.upper() # convert string to uppercase
    yourSeqLength = len(refAllele)
    fetchedGenomic = refGenomeRecord[chromosomeID][(genomicPosition-1) : (genomicPosition - 1 + yourSeqLength) ] # this is the fetched genomic sequence
    if fetchedGenomic == refAllele:
        validated = True

    return(validated) 


def repeatScan(chromosomeID, genomicPosition, repeatUnit):

    """Scan insertion or deletion site for tandem repeats in both directions. The genomic position should be such that the repeatUnit completely overlaps with any segment of the tandem repeat region otherwise the function will not find it. It requires a reference genome chromosome/contig identifier and the genomic position based on numbering system where the first nucleotide is 1. The repeatUnit is a short nt sequence (without context nucleotides)."""

    repeatUnit = repeatUnit.strip().upper()

    def scan5prime(chromosomeID, genomicPosition, repeatUnit):

        """Scan insertion or deletion site for tandem repeats in the 5-prime direction. It requires a reference genome chromosome/contig identifier and the genomic position based on numbering system where the first nucleotide is 1. The repeatUnit is a short nt sequence."""

        if vcfPosValidator(chromosomeID, genomicPosition, repeatUnit):
            genomicPosition -= 1
            return(scan5prime( chromosomeID, genomicPosition, repeatUnit ))
        else:
            return( genomicPosition - len(repeatUnit) ) # returns position of last nt before tandem repeat starts in 1-start based nucleotide numbering or repeat start in a 0-start based nt numbering system
            

    def scan3prime(chromosomeID, genomicPosition, repeatUnit):

        """Scan insertion or deletion site for tandem repeats in the 3-prime direction. It requires a reference genome chromosome/contig identifier and the genomic position based on numbering system where the first nucleotide is 1. The repeatUnit is a short nt sequence."""


        if vcfPosValidator(chromosomeID, genomicPosition, repeatUnit):
            genomicPosition += 1
            return(scan3prime( chromosomeID, genomicPosition, repeatUnit ))
        else:
            return( genomicPosition -1 ) # converting to 0-start based nucleotide numbering

    genomicPos1 = scan5prime(chromosomeID, genomicPosition, repeatUnit) + len(repeatUnit)
    genomicPos2 = scan3prime(chromosomeID, genomicPosition, repeatUnit) + len(repeatUnit) -1
    repeatSequence = refGenomeRecord[chromosomeID][genomicPos1 : genomicPos2]
#   print('repeatSequence:', repeatSequence, 'repeatUnit:', repeatUnit ) # for testing
#   print ('Tandem-test 4:\t', 'g.', chromosomeID, ':', genomicPos1, repeatUnit, '[', repeatSequence.count(repeatUnit), ']', sep='' )  # for testing
    print ( 'g.', chromosomeID, ':', genomicPos1, repeatUnit, '[', repeatSequence.count(repeatUnit), ']', sep='' ) # genomicPos1 + repeatUnit, '*', repeatCount should reconstruct the reference allele

######################################################################

def parseVCF( vcfFile ):

    """Read-in VCF file and parse it. It can handle uncompressed and gzip compressed VCF files as well."""

    try:
        vcf_reader = vcf.Reader(filename=vcfFile)
        print(vcfFile)
    except Exception as errVCF1:
        print('VCF error 1:', errVCF1)
        raise SystemExit

    
    try:
   #    print(next(vcf_reader))
#       print('samples:', vcf_reader.samples)
        print('metadata:', vcf_reader.metadata)
    #   print('filters:', vcf_reader.filters)
    #   print('infos:', vcf_reader.infos)
    except Exception as errVCF2:
        print('VCF error 2:', errVCF2)
        raise SystemExit


    
    try:
       for record in vcf_reader:
            for n in record.ALT:
                ref = str(record.REF).upper()
                if ref.find(',') > -1 :
                    raise ValueError ('ERROR: can not handle multiple reference alleles in VCFs', str(record) )
                alt = str(n).upper()
                pos = record.POS
                mutType = mutTyperVCF(ref, alt)
#               if mutation type
                print(pos, ref, alt, mutType)

    except Exception as errVCF3:
        print('VCF error 3:', errVCF3)
        raise SystemExit
    ### Trim common suffix
    # Trim a list of sequences by removing the longest common suffix while leaving all of them at least one character in length.





######################################################################
######    testing    #################################################

print('\n')
repeatScan( 'MT', 69, 'gg')
repeatScan( 'MT', 305, 'cc')
repeatScan( 'MT', 303, 'c')
repeatScan( 'MT', 66, 'g')
repeatScan( 'MT', 66, 'gg')
repeatScan( 'MT', 68, 'gg')
repeatScan( 'MT', 70, 'gg')
repeatScan( 'MT', 68, 'gg')
repeatScan( 'MT', 68, 'ggg')
repeatScan( 'MT', 65, 'gg')
repeatScan( 'MT', 66, 'gg')

######               #################################################
##################################################################################################################################


#################################################################
######################## Read in a VCF file(s).

print('\n####################### testing parseVCF #################################')

parseVCF( fname )


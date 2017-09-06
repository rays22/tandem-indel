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
######    testing    #################################################

#fname = '/Volumes/rs22/pancancer-R/example-files/24f887e1-ce82-40f3-9674-11102bd076c0.broad-snowman.20150918.somatic.indel.vcf.gz'
#fname = 'original.vcf.gz'
#fname = 'testdata/svcp-test.vcf.gz'
#fname = 'testdata/dkfz-test.vcf.gz'
#fname = 'testdata/consensus-test.vcf.gz'
##fname = 'testdata/vcf-testfile.vcf.gz' # binary garbage: BECAUSE out/out2.vcf:      POSIX tar archive
#fname = 'testdata/broad-test.vcf.gz' # end-of-line might be different than in other test files?
######               #################################################
######################################################################
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

    """Validate reference allele sequence and position. The refGenomeRecord argument should refer to a dictionary of chromosome or contig names as keys and the corresponding sequences as values (uppercase strings). Nucleotide numbering 1st base having position 1.."""

    validated = False # set default value

#   refGenomeRecord = readRefGenome( genomeFilePath, isFastaCompressed ) # use global variable

    refAllele = refAllele.upper() # convert string to uppercase
    yourSeqLength = len(refAllele)
    fetchedGenomic = refGenomeRecord[chromosomeID][(genomicPosition-1) : (genomicPosition - 1 + yourSeqLength) ] # this is the fetched genomic sequence
    if fetchedGenomic == refAllele:
        validated = True

    return(validated) 

def scan5prime(chromosomeID, genomicPosition, refAllele):

    """Scan insertion or deletion site for 5-prime repeats. It requires a 1-nt based genomic position as argument. """


    if vcfPosValidator(chromosomeID, genomicPosition, refAllele):
        genomicPosition -= 1
        return(scan5prime( chromosomeID, genomicPosition, refAllele ))
    else:
        return( genomicPosition - len(refAllele) ) # position of last nt before tandem repeat starts in 1-start based nucleotide numbering or repeat start in a 0-start based nt numbering system
        

def scan3prime(chromosomeID, genomicPosition, refAllele):

    """Scan insertion or deletion site for 3-prime repeats. It requires a 1-nt based genomic position as argument. """


    if vcfPosValidator(chromosomeID, genomicPosition, refAllele):
        genomicPosition += 1
        return(scan3prime( chromosomeID, genomicPosition, refAllele ))
    else:
        return( genomicPosition -1 ) # converting to 0-start based nucleotide numbering



def repeatScan ():

    """Scan sequence context of insertion or deletion site for tandem repeats. It requires consistently shifted (either 5' or 3' shifted) insertion/deletion positions in the input file to report back consistent boundary position. VCFs are OK as they are 5' shifted. """


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
#               print( str(record.REF).upper(), str(n).upper(), '\t-->' )
                ref = str(record.REF).upper()
                alt = str(n).upper()
#               print(mutTyperVCF(ref, alt))

######### for one record ##
#        rec1 = next(vcf_reader)
#        for n in rec1.ALT:
#            print( str(rec1.REF).upper(), str(n).upper() )
#            ref = str(rec1.REF).upper()
#            alt = str(n).upper()
#            if ref.startswith(alt):
#                ref = ref.lstrip(alt)
#                print(ref, '-', 'DEL')
#            if alt.startswith(ref):
#                ref = alt.lstrip(ref)
#                print('-', ref, 'INS')
###########################
#       for record in vcf_reader:
#           for n in record.ALT:
#               print( str(n).upper() )
#               print(type(record), record)
    except Exception as errVCF3:
        print('VCF error 3:', errVCF3)
        raise SystemExit
    ### Trim common suffix
    # Trim a list of sequences by removing the longest common suffix while leaving all of them at least one character in length.

######################################################################
######    testing    #################################################


#refGenome = readRefGenome( genomeFilePath, True )
refGenome = refGenomeRecord
print('\nTesting reference genome parsing. Chromosomes to follow:')
for keys in refGenome:
    print ('\t', keys)
print('\n')

toTest = ( 'MT', 66, 'g')
toTest = ( 'MT', 66, 'gg')
toTest = ( 'MT', 68, 'gg')
toTest = ( 'MT', 70, 'gg')
#toTest = ( 'MT', 69, 'gg')
toTest = ( 'MT', 305, 'cc')
toTest = ( 'MT', 303, 'c')

toTest = ( 'MT', 69, 'gg')
toTest = ( 'MT', 68, 'gg')
toTest = ( 'MT', 305, 'cc')
toTest = ( 'MT', 68, 'ggg')
toTest = ( 'MT', 69, 'gg')

if vcfPosValidator( toTest[0], toTest[1], toTest[2] ):
    print( 'validated', toTest )
else:
    print('failed validation', toTest )
### 
repeatUnit = toTest[2].upper()
startG = scan5prime(toTest[0], toTest[1], toTest[2]) + len(repeatUnit)
endG = scan3prime(toTest[0], toTest[1], toTest[2]) + len(repeatUnit) -1
#print('startG:', startG)
#print('endG:', endG)
repeatCount = int((endG-startG) / len(repeatUnit))
endP = startG + len(repeatUnit) * repeatCount 
print ('Tandem-test 1:', toTest[0], ':', startG, endP, endG, repeatUnit, '*', repeatCount ) # startG + repeatUnit, '*', repeatCount should reconstruct the reference allele
print ('Tandem-test 2:\t', 'g.', toTest[0], ':', startG, repeatUnit, '[', repeatCount, ']', sep='' ) # startG + repeatUnit, '*', repeatCount should reconstruct the reference allele

print(type(refGenomeRecord))
repeatCount2 = refGenomeRecord[toTest[0]][(startG) : (endG) ]
print( repeatCount2, repeatCount2.count(repeatUnit), repeatCount)

print ('Tandem-test 3:\t', 'g.', toTest[0], ':', startG, repeatUnit, '[', repeatCount2.count(repeatUnit), ']', sep='' ) # startG + repeatUnit, '*', repeatCount should reconstruct the reference allele

######               #################################################
##################################################################################################################################


#################################################################
######################## Read in a VCF file(s).

print('\n####################### testing parseVCF #################################')

parseVCF( fname )


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
# This module is to help create a unique filename for the output
from tempfile import mkstemp
# to get formatted date for simple timestamps
from datetime import date
# this module is to get the name of the current script omitting the directory part
from os import path
# Return a random element from the non-empty sequence.
from random import choice
# Module to display time
from time import strftime, gmtime, sleep
#################################################################

# parse user input

parser = ArgumentParser( description='Read-in CSV FILE that contains three columns in the following order:\n 1. chromosome identifier, 2. genomic position, and 3. DNA sequence of corresponding insertion alleles.\n\n %(prog)s allows you to check if the inserted sequences are related to the genomic sequences surrounding the insertion sites. You can do this with the original insertions in the input file or by simulating insertions at randomised insertion sites based on the ones in your input. Here an insertion allele is defined as being related to the sequence surrounding the insertion site if it creates or expands a tandem repeat. Reference genome file must be either a fasta or a gzip compressed fasta.gz file.', formatter_class=RawTextHelpFormatter )

# Limit the maximum number of allowed iterations:
imax = 32

# positional argument
parser.add_argument( 'REFERENCE_GENOME', help='Path to reference genome fasta file.' )
parser.add_argument( 'CSV_FILE', help='Path to input CSV FILE.' )
parser.add_argument( 'iterations', type=int, default=imax, help='Enter number of iterations to run (maximum %(default)s).' )

## Optional argument
#parser.add_argument('-i', '--inputScore', help='compute score for original input in addition to simulation scores', action='store_false')

#IMPLEMENT OPTION: compute score for original input or not? (default = no)

args = parser.parse_args()


######################################################################

genomeFilePath = args.REFERENCE_GENOME # get reference genome fasta file path from user input

######################################################################
### Read-in input file
# The input file is expected to contain the full path to the CSV file.
# The comma separated values should correspond to these comlumns:
#1. chromosome identifier, 2. genomic position, and 3. DNA sequence of corresponding insertion alleles
# PATH/INPUT_FILE_CSV, FILENAME AND PATH
######################################################################

iterationCount = args.iterations # How many iterations does the user requires for simulated insertions?

# Get the name of the present script for output
scriptname = path.basename(__file__) # get the name of the current script for the output filename

if iterationCount > imax:
    print('The number of requested iterations is too large. The maximum allowed value is:', imax ,  '.')
    raise SystemExit

if iterationCount < 1:
    print('Iteration must be a natural number.')
    raise SystemExit

fname = args.CSV_FILE

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

#### just for testing ###
# display time just to show that script is running and doing something
print(scriptname, strftime('is running at %Y-%m-%d %H:%M:%S', gmtime()), '(1)' )

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

######################################################################

def vcfPosValidator ( chromosomeID, genomicPosition, refAllele ):

    """Validate reference allele sequence and position. The refGenomeRecord argument should refer to a dictionary of chromosome or contig names as keys and the corresponding sequences as values (uppercase strings). Nucleotide numbering: 1st base having position 1."""

#   validated = False # set default value

    refAllele = refAllele.upper() # convert string to uppercase
    yourSeqLength = len(refAllele)
#   fetchedGenomic = refGenomeRecord[chromosomeID][(genomicPosition-1) : (genomicPosition - 1 + yourSeqLength) ] # this is the fetched genomic sequence
    fetchedGenomic = dict()
    try:
        fetchedGenomic = refGenomeRecord[chromosomeID][(genomicPosition-1) : (genomicPosition - 1 + yourSeqLength) ] # this is the fetched genomic sequence
    except KeyError:
        validated = False
    if fetchedGenomic == refAllele:
        validated = True
    else:
        validated = False # set default value

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
### results = ( repeatCount )
    return( results ) # genomicPos1 + repeatUnit, '*', repeatCount should reconstruct the reference allele

######################################################################

######################################################################
##############   works up to this point   ############################
######################################################################

def findTandemInsertions( chrom, pos, alt ):
#def findTandemInsertions( insertionData, randomiseData=False ):

    """Scan insertion data for tandem repeats from list of lists that contains insertions positions and the corresponding sequences of inserted nucleotides (in case of insertions)."""

    
    repeatChange = False # The effect of the variant on repeat-unit count ( True for gain or loss, False otherwise).

    ####################
####if mutType == 'insertion':
#   ####################
#   if alt.startswith(ref): # Trim common suffix.
#       alt = alt[len(ref):]
    ref = '-' # default reference allele for all insertions
    insertionStart = pos
#   insertionEnd = insertionStart + 1
    repeatScanResults = repeatScan(chrom, insertionStart + 1, alt)
    if repeatScanResults[4] == 0:
        repeatScanResults = repeatScan(chrom, insertionStart - len(ref), alt) # look upstream
    pos = repeatScanResults[2] # NOTE: do a 5-prime shift as for deletions.
    end = pos + 1
    repeats = repeatScanResults[4]
    if repeats > 0: # 'repeats == 1' means duplication
        repeatChange = True

#   cosmicFormatted = ( chrom, pos, end, ref, alt, mutType )
#   cosmicFormatted = map( str, cosmicFormatted ) # convert to strings for printing
#   repeatScanResults = map( str, repeatScanResults ) # convert to strings for printing
#   print( sampleID, delimiter.join(cosmicFormatted), validation, delimiter.join(repeatScanResults), repeatChange, sep=delimiter, file=output )
    ####################

    return( repeatChange )


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

# Load reference genome into memory - subsequent functions use it repeatedly
refGenomeRecord = readRefGenome( genomeFilePath, isFastaCompressed ) # get genomeFilePath from user input

###########################################################################
# get chromosome symbols for categories (a set of counters) from input data
# count the total number of elements in each category
###########################################################################


#### just for testing ###
# display time just to show that script is running and doing something
print(scriptname, strftime('is still running at %Y-%m-%d %H:%M:%S', gmtime()), '(2)' )


ctgrs = dict() # dictionary of categories
chromosomesPositions = dict() # original genomic insertion positions categorised by chromosomes

############################################################################## 
# calculate total insertion count and count of repeat altering insertions from original data
for i in inputdata:
    chrom = i[0]
    pos = int(i[1])
    alt = i[2]

    if i[0] not in ctgrs:
        ctgrs[ i[0] ] = [1] # initialise counter for 'total'
        chromosomesPositions[ i[0] ] = [pos] # initialise insertion site array for each chromosomes 
        # initialise counter for repeat altering insertions for original data
        if findTandemInsertions(chrom, pos, alt):
            ctgrs[ i[0] ].append(1) # initialise counter for repeat altering insertions for original data
        else:
            ctgrs[ i[0] ].append(0) # initialise counter for repeat altering insertions for original data
    else:
        ctgrs[ i[0] ][0] += 1 # increment count for total
        chromosomesPositions[ i[0] ].append(pos) # initialise insertion site array for each chromosomes 
        if findTandemInsertions(chrom, pos, alt):
            ctgrs[ i[0] ][1] += 1 # increment count for repeat altering insertions for original data
        else:
            continue # do not increment count for repeat altering insertions for original data


############################################################################## 
# create placeholders for values calculated from simulated insertion data
for i in range(iterationCount):
    for key in ctgrs:
        ctgrs[key].append(0)


# The first two positions in the lists are already populated with the total insertion count and repeat-altering count. Any new data from the simulations should be appended after these two.
indexShift = 2



# calculate total insertion count and count of repeat altering insertions from simulated data
for c in range(iterationCount):
    for i in inputdata:
        chrom = i[0]
        try:
            pos = choice( chromosomesPositions[ i[0] ] ) #shuffle insertion positions
        except KeyError:
            print('ERROR e1r: uninitialised chromosome dictionary!')
        alt = i[2]
        
        if findTandemInsertions(chrom, pos, alt):
            ctgrs[ i[0] ][indexShift + c] += 1 # increment count for repeat altering insertions for original data
        else:
            continue # do not increment count for repeat altering insertions for original data


############################################################################## 



# reformat dictionary keys for better sorting
formattedCtgs = dict()

for key in ctgrs:
    k = -1
    try:
        k = int(key)
    except ValueError:
#   except Exception as e1:
#       print('error:', e1)
        newKey = key
    if 9 < k < 23:
        newKey = ''.join(['chr', key ])
    elif 0 < k < 10:
        newKey = ''.join(['chr0', key ])
    else:
        newKey = key
    formattedCtgs[newKey] = ctgrs[key]



##############################################################################
# convert the dictionary of categories into a list of list to append more data
##############################################################################
# Sorting the results.
## convert dictionary into list of lists
chrmCategories = list( map( list, formattedCtgs.items() ) )
## sort list of lists by first item in each list
chrmCategories.sort(key= lambda x: x[0].lower())

# initialise list with placeholder for first putative column 
myFirstColumn = [ '#description', ['total', 'original_RA'] ]

# add the optional number of simulation headers from user input
for i in range(iterationCount):
    myFirstColumn[1].append('simulated_RA')

myColumns = [ myFirstColumn ]

# append chromosome categories and total counts
for i in chrmCategories:
    myColumns.append(i)

#################################################################
# Produce output for the results.
#################################################################
output= StringIO() # buffered stream as output

delimiter = '\t' # for output formatting (tsv)


# placeholder list for transposed rows
myRows = [ [] ] # a list of lists

################################################################
#### Get table dimensions.
## 1. Get the number of rows to create for the transposition. In more detail, gGet number of items in one of the 2 expected items from the list of lists to help create the required number of rows:
rowCount = len(myColumns[0][1]) # The item count should be the same for all list items and should have a minimum value of 2.
# 2. Get number of columns to create by the transposition:
columnCount = len(myColumns)

#### just for testing ###
#print('rowCount:', rowCount, 'columnCount:', columnCount)
print(scriptname, 'produced', rowCount, 'rows', columnCount, 'columns')
#print(myColumns[0])

# transpose table (which is a list of lists)
### make the first row
for i in range(columnCount):
    myRows[0].append(myColumns[i][0])


##### just for testing ###
#print('initial myRows', myRows , '\n')
##for i in myRows:
##    print('items in myRows \n', i)

# make placeholders for subsequent rows
for i in range(rowCount):
    myRows.append( [] )

##### just for testing ###
#print('empty myRows with placeholders', myRows , '\n')
##for i in myRows:
##    print('items in myRows \n', i)

for i in range(columnCount):
#   print('i', i)
    for j in range(rowCount):
        myRows[j+1].append(myColumns[i][1][j])


#for j in range(rowCount):
##   print('j', j)
#    singleRow = []
##   for i in range(columnCount):
#    for i in myColumns[j][1]
##       print('i', i, 'j', j)
#        singleRow.append(myColumns[i][1][i])
#    myRows.append(singleRow)

firstThreeRows = myRows

# print('transposed table')
for r in firstThreeRows:
    print(delimiter.join( map(str, r) ), file=output )



###################################
# write results into an output file
###################################
contents = output.getvalue()
output.close()

myprefix = '-'.join( ['out', scriptname, date.today().strftime('%Y-%m-%d-')]) # this is for the output filename
mysuffix = '.tsv' # this is for the output filename
tempfile, path = mkstemp(suffix=mysuffix, prefix=myprefix, dir='./', text=True) # creates a unique output filename that will not be overwritten next time you run the script

with open(tempfile, 'w') as f_out:
    f_out.write(contents)

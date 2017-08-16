#! /usr/bin/env python3
# -*- coding: utf-8 -*-
""" Get 3' shifted HGVS syntax for insertion type mutations using NCBI Variation Services
     GET /spdi/{spdi}/hgvs
 at https://api.ncbi.nlm.nih.gov/variation/v0/. """

__author__  = "Ray Stefancsik"
__version__ = "0.1"

#################################################################
import requests
import sys
# Import command-line parsing module from the Python standard library.
import argparse
from time import sleep
#################################################################

parser = argparse.ArgumentParser( description='Get correct HGVS syntax for mutation genomic coordinates from input file using NCBI Variation Services. INPUT FILE: GRCh37/38 coordinates in tabulated format (COSMIC format). FOR EXAMPLE:\n#chr\tCOLUMN_IGNORED\tpos1\tpos2\tref\talt\tCOLUMN_IGNORED\tCOLUMN_IGNORED\tCOLUMN_IGNORED\tsampleID(IGNORED) \n17\t37\t7577566\t7577567\t-\tA\tNULL\t0\t0\tS0123 \n22\t37\t20073728\t20073729\t-\tT\tNULL\t0\t0\tS0123 \n18\t37\t55990549\t55990550\t-\tCTCT\tNULL\t0\t0\tS0123 \n1\t37\t53925373\t53925374\t-\tCCGCCC\tNULL\t0\t0\tS0123 \nX\t37\t47030524\t47030525\t-\tCGACTAT\tNULL\t0\t0\tS0123 \n19\t37\t6230649\t6230650\t-\tCAGGTGGTT\tNULL\t0\t0\tS0123 \n9\t37\t90502205\t90502206\t-\tG\tNULL\t0\t0\tS0123 \n16\t37\t72992264\t72992265\t-\tTT\tNULL\t0\t0\tS0123 \n5\t37\t37024707\t37024708\t-\tA\tNULL\t0\t0\tS0123 \n10\t37\t23481740\t23481741\t-\tCGG\tNULL\t0\t0\tS0123 \n10\t37\t23481740\t23481741\t-\tCGG\tNULL\t0\t0\tS0123 \n7\t37\t128481288\t128481289\t-\tGAT\tNULL\t0\t0\tS0123 ', formatter_class=argparse.RawTextHelpFormatter )

parser.add_argument( 'input_file', help='Path to your input file.' )

# 1. Optional, but mutually exclusive user options:
group1 = parser.add_mutually_exclusive_group()
group1.add_argument('-c', '--csv', help='Use comma as column separator (default).', action='store_true')
group1.add_argument('-t', '--tab', help='Use tab as column separator.', action='store_true')

# 2. Optional, but mutually exclusive user options:
group2 = parser.add_mutually_exclusive_group()
group2.add_argument('-37', '--grch37', help='Use GRCh37 coordinates (default).', action='store_true')
group2.add_argument('-38', '--grch38', help='Use GRCh37 coordinates.', action='store_true')

# 3. Optional, but mutually exclusive user options:
group3 = parser.add_mutually_exclusive_group()
group3.add_argument('-i', '--cosmic', help='Use COSMIC format columns (default).', action='store_true')
group3.add_argument('-v', '--vcf', help='Use VCF format columns (NOT IMPLEMENTED).', action='store_true')

args = parser.parse_args()

######################################################################

# input filename
fname = args.input_file

#################################################################
# Read in a data file.
#################################################################

try:
    fhand = open( fname )
except:
    print('File cannot be opened:', fname)

### 1. Input file format:
if args.tab:
    delimiter = '\t' #tab separated values format
else:
    delimiter = ',' #csv file format
### 2. Input file format:
if args.vcf:
    vcf = 'y' #VCF-like format
    delimiter = '\t' #VCF assumes tab separated values format
else:
    vcf = 'n' #COSMIC-type format with extra column[1] for genome build


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

# Chose genome version here:
if args.grch38:
    RefSeqChr = GRCh38 #Use GRCh38 genome build coordinates. 
else:
    RefSeqChr = GRCh37  #Use GRCh37 genome build coordinates. 


# Read in a data line-by-line
for line in fhand:
    if line.startswith('#'): # skip header line
        continue
    line = line.rstrip()
    columns = line.split(delimiter)
    if vcf == 'y': #VCF-type columns
        try:
            chrm = columns[0]
            pos1 = columns[1]
            pos2 = '?'
            ref = columns[3]
            alt = columns[4]
            sample = columns[2] # sample identifier
        except:
            print( 'ERROR: incorrect format\t', line)
            sys.exit()
    else: #COSMIC-type columns
        try:
            chrm = columns[0]
            pos1 = columns[2]
            pos2 = columns[3]
            ref = columns[4] # use '-' for insertions
            alt = columns[5] # use '-' for deletions
#           sample = columns[9] # sample identifier
        except:
            print( 'ERROR: incorrect format\t', line)
            sys.exit()

    try:
        nc = RefSeqChr[chrm] # convert chromosome number to RefSeq NC accession number
    except:
        nc = 'ERROR: Not a valid chromosome.'

## Allele in SPDI syntax of SEQUENCE-ID:POSITION:DEL-SEQUENCE:INS-SEQUENCE
##  (or just deletion length: SEQUENCE-ID:POSITION:DEL-LEN:INS-SEQUENCE)
#
#### example input:
##   NC_000001.10:32741982:0:CCCAGAAGC
##   NC_000001.10%3A32741982%3A0%3ACCCAGAAGC/hgvs
##   NC_000007.13:g.128481289_128481291delGAT
##   NC_000007.13:128481287:CGATG:CG -> NC_000007.13:g.128481290_128481292delATG ##!! 3' shifted del example!!
##   7	128481287	CGATG	CG##!! 3' shifted del example!!
##   7:g.128481289_128481291delGA##!! 3' shifted del example!!T
#############################################################
# NOTE: The SPDI data model uses a 0-based coordinate system.
#############################################################

    if ref == '-': # insertions
        lngth = '0'
        mtype = 'ins'
        hgvs = chrm + ':g.' + pos1 + '_' + pos2 + mtype + alt
        spdi = '%3A'.join( [nc, pos1, lngth, alt] )
    else:
        try:
            pos01 = str( int(pos1) - 1 )
        except:
            print('ERROR: incorrect input format')
        try:
            pos02 = int(pos2) - 1
        except:
            pos02 = int(pos01) + len(ref) - 1
        if alt == '-': # deletions
            alt0 = ''
            lngth = str(len(ref))
            mtype = 'del'
            #pos2 = str( int(pos1) + len(ref) )
            hgvs = chrm + ':g.' + pos1 + '_' + pos2 + mtype + ref
            spdi = '%3A'.join( [nc, pos01, ref, alt0] )
        else: # substitutions and complex mutations
            #hgvs = chrm + ':g.' + pos1 + '_' + pos2 + ref + mtype + alt
            lngth = pos02 - int(pos01)
            if lngth == 0:
                mtype = 'sub'
                mtypeS = '>'
                hgvs = chrm + ':g.' + pos1 + ref + mtypeS + alt
                spdi = '%3A'.join( [nc, pos01, ref, alt] )
            else:
                mtype = 'complex'
                mtypeS = '>'
                hgvs = chrm + ':g.' + pos1 + ref + mtypeS + alt
                spdi = '%3A'.join( [nc, pos01, ref, alt] )
    if vcf == 'y': # for VCF-format files only
        data.append( [sample, hgvs, spdi, chrm, pos1, pos2, ref, alt, mtype.upper()] )
    else: # for COSMIC-type inputs
        data.append( [hgvs, spdi, chrm, pos1, pos2, ref, alt, mtype.upper()] )
# close file
fhand.close()

## test print of converted data
#for i in data:
#    print (i)

#################################################################
#################################################################
#################################################################

# GET /spdi/{spdi}/hgvs


def GetHGVSsyntax(SPDIsyntax):
    """GET HGVS syntax for spdi input from NCBI Variation Services """
    ## server for NCBI Variation Services
    server = 'https://api.ncbi.nlm.nih.gov/variation/v0/spdi/'
    headers={ "Accept" : "application/json"}
    options = '/hgvs'
    r = requests.get(server+SPDIsyntax+options, headers)
    if not r.ok:
        #return(r.raise_for_status())
        #return('ERROR:could not convert it!')
        return( r.json())
    decoded = r.json()
    #print(repr(decoded))
    #print( type(decoded), decoded )
    result = decoded['data']['hgvs']
    # convert to human readable chromosome number
    result[result.find('NC_00000')+8 : result.find('.')]
    chrm = result[result.find('NC_00000')+8 : result.find('.')]
    if chrm == '23': # convert chromosome number to letter
        chrm = 'X'
    if chrm == '24': # convert chromosome number to letter
        chrm = 'Y'
    syntax2 = result[result.find(':') : ]
    result = ''.join([ chrm, syntax2 ])
    return(result)
    #print(type(result))

# counter for sleep functionality
sleepcounter = 0

for i in data:
    if vcf == 'y': #yourMutType would be incorrect for VCF-type input
        yourMutType = i[8]
        yourInput = ','.join(i[3:8])
    else: # for COSMIC-type inputs
        yourMutType = i[7]
        yourInput = ','.join(i[2:7])
#####################
    sleepcounter += 1
    # Guidelines for Scripting Calls to NCBI Servers
    # https://www.ncbi.nlm.nih.gov/home/about/policies.shtml
    if sleepcounter % 3 == 0: # Adhering to NCBI guidelines for scripting calls to their servers.
        sleep( 1 )
#####################
    if vcf == 'y': #yourMutType would be incorrect for VCF-type input
        print( yourInput, GetHGVSsyntax(i[2]), sep='\t')
    else: # for COSMIC-type inputs
        print( yourInput, GetHGVSsyntax(i[1]), sep='\t')

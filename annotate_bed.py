#!bin/usr/python3
# -*- coding: utf-8 -*-.
import sys
from os import path
import argparse
from shutil import which
from subprocess import call
from time import sleep

parser=argparse.ArgumentParser(description='''Find overlapping genomic regions of a query file in a genome file using bedtools. Original regions of both files for each overlap will be also saved (options '-wa' '-wb' from bedtools). Gene IDs must identify genes in each genomic region of the genome file.''')
parser.add_argument('query_file',help='File of genomic regions to query.')
parser.add_argument('genome_file',help='Search overlappings in this file containing (annotated) genomic regions')
parser.add_argument('-o','--output',help='Save results into this file instead of the default.')

args = parser.parse_args()
if which('bedtools') is None:sys.exit("Couldn't locate 'bedtools'. Check it is installed in your system.")

if args.output is None:
 dirn,filen=path.split(args.query_file)
 if path.splitext(filen)[1]=='.gz': filen=path.splitext(args.query_file)[0]
 args.output=path.join(dirn,"annot_"+filen)
 if path.isfile(args.output):
  sys.exit(f"File '{args.output}' exists. Choose another location or force replacing file with '-o' option.")

#File B is loaded into memory (most of the time).
cmd_bedtools = "bedtools intersect -a {query_file} -b {genome_file} -wa -wb".format(**vars(args))
#Dprint(cmd_bedtools)
if True:
	with open(args.output,'w') as outf:
		retcode=call(cmd_bedtools.split(),stdout=outf)
	#Vprint("Genes from '{}' found in '{}' regions saved in: '{}'".format(args.genome_file,args.query_file,args.output))
	#m="File '{}' created. Find the enriched GO terms running:\npython3 {} -o {}"
	#print(m.format(args.output, 'goea_bedfile.py', path.join(path.dirname(args.query_file), 'results')), file=sys.stderr)
else:retcode=call(cmd_bedtools.split())
print("Annotated BED file saved into: '{}'".format(args.output))


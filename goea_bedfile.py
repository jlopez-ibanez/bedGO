#!bin/usr/python3
# -*- coding: utf-8 -*-.
import sys
from os import path,mkdir
import argparse
from subprocess import call
from tempfile import gettempdir,mkstemp
from time import strftime #,time
import gzip

pydir=path.curdir if path.dirname(sys.argv[0]) == '' else path.dirname(sys.argv[0])
pydir=path.relpath(pydir)
parser=argparse.ArgumentParser(description='''Run TopGO to find overrepresented GO terms in a file containing GeneIDs.''')
parser.add_argument('annotated_regions', help='An annoted .bed file (file of genomic regions with geneIDs associated).')
parser.add_argument('-m','--mapfile',default=path.join(pydir,'gene2go.gz'), help='Mapping of GeneIDs to GO IDs (gene2go.gz).')
parser.add_argument('-t','--taxid',default='9606', help="TaxID to search in 'gene2go.gz' for creating the map file.")
parser.add_argument('-T','--topgo',default=path.join(pydir,"TopGO.r"), help="Path to the script 'TopGO.r'")
parser.add_argument('-c','--geneIDcol',default=-1,type=int, help='GeneIDs are in this column of the input.')
parser.add_argument('-go','--go_ontology',default='BP',type=str, help='Comma-separated list of GO categories to include in the enrichment analysis: BP, MF or CC.')
parser.add_argument('-pv','--threshold', default=0.001, type=float, help='Threshold to filter out results.')
parser.add_argument('-o','--output', help='Results will be saved in this folder.')
args = parser.parse_args()


if args.output is None:
 args.output=path.dirname(args.annotated_regions)
elif not path.isdir(path.dirname(args.output)): mkdir(args.output)
elif path.isfile(args.output): sys.exit("ERROR: Using path to a filename instead of folder: '{}'".format(args.output))
valid_go_ontologies=['BP','MF','CC']
args.go_ontology=",".join([ont for ont in  args.go_ontology.split(",") if ont in valid_go_ontologies])
outputf=args.go_ontology.split(",")[0]+"_"+path.basename(args.annotated_regions)+'.tsv'
args.output=path.join(args.output,outputf)


if not path.isfile(args.mapfile):
 sys.exit("ERROR: Couldn't locate file mapping GeneIDs to GO terms (looking for: '{}')".format(args.mapfile))
with gzip.open(args.mapfile,'rt') as inpf:
	try: line=inpf.readline().split()
	except IOError: sys.exit(f"ERROR: File mapping Gene ID to GO IDs ({args.mapfile}) must be gzipped.")
print("Using '{}' to map genes to GO terms".format(args.mapfile),file=sys.stderr)
if len(line)>2:
	geneid2go={}
	with gzip.open(args.mapfile, "r") as inpf:
		for line in inpf:
			#if line.decode()[0]=="#":continue
			l=line.decode().strip().split("\t")[:3]
			#if has entered here, this shouldn't happen 
			if len(l[1].split(","))>1:sys.exit("ERROR: Found an unexpected number of columns in GeneID mapping file. Check '{}' and run again".format(args.mapfile))
			if not l[0]==args.taxid:continue
			if not l[1] in geneid2go: geneid2go[l[1]]=[l[2]]
			else:geneid2go[l[1]].append(l[2])
	#9606_+basename()
	newmapfile=path.join(path.dirname(args.mapfile),args.taxid+"_"+path.basename(args.mapfile))
	r="{}\t{}"
	print("Saving mappings of taxid '{}' into '{}'".format(args.taxid,newmapfile), file=sys.stderr)
	with gzip.open(newmapfile,'wt') as outf:
		lines=(r.format(k,",".join(set(v))) for k,v in geneid2go.items())
		outf.write(("\n".join(lines)+"\n"))
	args.mapfile=newmapfile

#Dprint(vars(args))
#args.topgo=path.join(pydir,"TopGO.r")
cmd_r="Rscript {topgo} {annotated_regions} {mapfile} {geneIDcol} {output} -go {go_ontology} -pv {threshold}".format(**vars(args))
print("Running:",cmd_r,file=sys.stderr)
#el nombre del archivo de salida sin el prefijo de ontologia
basef="_".join(path.basename(path.splitext(args.output)[0]).split("_")[1:])
logfile=path.join(gettempdir(),'topgo_'+strftime('%d%H%M')+basef+".log")

with open(logfile,'w') as outf: call(cmd_r.split(),stdout=sys.stdout,stderr=outf)
print("Details about topGO execution saved into: '{}'".format(logfile))


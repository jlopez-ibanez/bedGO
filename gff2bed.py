#!bin/usr/python3
# -*- coding: utf-8 -*-.
import sys
from os import path
import gzip
import argparse
import json
import re

pydir=path.curdir if path.dirname(sys.argv[0]) == '' else path.dirname(sys.argv[0])

parser=argparse.ArgumentParser(description='Convert gff file to bed format. Input may be gzipped or not but output will be always gzipped. In case genes in each region were annotated, their IDs will be included in the file.')
parser.add_argument('genome_file', help='A GFF file either compressed (gzip) or not.')
parser.add_argument('-m','--mappings', default='mappings_refseq2ucsc.json', metavar='SEQID_MAPPINGS',help="A file including conversions of the seq IDs used in 'genome_file' to UCSC seqID format.")
parser.add_argument('--regid',action='store_true',default=False, help='Each region is identified by a primary ID. Default is NOT to include that ID in the output.')
parser.add_argument('-o', '--output', help='This option overrides default output filename.')
args = parser.parse_args()

args.mappings=path.join(pydir,args.mappings)
f,ext=path.splitext(args.mappings)

if ext == '.json':
 with open(args.mappings, 'r') as inpf: seqid_ucsc = json.load(inpf)
else:
 with open(args.mappings, 'r') as inpf:
  seqid_ucsc = {}
  for i,line in enumerate(inpf):
   seqid,ucsc = line.strip().split()[:2]
   if seqid.lower().startswith('chr'):
    ucsc,seqid = seqid,ucsc
   elif not ucsc.lower().startswith('chr'):
    sys.exit("ERROR: UCSC seqID not found in line {i} ('{ucsc}'-'{sequid}') of the mapping file ({args.mappings})")
   seqid_ucsc[seqid.upper()] = ucsc.lower()

d,f = path.split(args.genome_file)
basef,ext = path.splitext(f)

if args.output is None:
 outfile = path.join(d,basef.split(".gff3")[0]+'.bed.gz')
 if path.isfile(outfile): sys.exit(f"File '{outfile}' already created. Choose another location or force replacing file with  '-o' option.")
else: outfile = args.output
if not path.splitext(outfile)[-1] == ".gz":outfile+=".gz"
#Dprint("Results will be saved into:", outfile,file=sys.stderr)

gid_pattern=re.compile('(GeneID)(?::|=)([\w\d]+)',re.IGNORECASE)
rid_pattern=re.compile('ID=([\w\d]+)')
#the file include genes identified in the region
ncols=5 if args.regid else 4
lf='\t'.join(['{}']*ncols)+'\n'

if ext=='.gz':
 inpf=gzip.open(args.genome_file,'rt')
else:
 inpf=open(args.genome_file,'rt')

#it could be faster checking first all seqIDs are in my dict then processing the file?
with gzip.open(outfile,'wt') as outf:
 for i,line in enumerate(inpf):
  #if chr(line[0])=='#':continue
  if line[:1]=='#':continue
  if not line[:3]=='NC_':continue
  l=line.strip().split("\t")
  if l[2]=='gene':
   seqid = l[0].split(".")[0].upper()
   if seqid.lower().startswith('chr'):
    eq_seqid = seqid
   elif seqid in seqid_ucsc:
    eq_seqid = seqid_ucsc[seqid]
   else: sys.exit(f"ERROR: Couldn't find a valid mapping of {seqid} to UCSC seqID. Please provide a file with mappings using '--mappings' argument.")
   gids=gid_pattern.findall(l[-1])
   catgid=gids[0][0]
   gene_id=",".join([gid for cat,gid in gids])
   if len(gids)>1:
    regid=rid_pattern.findall(l[-1])[0]
    print("Multiple {} for {}: {}. Check this is not an error...".format(catgid,regid, gene_id), file=sys.stderr)
   if args.regid:
    regid=rid_pattern.findall(l[-1])[0]
    newl=lf.format(eq_seqid,int(l[3])-1,int(l[4])-1,regid,gene_id)
   else:
    newl=lf.format(eq_seqid,int(l[3])-1,int(l[4])-1,gene_id)
   outf.write(newl)
 inpf.close()
print("BED format of input file saved into: {}".format(outfile))

# BedGO

Enrichment analysis of Gene Ontology (GO) Terms associated to the genomic regions of a BED file.

GO terms considered for the enrichment may correspond to any of these categories:
 - Biological Process (BP)
 - Cellular Component (CC)
 - Molecular Function (MF)

**Only GO terms from the Biological Process (BP) category are searched by default.**

The scripts to be run will depend on the files from which the user starts from. The most simple usage would be finding GO terms enriched in a list of genes given a reference genome in BED format. See _Usage_ section for detailed information about each script.

## Requirements
   
The number indicates the version with which were tested the commands in the _Usage_ section. 
 
- Python (3.7)
- R (3.6)
  - topGO (v. 2.37) from [Bioconductor](http://www.bioconductor.org/)
- Bedtools (2.30.0)

The recommended way of installing BedGO is using [conda](https://docs.conda.io/en/latest/miniconda.html):

	`conda create -n BedGO -c bioconda python bioconductor-topgo bedtools`

With conda you don't have to modify *system* files nor interfer with the current installation in your system. Instructions about how to install conda in your system: https://docs.conda.io/en/latest/miniconda.html
 
## Usage 
Commands in this sections assume you are working in the same folder as the python scripts.
The procedence and instructions for downloading the files used in each example it is described in _Data Sources_ section. 

###1. Generate a BED file from a reference genome in GFF3 format
#### _gff2bed.py_
 Extract from a _GFF3_ file the coordinates of the chromosomic regions and the GeneIDs included in them. _GFF3_ is a specification of the _General Feature Format_ ([GFF](http://gmod.org/wiki/GFF3)). There are other (older) specifications based on GFF that **may not be a valid input**.

 Sequence IDs would be automatically converted to [https://genome.ucsc.edu/FAQ/FAQformat.html](UCSC) format (chr1, chrX, chrM etc.). The execution will stop if the conversion of sequence IDs fails. Use '--mappings' option to indicate using a file instead of the JSON files for this conversion.
 The filename of the input will be used for saving the results but replacing '.gff3' with '.bed'. This can be overriden with  '--output' option.

		`python3 gff2bed.py --mappings <seqids-refseq2ucsc> <genome.version.gff3.gz>`

###2. Annotate a BED file
#### _annotate_bed.py_
 Creates a new BED file with the GeneIDs corresponding to each genomic region using a reference annoted genome.
 This script will fail if bedtools is not installed. See _Setup_ section for more information.
 A new file prefixed by 'annoted_' will be saved in the same location as the input BED file. This can be overriden with '-o' option.

		`python3 annotate_bed.py <genomic_regions.bed> <genome_assembly.bed.gz>`
		

###3. Find overrepresented GO terms in a BED file (annotated)
#### _goea_bedfile.py_
 This script needs to load `gene2go.gz` (see _Data Sources_ section) to create a mapping file with the format required by TopGO (see _Setup_ section for more information).
 Use '-m' option to indicate the location of `gene2go.gz` file.
 In case you want to run the script with another BED file, this will create a mapping file with the format used by topGO.
 Resulting files will be saved in the folder indicated using '-o' option. Otherwise they will be created in the current working directory.
 By default only enriched terms from Biological Process (BP) category are considered. GO terms from Cellular Component (CC) and Molecular Function (MF) categories can be selected or included in the same execution with '-go' option.
 	
		`python3 goea_bedfile.py <genomic_regions.bed> -m geneid-goids.map.gz -o results/`
		
## Example

Find the enriched GO terms associated with the genes found in the genomic regions of the file *example/example.bed*
Note that because of its size, **the reference genome used here is not included in the *example* folder**. You'll need to download it from [here](https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.109/ref_GRCh38.p12_top_level.gff3.gz)

1.
	`python3 gff2bed.py --mappings mappings_refseq2ucsc.json example/ref_GRCh38.p12_top_level.gff3.gz`
2. 
	`python3 annotate_bed.py example/ref_GRCh38.p12_top_level.bed.gz example/example.bed`
3. 
	`python3 goea_bedfile.py example/annot_example.bed -m gene2go.gz -go BP,MF,CC -o results/`

## Data Sources

The files included with BedGO and used in the example commands of _Usage_ section were downloaded from the NCBI ftp site.

### Reference Genome in '.gff3' format

	https://ftp.ncbi.nlm.nih.gov/genomes/archive/old_refseq/Homo_sapiens/ARCHIVE/ANNOTATION_RELEASE.109/ref_GRCh38.p12_top_level.gff3.gz
	
  >This can be used as input of `gff2bed.py` to generate a list of Gene IDs associated to genomic regions. This list need to be created only once and used many times as input of `annotate_bed.py` for annotating a given BED file.
  
### Mapping of Gene IDs to GO terms.

	https://ftp.ncbi.nih.gov/gene/DATA/gene2go.gz

 >This file contains mappings of Gene IDs from multiple organisms to GO IDs. Place the file in the same folder as the `goea_bedfile.py` script or use `-m` option to indicate another location. 
	Note that since it was downloaded, this file may have changed and differ from the one provided with BedGO.


#!/usr/bin/env Rscript

library(topGO)

orders <- commandArgs(trailingOnly = TRUE)

# Command-line args as a matrix of 2 columns
args_matrix <- matrix(orders, ncol = 2, byrow = T)
                                                  
# These args are mandatory
input_file <- args_matrix[1,1]
map_file <- args_matrix[1,2]
gene_col <- as.integer(args_matrix[2,1])
output_f <- args_matrix[2,2]

# Until END these args are not used... ALWAYS a list of genes to perform enrichment...
score <- args_matrix[match("-score", args_matrix),2]
if (!is.na(score)){
  scores <- unlist(strsplit(score, ",|-"))
  gene_score_col <- as.integer(scores[1])
  gene_score <- as.numeric(scores[2])
}
#### END ####
mode <- unlist(strsplit(args_matrix[match("-go", args_matrix),2], ","))
pval_thres <- as.numeric(args_matrix[match("-pv", args_matrix),2])

errors <- ""
map <- readMappings(map_file)
gene_universe <- read.table(file = map_file, header = F, comment.char="", sep="\t", quote="")[,1]
#this object is an instance of the class groupStats to be used as 2nd arg to getSigGroups
classic <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher_Test")
cat(sprintf("Processing '%s' ...\n", input_file))
# ERROR Empty File
fsize<-file.info(input_file)$size
if (is.na(fsize) | fsize == 0) {
	cat("ERROR! Non-existing or empty file: '",input_file,"' (",fsize,")",sep="",end='\n')
}
bed_file <- read.table(input_file, header = F, sep = "\t", quote = "")
if (is.na(score)){
	#colum of the input file with the genes to perform enrichment
	if (gene_col<0) {gene_col=gene_col+1+length(bed_file[0,])}
	genes_list <- bed_file[, gene_col]
} else {
	genes_list <- bed_file[bed_file[,gene_score_col] > gene_score, gene_col]
}
genes_list <- factor(as.integer(gene_universe %in% genes_list))
names(genes_list) <- gene_universe

# ERROR No gene-score above cut-off
if (length(genes_list) == 0)	{
	errors <- paste(errors,'ERROR! No score above cut-off',sep='\n')
}
# ERROR No mapping genes
if(!any(genes_list %in% gene_universe))	{
	errors <- paste(errors,'ERROR! Not found any mappings of GeneIDs to GO terms. Check you are selecting the column with the GeneIDs in the input file and the right file for mapping GeneIDs to GO terms.',sep='\n')
}
if (!errors=='') {
	cat(errors,'\n')
	q()
}
#FullBasicTopGoAnalysis
genes_query <- names(genes_list)[genes_list == 1]
#The geneList object (genes_list here) is a named factor that indicates which genes are interesting and which not. To compute such a named vector with the user own predefined list of interesting genes.
flush.console()
for (ont in mode)	{
	#print(sprintf("Generating %s Results", ont))
	cat("Searching enriched GO terms from:",ont,end='\n')
	if (ont!='BP') {
   	results_filename <- paste(dirname(output_f), sub("BP", ont,basename(output_f)), sep = "/")
  } else {	results_filename <- output_f}
  #cat(results_filename)
  graph_pdfname <- paste(dirname(output_f), "graph_",basename(output_f), sep = "")
  TopGOobject <- new("topGOdata", ontology = ont, allGenes = genes_list, annot = annFUN.gene2GO, gene2GO = map)
  #The main function for running the GO enrichment
 	resultsclassic <- getSigGroups(TopGOobject, classic)
 	sco<-score(resultsclassic)
 	tops <- length(sco[sco < pval_thres])
 	#topsgraph <- length(sco[sco < 0.00001])
 	if (tops == 0) {
	 	cat("Zero results above the threshold used\n")
 	} else {
 		results_table <- GenTable(TopGOobject, classic = resultsclassic, orderBy = "classic", ranksOf = "classic", topNodes = tops)
 		genes_in_goes <- genesInTerm(TopGOobject, results_table$GO.ID)
 		colnames(results_table)[6] <- "Fisher Test"
 		intersected <- lapply(genes_in_goes, function(x) intersect(x, genes_query))
 		genes_of_goes <- unlist(lapply(intersected, function(x) paste(x, collapse = ":")))
 		results_table <- cbind(results_table, genes_of_goes)
 		colnames(results_table)[7] <- "Significant Genes"
 		write.table(results_table, file = results_filename, sep = "\t", quote = FALSE, row.names = F, col.names = T)
 		cat("Results saved into: '",results_filename,"'\n",sep="")
	}
}
  #printGraph(TopGOobject, resultsclassic, firstSigNodes = topsgraph, fn.prefix = graph_pdfname, useInfo = "all", pdfSW = T)

# script for finding closest features using GRanges
library(GenomicRanges)
library(rtracklayer)

ref = import.gff3('~/pvulgaris/annotation/Pvulgaris_218_gene.gff3')

# extract genes
genes = subset(ref, ref$type == 'gene')
snps = import.bed('sign_snps.bed')


find_dist = distanceToNearest(genes, snps, ignore.strand = T)

# convert to db
find_dist = as.data.frame(find_dist)
max_dist = 50000

close_SNPs = find_dist[find_dist$distance<=max_dist,]

# get geneID

closestGenes = genes[close_SNPs$queryHits,]

final_genes = as.data.frame(closestGenes$ID)
genes_dist = cbind.data.frame(final_genes, close_SNPs$distance)
colnames(genes_dist)=c('ID', 'distance')

write.table(final_genes, 'sign_genes.txt', sep='\t', col.name= F, row.name = F, quote = F)
write.table(genes_dist, 'genes_info.txt', sep='\t', col.name= F, row.name = F, quote = F)


#### helper functions for LGapp

########################################
#### Data Conversion ###################
########################################
convertData =  function(pf, vcf){ # project folder, vcf file
		setwd(pf)
		dir.create('Data_conversion')
		setwd('Data_conversion')
		file.copy(as.character(vcf[1,4]), '.')
		vcf2lfmm(as.character(vcf[1,1]))
		}

########################################
#### Climatic Data #####################
########################################
#### Plot map of coordinates
plotMap=function(coord){
	n=read.csv(coord, sep=',', header=T)
### you have to modify the dataset
	coord=paste(n[,2], n[,3], sep=':')
	f=cbind.data.frame(n, coord)
	colnames(f)[1]= 'ID'
	l=gvisMap(f, 'coord', 'ID', options=list(showTip=TRUE, showLine=F, enableScrollWheel=TRUE, mapType='satellite', useMapTypeControl=TRUE, width=600,height=600))
	return(l)
	}


#### extract bioclimatic variables
extractBiovar=function(coord, climvar, folder, updateProgress = NULL){ ### need also the location of the
	for (i in 1:10){
                if (is.function(updateProgress)) {
                        text = ' ...'
                        updateProgress(detail = text)
                        }
                }
	n=read.csv(coord, sep=',', header=T)
	pos=n[,c(3,2)]
	biofiles=list.files(path=folder, pattern='bil', full.names=T)
	variables=stack(biofiles)
	vals=extract(variables, pos)
	vals=vals[,climvar]
	vals=cbind.data.frame(n[,1], vals)
	colnames(vals)[1]='Genotype'
	vals
	}

getPCA=function(climFile){
	n=read.table(climFile, sep=',', header=T, row.name=1)
	pcanalysis=PCA(scale(n))
	pcanalysis
	}
	
getID=function(climFile){
	n=read.table(climFile, sep=',', header=T, row.name=1)
	rownames(n)
	}

getTESS_struct = function(genofile, coordfile, k, ploidy, rep){
        genotype = as.matrix(read.table(genofile, sep=' '))
        genotype[genotype == 9] = NA
        coordinates = as.matrix(read.table(coordfile, sep = ',', header = T, row.name = 1))
        coordinates = coordinates[, c(2,1)]
	set.seed(101)
	tess3.obj <- tess3(X = genotype, coord = coordinates, K = 1:k, method = "projected.ls", ploidy = ploidy, rep = rep)
	pdf('Cross-Entropy_profile_by_Number_of_K.pdf')
	plot(tess3.obj, pch = 19, col = "blue", xlab = "Number of ancestral populations", ylab = "Cross-validation score")
	dev.off()
	saveRDS(tess3.obj, '../tess3.rds')
	tess3.obj
        }

exportTESS = function(tess_obj, k, coordfile, updateProgress = NULL){
        for (i in 1:10){
                if (is.function(updateProgress)) {
                        text = ' ...'
                        updateProgress(detail = text)
                        }
                }
	coordinates = as.matrix(read.table(coordfile, sep = ',', header = T, row.name = 1))
	coordinates = coordinates[, c(2,1)]
	q.matrix <- qmatrix(tess_obj, K = k)
	Qm=q.matrix
	colnames(Qm) = paste('Q', 1:k, sep='')
	write.table(Qm, paste('Qmatrix_K', k, '.csv', sep=''), sep=',', row.names=F, col.names=T, quote=F)
	my.colors = colorpanel(ncol(q.matrix), 'red', 'green', 'blue')
	my.palette <- CreatePalette(my.colors, 9)
	pdf(paste('TESS3_pop_struct_summary_K', k, '.pdf', sep=''))
	barplot(q.matrix, border = NA, space = 0, main = "Ancestry matrix", xlab = "Individuals", 
		ylab = "Ancestry proportions", col.palette = my.palette) -> bp
	plot(q.matrix, coordinates, method = "map.max", interpol = FieldsKrigModel(10), main = "Ancestry coefficients", 
		xlab = "Longitude", ylab = "Latitude", resolution = c(600,600), cex = .4, col.palette = my.palette)
	dev.off()
	}

createqqmanDF = function(pval_data){
	p_file = list.files('Data_conversion', pattern = 'vcfsnp', full.names=T)
	markerPos = read.table(p_file, sep=' ')
	final_dat = cbind.data.frame(paste(markerPos$V1, markerPos$V2, sep='_'), as.numeric(markerPos$V1), markerPos$V2, pval_data)
	colnames(final_dat) = c('SNP', 'CHR', 'BP', 'P')
	final_dat
	}

exportFst = function(dat, k, padj){
	pdf(paste('Fst_results_plot_K', k, '.pdf', sep=''), height=3.5, width=7)
	hist(dat$P, col='lightblue', main = 'Histogram of P values distribution', xlab = 'P values')
	manhattan(dat, ylim=c(0, -log10(min(dat$P)) + 0.5))
	dev.off()
	dat = cbind.data.frame(dat, p.adjust(dat$P, method=padj))
	colnames(dat)[5]='Padj'
	write.table(dat, paste('Fst_results_table_K', k, '_', padj, '.csv', sep=''), sep=',',  col.name=T, row.name=F, quote=F)
	}
	
createEnv = function(filein){
	allVar = read.table(filein, sep=',', header=T, row.name=1)
	for (i in 1:ncol(allVar)){
		singleVar=allVar[,i]
		envFileName=paste(colnames(allVar)[i], 'env', sep='.')
		write.env(singleVar, envFileName)
		}
	}

run_LFMM = function(envFile, k, rep, miss) {
# this function should return the matrix for manhatthan plot		
	lfmm_file = list.files('../Data_conversion', pattern = 'lfmm')
	file.copy(paste('../Data_conversion/', lfmm_file, sep=''), '.')
	obj.lfmm = lfmm(lfmm_file, envFile, K=k, missing.data=miss, repetitions=rep, project='new')
	zs = z.scores(obj.lfmm, K=k)
	zs.median = apply(zs, MARGIN = 1, median)
	lambda = median(zs.median^2)/qchisq(0.5, df = 1)
	adj.p.values = pchisq(zs.median^2/lambda, df = 1, lower = FALSE)
#	assoc_res = createqqmanDF(adj.p.values)
	p_file = list.files('../Data_conversion', pattern = 'vcfsnp', full.names=T)
        markerPos = read.table(p_file, sep=' ')
        assoc_res = cbind.data.frame(paste(markerPos$V1, markerPos$V2, sep='_'), as.numeric(markerPos$V1), markerPos$V2, adj.p.values)
	colnames(assoc_res) = c('SNP', 'CHR', 'BP', 'P')
	final_res = list(assoc_res, lambda)
	final_res
	}

exportLFMM = function(res_df, k, padj, phenoFile){
	suffix = strsplit(phenoFile, '.env')[[1]]
	pdf(paste('Association_results_K', k, '_', suffix, '.pdf', sep=''), height=3.5, width=7)
	hist(res_df$P, col='lightblue', main = 'Histogram of P values distribution', xlab = 'P values')
	manhattan(res_df, ylim=c(0, -log10(min(res_df$P)) + 0.5))
	dev.off()
	print(padj)
	res_final = cbind.data.frame(res_df, p.adjust(res_df$P, method = padj))
	colnames(res_final)[5] = 'Padj'
	write.table(res_final, paste('Association_results_table_K', k,'_', suffix, '_', padj, '.csv', sep=''), sep=',', col.name=T, row.name=F, quote=F)
	}

get_SNPs_Ranges=function(assoc_file, padj){
	pType = 'P'
	pos_file = list.files('Data_conversion', pattern='vcfsnp', full.names=T)
	pos_file = read.table(pos_file, sep=' ')
	SNPs_info = read.table(assoc_file, sep=',', header=T)
	if (padj){
		pType='Padj'
		}
	SNPs_data = GRanges(seqnames=pos_file$V1, 
		ranges=IRanges(start=pos_file$V2, end = pos_file$V2), 
		scores = SNPs_info[,pType])
	SNPs_data
	}

get_annot = function(genes, SNPs, pval, kb, updateProgress = NULL){
	for (i in 1:10){
		if (is.function(updateProgress)) {
               		text = ' ...'
         	      	updateProgress(detail = text)
               		}
		}
	find_closest = distanceToNearest(genes, SNPs, ignore.strand = T) # get distance
	genes_dist = as.data.frame(find_closest) # DF of distances
	sign_SNPs = as.data.frame(SNPs) # df of SNPs
	sign_SNPs = subset(sign_SNPs, sign_SNPs$scores <= pval) # subset for sign SNPs
	sign_SNPs_dist = subset(genes_dist, genes_dist$subjectHits %in% rownames(sign_SNPs)) # extract sign SNPs distance
	sign_SNPs_good = subset(sign_SNPs_dist, sign_SNPs_dist$distance <= kb*1000) # extract genes by distance
	if (dim(sign_SNPs_good)[1]> 0){
		genes_df= as.data.frame(genes) # df of annot
		genes_infos = genes_df[sign_SNPs_good$queryHits,] # subset for good genes
		SNPs_infos = data.frame() 
		for (i in sign_SNPs_good$subjectHits){
			SNPs_infos = rbind.data.frame(SNPs_infos, subset(sign_SNPs, rownames(sign_SNPs) == i))
			}
		final_data = cbind.data.frame(genes_infos[,c('ID', 'seqnames', 'start', 'end')], SNPs_infos$start, sign_SNPs_good$distance, SNPs_infos$scores)
		colnames(final_data)=c('GeneID', 'Chr', 'Gene_Start', 'Gene_End', 'SNPs', 'distance', 'P')
	#	final_data=unique(final_data)
		final_data
		}
	else {
		final_data = rbind.data.frame(rep('NA', 7))
		colnames(final_data)=c('GeneID', 'Chr', 'Gene_Start', 'Gene_End', 'SNPs', 'distance', 'P')
		final_data
		}
	}
	






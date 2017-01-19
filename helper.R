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
		SNP_pos=read.table(as.character(vcf[1,1]), sep='\t')
		write.table(SNP_pos[,1:2], file='SNPs_positions.txt', sep='\t', row.name=F, col.name=F, quote=F)
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
extractBiovar=function(coord, climvar, folder){ ### need also the location of the
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
	tess3.obj <- tess3(X = genotype, coord = coordinates, K = 1:k, method = "projected.ls", ploidy = ploidy, rep = rep)
	pdf('Cross-Entropy_profile_by_Number_of_K.pdf')
	plot(tess3.obj, pch = 19, col = "blue", xlab = "Number of ancestral populations", ylab = "Cross-validation score")
	dev.off()
	tess3.obj
        }

exportTESS = function(tess_obj, k, coordfile){
	coordinates = as.matrix(read.table(coordfile, sep = ',', header = T, row.name = 1))
	coordinates = coordinates[, c(2,1)]
	q.matrix <- qmatrix(tess_obj, K = k)
	Qm=q.matrix
	colnames(Qm) = paste('Q', 1:k, sep='')
	write.table(Qm, 'Qmatrix.csv', sep=',', row.names=F, col.names=T, quote=F)
	my.colors = rainbow(ncol(q.matrix))
	my.palette <- CreatePalette(my.colors, 9)
	pdf('TESS3_pop_struct_summary.pdf')
	barplot(q.matrix, border = NA, space = 0, main = "Ancestry matrix", xlab = "Individuals", 
		ylab = "Ancestry proportions", col.palette = my.palette) -> bp
	plot(q.matrix, coordinates, method = "map.max", interpol = FieldsKrigModel(10), main = "Ancestry coefficients", 
		xlab = "Longitude", ylab = "Latitude", resolution = c(600,600), cex = .4, col.palette = my.palette)
	dev.off()
	}


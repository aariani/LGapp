#### helper functions

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



	



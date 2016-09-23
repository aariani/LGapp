#### helper functions

plotMap=function(coord){
	n=read.csv(coord, sep=',', header=T)
### you have to modify the dataset
	coord=paste(n[,2], n[,3], sep=':')
	f=cbind.data.frame(n, coord)
	colnames(f)[1]= 'ID'
	l=gvisMap(f, 'coord', 'ID', options=list(showTip=TRUE, showLine=F, enableScrollWheel=TRUE, mapType='satellite', useMapTypeControl=TRUE, width=600,height=600))
	return(l)
	}



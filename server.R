# server.R
library(shinythemes)
library(LEA)
library(dismo)
#library(rasterVis)
#library(maptools)
#library(rgeos)
library(googleVis)
library(ChemometricsWithR)

#source('plotData.R') ### script with all the plot data informations look at https://pakillo.github.io/R-GIS-tutorial/#mapping
source('helper.R')


shinyServer(function(input, output){
	shinyFileChoose(input, 'vcf',  root=c(home=path.expand('~')))
	shinyDirChoose(input, 'outdir', root=c(home=path.expand('~')))
	ProjFolder=reactive({parseDirPath(c(home=path.expand('~')), input$outdir)}) # THIS IS THE MAIN PROJECT FOLDER
#########################################
#### Conversion Tab START ###############
	vcf=reactive({parseFilePaths(c(home=path.expand('~')), input$vcf)}) # get VCF file
	prefix=reactive({strsplit(as.character(vcf()[1,1]), 'vcf')[[1]][1]}) # get file prefix
# Copy data in local folder
	getData=reactive({
		file.copy(as.character(vcf()[1,4]), '.')
		vcf2lfmm(as.character(vcf()[1,1]))
		SNP_pos=read.table(as.character(vcf()[1,1]), sep='\t')
		write.table(SNP_pos[,1:2], file=paste(prefix(), 'positions.txt', sep=''), sep='\t', row.name=F, col.name=F, quote=F)
		})
# Transfer and clean files
	cleanData=reactive({
		allfiles=list.files('.', pattern=prefix())
		for (i in allfiles)
			file.copy(i, ProjFolder())		
			file.remove(allfiles)
		})
## Output data	
	output$text1= renderPrint({
		if (!is.null(input$vcf))
			getData()
		})
	output$text2=renderText({
		if (is.null(input$outdir))
			return(NULL)	
		tryCatch({cleanData()}, error = function(e){})
		paste('Selected', ProjFolder(), 'as Project Folder')
		})

########################################
#### Climatic Data Start ###############
#############################
#### Download ###############
## you can use renderGvis for rendering an html data type with the coordinates
## but you need to create a function in the main body of the server.R
	shinyDirChoose(input, 'climdir', root=c(home=path.expand('~'))) ## set up directory
	getGvismap=function(){
		coordFile=input$coord
		plotMap(coordFile$datapath)
		}
	output$bioclim=renderGvis({
		coordFile=input$coord
		if (is.null(coordFile)) return(NULL)
		getGvismap()
		})	
## get bioclimatic data
	getBioclimData=reactive({
		coordFile=input$coord
		climFolder=parseDirPath(c(home=path.expand('~')), input$climdir)
		extractBiovar(coordFile$datapath, input$climvar, climFolder) ### need also the location of the folder
		})
## Output table
	output$climTable=renderDataTable({
		if (is.null(input$climdir))
			return(NULL)
		getBioclimData()
#		if (input$download_clim == T){
#			write.table(getBioclimData(), paste(ProjFolder(), '/climatic-data', Sys.Date(),'.csv', sep=''), sep=',', col.names=T, quote=F, row.names=F)
#			}
		})
## Download buttons
#	output$download_clim=downloadHandler(
#		filename=function(){
#			paste('climatic_data-', Sys.Date(), '.csv', sep='')
#			},
#		content=function(file){
#			write.table(getBioclimData(), file,  sep=',', col.names=T, quote=F, row.names=F)
#			}
#		)

##############################
#### PCA #####################
	getPCAdata=reactive({
		climFile=input$climDat
		getPCA(climFile$datapath)
		})

	genoID=reactive({
		climFile=input$climDat
		getID(climFile$datapath)	
		})
	
	output$PCAplot=renderPlot({
		if (!is.null(input$climDat))
			biplot(getPCAdata(), show.names='loadings',loading.col='black')
		})
	output$PCAsummary=renderPrint({
		if (!is.null(input$climDat))
			summary(getPCAdata())
		})

	output$pc_coord=downloadHandler(
		filename=function(){
			paste('Coordinates_', input$n_PCs, '_PC_', Sys.Date(), '.csv', sep='')
			},
		content=function(file){
			PCcord=getPCAdata()$scores[,1:as.numeric(input$n_PCs)]
			rownames(PCcord)=genoID()
			write.table(PCcord, file, sep=',', col.names=T, quote=F, row.name=T)
			}
		)
	output$pc_load=downloadHandler(
	#	if (!is.null(input$n_PCs))	
		filename=function(){
			paste('Loadings_', input$n_PCs, '_PC_', Sys.Date(), '.csv', sep='')
			},
		content=function(file){
			write.table(getPCAdata()$loadings[,1:as.numeric(input$n_PCs)], file, sep=',', col.names=T, quote=F, row.names=T)
			}
		)
################################
### Population Structure #######
	shinyFileChoose(input, 'geno', root=c(home=path.expand('~')))
### get geno file
	geno=reactive({parseFilePaths(root=c(home=path.expand('~')), input$geno)})
### do sNMF analysis
	sNMF_analysis=reactive({
	if (!is.null(geno()))
		file.copy(as.character(geno()[1,4]), '.')
		snmf(as.character(geno()[1,1]), K=input$K_range[1]:input$K_range[2], entropy=T, rep=input$rep, project='new')
### maybe also remove the geno file. Add a name to the snmf. Remoive the file and then call the project
		})
### get best run
	best_run=reactive({
		ce=cross.entropy(sNMF_analysis(), K=input$n_K)
		best=which.min(ce)
		best
		})
### plot CE for choosing best K
	output$CEplot=renderPlot({
		if (!is.null(input$geno))
			plot(sNMF_analysis(), lwd = 5, col = "red", pch=1)
			})

})






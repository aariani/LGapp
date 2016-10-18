# server.R
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
#########################################
#### Conversion Tab START ###############
## you need to use a reactive expression for getting the file name with the parseFilePaths
## function and do ti outside of the output$text1 thing
	vcf=reactive({parseFilePaths(c(home=path.expand('~')), input$vcf)})
## get the prefix of the file	
	prefix=reactive({strsplit(as.character(vcf()[1,1]), 'vcf')[[1]][1]})
## Get output folder
	outfold=reactive({parseDirPath(c(home=path.expand('~')), input$outdir)})
## Get data file in a local folder
	getData=reactive({
## Copy the file in the local folder
		file.copy(as.character(vcf()[1,4]), '.')
		vcf2lfmm(as.character(vcf()[1,1]))
		})
## Transfer and clean files
	cleanData=reactive({
		allfiles=list.files('.', pattern=prefix())
		for (i in allfiles)
			file.copy(i, outfold())		
			file.remove(allfiles)
		})
## Output data	
	output$text1= renderPrint({
		if (!is.null(input$vcf))
			##paste('Conversion Summary')
			getData()	
		})
	output$text2=renderText({
		if (is.null(input$outdir))
			return(NULL)	
		cleanData()
		paste('Converting Files into', outfold(), 'folder')
		})
#### Conversion Tab chunk END ##########
########################################

########################################
#### Climatic Data Start ###############

#############################
#### Download ###############
#############################

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
#		if (!is.null(input$climdir))
#			return(NULL)
		coordFile=input$coord
		climFolder=parseDirPath(c(home=path.expand('~')), input$climdir)
		extractBiovar(coordFile$datapath, input$climvar, climFolder) ### need also the location of the folder
		})
## Output table
	output$climTable=renderDataTable({
		if (is.null(input$climdir))
			return(NULL)
		getBioclimData()
		})
## Download buttons
	output$download_clim=downloadHandler(
		filename=function(){
			paste('climatic_data-', Sys.Date(), '.csv', sep='')
			},
		content=function(file){
			write.table(getBioclimData(), file,  sep=',', col.names=T, quote=F, row.names=F)
			}
		)

##############################
#### PCA #####################
##############################
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
})






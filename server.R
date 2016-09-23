# server.R
library(LEA)
library(dismo)
library(rasterVis)
library(maptools)
library(rgeos)
library(googleVis)
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
		extractBiovar(coordFile$datapath, input$climdata)
		})
## Output table

## Download buttons


})






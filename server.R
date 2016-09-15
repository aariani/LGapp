# server.R
library(LEA)
library(shinyFiles)


shinyServer(function(input, output){
	shinyFileChoose(input, 'vcf',  root=c(home='~'))
	shinyDirChoose(input, 'outdir', root=c(home='~'))

#### Conversion Tab chunk START
## you need to use a reactive expression for getting the file name with the parseFilePaths
## function and do ti outside of the output$text1 thing
	vcf=reactive({parseFilePaths(c(home='~'), input$vcf)})
## get the prefix of the file	
	prefix=reactive({strsplit(as.character(vcf()[1,1]), 'vcf')[[1]][1]})
## Get output folder
	outfold=reactive({parseDirPath(c(home='~'), input$outdir)})
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
#### Conversion Tab chunk END
	})

# server.R
library(LEA)
library(shinyFiles)


shinyServer(function(input, output){

	shinyFileChoose(input, 'vcf',  root=c(home='~'))
	shinyDirChoose(input, 'outdir', root=c(home='~'))

#### you need to use a reactive expression for getting the file name with the parseFilePaths
### function and do ti outside of the output$text1 thing
	vcf=reactive({parseFilePaths(c(home='~'), input$vcf)})

### get the prefix of the file	
	prefix=reactive({strsplit(as.character(vcf()[1,1]), 'vcf')[[1]][1]})

## Get output folder
	outfold=reactive({parseDirPath(c(home='~'), input$outdir)})

### Get data file
	getData=reactive({
		if (is.null(input$vcf))
			return(NULL)
#### Copy the file in the local folder
		file.copy(as.character(vcf()[1,4]), '.')
		vcf2lfmm(as.character(vcf()[1,1]))
		})

### clean data
	cleanData=reactive({
		if (is.null(input$outdir))
			return(NULL)	
		allfiles=list.files('.', pattern=prefix())
		for (i in allfiles)
			file.copy(i, outfold())		
		
		file.remove(allfiles)
		})

	
	output$text1= renderPrint({
		if (is.null(input$vcf))
			return(NULL)
		getData()

		cleanData()
		paste('Converting', vcf()[1,1], 'into', outfold())
### output data
		})
	})

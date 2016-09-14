# server.R
library(LEA)
library(shinyFiles)


shinyServer(function(input, output){

	shinyFileChoose(input, 'vcf',  root=c(home='~'))
	shinyDirChoose(input, 'outdir', root=c(home='~'))

#### you need to use a reactive expression for getting the file name with the parseFilePaths
### function and do ti outside of the output$text1 thing
	vcf=reactive({parseFilePaths(c(home='~'), input$vcf)})

	getData=reactive({
		if (is.null(input$vcf))
			return(NULL)
		### I dunno why I cannot make this function works
		vcf2lfmm(as.character(vcf()[1,4]))
		})

	output$text1= renderText({
		if (is.null(input$vcf))
			return(NULL)
##### You need tio get the parseFilePath function for getting the path of the file		
#		vcf=parseFilePaths(c(home='~'), input$vcf)
#		inputFilePath=input$vcf
#		inputFilePath=inputFilePath$'0'
#		inputFilePath=paste(inputFilePath[2:length(inputFilePath)], collapse='/')
#		vcf2lfmm(as.character(s[1,4]))
		getData()
		paste('Converting VCF file named',  vcf()[1,4])
### output data
		})
	})

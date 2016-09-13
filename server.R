# server.R
library(LEA)
library(shinyFiles)


shinyServer(function(input, output){

	shinyFileChoose(input, 'vcf',  root=c(home='~'))
	shinyDirChoose(input, 'outdir', root=c(home='~'))
#	outform=reactive({
#		switch(input$outformat,
#		'lfmm'='Landscape Genomic analysis (LFMM format)',
#		'geno'='Population clustering analysis (geno format)'
#		)
#	})
	output$text1= renderText({
		formdat=switch(input$outformat,
			'lfmm'='Landscape Genomic analysis (LFMM format)',
			'geno'='Population clustering analysis (geno format)'
			)
		print(paste('Converting VCF file for', formdat, 'in', getwd(), 'AND', input$vcf))
### output data
		})

	

	})

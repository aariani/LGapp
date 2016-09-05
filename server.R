# server.R
library(LEA)
library(shinyFiles)


shinyServer(function(input, output){

	shinyFileChoose(input, 'vcf',  root=c(wd=getwd()))
	shinyDirChoose(input, 'outdir', root=c(wd=getwd()))
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
		paste('Converting VCF file for', formdat, 'in', input$outdir)
### output data
		})

	

	})

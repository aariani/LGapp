# server.R
library(LEA)

shinyServer(function(input, output){
	
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
		paste('Converting VCF file for', formdat, 'in', getwd())
### output data
		})

	

	})

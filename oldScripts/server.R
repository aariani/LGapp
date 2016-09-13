# server.R
library(LEA)


shinyServer(function(input, output){

	output$text1= renderText({
		formdat=switch(input$outformat,
			'lfmm'='Landscape Genomic analysis (LFMM format)',
			'geno'='Population clustering analysis (geno format)'
			)
		inputVCF=input$vcf
		if (is.null(inputVCF))
			return(NULL)
		file.rename(inputVCF[1,4], paste(inputVCF[1,4], 'vcf', sep='.'))
		vcf2lfmm(paste(inputVCF[1,4], 'vcf', sep='.'), output.file='conversion')
		print(paste('Converting VCF file', inputVCF[1,4], 'in', formdat, 'in', class(inputVCF[1,4])))
		})
	

})

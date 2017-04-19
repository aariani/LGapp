# server.R
library(LEA)
library(dismo)
library(googleVis)
library(ChemometricsWithR)
library(tess3r)
library(maps)
library(qqman)
library(gplots)
library(GenomicRanges)
library(rtracklayer)

source('helper.R')

shinyServer(function(input, output, session){
	shinyFileChoose(input, 'vcf',  root=c(home=path.expand('~')))
	shinyDirChoose(input, 'outdir', root=c(home=path.expand('~')))
	ProjFolder=reactive({parseDirPath(c(home=path.expand('~')), input$outdir)}) # THIS IS THE MAIN PROJECT FOLDER
#########################################
#### Conversion Tab START ###############
	vcf=reactive({parseFilePaths(c(home=path.expand('~')), input$vcf)}) # get VCF file
# Convert Data
	getData=reactive({
		convertData(ProjFolder(), vcf())
		})
# Output data	
	output$text1= renderPrint({
		if (!is.null(input$vcf))
			getData()
		})
	output$text2=renderText({
		if (is.null(input$outdir))
			return(NULL)	
		paste('Selected', ProjFolder(), 'as Project Folder')
		})

################################
#### Climatic Data Start #######
################################
#### Download ##################
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
                progress <- shiny::Progress$new()
                progress$set(message = "Extracting data", value = 0)
                on.exit(progress$close())
                updateProgress <- function(value = NULL, detail = NULL) {
                        if (is.null(value)) {
                                value <- progress$getValue()
                                value <- value + (progress$getMax() - value) / 5
                                }
                progress$set(value = value, detail = detail)
                	}
		coordFile=input$coord
		climFolder=parseDirPath(c(home=path.expand('~')), input$climdir)
		extractBiovar(coordFile$datapath, input$climvar, climFolder, updateProgress) ### need also the location of the folder
		})
## Output table
	output$climTable=renderDataTable({
		if (is.null(input$climdir))
			return(NULL)
		getBioclimData()
		})

	observeEvent(input$download_clim, {
			setwd(ProjFolder())
			dir.create('Bioclimatic_data')
			write.table(getBioclimData(), paste('Bioclimatic_data/Climatic_data','.csv', sep=''), 
			sep=',', col.names=T, quote=F, row.names=F)
			session$sendCustomMessage(type = 'testmessage',
      				message = 'Data has been downloaded in the Bioclimatic_data folder')
			})	
##############################
#### PCA #####################
	getPCAdata=reactive({
		climFile=input$climDat
		getPCA(climFile$datapath)
		})
	output$PCAplot=renderPlot({
		if (!is.null(input$climDat))
			biplot(getPCAdata(), show.names='loadings',loading.col='blue')
		})
	output$PCAsummary=renderPrint({
		if (!is.null(input$climDat))
			summary(getPCAdata())
		})
	observeEvent(input$download_PCA, {
		if (input$n_PCs > 0)
			write.table(getPCAdata()$scores[,1:input$n_PCs], paste('Bioclimatic_data/Climatic_PCA_coordinates_', 
				input$n_PCs, '_PCs', '.csv', sep=''), sep=',', col.names=T, quote=F, row.name=T)
			write.table(getPCAdata()$loadings[,1:input$n_PCs], paste('Bioclimatic_data/Climatic_PCA_loadings_', 
				input$n_PCs, '_PCs', '.csv', sep=''), sep=',', col.names=T, quote=F, row.names=T)
			pdf('Bioclimatic_data/PCA_biplot.pdf')
			biplot(getPCAdata(), show.names='loadings',loading.col='blue')
			dev.off()
			session$sendCustomMessage(type = 'testmessage',
                                message = 'Data has been downloaded in the Bioclimatic_data folder')
			})						
		
################################
### Population Structure #######
	TESS_analysis = eventReactive(input$analyze_geno, {
		setwd(ProjFolder())
		dir.create('Population_Structure')
		setwd('Population_Structure')
		fname = list.files('../Data_conversion', pattern = 'lfmm', full.names = T)
		coordFile2=input$coord2
		getTESS_struct(fname, coordFile2$datapath, input$K, input$ploidy, input$rep)
		})

### plot CE for choosing best K
	output$CEplot=renderPlot({
		if (input$analyze_geno > 0)
			plot(TESS_analysis(), pch = 19, col = "blue", xlab = "Number of ancestral populations", ylab = "Cross-validation score")
		})

	observeEvent(input$res, {
		if (input$n_K > 0)
			progress <- shiny::Progress$new()
			progress$set(message = "Exporting data", value = 0)
                	on.exit(progress$close())
                	updateProgress <- function(value = NULL, detail = NULL) {
                        if (is.null(value)) {
                                value <- progress$getValue()
                                value <- value + (progress$getMax() - value) / 5
                                }
                	progress$set(value = value, detail = detail)
                        }
			coordFile2=input$coord2
			exportTESS(TESS_analysis(), input$n_K, coordFile2$datapath, updateProgress) 
			session$sendCustomMessage(type = 'testmessage',
                                message = 'Data has been downloaded in the Population_Structure folder')
			})

##################################
#### Selection scan analysis######
##################################
#### Population differentiation ##
	Fst = eventReactive(input$fst_analysis, {
		setwd(ProjFolder())
		tess3struct = readRDS('tess3.rds')
		p.values = pvalue(tess3struct, K = input$Kfst)
		createqqmanDF(p.values)
		})

	output$fst_pvals = renderPlot({
		hist(Fst()[,4], col='lightblue', main = 'Histogram of P values distribution', xlab = 'P values')
#		plot(1,1)
		})

	output$fst_manhattan = renderPlot({
			manhattan(Fst(), ylim=c(0, -log10(min(Fst()[,4])) + 0.5))
		})

	observeEvent(input$fst_res, {
		dir.create('Population_differentiation')
		setwd('Population_differentiation')
		exportFst(Fst(), input$Kfst, input$padj)
		session$sendCustomMessage(type = 'testmessage',
                                message = 'Data has been downloaded in the Population_differentiation folder')
		})

###################################
#### Association analysis #########
	observeEvent(input$writeEnv, {
		setwd(ProjFolder())
		phenofile = input$pheno_all
		dir.create('Association_analysis')
		setwd('Association_analysis')
		createEnv(phenofile$datapath)
		session$sendCustomMessage(type = 'testmessage',
                                message = 'ENV file has been downloaded in the Association_analysis folder')
		})
### LFMM analysis
	LFMM_analysis = eventReactive(input$lfmm_analysis, {
		setwd(paste(ProjFolder(), 'Association_analysis', sep='/'))
		envFile = input$pheno
		run_LFMM(envFile$name, input$LF, input$rep_lfmm, input$miss)
		})

	output$lfmm_pvals=renderPlot({
		if (input$lfmm_analysis > 0)
			hist(LFMM_analysis()[[1]][,4], col='lightblue', main = 'Histogram of P values distribution', xlab = 'P values')
		})

	output$lfmm_manhatthan=renderPlot({
		manhattan(LFMM_analysis()[[1]], ylim=c(0, -log10(min(LFMM_analysis()[[1]][,4])) + 0.5))
		})

	observeEvent(input$lfmm_res, {
		setwd(paste(ProjFolder(), 'Association_analysis', sep='/'))
		envFile = input$pheno
		exportLFMM(LFMM_analysis()[[1]], input$LF, input$padj2, envFile$name)
		session$sendCustomMessage(type = 'testmessage',
                                message = 'Data has been downloaded in the Association_analysis folder')		
		})


###################################
#### SNPs Annotation ##############
	shinyFileChoose(input, 'annot',  root=c(home=path.expand('~')))
	GFF3_genes=reactive({
		f=parseFilePaths(c(home=path.expand('~')), input$annot)
		f=import.gff3(as.character(f$datapath))
		genes=subset(f, f$type == 'gene') 
		genes
		})
	
	output$gff3infos=renderText({
		if (is.null(input$annot))
                        return(NULL)
		paste('GFF3 file loaded')
		})

	SNPs=reactive({
		setwd(ProjFolder())
		sign_file = input$assoc_res
		snps = get_SNPs_Ranges(sign_file$datapath, input$padj_type)
		})
	
	output$annot_res = renderDataTable({
		if(is.null(input$assoc_res))
                        return(NULL)
		progress <- shiny::Progress$new()
                progress$set(message = "Computing data", value = 0)
                on.exit(progress$close())
                updateProgress <- function(value = NULL, detail = NULL) {
                        if (is.null(value)) {
                                value <- progress$getValue()
                                value <- value + (progress$getMax() - value) / 5
                                }
                	progress$set(value = value, detail = detail)
                	}
                get_annot(GFF3_genes(), SNPs(), input$pval_max, input$dist_kb, updateProgress)
		})

	observeEvent(input$annot_res, {
		setwd(ProjFolder())
		dir.create('SNPs_annotation')
		snps_file = input$assoc_res
		annot_res_index = strsplit(snps_file$name, '.csv')[[1]]
		progress <- shiny::Progress$new()
                progress$set(message = "Extract data", value = 0)
                on.exit(progress$close())
                updateProgress <- function(value = NULL, detail = NULL) {
                        if (is.null(value)) {
                                value <- progress$getValue()
                                value <- value + (progress$getMax() - value) / 5
                                }
                        progress$set(value = value, detail = detail)
                        }
		write.table(get_annot(GFF3_genes(), SNPs(), input$pval_max, input$dist_kb, updateProgress),
			paste('SNPs_annotation/Annotation_', annot_res_index, '.csv', sep=''),
			sep=',', row.names=F, col.names=T, quote=F)
		session$sendCustomMessage(type = 'testmessage',
                                message = 'Data has been downloaded in the SNPs_annotation folder')
		})		
})

# server.R
library(LEA)
library(dismo)
library(googleVis)
library(ChemometricsWithR)
library(tess3r)
library(maps)
library(qqman)

source('helper.R')

shinyServer(function(input, output){
	shinyFileChoose(input, 'vcf',  root=c(home=path.expand('~')))
	shinyDirChoose(input, 'outdir', root=c(home=path.expand('~')))
	ProjFolder=reactive({parseDirPath(c(home=path.expand('~')), input$outdir)}) # THIS IS THE MAIN PROJECT FOLDER
#	setwd(ProjFolder)
#	s=list.files(pattern='RData')
#	if (length(s)==1)
#		load(s)
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
	observeEvent(input$download_clim, {
			setwd(ProjFolder())
			dir.create('Bioclimatic_data')
			write.table(getBioclimData(), paste('Bioclimatic_data/Climatic_data','.csv', sep=''), 
			sep=',', col.names=T, quote=F, row.names=F)
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
			coordFile2=input$coord2
			#save.image('../log_infos.RData')
			exportTESS(TESS_analysis(), input$n_K, coordFile2$datapath) 
			})

##################################
#### Selection scan analysis######
##################################
#### Population differentiation ##
	Fst = eventReactive(input$fst_analysis, {
		setwd(ProjFolder())
		tess3struct = readRDS('tess3.rds')
		dir.create('Population_differentiation')
		setwd('Population_differentiation')
		p.values = pvalue(tess3struct, K = input$Kfst)
		p.values
		})

	output$fst_pvals = renderPlot({
		hist(Fst(), col='lightblue', main = 'Histogram of P values distribution', xlab = 'P values')
		})
	
	output$fst_manhattan = renderPlot({
		if (input$fst_analysis > 0)	
			dat = createqqmanDF(Fst())
			manhattan(dat)
		})

})






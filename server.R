# server.R
library(LEA)
library(dismo)
library(googleVis)
library(ChemometricsWithR)

source('helper.R')

shinyServer(function(input, output){
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
## Output data	
	output$text1= renderPrint({
		if (!is.null(input$vcf))
			getData()
		})
	output$text2=renderText({
		if (is.null(input$outdir))
			return(NULL)	
		paste('Selected', ProjFolder(), 'as Project Folder')
		})

########################################
#### Climatic Data Start ###############
#############################
#### Download ###############
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
			biplot(getPCAdata(), show.names='loadings',loading.col='black')
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
			})						
		
################################
### Population Structure #######
	sNMF_analysis = eventReactive(input$analyze_geno, {
		setwd(ProjFolder())
		dir.create('Population_Structure')
		setwd('Population_Structure')
		fname = list.files('../Data_conversion', pattern = 'geno', full.names = T)
		file.copy(fname, '.')
		snmf(list.files('.', pattern ='geno'), K=input$K_range[1]:input$K_range[2], entropy=T, rep=input$rep, project='new')
		})
### get best run
	best_run=reactive({
		ce=cross.entropy(sNMF_analysis(), K=input$n_K)
		best=which.min(ce)
		best
		})
### plot CE for choosing best K
	output$CEplot=renderPlot({
		if (input$analyze_geno > 0)
			plot(sNMF_analysis(), lwd = 5, col = "red", pch=1)
		})
	observeEvent(input$K_matrix, {
		if (input$n_K > 0) ## i think in this step you will need also the G matrix to download for subsequent steps
			Qm = Q(sNMF_analysis(), input$n_K, best_run())
			colnames(Qm) = paste('Q', 1:input$n_K, sep='') 
			write.table(Qm, 'Qmatrix.csv', sep=',', row.names=F, col.names=T, quote=F)
			})
})






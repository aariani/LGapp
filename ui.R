## ui.R
mycss <- "
#plot-container {
  position: relative;
}
#loading-spinner {
  position: absolute;
  left: 50%;
  top: 50%;
  z-index: -1;
  margin-top: -33px;  /* half of the spinner's height */
  margin-left: -33px; /* half of the spinner's width */
}
#plot.recalculating {
  z-index: -2;
}
"

shinyUI(navbarPage('LGapp', theme=shinytheme('cosmo'),
######################################
#### Data conversion #################
	tabPanel('Project Home', icon=icon("home"),
		sidebarPanel(
			h2('Start your project'),
			shinyDirButton('outdir', 'Choose Project Folder', 'Please select output folder'),
			helpText('Select output directory where exporting the files produced during the analysis'),
			h3(textOutput('text2')),
			tags$hr(),
			shinyFilesButton('vcf', 'Choose VCF file', 'Please select a file', multiple=F),
			helpText('Convert VCF files for downstream analysis'),
			verbatimTextOutput('text1')
			),
		mainPanel(
			h1(strong('Welcome to LGapp!!!'), align='center'), 
			h2('a Shiny web-app for landscape genomics analysis within R', align='center'),
			h3('For starting your landscape genomics project please select a ', strong('Project Folder'), ' and your SNPs file (in VCF format)', align='center')
			)
		),

#######################################
#### Climatic Data ###################
	navbarMenu('Bioclimatic Data', icon=icon('globe'),
		tabPanel('Download', icon=icon('download'),
		sidebarPanel(
			h2('Download bioclimatic data'),
			tags$hr(),
			checkboxGroupInput('climvar',
				label='Select Climatic Data', 
				choices=c('BIO1 = Annual Mean Temperature'='bio_1',
					'BIO2 = Mean Diurnal Range'='bio_2',
					'BIO3  =  Isothermality (BIO2/BIO7) (* 100)'='bio_3',
					'BIO4  =  Temperature Seasonality (standard deviation *100)'='bio_4',
					'BIO5  =  Max Temperature of Warmest Month'='bio_5',
					'BIO6  =  Min Temperature of Coldest Month'='bio_6',
					'BIO7  =  Temperature Annual Range (BIO5-BIO6)'='bio_7',
					'BIO8  =  Mean Temperature of Wettest Quarter'='bio_8', 
					'BIO9  =  Mean Temperature of Driest Quarter'='bio_9',
					'BIO10  =  Mean Temperature of Warmest Quarter'='bio_10',
					'BIO11  =  Mean Temperature of Coldest Quarter'='bio_11',
					'BIO12  =  Annual Precipitation'='bio_12',
					'BIO13  =  Precipitation of Wettest Month'='bio_13',
					'BIO14  =  Precipitation of Driest Month'='bio_14',
					'BIO15  =  Precipitation Seasonality (Coefficient of Variation)'='bio_15',
					'BIO16  =  Precipitation of Wettest Quarter'='bio_16',
					'BIO17  =  Precipitation of Driest Quarter'='bio_17',
					'BIO18  =  Precipitation of Warmest Quarter'='bio_18',
					'BIO19  =  Precipitation of Coldest Quarter'='bio_19'),
				selected=c('bio_1', 'bio_2', 'bio_3', 'bio_4', 'bio_5', 'bio_6', 'bio_7',
					'bio_8', 'bio_9', 'bio_10', 'bio_11', 'bio_12', 'bio_13',
					'bio_14', 'bio_15', 'bio_16', 'bio_17', 'bio_18', 'bio_19'),
				width='100%'
				),
			tags$hr(),
			fileInput('coord', 'Select file with coordinates', accept='.csv'),
			helpText('The conversion requires a comma separated csv file with 3 columns. The First colum should contains the Genotype ID, the second the Latitude and the third the Longitude'),
			tags$hr(),
			shinyDirButton('climdir', 'Select Folder with bioclimatic data (bil files)', 'Please select folder with climatic data'),
			helpText('Download climatic data to your computer from the ', a('WorldClim database', href='http://www.worldclim.org/current')),
			tags$hr(),
			actionButton('download_clim', label = 'Download Bioclimatic Variables')),
			mainPanel(
				htmlOutput('bioclim'),
				dataTableOutput('climTable')
				)
			),
		tabPanel('PCA', icon=icon('object-group'),
			sidebarPanel(
				h2('PCA on bioclimatic variables'),
				tags$hr(),
				p('PCA analysis for identifying the most representative bioclimatic variables in your dataset (to be used for Landscape Genomic analysis).
				You can also use the PC data as phenotypes'),
				tags$hr(),
				fileInput('climDat', 'Select file with climatic data', accept='.csv'),
				helpText('You can use directly the output of the Bioclimatic Data download step.
					Otherwise you can upload a comma separated csv file having the genotype ID in the first column'),
				tags$hr(),
				h3('Download PCA coordinates and loadings'),
				numericInput('n_PCs', label=p('Select number of PCs to download'), value = 0, min = 0),
				tags$hr(),
				actionButton('download_PCA', 'Download PCA data')
				),
### output biplot in an image, and table with loadings
			mainPanel(
				plotOutput('PCAplot', width='100%',height='800px'),
				verbatimTextOutput('PCAsummary')
				)
			)),
#################################
#### Population structure #######
	tabPanel('Population structure',icon=icon('bar-chart'),
		sidebarPanel(
			h2('Population structure analysis with TESS3'),
			tags$hr(),
			fileInput('coord2', 'Select file with coordinates', accept='.csv'),
			tags$hr(),
			numericInput('ploidy', label=p('Select ploidy'), value = 2, min = 1),
			tags$hr(),
			sliderInput('K', label='Select K range', min=1, max=20, value=6),
			tags$hr(),
			sliderInput('rep', label='Select number of repetition for each K', min=1, max=100, value=10),
			tags$hr(),
			actionButton('analyze_geno', 'Start Population Clustering analysis'),#, 'Please select a file', multiple=F),
			helpText('This step will analyze the genotype file created during the Data Conversion Step'),
			tags$hr(),
			numericInput('n_K', label='Select best number of K', value = 0, min = 0),
			tags$hr(),
			actionButton('res', 'Download Results')
			),
		mainPanel(
			tags$head(tags$style(HTML(mycss))),
			div(id = "plot-container",
			tags$img(src = "spinner.gif",
			id = "loading-spinner"),
			plotOutput('CEplot', width='100%',height='800px'))
			)
		),	
##############################
#### Selection scan ##########
	navbarMenu('Selection Scan', icon=icon('cogs'),
		tabPanel('Population differentiation', icon=icon('arrows-alt'),
			sidebarPanel(
				h2('Population differentiation analysis based on Fst statistic'),
				tags$hr(),
				sliderInput('Kfst', label='Select number of K for the analysis', value=2, min=0, max=20),
				tags$hr(),
				actionButton('fst_analysis', 'Start Analysis'),
				tags$hr(),
				selectInput('padj', label='Select P values correction method', c('None' = 'none', 
					'Bonferroni'='bonferroni', 'Holm'='holm', 'Hochberg'='hochberg', 'Hommel'='hommel', 
					'Benjamini & Hochberg'='BH', 'Benjamini & Yekutieli'='BY', 'FDR'='fdr')),
				tags$hr(),
				actionButton('fst_res', 'Export Results')
				),
			mainPanel(
				plotOutput('fst_pvals'),
				plotOutput('fst_manhattan')
				)
		),
		tabPanel('Association analysis', icon=icon('filter'),
			sidebarPanel(
				h2('Association analysis with LFMM'),
				tags$hr(),
				h3('Convert Phenotype file'),
				fileInput('pheno_all', 'Select phenotype file', accept='.csv'),
				actionButton('writeEnv', 'Convert Files'),
				helpText('This step convert the Phenotype file in single variable files'),
				tags$hr(),
				fileInput('pheno', 'Select Variable file', accept = '.env'),
				helpText('Select one the variable file converted in the previous step'),
				#tags$hr(),
				sliderInput('LF', label='Select number of K for the association analysis', value=2, min=0, max=20),
				sliderInput('rep_lfmm', label='Select number of repetition for each K', min=1, max=100, value=10),
				checkboxInput('miss', 'Missing data', value = F),
				tags$hr(),
				actionButton('lfmm_analysis', 'Run association analysis')
				
		)	
		)
	)
### after this keep the parenthesis
	)
)

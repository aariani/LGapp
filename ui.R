## ui.R

shinyUI(navbarPage('LGapp',
######################################
#### Data conversion #################
	tabPanel('Data Conversion',
		sidebarPanel(
			h2('Convert VCF files'),
			tags$hr(),
			helpText('Convert VCF files for downstream analysis'),
			shinyFilesButton('vcf', 'Choose VCF file', 'Please select a file', multiple=F),
			tags$hr(),
			helpText('Select output directory where export the converted files'),
			shinyDirButton('outdir', 'Choose output folder', 'Please select output folder')
			),
		mainPanel(
			h2('Conversion Summary'),
			verbatimTextOutput('text1'),
			h3(textOutput('text2'))
			)
		),

#######################################
#### Climatic Data ###################
	navbarMenu('Bioclimatic Data',
		tabPanel('Download',
		sidebarPanel(
			h2('Download bioclimatic data'),
			tags$hr(),
			checkboxGroupInput('climdata',
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
			helpText('The conversion requires a comma separated csv file with 3 columns. The First colum should contains the Genotype ID, the second the Latitude and the third the longitude'),
			tags$hr(),
			helpText('Data downloaded from the ', a('WorldClim database', href='http://www.worldclim.org/')),
			helpText('This step will plot the positions of the coordinates in the input file on a map, 
				and will allow to download the bioclimatic variables for each of the coordinates in the input file')),
			mainPanel(
				htmlOutput('bioclim'),
				downloadButton('downloadData', 'Download')
				)
			),
		tabPanel('PCA',
			sidebarPanel(
				h2('PCA on bioclimatic variables')
				)
			)
		)
	)
)

## ui.R

shinyUI(navbarPage('LGapp',
	tabPanel('Data Conversion',
		sidebarPanel(
			h2('Convert VCF files'),
			tags$hr(),
			helpText('Convert VCF files to geno or LFMM format for further analysis'),
			helpText('This will create a geno and an lfmm file in the directory selected as output'),
			shinyFilesButton('vcf', 'Choose VCF file', 'Please select a file', multiple=F),
	#		selectInput('outformat',
	#			'Select Output Format',
	#			choices=c('lfmm', 'geno'),
	#			),
			tags$hr(),
			helpText('Select output directory'),
			shinyDirButton('outdir', 'Choose output folder', 'Please select output folder')
			),
		mainPanel(textOutput('text1'))
		),
	navbarMenu('Bioclimatic Data',
		tabPanel('Download',
		sidebarPanel(
			h2('Download bioclimatic data'),
			tags$hr(),
			selectInput('climdata',
				'Select Climatic Data', 
				choices=c('All',
					'BIO1 = Annual Mean Temperature',
					'BIO2 = Mean Diurnal Range')
				),
			tags$hr(),
			p('Data downloaded from the ', a('WorldClim database', href='http://www.worldclim.org/'))
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

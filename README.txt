Order of usage of scripts for the workflow to detect AI/LOH from BAF data and confirm CNV calls: 1 - extractBAF.py, 2 - detectAILOH_window.py or detectAILOH_cbs.py, 3 - confirmCNV.py

extractBAF.py
	Extracts BAF values from a gVCF file
	Outputs: 	Two tab-separated files containing BAF values of all SNPs and BAF values of filtered SNPs for each chromosome
			A tab-separated file containing regions of poor coverage (depth below 10) larger than 100bp for each chromosome
					
	Usage:		python extractBAF.py 		

		Positional arguments:
			-i 	path to input gVCF Files 
                	-o 	output directory for BAF files
				
		Optional arguments:
			-n 	sample name, default: 'sample' 
			--minDp 	minimum depth for filtering, default: 0
                	--mingQual 	minimum genotyping quality for filtering, default: 0
			-v	visualize extracted BAF data
				
detectAILOH_window.py
	Detects AI/LOH events from BAF data using a fixed size sliding window
	Outputs:	BED-formatted file containing detected AI/LOH calls
	
	Usage:		python detectAILOH_window.py
	
		Positional arguments:
			-n 	sample name (sample name must match one used for BAF extraction)		
			-b 	path to directory containing BAF files for all chromosomes
			-g	path to BED file containing coordinates of chromosomes under investigation
                	-o 	output directory
			-n 	sample name (sample name must match one used for BAF extraction)
				
		Optional arguments:
			--windowSize	size of the sliding window, default: 50000 
                	--windowStep 	sliding step of the window, default: 50000
			--ref	path to directory containing SNP reference set
			--ncores 	number of cores to be used parallelly
				
detectAILOH_cbs.py
	Detects AI/LOH events from BAF data using Circular Binary Segmentation (CBS)
	Needs to be run from the same directory where script cbsBAF.R is located
	Outputs:	BED-formatted file containing detected AI/LOH calls
	
	Usage:		python detectAILOH_cbs.py
	
		Positional arguments:
			-n 	sample name (sample name must match one used for BAF extraction)
			-b 	path to directory containing BAF files for all chromosomes
                	-o 	output directory
							
		Optional arguments:
			--hom	BAF cutoff for homozygous SNPs, default: 0.9, i.e. SNPs with BAF > 0.9 are homozygous
			--ai	mBAF cutoff for calling AI, default: 0.15
			--loh	mBAF cutoff for calling LOH, default: 0.4
				
confirmCNV.py
	Identify TP CNV calls by comparing CNV calls against AI/LOH calls
	Needs bedtools added in PATH
	Outputs:	BED-formatted file containint TP CNV calls
	
	Usage:		python confirmCNV.py
	
		Positional arguments:
			-n 	sample name
			-c	path to BED-formatted file containing CNV calls
			-s	path to BED-formatted file containing AI/LOH calls
			-o 	output directory
				
				-
				
				



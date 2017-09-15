/*
 *  Parse the input parameters
 */
params.genome     = "data/genome.fa"
params.variants   = "data/known_variants.vcf.gz"
params.blacklist  = "data/blacklist.bed"
params.reads      = "data/reads/ENCSR000COQ1_{1,2}.fastq.gz" 
params.results    = "results" 
params.gatk       = "/opt/broad/GenomeAnalysisTK.jar" 

genome_file     = file(params.genome)
variants_file   = file(params.variants)
blacklist_file  = file(params.blacklist)
reads_ch        = Channel .fromFilePairs( params.reads )
GATK            = params.gatk
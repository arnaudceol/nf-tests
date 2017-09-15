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

reads_ch.println()

/*
 * Process 1A: Create a FASTA genome index with samtools
 */
process '1A_prepare_genome_samtools' {
	input: 
		file genome from genome_file
	
	output:
		file "${genome}.fai" into genome_index_ch
		
		
	script:
	"""
		samtools faidx ${genome} 	
	"""
		
}


/*
 * Process 1B: Create a FASTA genome sequence dictionary with Picard for GATK
 */

process '1B_prepare_genome_picard' {

  input:
      file genome from genome_file

  output:
      file "${genome.baseName}.dict" into genome_dictionary_ch

  script:
  """
  PICARD=`which picard.jar`
  java -jar \$PICARD CreateSequenceDictionary R= $genome O= ${genome.baseName}.dict
  """
}


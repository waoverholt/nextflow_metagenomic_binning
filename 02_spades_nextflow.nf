#!/usr/bin/env nextflow

params.reads = "$HOME/work/test_metaGs/02_qaqc/*_{1,2}.fastq.gz"
params.output_dir = "$HOME/work/test_metaGs/03_Assembly_noerror-correct"

Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs }

process metaspades {
	conda '/home/overholt/.conda/envs/metagenomics'
    
	tag {pair_id}
    publishDir "${params.output_dir}", mode: 'copy'
	
	cpus = '4' 
	clusterOptions = '--mem-per-cpu=10G'
	time = '4h'
 
    input:
    set pair_id, file(reads) from read_pairs

    output:  
    set pair_id, file("${pair_id}/*") into spades_assembly

    script:
    """
    spades.py --meta -o ${pair_id} -1 ${reads[0]} -2 ${reads[1]} -t ${task.cpus} \
    -m 50 --tmp-dir $TMP --only-error-correction
    """
}

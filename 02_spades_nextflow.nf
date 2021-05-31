#!/usr/bin/env nextflow

params.reads = "/home/overholt/kuesel_data/Siberia/02_qaqc/deplete_euks/error_correct/*_R{1,2}.fastq"
params.output_dir = "/home/overholt/kuesel_data/Siberia/03_Assembly/"

Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs }

process metaspades {
	conda '/home/overholt/.conda/envs/metagenomics'
    
	tag {pair_id}
    publishDir "${params.output_dir}", mode: 'copy'
	
	cpus = '4' 
	clusterOptions = '--mem-per-cpu=80G'
	time = '96h'
 
    input:
    set pair_id, file(reads) from read_pairs

    output:  
    set pair_id, file("${pair_id}/*") into spades_assembly

    script:
    """
    export OMP_NUM_THREADS=${task.cpus}
    spades.py --meta -o ${pair_id} -1 ${reads[0]} -2 ${reads[1]} -t ${task.cpus} \
    -m 300 --tmp-dir $TMP --only-assembler 
    """
}

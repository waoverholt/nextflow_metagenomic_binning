#!/usr/bin/env nextflow

params.reads = "/home/overholt/kuesel_data/Siberia/02_qaqc/deplete_euks/*_R{1,2}.fastq"
params.output_dir = "/home/overholt/kuesel_data/Siberia/02_qaqc/deplete_euks/error_correct/"

Channel
    .fromFilePairs( params.reads )
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .set { read_pairs }

process tadpole {
	tag {pair_id}
    publishDir "${params.output_dir}", mode: 'copy'

    executor = 'slurm'	
	cpus = '5' 
	clusterOptions = '--mem-per-cpu=11G'
	time = '24h'
 
    input:
    set pair_id, file(reads) from read_pairs

    output:  
    set pair_id, file("${pair_id}/*") into spades_assembly

    script:
    """
    #!/usr/bin/env bash
    CONDA_BASE=\$(conda info --base)
    source \$CONDA_BASE/etc/profile.d/conda.sh

    conda activate metagenomics

    

    #tadpole.sh in=${reads[0]} in2=${reads[1]} \
    #    out=${pair_id}_err_correct_R1.fastq out2=${pair_id}_err_correct_R2.fastq \
    #    mode=correct k=50 threads=${task.cpus} \
    #    prealloc=t -Xmx50g 
    """
}

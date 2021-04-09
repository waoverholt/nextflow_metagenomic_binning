#!/usr/bin/env nextflow

params.input_dir = "$HOME/data/Projects/example_metaGs/test_metaGs/03_Assembly/"
params.output_dir = "$HOME/data/Projects/example_metaGs/test_metaGs/04_Binning/"
params.qaqc = "$HOME/data/Projects/example_metaGs/test_metaGs/02_qaqc/"

sample_dir = Channel
    .fromPath("${params.input_dir}/H*/", type: 'dir')

process trim_scaffolds {
    tag {sample_id.baseName}

    executor 'local'
    cpus '1'

    input:
    file(sample_id) from sample_dir
    //file(contigs) from contigs_ch

    output:
    tuple file(sample_id), file("${sample_id}/scaffolds_1000.fasta") into trimmed_1000bp_scaffold
    file("${sample_id}/scaffolds_1000_ids.txt") into trimmed_1000bp_ids
    tuple file(sample_id), file("${sample_id}/scaffolds_3000.fasta") into trimmed_3000bp_scaffold
    file("${sample_id}/scaffolds_3000_ids.txt") into trimmed_3000bp_ids

    script:
    """
    #remove spaces in the header & concatenate multiline fasta files into 1 line (header1 /n sequence1 /n header2...)
    sed -e 's/ /_/g' ${sample_id}/scaffolds.fasta | awk '/^[>;]/ { if (seq) { print seq }; seq=""; print } /^[^>;]/ { seq = seq \$0 } END { print seq }' > temp_scaffolds.fasta

    #get contigs/scaffolds >1kb & >3kb
    cat temp_scaffolds.fasta | paste -d ";" - -  | perl -F\\; -ane 'if(length(\$F[1]) > 999) {print \$F[0]."\n".\$F[1]}' > scaffolds_1000.fasta
    cat temp_scaffolds.fasta | paste -d ";" - -  | perl -F\\; -ane 'if(length(\$F[1]) > 2999) {print \$F[0]."\n".\$F[1]}' > scaffolds_3000.fasta

    #get the header file that binsanity needs
    grep "^>" scaffolds_1000.fasta | sed -e 's/>//g' > scaffolds_1000_ids.txt
    grep "^>" scaffolds_3000.fasta | sed -e 's/>//g' > scaffolds_3000_ids.txt

    #make sure these get copied back to the assembly directory, in case something goes wrong (this is a bit lazy, but I couldn't get nextflow to do this easily)
    cp scaffolds_1000.fasta ${params.input_dir}/${sample_id.baseName}/scaffolds_1000.fasta
    cp scaffolds_1000_ids.txt ${params.input_dir}/${sample_id.baseName}/scaffolds_1000_ids.txt
    cp scaffolds_3000.fasta ${params.input_dir}/${sample_id.baseName}/scaffolds_3000.fasta
    cp scaffolds_3000_ids.txt ${params.input_dir}/${sample_id.baseName}/scaffolds_3000_ids.txt
    """
}

qaqc_path = Channel
    .fromPath( "${params.qaqc}")
qaqc_file = qaqc_path.first() 

process run_metawrap {
    tag {sample_id.baseName}

    cpus '4'

    input:
    tuple file(sample_id), file(trimmed) from trimmed_1000bp_scaffold
    file(qaqc_dir) from qaqc_file

    output:
    file("metawrap_binning/work_files/") into metawrap_bams 
    tuple file(sample_id), file("metawrap_binning/maxbin2_bins/") into maxbins 
    file("metawrap_binning/metabat2_bins") into metabat_bins 
    
    shell:
    """
    #!/usr/bin/env bash
    CONDA_BASE=\$(conda info --base)
    source \$CONDA_BASE/etc/profile.d/conda.sh

    conda activate metawrap
    metawrap binning -o metawrap_binning/ -t ${task.cpus} --universal -a ${trimmed} --maxbin2 --metabat2 ${qaqc_dir}/*fastq

    cp -r metawrap_binning ${params.output_dir}/${sample_id.baseName}/
    
    mkdir ${params.output_dir}/${sample_id.baseName}/logs/
    cp .command.sh ${params.output_dir}/${sample_id.baseName}/logs/metawrap_binning_command.sh
    cp .command.log ${params.output_dir}/${sample_id.baseName}/logs/metawrap_binning.log
    """
}

process binsanity_profiler {
    tag {sample_id.baseName}

    cpus 10

    input:
    path(bams) from metawrap_bams
    tuple file(sample_id), file(trimmed) from trimmed_3000bp_scaffold
    file(trimmed_ids) from trimmed_3000bp_ids

    output:
    file("binsanity.cov") into binsanity_cov
    tuple file(sample_id), file("binsanity.cov.x100.lognorm") into binsanity_lognorm
    file("binsanity_binning/*") into binsanity_profiler_out

    shell:
    """
    #!/usr/bin/env bash
    CONDA_BASE=\$(conda info --base)
    source \$CONDA_BASE/etc/profile.d/conda.sh

    conda activate metagenomics
    Binsanity-profile -i ${trimmed} -s ${bams} -T ${task.cpus} -o binsanity_binning/ --ids ${trimmed_ids} -c binsanity
    cp binsanity.c* ${params.output_dir}/${sample_id.baseName}/
    
    cp .command.sh ${params.output_dir}/${sample_id.baseName}/logs/binsanity_profiler_command.sh
    cp .command.log ${params.output_dir}/${sample_id.baseName}/logs/binsnity_profiler.log
 
"""
}

process binsanity_binning {
    tag {sample_id.baseName}
   
    cpus 10

    input:
    tuple file(sample_id), file(cov) from binsanity_lognorm

    output:
    path "binsanity_binning/*" into binsanity_bins_out
    tuple file(sample_id), path("binsanity_binning/BinSanity-Final-bins") into final_binsanity_bins 
    tuple file(sample_id), file("binsanity_binning/BinSanity-Final-bins/*fna") into list_binsanity_bins

    shell:
    """
    #!/usr/bin/env bash
    CONDA_BASE=\$(conda info --base)
    source \$CONDA_BASE/etc/profile.d/conda.sh

    conda activate metagenomics
    Binsanity-wf -f ${sample_id} -l scaffolds_3000.fasta -c ${cov} -o binsanity_binning --threads ${task.cpus}

    cp .command.sh ${params.output_dir}/${sample_id.baseName}/logs/binsanity_binning_command.sh
    cp .command.log ${params.output_dir}/${sample_id.baseName}/logs/binsnity_binning.log
     
    """
}

process rename_binsanity_bins {
    tag {sample_id.baseName}
    
    input:
    tuple file(sample_id), path(binsanity) from final_binsanity_bins
    
    output:
    path "binsanity_binning/final-renamed/" into binsanity_renamed_bins
    file "binsanity_binning/renaming_bins.txt" into renamed_bins
    file "binsanity_binning/final-renamed/*fa" into renamed_bin_files

    shell:
    """
    mkdir -p binsanity_binning/final-renamed/
    i=0
    for file in \$(find ${binsanity}/ -name "*.fna"); do 
        i=\$((i+1)) 
        echo \$file final-renamed/bin\$i.fa >> binsanity_binning/renaming_bins.txt
        mv \$file binsanity_binning/final-renamed/bin\$i.fa
    done

    cp -r binsanity_binning ${params.output_dir}/${sample_id}/
    """
}


process clean_up {
    tag {sample_id.baseName}
    
    input:
    tuple file(sample_id), file(bin) from list_binsanity_bins

    shell:
    """
    # The || will allow the script to finish in case you are "resuming" a
    # a previous run
    rm -r ${params.output_dir}/${sample_id.baseName}/work_files/maxbin2_out/ || true
    rm ${params.output_dir}/${sample_id.baseName}/work_files/assembly* || true
    rm ${params.output_dir}/${sample_id.baseName}/work_files/*sam || true
    """
}

process refine_bins {
    tag {sample_id.baseName}
    cpus 10
    memory '200 GB'

    input:
    tuple file(sample_id), path(max) from maxbins
    path meta from metabat_bins
    path binsanity from binsanity_renamed_bins

    output:
    file("metawrap_refinement_50_10/*") into metawrap_bins

    shell:
    """
    #!/usr/bin/env bash
    CONDA_BASE=\$(conda info --base)
    source \$CONDA_BASE/etc/profile.d/conda.sh

    conda activate metawrap
    metawrap bin_refinement -o metawrap_refinement_50_10 -t ${task.cpus} \
        -c 50 -x 10 \
        -A ${max} -B ${meta} -C ${binsanity}
    cp -r metawrap_refinement_50_10 ${params.output_dir}/${sample_id.baseName}/

    cp .command.sh ${params.output_dir}/${sample_id.baseName}/logs/metawrap_refine_command.sh
    cp .command.log ${params.output_dir}/${sample_id.baseName}/logs/metawrap_refine.log
 
    """
}

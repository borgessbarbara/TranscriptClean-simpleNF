// TranscriptClean module

process EXTRACT_SPLICE_JUNCTIONS {
    tag "Extract splice junctions from GTF file"
    label 'process_high'

    input:
    tuple val(meta), path(gtf), path(fasta)

    output:
    tuple val(meta), path("${prefix}_spliceJns.txt"), emit: txt
    path  "versions.yml",                            emit: versions

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    // test if it is actually transcriptclean --version
    """
    python home/biodocker/TranscriptClean-2.0.2/accessory_scripts/get_SJs_from_gtf.py \\
        --f ${gtf} \\
        --g ${fasta} \\
        --o ${prefix}_spliceJns.txt \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transcriptclean: \$(echo \$(transcriptclean --version 2>&1) | sed 's/^.*transcriptclean //; s/Using.*\$//' ))
    END_VERSIONS
    """
    stub:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}"_spliceJns.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transcriptclean: \$(echo \$(transcriptclean --version 2>&1) | sed 's/^.*transcriptclean //; s/Using.*\$//' ))
    END_VERSIONS
    """
}

process TRANSCRIPTCLEAN {
    tag "TranscriptClean"
    label 'process_high'

    input:
    tuple val(meta), path(sam), path(fasta), path(vcf), path(splice_junctions)

    output:
    tuple val(meta), path("*.sam"), emit: sam
    tuple val(meta), path("*.fasta"), emit: fasta
    tuple val(meta), path("*.TE.log"), path("*.log"), emit: log
    path  "versions.yml",           emit: versions

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    def variant_flag = params.variant_aware && vcf ? "--variants ${vcf}" : ''
    def splice_flag = params.sj_correction && splice_junctions ? "--spliceJns ${splice_junctions}" : ''
    // def other_params = params.other_params ? params.other_params : ''
    // test if it is actually transcriptclean --version
    """
    python home/biodocker/TranscriptClean-2.0.2/TranscriptClean.py \\
        --SAM ${sam} \\
        --genome ${fasta} \\
        ${variant_flag} \\
        ${splice_flag} \\
        --outprefix "${prefix}" \\
        $args 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transcriptclean: \$(echo \$(transcriptclean --version 2>&1) | sed 's/^.*transcriptclean //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}"_transcriptClean.sam
    touch "${prefix}"_transcriptClean.fasta
    touch "${prefix}"_transcriptClean.TE.log
    touch "${prefix}"_transcriptClean.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transcriptclean: \$(echo \$(transcriptclean --version 2>&1) | sed 's/^.*transcriptclean //; s/Using.*\$//' ))
    END_VERSIONS
    """
}

process GENERATE_REPORT {
    tag "Generate TranscriptClean report"
    label 'process_high'

    input:
    tuple val(meta), path(transcriptclean_logs)

    output:
    tuple val(meta), path("*_report.pdf"), emit: pdf
    path  "versions.yml",          emit: versions

    script:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    Rscript home/biodocker/TranscriptClean-2.0.2/generate_report.R \\
    ../results/${prefix}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transcriptclean: \$(echo \$(transcriptclean --version 2>&1) | sed 's/^.*transcriptclean //; s/Using.*\$//' ))
    END_VERSIONS
    """

    stub:
    args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}"_report.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transcriptclean: \$(echo \$(transcriptclean --version 2>&1) | sed 's/^.*transcriptclean //; s/Using.*\$//' ))
    END_VERSIONS
    """
} 

workflow RUN_TRANSCRIPTCLEAN {

    ch_versions = Channel.empty()

    Channel.fromTuple(params.sam, params.fasta).set { ch_sam_fasta }

    if (params.gtf) {
        Channel.fromPath(params.gtf).ifEmpty(null).set  { ch_gtf } 
    }
    
    if (params.splice_junctions) {
        Channel.fromPath(params.splice_junctions).ifEmpty(null).set { ch_sjs } 
    }

    if (params.vcf) {
        Channel.fromPath(params.vcf).ifEmpty(null).set  { ch_vcf }
    }

    if (ch_sam_fasta.ifEmpty(null) {
        error("Required SAM and FASTA files not inputed.")
    } else if (params.extract_sjs && params.gtf && params.sj_correction){
        EXTRACT_SPLICE_JUNCTIONS(ch_gtf.combine(ch_sam_fasta.map{it[1]})).set { ch_ex_sjs }
        TRANSCRIPTCLEAN(ch_sam_fasta.combine(ch_ex_sjs.out.map{it})).set { ch_transcriptclean_res }
        GENERATE_REPORT(ch_transcriptclean_res.collect()).set { ch_report }
    } else if (params.extract_sjs && params.gtf && params.sj_correction && params.variant_aware && params.vcf) {
        EXTRACT_SPLICE_JUNCTIONS(ch_gtf.combine(ch_sam_fasta.map{it[1]})).set { ch_ex_sjs }
        TRANSCRIPTCLEAN(ch_sam_fasta.combine(ch_ex_sjs.out.map{it}, ch_vcf)).set { ch_transcriptclean_res }
        GENERATE_REPORT(ch_transcriptclean_res.collect()).set { ch_report }
    } else if (params.sj_correction && params.splice_junctions){
        TRANSCRIPTCLEAN(ch_sam_fasta.combine(ch_sjs)).set { ch_transcriptclean_res }
        GENERATE_REPORT(ch_transcriptclean_res.collect()).set { ch_report }
    } else if (params.sj_correction && params.variant_aware && params.splice_junctions && params.vcf){
        TRANSCRIPTCLEAN(ch_sam_fasta.combine(ch_sjs, ch_vcf)).set { ch_transcriptclean_res }
        GENERATE_REPORT(ch_transcriptclean_res.collect()).set { ch_report }
    } else if (params.variant_aware && params.vcf) {
        TRANSCRIPTCLEAN(ch_sam_fasta.combine(ch_vcf)).set { ch_transcriptclean_res }
        GENERATE_REPORT(ch_transcriptclean_res.collect()).set { ch_report }
    } else {
        TRANSCRIPTCLEAN(ch_sam_fasta).set { ch_transcriptclean_res }
        GENERATE_REPORT(ch_transcriptclean_res.collect()).set { ch_report }
    }

    ch_versions = Channel.fromPath("/versions.yml")

    emit:
    extracted_sjs = ch_ex_sjs.out.ifEmpty(null)
    transcriptclean_results = ch_transcriptclean_res.out
    transcriptclean_report = ch_report.out
    versions       = ch_versions                

}

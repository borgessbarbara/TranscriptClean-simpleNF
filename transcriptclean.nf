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
    """
    python /home/biodocker/TranscriptClean-2.0.2/accessory_scripts/get_SJs_from_gtf.py \\
        --f ${gtf} \\
        --g ${fasta} \\
        --o ${prefix}_spliceJns.txt \\
        $args

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
    def variant_flag = params.variant_aware && vcf.name != 'null' ? "--variants ${vcf}" : ''
    def splice_flag = params.sj_correction && splice_junctions.name != 'null' ? "--spliceJns ${splice_junctions}" : ''
    """
    python /home/biodocker/TranscriptClean-2.0.2/TranscriptClean.py \\
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
    Rscript /home/biodocker/TranscriptClean-2.0.2/generate_report.R \\
    ../results/${prefix} \\
    ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        transcriptclean: \$(echo \$(transcriptclean --version 2>&1) | sed 's/^.*transcriptclean //; s/Using.*\$//' ))
    END_VERSIONS
    """
} 

workflow RUN_TRANSCRIPTCLEAN {
    main:
    ch_versions = Channel.empty()

    // Create input channels
    ch_sam_fasta = Channel.fromFilePairs(params.sam, flat: true, checkIfExists: true)
        .combine(Channel.fromPath(params.fasta, checkIfExists: true))
        .map { meta, sam, fasta -> [meta, sam, fasta] }

    if (params.gtf) {
        ch_gtf = Channel.fromPath(params.gtf, checkIfExists: true)
    } else {
        ch_gtf = Channel.empty()
    }
    
    if (params.splice_junctions) {
        ch_sjs = Channel.fromPath(params.splice_junctions, checkIfExists: true)
    } else {
        ch_sjs = Channel.empty()
    }

    if (params.vcf) {
        ch_vcf = Channel.fromPath(params.vcf, checkIfExists: true)
    } else {
        ch_vcf = Channel.empty()
    }

    if (params.extract_sjs && params.gtf && params.sj_correction) {
        ch_ex_sjs = EXTRACT_SPLICE_JUNCTIONS(ch_gtf.combine(ch_sam_fasta.map{ [it[0], it[2]] }))
        ch_transcriptclean_res = TRANSCRIPTCLEAN(ch_sam_fasta.combine(ch_vcf, ch_ex_sjs.out.txt))
    } 
    else if (params.sj_correction && params.splice_junctions) {
        ch_transcriptclean_res = TRANSCRIPTCLEAN(ch_sam_fasta.combine(ch_vcf, ch_sjs))
    }
    else if (params.variant_aware && params.vcf) {
        ch_transcriptclean_res = TRANSCRIPTCLEAN(ch_sam_fasta.combine(ch_vcf))
    }
    else {
        ch_transcriptclean_res = TRANSCRIPTCLEAN(ch_sam_fasta)
    }

    ch_report = GENERATE_REPORT(ch_transcriptclean_res.out.log.collect())
    ch_versions = ch_report.out.versions.mix(ch_transcriptclean_res.out.versions)

    emit:
    extracted_sjs = ch_ex_sjs ? ch_ex_sjs.out.txt : Channel.empty()
    transcriptclean_results = ch_transcriptclean_res.out
    transcriptclean_report = ch_report.out
    versions = ch_versions
}

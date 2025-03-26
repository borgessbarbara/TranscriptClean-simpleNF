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

workflow {
    main:
    // Validate mandatory inputs
    // Create input channels

    if (!params.sam) {
        error "Mandatory parameter --sam not specified" 
    } else {
        ch_sam = Channel.fromFilePath(params.sam, checkIfExists: true)
    }
    if (!params.fasta) {
        error "Mandatory parameter --fasta not specified"
    } else {
        ch_fasta = Channel.fromPath(params.fasta, checkIfExists: true)
    }
    
    // Combine with metadata
    ch_sam_fasta = ch_sam.combine(ch_fasta)

    // Handle optional inputs
    ch_gtf = params.gtf ? Channel.fromPath(params.gtf, checkIfExists: true) : Channel.empty()
    ch_sjs = params.splice_junctions ? Channel.fromPath(params.splice_junctions, checkIfExists: true) : Channel.empty()
    ch_vcf = params.vcf ? Channel.fromPath(params.vcf, checkIfExists: true) : Channel.empty()

    // Main processing with proper output handling
    if (params.extract_sjs && params.gtf && params.sj_correction) {
        EXTRACT_SPLICE_JUNCTIONS(ch_gtf.combine(ch_fasta))
        TRANSCRIPTCLEAN(ch_sam_fasta.combine(ch_vcf, EXTRACT_SPLICE_JUNCTIONS.out.txt))
        GENERATE_REPORT(TRANSCRIPTCLEAN.out.log.collect())
        
        // Define outputs only for this branch
        ext_sjs_results = EXTRACT_SPLICE_JUNCTIONS.out.txt
        transcriptclean_results = TRANSCRIPTCLEAN.out
        report_output = GENERATE_REPORT.out
    } 
    else if (params.sj_correction && params.splice_junctions) {
        TRANSCRIPTCLEAN(ch_sam_fasta.combine(ch_vcf, ch_sjs))
        GENERATE_REPORT(TRANSCRIPTCLEAN.out.log.collect())
        
        transcriptclean_results = TRANSCRIPTCLEAN.out
        report_output = GENERATE_REPORT.out
    }
    else if (params.variant_aware && params.vcf) {
        TRANSCRIPTCLEAN(ch_sam_fasta.combine(ch_vcf))
        GENERATE_REPORT(TRANSCRIPTCLEAN.out.log.collect())
        
        transcriptclean_results = TRANSCRIPTCLEAN.out
        report_output = GENERATE_REPORT.out
    }
    else {
        TRANSCRIPTCLEAN(ch_sam_fasta)
        GENERATE_REPORT(TRANSCRIPTCLEAN.out.log.collect())
        
        transcriptclean_results = TRANSCRIPTCLEAN.out
        report_output = GENERATE_REPORT.out
    }
}

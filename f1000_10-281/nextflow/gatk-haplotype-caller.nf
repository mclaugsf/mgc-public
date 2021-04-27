#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    ============================================================================
     :  Git version: ${version}
    ============================================================================
    Usage:

       Mandatory arguments:
           --out_dir        Output Directory (for things that aren't intermediate files to keep)
           --vcf            VCF to annotate w/snpEff
    """.stripIndent()
}

workingDir = new File(".").getCanonicalPath()

baseDir = file("$baseDir")

def proc_git = "git -C $baseDir rev-parse HEAD".execute()
//version = proc.text.trim().substring(0,6)
version = proc_git.text.trim()
//println params.g


params.help = false
if (params.help){
    println(params)publishDir params."out_dir", mode: 'copy', overwrite: true, pattern: '*vcf.gz*'
    helpMessage()
    exit 0
}

gatk_haplotype_caller_input = Channel.fromPath(params."bam")
  .map {it -> [it.baseName.replaceFirst(/.bam$/, ""), it]}
  .map {it -> [it[0], it[1], file(it[1] + '.bai')]}

process gatk_haplotype_caller {
    tag {'gatk_haplotype_caller' + '-' + sample}
    publishDir params."out_dir", mode: 'copy', overwrite: true, pattern: '*vcf.gz*'

	container 'broadinstitute/gatk:4.1.6.0'
    memory '16 G'
    cpus 4

    input:
    set val(rsp), file(bam), file(bam_index) from gatk_haplotype_caller_input
    file(ref) from file(params."ref")
    file(fai) from file(params."ref" + '.fai')
    file(dict) from file(params."ref".replaceFirst(/\.fa$|\.fasta$/, ".dict"))

    script:
    cpu    = task.cpus
    memory = task.memory.toGiga()
    """
    /gatk/gatk \
        --java-options \
        "-Xmx${memory}G -Xms${memory}G" \
    HaplotypeCaller \
        --native-pair-hmm-threads ${cpu} \
        --input ${bam} \
        --output ${rsp}.vcf.gz \
        --reference $ref
    """

}

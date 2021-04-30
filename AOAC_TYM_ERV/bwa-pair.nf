#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    ============================================================================
     :  Git version: ${version}
    ============================================================================
    Usage:

       Mandatory arguments:
           --fq             FastQ1 Pair (i.e. JL.R{1,2}.fastq.gz)
           --r              Fasta Reference
           --out_dir        Output Directory (for things that aren't intermediate files to keep)

    """.stripIndent()
}

//paramater parsing and error checking:

workingDir = new File(".").getCanonicalPath()

baseDir = file("$baseDir")

def proc_git = "git -C $baseDir rev-parse HEAD".execute()
//version = proc.text.trim().substring(0,6)
version = proc_git.text.trim()
//println params.g

params.help = false
if (params.help){
    println(params)
    helpMessage()
    exit 0
}

bwa_index = Channel
              .fromPath(params."r") // + '.{amb,ann,bwt,pac,sa}')
              .map{it -> [it.simpleName.replaceFirst(/_.*/, ""), it, file(it + '.{amb,ann,bwt,pac,sa}')]}
              //.map{it -> [file(it[0]), file(it[0] + '.{amb,ann,bwt,pac,sa}')]}

Channel
    .fromFilePairs(params."fq")
    //take everything before the underscore
    //.map{it -> [it[0].replaceFirst(/_.*/, ""), it[1], file(params."r"), file(params."r" + '.{amb,ann,bwt,pac,sa}')]}
    .map{it -> [it[0].replaceFirst(/_.*/, ""), it[1]]}
    .phase(bwa_index)
    .map{it -> [it[0][0], it[0][1], it[1][1], it[1][2]]}
    .set {bwa_input}

process bwa {
    publishDir params."out_dir", mode: 'copy', overwrite: true

    tag {'bwa' + '-' + id_run}

    container 'medicinalgenomics/bwa_sambamba_picard'

    memory '30 G'
    cpus 8

    input:
    set val(id_run), file(fq), file(fasta_ref), file(bwa_index) from bwa_input

    output:
    file("${id_run}.bam") into bwa_bam
    file("${id_run}.bam.bai") into bwa_bam_index

    script:
    cpu    = task.cpus
    id_sample = id_run
    """
    bwa mem -R "@RG\\tID:$id_run\\tPU:$id_run\\tSM:$id_run\\tLB:$id_run\\tPL:illumina" \
            -t $cpu -M $fasta_ref ${fq[0]} ${fq[1]} | \
    samtools view -hu - \
            | sambamba sort --tmpdir=. /dev/stdin -o ${id_sample}.bam
    """

}

#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    ============================================================================
     :  Git version: ${version}
    ============================================================================
    Usage:

       Mandatory arguments:
           --out_dir        Output Directory (for things that aren't intermediate files to keep)
           --f              FastQ Path (i.e. JL.R1.fastq.gz)
    """.stripIndent()
}

def proc_git = "git -C $baseDir rev-parse HEAD".execute()
version = proc_git.text.trim()

params.help = false
if (params.help){
    println(params)
    helpMessage()
    exit 0
}

Channel
	.fromFilePairs(params."f")
	.set {megahit_input}

process megahit {
    errorStrategy 'finish' //don't kill all the jobs if one dies
	publishDir params."out_dir" + '/' + 'megahit/', mode: 'copy', overwrite: true, pattern: '*'

	tag {'megahit' + '-' + id_run}

    container 'medicinalgenomics/megahit'
    memory '8 G'

    cpus 4

    input:
    set val(id_run), file(fq) from megahit_input

    output:
    file("${id_run}") into megahit_output

    script:
    cpu = task.cpus
    memory = "${task.memory.toGiga()}"
    """
    mkdir tmp

    megahit \
    	-1 ${fq[0]} \
    	-2 ${fq[1]} \
    	-o ${id_run} \
    	-t ${cpu} \
    	-m ${memory} \
        --tmp-dir tmp
    """
}

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
	.set {spades_input}

//spades_input
//    .view()


process megahit {
    errorStrategy 'finish' //don't kill all the jobs if one dies
	publishDir params."out_dir" + '/' + 'megahit/', mode: 'copy', overwrite: true, pattern: '*'

	tag {'megahit' + '-' + id_run}

    container 'medicinalgenomics/megahit'
    memory '30 G'

    cpus 8

    input:
    set val(id_run), file(fq) from spades_input

    output:
    file("${id_run}") into spades_output

    script:
    cpu = task.cpus
    memory = "${task.memory.toGiga()}"
    ////using /home/tools/SPAdes-3.14.1-Linux/bin/spades.py -1 /NGS/Data/PathoSEEK/AOAC_PLATES/BF3/BF3_R1_001.fastq.gz -2 /NGS/Data/PathoSEEK/AOAC_PLATES/BF3/BF3_R2_001.fastq.gz -o /NGS/Data/PathoSEEK/AOAC_PLATES/BF3/BF3_spades2 -k33 -t8 -m192
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

#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    Usage:

       Mandatory arguments:
       	   --r              Fasta Reference
       	   --out_dir        Output Directory (for things that aren't intermediate files to keep)
       Optiona:
           --blast_only     Only make BLAST databases (no bwa)
    """.stripIndent()
}

params."blast_only" = 'false'

params.help = false
if (params.help){
    println(params)
    helpMessage()
    exit 0
}

ref_input = Channel.fromPath(params.r)
			  .map {it -> [it.baseName, it]}

ref_input.into {
	ref_input1
	ref_input2
	ref_input3
}

process dict_and_fai {

	publishDir params."out_dir" + '/', mode: 'copy', overwrite: true, pattern: '*.dict'
	publishDir params."out_dir" + '/', mode: 'copy', overwrite: true, pattern: '*.fai'

	tag {'dict_and_fai' + "-" + prefix}

	container 'medicinalgenomics/bwa_sambamba_picard'
    memory '30G'
    cpus 1

    input:
    set val(prefix), file(fasta_ref) from ref_input1
    val(blast_only) from params."blast_only"

    output:
    set val(prefix), file("${fasta_ref}.fai"), file("${prefix}.dict") into dict_file_fai_out

    script:
    """
    if (${blast_only} == 'true')
    then
        echo "only making BLAST database...";
        touch ${fasta_ref}.fai
        touch ${prefix}.dict
    else
        samtools faidx ${fasta_ref}

        java -jar /opt/picard.jar \
    	   CreateSequenceDictionary \
    	   REFERENCE=${fasta_ref}
    	   OUTPUT=${prefix}.dict
    fi
    """
}

process bwa_index {

	publishDir params."out_dir" + '/', mode: 'copy', overwrite: true, pattern: '*.{amb,ann,bwt,pac,sa}'
	publishDir params."out_dir" + '/', mode: 'copy', overwrite: true, pattern: '*.{fasta,fa}'

	tag {'bwa_index' + "-" + prefix}

	container 'medicinalgenomics/bwa_sambamba_picard'
    memory '30G'
    cpus 1

    input:
    set val(prefix), file(fasta_ref) from ref_input2
    val(blast_only) from params."blast_only"

    output:
    set val(prefix), file("*.{amb,ann,bwt,pac,sa}"), file(fasta_ref) into bwa_index_out

    script:
    """
    if (${blast_only} == 'true')
    then
        echo "only making BLAST database...";
        touch nothing.amb
        touch nothing.ann
        touch nothing.bwt
        touch nothing.pac
        touch nothing.sa
    else
        bwa index ${fasta_ref}
    fi
    """
}

process blast_index {
    errorStrategy 'ignore' //don't kill all the jobs if one dies
	publishDir params."out_dir" + '/', mode: 'copy', overwrite: true, pattern: '*.{ndb,nhr,nin,not,nsq,ntf,nto}'

	tag {'blast_index' + '-' + prefix}

	container 'ncbi/blast'
	memory '4G'
	cpus 1

	input:
    set val(prefix), file(fasta_ref) from ref_input3

    output:
    set val(prefix), file("*.{ndb,nhr,nin,not,nsq,ntf,nto}") into blast_index_out

    script:
    """
    makeblastdb -dbtype nucl \
    			-in ${fasta_ref}
    """
}

#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    ============================================================================
     :  Git version: ${version}
    ============================================================================
    Usage:

       Mandatory arguments:
       	   --fasta          denovo-assembled FASTA reference file from the same BAM file used in --bam_file
           --bam            BAM filed mapped to denovo assembled reference
           --blast_tsv      Results of blasting the assembly to blast nt db (see blast-fungus-blobtools.nf)
           --out_dir        Output Directory (for things that aren't intermediate files to keep)
           --json           JSON file with strain names to put into the plots
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

bam_file = Channel.fromPath(params.bam)
			      .map {it -> [it.simpleName, it]}

fasta_file = Channel.fromPath(params.fasta)
					.map {it -> [it.simpleName.replaceFirst(/-megahit/, ""), it]}

blast_file = Channel.fromPath(params.blast_tsv)
                    .map {it -> [it.simpleName, it]}

process read_priority_file {
    publishDir params."out_dir" + '/' + 'priority-list.txt', mode: 'copy', overwrite: true, pattern: 'priority-list.txt'
    container 'medicinalgenomics/r-with-bam-vcf'
    memory '500 MB'
    cpus 1

    input:
    file(rsp_file) from file(params."rsp_file")
    file(rsps_txt) from file('priority-list.txt')

    output:
    stdout into rsps

    """
    cat ${rsp_file}
    cp ${rsp_file} rsps.txt
    """
}

rsps
   .flatMap {n -> n.split(/\n/).collect()}
   .set{rsps}

bam_file
	.phase(fasta_file)
    .map {it -> it.flatten()}
    .map {it -> [it[0], it[1], it[1] + '.bai', it[3]]}
    .phase(blast_file)
    .map {it -> it.flatten()}
    .map {it -> [it[0], it[1], it[1] + '.bai', it[3], it[5]]}
	.set {blobtools_create_input}

process blobtools_create {
	publishDir params."out_dir" + '/' + 'coverage_json/', mode: 'copy', overwrite: true, pattern: '*.json'

    container 'medicinalgenomics/blobtools'
    memory '4 GB'
    cpus 1

    input:
    set val(sample), file(bam), file(bam_index), file(fasta), file(blast_tsv) from blobtools_create_input

    output:
    set val(sample), file("*.json") into blobtools_create_output

    script:
    """
    mkdir ${sample}-blobtools-create

    /blobtools/blobtools create -i ${fasta} \
    				   -b ${bam} \
    				   -o ${sample}-blobtools-create/${sample} \
                       --hitsfile ${blast_tsv}

    cp ${sample}-blobtools-create/*.json .
    """
}

process blobtools_plot {
	publishDir params."out_dir" + '/' + 'coverage_plots/', mode: 'copy', overwrite: true, pattern: '*read_cov.png'
	publishDir params."out_dir" + '/' + 'barplot/', mode: 'copy', overwrite: true, pattern: '*barplot.png'
	publishDir params."out_dir" + '/' + 'stats/', mode: 'copy', overwrite: true, pattern: '*stats.txt'

	container 'medicinalgenomics/blobtools'
    memory '4 GB'
    cpus 1

    input:
    set val(sample), file(sample_json) from blobtools_create_output
    file(psp_json) from file(params."json")

    output:
    set val(sample), file("${sample}*barplot.png"), file("${sample}*read_cov.png"), file("${sample}*stats.txt") into blobtools_plot_out

    script:
    """
    strain_name=`jq -r .${sample} ${psp_json}`
    /blobtools/blobtools plot -i ${sample_json} --out \$strain_name

    mv *.blobDB.json.bestsum.phylum.p8.span.100.blobplot.bam0.png ${sample}-\${strain_name}-barplot.png

    mv *.blobDB.json.bestsum.phylum.p8.span.100.blobplot.read_cov.bam0.png ${sample}-\${strain_name}-read_cov.png

    mv *.blobDB.json.bestsum.phylum.p8.span.100.blobplot.stats.txt ${sample}-\${strain_name}-stats.txt
    """
}

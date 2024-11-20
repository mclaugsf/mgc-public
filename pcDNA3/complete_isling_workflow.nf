#!/usr/bin/env nextflow

def helpMessage() {
    log.info"""
    ============================================================================
     :  Git version: ${version}
    ============================================================================
    Usage:

       Mandatory arguments:
       	   --vector_ref     Vector reference bwa indexed FASTA
       	   --fq             Paired FASTQ data for processing

       Optional arguments:
       	   --host_ref       Host reference bwa indexed FASTA (DEFAULT: ${params."host_ref"})
           --out_dir        Output Directory (for things that aren't intermediate files to keep)

    """.stripIndent()
}

//defaults:
params."out_dir" = 'out'
params."isling_dir" = '/home/ubuntu/software/isling/'
params."host_ref" = "/Juicer/isling/ref/hs1.fa"
//params."vector_ref" = "/Juicer/isling/ref/pCMV-Spike-sequence-7286-bps.fasta"
params."isling_params" = "-a -Y -A 1 -B 2 -O 6,6 -E 1,1 -L 0,0 -T 10 -h 200"
params."isling" = true

def proc_git = "git -C $baseDir rev-parse HEAD".execute()

version = proc_git.text.trim()

params.help = false
if (params.help){
    println(params)
    helpMessage()
    exit 0
}

Channel
    .fromFilePairs(params."fq")
    .map{it -> [it[0], it[1], file(params."vector_ref"), file(params."vector_ref" + '.{amb,ann,bwt,pac,sa}'), params."isling", params."isling_params"]}
    .set {bwa_map_vector_input}

process bwa_map_vector {
    errorStrategy 'finish' //don't kill all the jobs if one dies
    //publishDir params."out_dir" + "/bwa_map_vector/", mode: 'copy', overwrite: true

    tag {'bwa' + '-' + id_run}

    container 'medicinalgenomics/bwa_sambamba_picard'
    memory '30 G'
    cpus 6

    input:
    set val(id_run), file(fq), file(fasta_ref), file(bwa_index), val(isling), val(isling_params) from bwa_map_vector_input

    output:

    set val(id_run), file("${id_run}*.bam"), file("${id_run}*bam.bai") into bam_input

    script:
    cpu    = task.cpus
    id_sample = id_run
    """
    if ("${isling}" == "true")
    then
        bwa mem -R "@RG\\tID:$id_run\\tPU:$id_run\\tSM:$id_run\\tLB:$id_run\\tPL:illumina" \
                -t $cpu ${isling_params} $fasta_ref ${fq[0]} ${fq[1]} | \
        samtools view -hu - \
            | sambamba sort --tmpdir=. /dev/stdin -o ${id_sample}-vector-mapped.bam
    else
        bwa mem -R "@RG\\tID:$id_run\\tPU:$id_run\\tSM:$id_run\\tLB:$id_run\\tPL:illumina" \
                -t $cpu -M $fasta_ref ${fq[0]} ${fq[1]} | \
        samtools view -hu - \
            | sambamba sort --tmpdir=. /dev/stdin -o ${id_sample}-vector-mappped.bam
    fi
    """
}

bam_input.into {
	bam_input1
	bam_input2
	bam_input3
}

process bam_to_fastq {

	tag {'bam_to_fastq' + "-" + rsp}
	publishDir params."out_dir" + '/fastq-from-vector-mapping/', mode: 'copy', overwrite: true, pattern: "*.fastq.gz"

	container 'medicinalgenomics/samtools:latest'
	memory '4 G'
	cpus 4

	input:
	set val(rsp), file(bam), file(bam_index) from bam_input1

	output:
	set val(rsp), file("${rsp}*_R1_001.fastq.gz"), file("${rsp}*_R2_001.fastq.gz") into bam_to_fastq_out

	script:
	"""
	all_contigs=`samtools idxstats ${bam} | cut -f 1 | grep -v "^*\$" | perl -pe "s/\\n/ /g" | perl -pe "s/ \$//"`

	samtools view -u ${bam} \${all_contigs} | \
	samtools collate -O - | \
	samtools fastq -1 ${rsp}-vector-mapped_R1_001.fastq.gz -2 ${rsp}-vector-mapped_R2_001.fastq.gz 
	"""
}

bam_to_fastq_out
	.map{it -> [it[0], [it[1], it[2]], file(params."host_ref"), file(params."host_ref" + '.{amb,ann,bwt,pac,sa}'), params."isling", params."isling_params"]}
	.set {bwa_input}

process bwa_map_host {
    errorStrategy 'finish'
    publishDir params."out_dir" + "/bwa_map_host-vector-mapped-reads/", mode: 'copy', overwrite: true

    tag {'bwa' + '-' + id_run}

    container 'medicinalgenomics/bwa_sambamba_picard'
    memory '30 G'
    cpus 8

    input:
    set val(id_run), file(fq), file(fasta_ref), file(bwa_index), val(isling), val(isling_params) from bwa_input

    output:

    set val(id_run), file("${id_run}*.bam"), file("${id_run}*.bam.bai") into bam_map_host_output

    script:
    cpu    = task.cpus
    id_sample = id_run
    """
    if ("${isling}" == "true")
    then
        bwa mem -R "@RG\\tID:$id_run\\tPU:$id_run\\tSM:$id_run\\tLB:$id_run\\tPL:illumina" \
                -t $cpu ${isling_params} $fasta_ref ${fq[0]} ${fq[1]} | \
        samtools view -hu - \
            | sambamba sort --tmpdir=. /dev/stdin -o ${id_sample}-host-mapped.bam
    else
        bwa mem -R "@RG\\tID:$id_run\\tPU:$id_run\\tSM:$id_run\\tLB:$id_run\\tPL:illumina" \
                -t $cpu -M $fasta_ref ${fq[0]} ${fq[1]} | \
        samtools view -hu - \
            | sambamba sort --tmpdir=. /dev/stdin -o ${id_sample}-host-mapped.bam
    fi
    """
}

process host_sort_bam_by_read_name {
	tag {"host_sort_bam_by_read_name" + "-" + rsp}
	publishDir params."out_dir" + '/bam-sorted-by-read-name/host', mode: 'copy', overwrite: true, pattern: "*.bam"

	container 'medicinalgenomics/bwa_sambamba_picard'

	memory '12 G'
	cpus 4

	input:
	set val(id_run), file(input_bam), file(input_bam_index) from bam_map_host_output

	output:
	set val(id_run), file("${id_run}-host-mapped-sorted-by-read-name.bam") into host_bam_sorted_by_read_name

	script:
	"""
	mkdir tmp
	samtools sort -n ${input_bam} -T tmp/tmp > ${id_run}-host-mapped-sorted-by-read-name.bam
	"""
}

process vector_mapped_reads_only {
	tag {"vector_mapped_reads_only" + "-" + rsp}
	publishDir params."out_dir" + '/vector_mapped_reads_only/', mode: 'copy', overwrite: true, pattern: "*.bam*"

	container 'medicinalgenomics/samtools:latest'

	memory '12 G'
	cpus 4

	input:
	set val(id_run), file(input_bam), file(input_bam_index) from bam_input3

	output:
	set val(id_run), file("${id_run}-vector-mapped-reads-only.bam"), file("${id_run}-vector-mapped-reads-only.bam.bai") into vector_mapped_reads_only_output

	script:
	"""
	all_contigs=`samtools idxstats ${input_bam} | cut -f 1 | grep -v "^*\$" | perl -pe "s/\\n/ /g" | perl -pe "s/ \$//"`
	samtools view --bam ${input_bam} \${all_contigs} > ${id_run}-vector-mapped-reads-only.bam
	samtools index ${id_run}-vector-mapped-reads-only.bam
	"""
}

process vector_sort_bam_by_read_name {
	tag {"vector_sort_bam_by_read_name" + "-" + rsp}
	publishDir params."out_dir" + '/bam-sorted-by-read-name/vector/', mode: 'copy', overwrite: true, pattern: "*.bam"

	container 'medicinalgenomics/bwa_sambamba_picard'

	memory '12 G'
	cpus 4

	input:
	set val(id_run), file(input_bam), file(input_bam_index) from vector_mapped_reads_only_output

	output:
	set val(id_run), file("${id_run}-vector-mapped-sorted-by-read-name.bam") into vector_bam_sorted_by_read_name

	script:
	"""
	mkdir tmp
	samtools sort -n ${input_bam} -T tmp/tmp > ${id_run}-vector-mapped-sorted-by-read-name.bam
	"""
}

host_bam_sorted_by_read_name
	.phase(vector_bam_sorted_by_read_name)
	.map{it -> [it[0][0], it[0][1], it[1][1], file(params."isling_dir")]}
	.set{find_ints_input}

process find_ints {
	tag {"isling_find_ints" + "-" + rsp}
	publishDir params."out_dir" + '/find_ints_unfiltered/', mode: 'copy', overwrite: true, pattern: "*.tsv"

	container 'szsctt/isling:latest'

	input:
	set val(id_run), file(host_bam), file(vector_bam), file(isling_dir) from find_ints_input

	output:
	set val(id_run), file("*UNFILTERED.tsv"), file(isling_dir) into find_ints_output

	script:
	"""
	python ${isling_dir}/scripts/find_ints.py \
		--host  ${host_bam} \
		--virus ${vector_bam} \
		--mean-template-length 80 \
		--integrations ${id_run}-find-ints-UNFILTERED.tsv \
		--map-thresh 20 \
		--tolerance 3 \
		--nm-pc 0.6 \
		--nm-diff 2
	"""
}

process filter_find_ints {
	tag {"filter_find_ints" + "-" + rsp}

	publishDir params."out_dir" + '/find_ints_filtered/', mode: 'copy', overwrite: true, pattern: "*.tsv"

	container 'szsctt/isling:latest'

	input:
	set val(id_run), file(unfiltered_ints_input), file(isling_dir) from find_ints_output

	output:
	set val(id_run), file("*FILTERED.tsv") into filter_find_ints_output

	script:
	"""
	python3 ${isling_dir}/scripts/filter.py \
		-i ${unfiltered_ints_input} \
		-k ${id_run}-FILTERED.tsv \
		-e ${id_run}-REMOVED.tsv \
		-c '(HostEditDist <= 5) and \
			(ViralEditDist <= 5) and \
			(NoAmbiguousBases < 20 or Type == discordant) \
			and (PossibleVectorRearrangement == False) \
			and (PossibleHostTranslocation == False)'
	"""
}


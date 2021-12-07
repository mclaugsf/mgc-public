def helpMessage() {
    log.info"""
    ============================================================================
     :  Git version: ${version}
    ============================================================================
    Usage:

       Mandatory arguments:
           --assembly_ref   assembled references to run through BLAST workflow
           --out_dir        Output Directory (for things that aren't intermediate files to keep)
    Optional:
           --nt_blast_db    Change location of nt BLAST db (DEFAULT: ${params."nt_blast_db"})
           --txids          Filter to this set of taxonomy IDs (DEFAULT ${params."txids"})
           --filter         filter the txids? true/false (DEFAULT ${params."filter"})
           --priority       list of PSP accessions to prioritize

    """.stripIndent()
}

params."nt_blast_db_dir" = '/NGS/blast-nt/'
params."blast_db_name"   = 'nt'

params."assembly_ref" = "/Juicer/shroomapedia/megahit/fasta/*.fasta"
params."txids" = "/home/ubuntu/software/mgc/txids/bacteria-and-fungi.txids"
params."filter" = 'true'

process read_priority_file {
    container 'medicinalgenomics/r-with-bam-vcf'
    memory '500 MB'
    cpus 1

    input:
    file(priority_list) from file(params."priority")

    output:
    stdout into priority_list

    """
    cat ${priority_list}
    """
}

priority_list
   .flatMap {n -> n.split(/\n/).collect()}
   .set{priority_list}

blast_nt_input = Channel.fromPath(params."assembly_ref")
						.map {it -> [it.simpleName, it, params."blast_db_name", file(params."nt_blast_db_dir"), file(params."txids"), params."filter"]}
						.phase(priority_list)
						.map {it -> it.flatten()}
						.map {it -> [it[0], it[1], it[2], it[3], it[4], it[5]]}

process blastn_nt {
    publishDir params."out_dir" + '/' + 'blast-nt/', mode: 'copy', overwrite: true, pattern: "*blast-out.txt"
    container 'ncbi/blast'
    memory '30 G'
    cpus 8

    input:

    set val(sample_name), file(query_fasta), val(blast_db_name), file(blast_db), file(txids), val(filter) from blast_nt_input
    //set val(primer_ref_name), file(query_fasta) from primer_ref1

    output:
    set val(sample_name), file("*blast-out.txt") into blastn_nt_out

    script:
    cpu = task.cpus
    """
    BLASTDB=\"${blast_db}/\";
    export BLASTDB;

    if (${filter} == 'true')
	then
		echo "BLASTing with txids filter using ${txids}...";
	    blastn \
	 		-query ${query_fasta} \
	 		-db ${blast_db}/${blast_db_name} \
	 		-taxidlist ${txids} \
	 		-outfmt "6 qseqid staxids bitscore std" \
	 		-max_target_seqs 1 \
	 		-max_hsps 1 \
	 		-evalue 1e-25 \
	 		-num_threads ${cpu} \
	 		> ${sample_name}.blast-out.txt
	else
		echo "BLASTing to nr-database (no txid filter)...";
		blastn \
	 		-query ${query_fasta} \
	 		-db ${blast_db}/${blast_db_name} \
	 		-outfmt "6 qseqid staxids bitscore std" \
	 		-max_target_seqs 1 \
	 		-max_hsps 1 \
	 		-evalue 1e-25 \
	 		> ${sample_name}.blast-out.txt
	fi
    """
}

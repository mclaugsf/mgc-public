<h2>blobplots using blobtools</h2>

The following code describes how reads were assembled with megahit, mapped back to the megahit assembly, and then used to generate blobplots https://github.com/DRL/blobtools.

Dockerized nextflow workflows are included that should enable anyone to reproduce these results the same way Medicinal Genomics does.
For details on nextflow see:
https://www.nextflow.io/

<h3>Creating reference asssemblies</h3>

Reference assemblies were create dusing megahit https://github.com/voutcn/megahit

the DockerFile used to create the container is here: https://github.com/mclaugsf/mgc-public/blob/master/fungus/Docker/megahit/Dockerfile

The megaghit nextflow workflow can be run like so:

<pre>
nextflow nextflow/megahit.nf --f "*_R{1,2}_001.fastq.gz" --out_dir out
</pre>

This snippet of code (which can also be viewed in the nextflow workflow itself) that calls megahit is here:

<pre>
megahit \
    	-1 ${fq[0]} \
    	-2 ${fq[1]} \
    	-o ${id_run} \
    	-t ${cpu} \
    	-m ${memory} \
      --tmp-dir tmp
</pre>

<h3>Mapping reads back to megahit assembly</h3>

The reads are mapped back to the megahit assembly using bwa http://bio-bwa.sourceforge.net/, samtools http://www.htslib.org/ and sambamba https://lomereiter.github.io/sambamba/

<pre>
bwa mem -R "@RG\\tID:$id_run\\tPU:$id_run\\tSM:$id_run\\tLB:$id_run\\tPL:illumina" \
            -t $cpu -M $fasta_ref ${fq[0]} ${fq[1]} | \
    samtools view -hu - \
            | sambamba sort --tmpdir=. /dev/stdin -o ${id_sample}.bam
</pre>

<h3>Mapping megahit assembly to BLAST database</h3>

The non-redundant BLAST nucleotide database was downloaded:

https://ftp.ncbi.nlm.nih.gov/blast/db/

Specifically from this location:
<pre>
nt*
taxdb*
</pre>

The BLAST command used (which is also included in the nextflow workflow) is the following:

<pre>
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
</pre>

The included nextflow workflow can be run like so:

<pre>
nextflow nextflow/blast-fungus-blobtools-priority-list.nf --assembly_ref PSP*.fasta --out_dir out --nt_blast_db_dir /NGS/blast-nt/ --blast_db_name nt --txids txids/bacteria-and-fungi.txids --filter true --priority accessions-to-process.txt

--nt_blast_db_dir   local path to where nt BLAST database is downloaded to

--assembly_ref would be a list of megahit assemblies named as PSP*.fasta
--priority would be a list of the accesions to process 

</pre>
For example, if you have the following 5 megahit assemblies:
<pre>
PSP10000.fasta
PSP10001.fasta
PSP10002.fasta
PSP10003.fasta
PSP10004.fasta
</pre>

to process all of them, the priority file should have the following lines in it:

<pre>
PSP10000
PSP10001
PSP10002
PSP10003
PSP10004
</pre>

Or leave out ones you don't want to process that might have a corresponding fasta file associated with them.

<h3>Ceating the blobplots</h3>

Finally, to create the blobplots we use these 3 inputs from the previous step:

The nextflow workflow was used to create the blobplots:

<pre>
nextflow nextflow/blobtools-denovo-with-blast-hits.nf --fasta '../megahit/fasta/*.fasta' --bam '../megahit/bam/*.bam' --blast_tsv 'blast-bacteria-fungi/out/blast-nt/*.blast-out.txt' --json psp-all.json --out_dir out
</pre>

The inputs are those generated from the previous steps:

1. --fasta = megahit assemblies from input FASTQ files
2. --bam = BAM files that result from mapping same input FASTQ back to assemblies
3. --blast_tsv = RESULT of the BLAST desribed above
4. --json psp-all.json from this repo - this just matches accessions up to sample names


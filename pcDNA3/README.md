Workflow for detecting integration events using isling: https://github.com/aehrc/isling/tree/master

<h3>Instructions:</h3>
This Nextflow workflow comes with a test data set and can be run as follows:

1. Checkout this version of the code in this repository under the master branch: `git checkout https://github.com/mclaugsf/mgc-public.git`
2. Checkout isling: `git checkout https://github.com/aehrc/isling.git` we used git hash `66f983fe7c8f5f2cb26ccf2cd23c6b3603adcb2f` for processing these data.  The parts of isling that get called here are the two python scripts: `${isling_dir}/scripts/find_ints.py` and `${isling_dir}/scripts/filter.py` and we use their Docker container in the workflow `szsctt/isling:latest` but call the code outside of the docker container.
3. Download hsa1 from UCSC: https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/ and create bwa indexes on it; this can be done w/the accompanying `prepare_reference.nf` Nextflow workflow: `nextflow prepare_reference.nf --r ref.fasta --out_dir .
4. This is written in NextFlow DSL 1.0 and uses an older version of nextflow.  The version we used specifically is `version 21.04.1 build 5556` which is available here: https://github.com/nextflow-io/nextflow/releases/tag/v21.04.1

After hs1 has been downloaded and indexed, the accompanying test dataset can be run through the workflow and takes approximately 1 minute to complete.  NOTE: you will need to have a computer with the necessary resources to run bwa, but you can edit the accompanying nextflow workflow to reduce these resources and give it a try still.

<h4>Running The Test</h4>
```
(base) run-test$ nextflow /home/ubuntu/software/mgc-public/pcDNA3/complete_isling_workflow.nf --vector_ref vector-fasta/pCMV-Spike-sequence-7286-bps.fasta --host_ref ref/hs1.fa --fq '../Test-R{1,2}*.fastq.gz'
N E X T F L O W  ~  version 21.04.1
Launching `/home/ubuntu/software/mgc-public/pcDNA3/complete_isling_workflow.nf` [prickly_mendel] - revision: 9bd2329b56
executor >  local (8)
[f3/23524a] process > bwa_map_vector (bwa-Test-R)                                      [100%] 1 of 1 ✔
[c0/31f622] process > bam_to_fastq (bam_to_fastq-Test-R)                               [100%] 1 of 1 ✔
[6c/c36861] process > bwa_map_host (bwa-Test-R)                                        [100%] 1 of 1 ✔
[01/1e2929] process > host_sort_bam_by_read_name (host_sort_bam_by_read_name-null)     [100%] 1 of 1 ✔
[f0/1eb7bb] process > vector_mapped_reads_only (vector_mapped_reads_only-null)         [100%] 1 of 1 ✔
[fc/cb610d] process > vector_sort_bam_by_read_name (vector_sort_bam_by_read_name-null) [100%] 1 of 1 ✔
[31/f1c66f] process > find_ints (isling_find_ints-null)                                [100%] 1 of 1 ✔
[67/478931] process > filter_find_ints (filter_find_ints-null)                         [100%] 1 of 1 ✔
Completed at: 20-Nov-2024 22:54:09
Duration    : 1m 4s
CPU hours   : 0.1
Succeeded   : 8
```

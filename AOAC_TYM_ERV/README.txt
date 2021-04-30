ubuntu@mgc_pacbio:~/software/mgc-public/AOAC_TYM_ERV$ nextflow bwa-pair.nf --help
N E X T F L O W  ~  version 20.04.1
Launching `bwa-pair.nf` [drunk_darwin] - revision: fd7a960a4d
[help:true]
============================================================================
 :  Git version: 3b574dc37485741a3e0e8821f05c35c0257bbb9c
============================================================================
Usage:

   Mandatory arguments:
       --fq             FastQ1 Pair (i.e. JL.R{1,2}.fastq.gz)
       --r              Fasta Reference
       --out_dir        Output Directory (for things that aren't intermediate files to keep)




ubuntu@mgc_pacbio:~/software/mgc-public/AOAC_TYM_ERV$ nextflow megahit.nf --help
N E X T F L O W  ~  version 20.04.1
Launching `megahit.nf` [goofy_meucci] - revision: e285302108
[help:true]
============================================================================
 :  Git version: 3b574dc37485741a3e0e8821f05c35c0257bbb9c
============================================================================
Usage:

   Mandatory arguments:
       --out_dir        Output Directory (for things that aren't intermediate files to keep)
       --f              FastQ Path (i.e. JL.R1.fastq.gz)

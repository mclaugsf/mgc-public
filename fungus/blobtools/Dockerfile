FROM continuumio/miniconda3
MAINTAINER Steve McLaughlin steve.mclaughlin@medicinalgenomics.com>
RUN conda install -c anaconda matplotlib docopt tqdm wget pyyaml git
RUN conda config --add channels bioconda
RUN conda config --add channels r
RUN conda config --add channels conda-forge
RUN conda install -c bioconda pysam --update-deps
RUN git clone https://github.com/DRL/blobtools.git
WORKDIR blobtools

RUN ./blobtools -h
# RUN ./blobtools create -i example/assembly.fna -b example/mapping_1.bam -t example/blast.out -o example/test
RUN wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -P data/
RUN tar zxf data/taxdump.tar.gz -C data/ nodes.dmp names.dmp
RUN ./blobtools nodesdb --nodes data/nodes.dmp --names data/names.dmp

RUN apt-get update && apt-get install -y --no-install-recommends jq #for passing off a JSON with accessions and sample names to make it easier to id the final plots.

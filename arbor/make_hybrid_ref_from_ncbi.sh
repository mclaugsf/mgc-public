#!/bin/bash

set -euo pipefail

wd=`pwd`

dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

[ -e JL_Mother ] && rm -rf JL_Mother
[ -e JL_Father ] && rm -rf JL_Father

#step1: Download Jamaican Lion Mother from NCBI
echo "Now Downloading JL Mother from NCBI...";
mkdir JL_Mother
cd JL_Mother
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs01/wgs_aux/JA/AT/IP/JAATIP01/JAATIP01.1.fsa_nt.gz
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs01/wgs_aux/JA/AT/IP/JAATIP01/JAATIP01.2.fsa_nt.gz
for i in `ls *gz`; do gunzip -c $i | grep ">"; done  | cut -f 1,10 -d " " | perl -pe "s/\,//" | perl -pe "s/>//" | perl -pe "s/ /\t/" > all_contigs.txt
cd ..

#step2: Download Jamaican Lion Father from NCBI
echo "Now Downloading JL Father from NCBI...";
mkdir JL_Father
cd JL_Father
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs01/wgs_aux/JA/AT/IQ/JAATIQ01/JAATIQ01.1.fsa_nt.gz
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs01/wgs_aux/JA/AT/IQ/JAATIQ01/JAATIQ01.2.fsa_nt.gz
wget https://sra-download.ncbi.nlm.nih.gov/traces/wgs01/wgs_aux/JA/AT/IQ/JAATIQ01/JAATIQ01.3.fsa_nt.gz
for i in `ls *gz`; do gunzip -c $i | grep ">"; done  | cut -f 1,10 -d " " | perl -pe "s/\,//" | perl -pe "s/_/\|/" | perl -pe "s/>//" | perl -pe "s/ /\t/" > all_contigs.txt
cd ..

echo "Now writing JL_Mother reference sequence w/Y Contigs..."
perl ${dir}/make_hybrid_ref_from_ncbi.pl JL_Mother JL_Father > JL_Mother_with_Y_Contigs-NCBI.fasta
echo "Done!"
echo ""
ref=`readlink -e JL_Mother_with_Y_Contigs.fasta.gz`
echo "Output reference: ${ref}"

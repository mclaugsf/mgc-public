#!/usr/bin/perl

use strict;
use warnings;
use FileHandle;

#my $cov_file = '/data/strainseek/Arbor-JL_female+Y-arrow-naming/remap-JL_Mother_with_Y_Contigs-NCBI/out.JL_Mother_with_Y_Contigs-NCBI/Arbor-JL_Mother_with_Y_Contigs-NCBI.bgz';

my $cov_file = $ARGV[0];
my $how_many_samples = $ARGV[1];

my $cmd = 'gunzip -c ' . $cov_file;
my $fh_IN = new FileHandle;
$fh_IN->open("$cmd|");

my $cov_th = 10;

while(<$fh_IN>) {
	chomp $_;
	my $i = $_;
	my @i = split(/\t/, $i);
	my $contig = $i[0];
	my $contig_pos = $i[1];
	my $keep = 0;
	my $total_seen = 0;

	for (my $j = 2; $j < scalar @i; $j++) {
		if ($i[$j] >= $cov_th) {
			$total_seen++;
		}
	}
	if ($total_seen >= $how_many_samples) {
		#print in BED format:
		my $s = $contig_pos - 1;
		my $e = $contig_pos;
		print "$contig\t$s\t$e\n";
	}
}

$fh_IN->close();
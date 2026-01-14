#!/usr/bin/perl

use strict;
use warnings;
use FileHandle;
use Data::Dumper;
use JSON;
use Getopt::Long;
use Pod::Usage;

MAIN_CODE: {
	my %mutation_type = (
    	'A' => { 'G' => 'Ts', 'C' => 'Tv', 'T' => 'Tv' },
    	'G' => { 'A' => 'Ts', 'C' => 'Tv', 'T' => 'Tv' },
    	'C' => { 'T' => 'Ts', 'A' => 'Tv', 'G' => 'Tv' },
    	'T' => { 'C' => 'Ts', 'A' => 'Tv', 'G' => 'Tv' },
	);
	#my $input_vcf = $ARGV[0];
	#my $exit_after = $ARGV[1];
	if(scalar @ARGV == 0) { pod2usage(); }

    my $opt = {'input_vcf'  => undef,
			   'exit_after' => undef,
			   'pass_only'  => undef,
			   'help' => undef};

    GetOptions($opt,
	           'input_vcf=s',
	           'exit_after=i',
	           'pass_only+',
	           'help+') ||  pod2usage();

	if (defined $opt->{'help'}) {
		pod2usage();
	}

	my $input_vcf  = $opt->{'input_vcf'};
	my $exit_after = $opt->{'exit_after'};

	my $bcftools_bin = '/home/ubuntu/software/bcftools-1.18/bcftools-1.18/bcftools';
	my $total_sample = `$bcftools_bin query -l $input_vcf | wc -l`;
	chomp $total_sample;
	if ($total_sample != 1) {
		die "ERROR: $total_sample samples detected this script only works for 1 sample\n"; 
	}
	my $cmd = $bcftools_bin . ' query -f "%CHROM\t%POS\t%REF\t%ALT\t%FILTER\t[%GT\t%AD]\n" ' . $input_vcf;
	my $fh_IN = new FileHandle;
	$fh_IN->open("$cmd|");
	my %counts;
	my $counter = 0;
	my $filter_count = 0;
	my @header = ("CHROM", "POS", "REF", "ALT", "FILTER", "GT", "AD");
	while(<$fh_IN>) {
		$counter++;
		if (defined $exit_after && $counter > $exit_after) {
			print STDERR "Exiting because we processed $exit_after already!\n";
			last;
		}
		chomp $_;
		my $i = $_;
		my @i = split(/\t/, $i);
		my %s;
		@s{@header} = @i;
		if (defined $opt->{'pass_only'}) {
			my $filter = uc $s{'FILTER'};
			if ($filter ne 'PASS') {
				$filter_count++;
				next;
			}
		} 
		my $ref_allele = $s{'REF'};
		my $alt_allele = $s{'ALT'};
		if ($s{'GT'} eq '0/0' || $s{'GT'} eq '0|0') {
			$counts{"skip"}{'ref'}{$s{'GT'}}++;
			next;
		}
		if ($s{'GT'} eq './.' || $s{'GT'} eq '.|.') {
			$counts{"skip"}{'no_call'}{$s{'GT'}}++;
			next;
		}
		$ref_allele = uc $ref_allele;
		$alt_allele = uc $alt_allele;
		my $gt = $s{'GT'};
		$gt =~ s/\//\|/;
		my @gt = split(/\|/, $gt);
		my @alt = split(/\,/, $alt_allele);

		my %seen;
		my $alt_gt_code;
		for my $gt (@gt) {
			if ($gt == 0) {
				next;
			}
			$alt_gt_code = $gt; #get this while we're here - this is the code for the alt allele we need to grab later.
			$seen{$gt}++;
		}
		my $count = scalar keys %seen;
		if ($count != 1) {
			$counts{'skip'}{'het-non-ref'}{$s{'GT'}}++;
			next;
		}
		my $alt = $alt[$alt_gt_code - 1];
		if (! defined $mutation_type{$ref_allele}{$alt}) {
			$counts{'skip'}{'indels'}++;
			next;
		}
		$counts{'Ts_Tv_count'}{$mutation_type{$ref_allele}{$alt}}++;
	}
	$counts{'Ts/Tv'} = $counts{'Ts_Tv_count'}{'Ts'}/$counts{'Ts_Tv_count'}{'Tv'};
	if (defined $exit_after) {
		$counter--;
	}
	$counts{'total_variants_processed'} = $counter;
	$counts{'filter_count'} = $filter_count;
	my $json = encode_json \%counts;
	print($json);
}

=head1 SYNOPSIS

ts_tv_sanity_check-v2.pl -input_vcf 
						 -pass_only *PASS only* 
						 -exit_after *optional i.e. stop after 10000 for testing*

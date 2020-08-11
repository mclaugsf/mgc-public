#!/usr/bin/perl

use strict;
use warnings;
use FileHandle;
use Getopt::Long;
use Pod::Usage;
use Path::Class qw / dir file /;

MAIN_CODE: {

	if(scalar @ARGV == 0) { pod2usage(); }

	my $opt = {'long_format' => undef,
			   'out_dir' => undef,
			   'ref' => undef,
			   'design_vcf' => undef,
			   'strand_flip_tsv' => undef};

	GetOptions($opt,
			   'long_format=s',
			   'out_dir=s',
			   'ref=s',
			   'design_vcf=s',
			   'strand_flip_tsv=s') || pod2usage();

	my $out_dir = dir($opt->{'out_dir'})->absolute();
	$out_dir->mkpath();

	my $rh_design_gt = load_in_design_vcf_genotypes($opt->{'design_vcf'});
	my $rh_strand_flip = load_in_strand_flip_tsv($opt->{'strand_flip_tsv'});
	parse_long_format_file($opt, $rh_strand_flip, $rh_design_gt);
}

sub parse_long_format_file {
	my $opt = shift;
	my $rh_strand_flip = shift;
	my $rh_design_gt = shift;
	my $fh_IN = new FileHandle;
	$fh_IN->open($opt->{'long_format'}) || die "ERROR: Cannot open input Long format file, $opt->{'long'}, for reading\n";
	my $parse_on = 0;
	my @header;
	my $parse_counter = 0;
	my %fh_vcf_OUT;
	while(<$fh_IN>) {
		$_ =~ s/\r\n/\n/g; #strip DOS carriage return
		chomp $_;
		my $i = $_;
		if ($i =~ /\[Data\]/) {
			$parse_on = 1;
			next;
		}
		if ($parse_on == 0) {
			next;
		}
		$parse_counter++;
		if ($parse_counter == 1) {
			@header = split(/\t/, $i);
			next;
		}
		my %gt_call;
		my @i = split(/\t/, $i);
		@gt_call{@header} = @i;
		if (defined $rh_design_gt->{$i[0]}) {
			my $snp_name = $gt_call{'SNP Name'};
			my $design_ref = $rh_design_gt->{$snp_name}{'ref'};
			my $design_alt = $rh_design_gt->{$snp_name}{'alt'};
			my $a1 = $gt_call{'Allele1 - Top'};
			my $a2 = $gt_call{'Allele2 - Top'};
			my $gc_score = $gt_call{'GC Score'};
			my $sample_id = $gt_call{'Sample ID'};

			if (! defined $fh_vcf_OUT{$sample_id}) {
				my $fh_OUT = new FileHandle;
				print STDERR "Now on $sample_id...\n";
				$fh_OUT->open(">$opt->{'out_dir'}/$sample_id-CannSNP.vcf");
				$fh_vcf_OUT{$sample_id} = $fh_OUT;
				write_vcf_header($fh_OUT, $sample_id, $opt->{'ref'});
			}

			my $probe_name = $rh_strand_flip->{$i[0]};
			my @tmp = split(/\-/, $probe_name);
			my $probe_type_id = $tmp[1];

			$probe_type_id =~ s/_(\d+)$//;
			if ($a1 eq '-' && $a2 eq '-') {
				my $gt = './.';
				print_vcf_line(\@i, $gt, $design_ref, $design_alt, $fh_vcf_OUT{$sample_id}, $gc_score);
				next; #no call
			}
			my ($a1_flipped, $a2_flipped) = ($a1, $a2);
			if ($probe_type_id =~ /B_F/ || $probe_type_id =~ /T_R/) {
				$a1_flipped =~ tr/ACGTacgt/TGCAtgca/;
				$a2_flipped =~ tr/ACGTacgt/TGCAtgca/;
			} elsif ($probe_type_id =~ /[MP]_/) {
				my $gt = indel_gt($design_ref, $design_alt, $a1, $a2);
				print_vcf_line(\@i, $gt, $design_ref, $design_alt, $fh_vcf_OUT{$sample_id}, $gc_score);
				next;
			}
			my $gt = snp_gt($a1_flipped, $a2_flipped, $design_ref, $design_alt);
			print_vcf_line(\@i, $gt, $design_ref, $design_alt, $fh_vcf_OUT{$sample_id}, $gc_score);
		}
	}
	$fh_IN->close();
	for my $sample_id (keys %fh_vcf_OUT) {
		my $fh_OUT = $fh_vcf_OUT{$sample_id};
		$fh_OUT->close();
	}
	$fh_IN->close();
}

sub load_in_strand_flip_tsv {
	my $strand_flip_tsv = shift;
	my $fh_IN = new FileHandle;
	$fh_IN->open($strand_flip_tsv) || die "ERROR: Cannot open $strand_flip_tsv as input for reading.\n";
	my $counter = 0;
	my %return_me;
	while(<$fh_IN>) {
		$counter++;
		if ($counter == 1) {
			next;
		}
		chomp $_;
		my $i = $_;
		my @i = split(/\t/, $i);
		my $locus_name = $i[0];
		my $ilmn_id = $i[1];
		$return_me{$locus_name} = $ilmn_id;
	}
	$fh_IN->close();
	return \%return_me;
}

sub print_vcf_line {
	my $ra_i   = shift;
	my $gt     = shift;
	my $ref    = shift;
	my $var    = shift;
	my $fh_OUT = shift;
	my $gc_score = shift;

	my $variant_id = $ra_i->[0];
	$variant_id =~ s/_arrow/\|arrow/;
	my @variant_id = split(/\_/, $variant_id);
	my $vcf_line = $variant_id[0] . "\t" . $variant_id[1] . "\t" . '.' . "\t" . $ref . "\t" . $var . "\t" . '.' . "\t" . 'PASS' . "\t" . '.' . "\t" . 'GT:GC' . "\t" . $gt . ':' . $gc_score;

	print $fh_OUT $vcf_line, "\n";
}

sub load_in_design_vcf_genotypes {
	#---------------------#
	my $design_vcf = shift;
	#---------------------#

	my $cmd = 'bcftools query -f \'%CHROM\_%POS\t%REF\t%ALT\n\'' . ' ' . $design_vcf;

    my $fh_IN = new FileHandle;
    $fh_IN->open("$cmd|");
    my %return_me;
    while(<$fh_IN>) {
    	chomp $_;
    	my $i = $_;
    	my @i = split(/\t/, $i);
    	my $chip_variant_id = $i[0];
    	$chip_variant_id =~ s/\|arrow/_arrow/;
    	$chip_variant_id = $chip_variant_id;
    	$return_me{$chip_variant_id}{'ref'} = $i[1];
    	$return_me{$chip_variant_id}{'alt'} = $i[2];
    }
    $fh_IN->close();
    return \%return_me;
}

sub indel_gt {
	#--------------#
	my $ref = shift;
	my $var = shift;
	my $a1  = shift;
	my $a2  = shift;
	#--------------#
	my $shorter;
	if (length($ref) > length($var)) {
		$shorter = 'var';
	} elsif (length($ref) < length($var)) {
		$shorter = 'ref';
	} else {
		die "Should not see this (it wouild be a SNP).\n";
	}
	my $gt;
	if ($a1 eq $a2) {
		if ($a1 eq 'D') {
			if ($shorter eq 'ref') {
				$gt = '0/0'; #homozygous, reference allele
			} elsif ($shorter eq 'var') {
				$gt = '1/1'; #homozygous, variant allele
			} else {
				die "should not see this\n";
			}
		} elsif ($a1 eq 'I') {
			if ($shorter eq 'ref') {
				$gt = '1/1'; #homozygous, variant allele
			} elsif ($shorter eq 'var') {
				$gt = '0/0'; #homozygous, reference allele
			} else {
				die "should not see this\n";
			}
		}
	} else {
		$gt = '0/1'; #heterozygous
	}
	return $gt;
}

sub snp_gt {
	#---------------------#
	my $a1_flipped = shift;
	my $a2_flipped = shift;
	my $ref = shift;
	my $var = shift;
	my $i = shift;
	#---------------------#

	my $a1_ref_or_var;
	my $a2_ref_or_var;
	#a1:
	my $total_ref = 0;
	my $total_var = 0;

	if ($a1_flipped eq $ref) {
		$a1_ref_or_var = 'ref';
		$total_ref++;
	} elsif ($a1_flipped eq $var) {
		$a1_ref_or_var = 'var';
		$total_var++;
	} else {
		die "ERROR: can't determine if a1 is ref or var\n$i\n";
	}
	#a2:
	if ($a2_flipped eq $ref) {
		$a2_ref_or_var = 'ref';
		$total_ref++;
	} elsif ($a2_flipped eq $var) {
		$a2_ref_or_var = 'var';
		$total_var++;
	} else {
		die "ERROR: can't determine if a2 is ref or var\n$i\n";
	}
	my $gt;
	if ($a1_ref_or_var eq $a2_ref_or_var) {
		if ($a1_ref_or_var eq 'ref') {
			$gt = '0/0';
		} elsif ($a1_ref_or_var eq 'var') {
			$gt = '1/1';
		} else {
			die "should not see this.\n";
		}
	} else {
		$gt = '0/1';
	}
	return $gt;
}

sub write_vcf_header {
	my $fh_OUT = shift;
	my $rsp = shift;
	my $ref = shift;

	my $ref_fai = $ref . '.fai';
	print $fh_OUT "##fileformat=VCFv4.2\n";
	print $fh_OUT '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">', "\n";
	print $fh_OUT '##FORMAT=<ID=GC,Number=G,Type=Float,Description="GenCall score">', "\n";
	my $fh_IN = new FileHandle;
	$fh_IN->open($ref_fai) || die "ERROR: Cannot open input, $ref_fai, for reading.\n";
	while(<$fh_IN>) {
		chomp $_;
		my $i = $_;
		my @i = split(/\t/, $i);
		my $contig_name = $i[0];
		my $contig_length = $i[1];
		printf $fh_OUT "##contig=<ID=$contig_name,length=$contig_length>\n";
	}
    print $fh_OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$rsp\n";
}

=head1 SYNOPSIS

cannsnp90_longformat_to_vcf -long_format -out_dir -ref -design_vcf -strand_flip_tsv

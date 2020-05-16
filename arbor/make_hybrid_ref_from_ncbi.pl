#!/usr/bin/perl

use strict;
use warnings;
use FileHandle;

my $mother_dir = $ARGV[0];
my $father_dir = $ARGV[1];

my @y_contigs = qw(000054F|arrow 000078F|arrow 000113F|arrow 000120F|arrow 000129F|arrow 000134F|arrow 000141F|arrow 000187F|arrow 000188F|arrow 000189F|arrow 000197F|arrow 000219F|arrow 000227F|arrow 000260F|arrow 000263F|arrow 000265F|arrow 000266F|arrow 000271F|arrow 000273F|arrow 000278F|arrow 000279F|arrow 000280F|arrow 000295F|arrow 000306F|arrow 000312F|arrow 000316F|arrow 000339F|arrow 000359F|arrow 000369F|arrow 000373F|arrow 000379F|arrow 000387F|arrow 000396F|arrow 000401F|arrow 000403F|arrow 000407F|arrow 000408F|arrow 000431F|arrow 000435F|arrow 000443F|arrow 000444F|arrow 000447F|arrow 000460F|arrow 000462F|arrow 000469F|arrow 000471F|arrow 000478F|arrow 000481F|arrow 000487F|arrow 000495F|arrow 000507F|arrow 000512F|arrow 000522F|arrow 000527F|arrow 000532F|arrow 000538F|arrow 000543F|arrow 000567F|arrow 000586F|arrow 000588F|arrow 000590F|arrow 000592F|arrow 000596F|arrow 000601F|arrow 000614F|arrow 000645F|arrow 000650F|arrow 000659F|arrow 000661F|arrow 000665F|arrow 000675F|arrow 000677F|arrow 000680F|arrow 000684F|arrow 000685F|arrow 000692F|arrow 000695F|arrow 000710F|arrow 000711F|arrow 000714F|arrow 000719F|arrow 000720F|arrow 000723F|arrow 000727F|arrow 000728F|arrow 000731F|arrow 000733F|arrow 000740F|arrow 000765F|arrow 000769F|arrow 000775F|arrow 000781F|arrow 000785F|arrow 000786F|arrow 000792F|arrow 000795F|arrow 000806F|arrow 000807F|arrow 000808F|arrow 000819F|arrow 000826F|arrow 000829F|arrow 000841F|arrow 000851F|arrow 000861F|arrow 000866F|arrow 000877F|arrow 000879F|arrow 000893F|arrow 000898F|arrow 000916F|arrow 000928F|arrow 000941F|arrow 000942F|arrow 000943F|arrow 000948F|arrow 000953F|arrow 000956F|arrow 000965F|arrow 000966F|arrow 000967F|arrow 000973F|arrow 000976F|arrow 000984F|arrow 000985F|arrow 000990F|arrow 000996F|arrow 001014F|arrow 001033F|arrow 001043F|arrow 001046F|arrow 001049F|arrow 001050F|arrow 001063F|arrow 001079F|arrow 001085F|arrow 001088F|arrow 001090F|arrow 001102F|arrow 001103F|arrow 001116F|arrow 001118F|arrow 001120F|arrow 001122F|arrow 001126F|arrow 001127F|arrow 001130F|arrow 001132F|arrow 001135F|arrow 001136F|arrow 001152F|arrow 001159F|arrow 001161F|arrow 001164F|arrow 001174F|arrow 001184F|arrow 001185F|arrow 001196F|arrow 001198F|arrow 001213F|arrow 001221F|arrow 001223F|arrow 001229F|arrow 001230F|arrow 001231F|arrow 001235F|arrow 001246F|arrow 001252F|arrow 001257F|arrow 001262F|arrow 001263F|arrow 001267F|arrow);

my @x_fasta = glob("$mother_dir/*fsa_nt.gz");
my @y_fasta = glob("$father_dir/*fsa_nt.gz");

#for Mother, keep all:
foreach my $i (@x_fasta) {
	print STDERR "Now reading $i...\n";
	system("gunzip -c $i");
}

my %y_contigs;
foreach my $i (@y_contigs) {
	$i =~ s/^(0+)//;
	$y_contigs{$i}++;
}

my %y_contig_to_ncbi;
my $fh_IN = new FileHandle;
$fh_IN->open("$father_dir/all_contigs.txt") || die "ERROR: Can't open file\n";

while(<$fh_IN>) {
	chomp $_;
	my $i = $_;
	my @i = split(/\t/, $i);
	my $ncbi_accession = $i[0];
	my $y_contig = $i[1];
	$y_contig_to_ncbi{$y_contig} = $ncbi_accession;
}
$fh_IN->close();

#for Father, only keep Y contigs:
my $t = scalar @y_contigs;
my $t_f = 0;
my %y_contig_ncbi;
foreach my $y_contig (@y_contigs) {
	if (defined $y_contig_to_ncbi{$y_contig}) {
		$t_f++;
		$y_contig_ncbi{$y_contig_to_ncbi{$y_contig}}++;
	}
}
if ($t != $t_f) { #they should all be there
	die "ERROR: $t,$t_f\n";
}

foreach my $y_fasta (@y_fasta) {
	print STDERR "Now reading $y_fasta...\n";
	my $cmd = 'gunzip -c ' . $y_fasta;
	my $fh_IN = new FileHandle;
	$fh_IN->open("$cmd|");
	my $print_me = 0;
	while(<$fh_IN>) {
		chomp $_;
		my $i = $_;
		if ($i =~ />/) {
			my @i = split(/ /, $i);
			my $contig_accession = $i[0];
			$contig_accession =~ s/^>//;
			if (defined $y_contig_ncbi{$contig_accession}) {
				$print_me = 1;
			} else {
				$print_me = 0;
			}
		}
		if ($print_me == 1) {
			print "$i\n";
		}
	}
	$fh_IN->close();
}

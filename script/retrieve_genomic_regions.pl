#!/usr/bin/perl

=head1

 Retrieves genomic regions containing matches identifed by run_patscan.pl, with specified number of flanking nucleotides on each side.

=cut

#use strict;
#use warnings;
use FindBin;
use lib ${FindBin::RealBin};
use miRcheck;

my $usage = qq'
perl $0 window_size sample_miRNA_matches  sample_genome.fa sample_miRNA_matches_genomic

* window_size -- number of flanking nucleotides to add to each side of putative miRNA
* sample_miRNA_matches -- list of sequences and genomic coorinates
* sample_genome.fa -- Fasta file of genomic sequence
* sample_miRNA_matches_genomic -- name of outputfile 
';

my ($WIN, $hitfile, $genome, $outfile);

$WIN = $ARGV[0]; # number of flanking nucleotides to add to each side of putative miRNA  
$hitfile =$ARGV[1]; # list of sequences and genomic coorinates
$genome = $ARGV[2]; #Fasta file of genomic sequence
$outfile = $ARGV [3]; #name of outputfile

if (not($outfile)) { print $usage; exit; }

$rseqs = hash_FASTA($genome);	# read fasta sequence file to hash, [for genome]
my %seqs = %{$rseqs};		# key: seq_id; value: sequence

open(O,">$outfile") || die "can not open output file $outfile $!\n";
open(F,$hitfile) || die "can not open sample_miRNA_matches file $hit_file $!\n";
while(<F>)
{
	chomp;
	$line = $_;
    	if ($line =~ m/query:(\S+)\s+qseq:(\S+) hit:(\S+) sense:(\S+) beg:(\d+) end:(\d+) hseq:(\S+)/)
	{
        	$id = $1;
        	$chr = $3;
        	$sense = $4;
        	$x = $5;
        	$y = $6;
        	$hseq = $7;
        	$chr =~ s/>//g;

		if ($sense eq 'sense')
		{
	    		$beg = max(1,$x - $WIN);
	    		$end = min(length($seqs{$chr}),$y + $WIN);
            		$s = '+';
            		$POS = min($WIN,$x-1);
        	}
		else
		{
            		$beg = min(length($seqs{$chr}),$x + $WIN);
            		$end = max(1,$y - $WIN);
            		$s = '-';
            		$POS = min($WIN,(length($seqs{$chr}) - $x)-1);
        	}

        	#retrieve flanking sequence 
		$name = "$id\_$chr\_$x$s";
        	$subseq = get_subseq($seqs{$chr},$beg,$end);
		$start = $POS ;
		$stop = $start + abs($x - $y);
		print O ">$name\t$chr\t$beg\t$end\t$hseq\t$start\t$stop\n$subseq\n";
    	}
	else
	{
		die "Error in line : $line\n";
	}
}
close(F);
close(O);

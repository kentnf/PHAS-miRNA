#!/usr/bin/perl

=head1
 01-14-2014: use strict and warnings, the fix several bugs
 11-10-2013: update mature miRNA to version 20, add function for search against hairpin sequence
 04-14-2012: Change the style for me
 04-13-2012: Add annotation of the program 
=cut

use strict;
use warnings;
use FindBin;
use IO::File;
use File::Path;
use Getopt::Long;

my $usage = qq'
usage: perl $0  [options]

  -s|small-RNA  	small RNA sequence (required)
  -g|genome		genome sequnece & genome index created by bowtie-build (required)
  -m|multi-hit		maximum of multiple hits for one small RNAs (default = 20)
  -b|sql		generate miRNA sql file (default = off)
  -t|target-predict	predict target from mRNA (default = off)
  -u|mRNA		mRNA sequence (required if -t|target-predict is on)
  -d|mRNA-index		mRNA index created by formatdb
  -r|miRNA-star		predict miRNA star (default = off)
  -e|sRNA-expr		expression of small RNA (required if -r|miRNA-star is on)
  -o|output-prefix	the prefix of output file (default = output)
  -h|?|help		print help info

* the input small RNA is : cleaned, uniq, remove rRNA and tRNA
* genome is the genome fasta file, as well as index build by bowtie-build

';

my ($help, $sRNA_seq, $genome_seq, $output, $max_multi_hit, $miRNA_sql_switch, $miRNA_target_switch, $mRNA_seq, $mRNA_seq_index,
    $miRNA_star_switch, $sRNA_expr);

GetOptions(
        "h|?|help"                      => \$help,
        "s|small-RNA=s"                 => \$sRNA_seq,
        "g|genome=s"         		=> \$genome_seq,
        "m|milti-hit=i"              	=> \$max_multi_hit,

        "b|sql"          		=> \$miRNA_sql_switch,
        "t|target-predict"		=> \$miRNA_target_switch,
        "u|mRNA=s"             		=> \$mRNA_seq,
        "d|mRNA-index=s"            	=> \$mRNA_seq_index,

	"r|miRNA-star"			=> \$miRNA_star_switch,
	"e|sRNA-expr=s"			=> \$sRNA_expr,

        "o|output=s"			=> \$output
);

my ($second, $minute, $hour) = localtime();
my $theTime = "$hour:$minute:$second";
print "Start at $theTime\n";

die $usage if $help;
die $usage unless $sRNA_seq;
die $usage unless $genome_seq;
if ($miRNA_target_switch) { die $usage unless $mRNA_seq; die $usage unless $mRNA_seq_index; }
if ($miRNA_star_switch) { die $usage unless $sRNA_expr; $miRNA_sql_switch = 1; }

$max_multi_hit ||= 20;
$output ||= 'output';

#################################################################
# checking input file and scripts				#
#################################################################

# check input file
unless (-s $sRNA_seq) { die "Error, can not find input small RNA sequence file $sRNA_seq \n"; }
unless (-s $genome_seq) { die "Error, can not find genome sequence file $genome_seq\n"; }
my @genome_index = ($genome_seq.".1.ebwt", $genome_seq.".2.ebwt", $genome_seq.".3.ebwt", $genome_seq.".4.ebwt", 
		    $genome_seq.".rev.1.ebwt", $genome_seq.".rev.2.ebwt");

foreach my $idx ( @genome_index ) {
	unless (-s $idx) { die "Error, can not find genome index file $idx\n"; }
}

if ($miRNA_target_switch) {
	unless (-s $mRNA_seq) { die "Error, can not find mRNA sequence file $mRNA_seq\n"; }
	my @mRNA_index = ($mRNA_seq_index.".nhr", $mRNA_seq_index.".nin", $mRNA_seq_index.".nsq");
	foreach my $mRNA_idx ( @mRNA_index ) {
		unless(-s $mRNA_idx) { die "Error, can not find mRNA index file $mRNA_idx\n"; }
	}
}

if ($miRNA_star_switch) {
	unless (-s $sRNA_expr) { die "Error, can not find small RNA expression file $sRNA_expr\n"; }
}

# check the mirBase
my $mirBase = ${FindBin::RealBin}."/mirBase/miRNA";	# mirBase fasta file name.
my @mirBase_index = ($mirBase.".1.ebwt", $mirBase.".2.ebwt", $mirBase.".3.ebwt", $mirBase.".4.ebwt",
		     $mirBase.".rev.1.ebwt", $mirBase.".rev.2.ebwt");

unless (-s $mirBase) { die "Error, can not find the mirBase sequence file $mirBase\n"; }
foreach my $mi_idx (@mirBase_index)
{
	unless (-s $mi_idx) { die "Error, can not find mirBase index file $mi_idx\n"; }
}

my $miRNA_mapping = ${FindBin::RealBin}."/mirBase/miRNA_mapping"; # miRNA mapping file
unless (-s $miRNA_mapping) { die "Error, can not find miRNA mapping file $miRNA_mapping\n"; }

# check output folder and set the output file name
my $output_folder = $output;
unless (-d $output_folder) { mkpath($output_folder, 0, 0777); }

=head1 Description of miRNA output file
                    compare
  small RNA[input] --------- genome[input]
   |                   | (bowtie)
   |             genome_match
   |                   | (filter milti-hit)
   |          genome_match_format
   |                   | (retireve_genomic_regions.pl)
   |               precursor
   |                   | (RNAfold)
   |                 fold
   |                   | (evaluate_miRNA_candidates.pl)
   |           hairpin & mature  ---------
   |              | (parse, RNAfold)     |
   |  hairpin_seq & hairpin_sql          |
   |         | (RNAfold) |               |             compare
   |  hairpin_fold       |    miRNA_temp & miRNA_seq ----------- mirBase
   |                     |         |        |            | (bowtie)
   |                     |         |        |       miRNA_match 
   |                     |         |        |            |
   |                     |         |        |   miRNA_conserved_sql     miRNA_mapping 
   |                     |         |--------+------------|--------------------|
   |                     |                  |            | (miRNA_sql.pl)
   |                     |    mRNA          |        output_sql
   |                     |      |-----------|            |
   |                     |            | (miRNA_target_pred.pl)
   |    		 |          target               |        sRNA_expr
   |                     |                               |          |
   ------------------------------------------------------------------
                                     | (miR_star.pl)
                                    star
=cut
my ($genome_match, $match, $precursor, $uniq_precursor,
    $fold, $eval, $hairpin, $mature,
    $sRNA_anno, $miRNA_match,
    $hairpin_sql, $hairpin_seq , $hairpin_fold,
    $miRNA_temp, $miRNA_seq,
    $miRNA_conserved_sql,
    $miRNA_nomatch, $miRNA_sql_output,
    $miRNA_target_output,
    $miRNA_star_output,
);
my $fd = $output_folder;

$genome_match = $fd."/".$output."_genome_match";	# small RNA aligned to genome
$match = $fd."/".$output."_genome_match_format";	# small RNA aligned to genome with multi-hit < max-multi-hit
$precursor = $fd."/".$output."_precursor";		# precursor candicate sequence retrieved base on match
#$uniq_precursor = $output."_uniq_precursor";		# ??
$fold = $fd."/".$output."_fold";			# precursor fold by RNAfold
$eval =	$fd."/".$output;				# evaluated result
$hairpin = $fd."/".$output."_hairpin";			# evaluated hairpin from fold
$mature = $fd."/".$output."_mature";			# evaluated mature from fold

$hairpin_sql = $fd."/".$output."_miRNA_hairpin_sql";	# hairpin info for sql
$hairpin_seq = $fd."/".$output."_miRNA_hairpin_seq";	# hairpin sequence
$hairpin_fold= $fd."/".$output."_miRNA_hairpin_fold";	# hairpin fold using hairpin sequence by RNAfold

$miRNA_temp = $fd."/".$output."_miRNA_temp";		# identified miRNA info
$miRNA_seq  = $fd."/".$output."_miRNA_seq";		# identified miRNA seq

my $sRNA_sql  = $fd."/".$output."_sRNA_sql";		# small RNA sql info
$sRNA_anno = $fd."/".$output."_sRNA_anno";              # small RNA annotation
                                                		
$miRNA_match = $fd."/".$output."_miRNA_match";		# miRNA align to mirBase using bowtie
$miRNA_conserved_sql = $fd."/".$output."_miRNA_conserved_seq";

$miRNA_nomatch = $fd."/".$output."_miRNA_nomatch";	# no match miRNA to mirBase
$miRNA_sql_output = $fd."/".$output."_miRNA_sql_output";# miRNA sql info
$miRNA_target_output = $fd."/".$output."_miRNA_target_output"; # miRNA target sequence
$miRNA_star_output = $fd."/".$output."_miRNA_star_output";     # miRNA star output

# check program and script
my $bin_dir = ${FindBin::RealBin}."/script";
my ($bowtie_program, $retrieve_genomic_regions_pl, $evaluate_miRNA_candidates_pl, $rnafold_program, 
    $miRNA_sql_pl, $miRNA_target_pred_pl, $miR_star_pl);

$bowtie_program = "$bin_dir/bowtie";
$retrieve_genomic_regions_pl = "$bin_dir/retrieve_genomic_regions.pl";
$evaluate_miRNA_candidates_pl = "$bin_dir/evaluate_miRNA_candidates.pl";
$rnafold_program = "$bin_dir/RNAfold";
$miRNA_sql_pl = "$bin_dir/miRNA_sql.pl";
$miRNA_target_pred_pl = "$bin_dir/miRNA_target_pred.pl";
$miR_star_pl = "$bin_dir/miR_star.pl";

my @script = ($bowtie_program, $retrieve_genomic_regions_pl, $evaluate_miRNA_candidates_pl, $rnafold_program, 
	      $miRNA_sql_pl, $miRNA_target_pred_pl, $miR_star_pl);
foreach my $s ( @script ) { unless (-s $s ) { die "Error, can not find scripts $s\n"; } }

#################################################################
# kentnf: main							#
#################################################################

my %seq; # key : small RNA ID; value: sequence
%seq = seq2hash($sRNA_seq);
compare2genome($sRNA_seq, $genome_seq, $genome_match);

my %num_hit; # key: small RNA ID; value: number of hit for this small RNA
%num_hit = filter_genome_match($genome_match, $match, $max_multi_hit);

#################################################################
# kentnf: subroutine						#
#################################################################
=head1 seq2hash

 function: small RNA sequ to hash

=cut
sub seq2hash
{
	my $sRNA_seq = shift;
	my %seq;
	my ($sID, $sSeq);
	my $sRNA_seq_fh = IO::File->new($sRNA_seq) || die "can't open small RNA sequence file: $sRNA_seq $!\n";
	while(<$sRNA_seq_fh>) {
		chomp;
		if (/>/) { s/>//; $sID = $_; } else { die "Error in read small RNA sequence file: $sRNA_seq\n"; }
		$sSeq = <$sRNA_seq_fh>;
		chomp($sSeq);
		$seq{$sID} = $sSeq;
	}
	$sRNA_seq_fh->close;

	return %seq;
}

=head1 compare2genome

 function: compare miRNA sequence with genome sequence 

=cut
sub compare2genome
{
	print "Compare to genome...";
	my $cmd_compare = "$bowtie_program -v 0 -k 21 -f $genome_seq $sRNA_seq $genome_match";
	print $cmd_compare."\n";
	system($cmd_compare) && die "Error in command $cmd_compare\n";
	print "Done\n";
}

=head1 filter_genome_match

 filter the bowtie align file base on the number of hit

=cut
sub filter_genome_match
{
	my ($genome_match, $match, $cutoff) = @_;
	
	my %num_hit;

	my $genomeB = IO::File->new("$genome_match") || die "can't open bowtie align file: $genome_match $!\n";
	while(<$genomeB>) {
        	my @a = split(/\t/, $_);
        	$num_hit{$a[0]}++;
	}
	$genomeB->close;

	$genomeB = IO::File->new("$genome_match") || die "can't open bowtie align file: $genome_match $!\n";
	my $match_fh = IO::File->new(">$match") || die "can't open filtered bowtie align file: $match $!\n";

	my ($length, $sense, $start, $end);
	while(<$genomeB>) 
	{
		my @a = split "\t";
		if ($num_hit{$a[0]} < $cutoff)  # only keep those with less than or equal to 20 hits
		{
			$length = length($seq{$a[0]});
			if ($a[1] eq "-") {
				$sense = "antisense";
				$end = $a[3] + 1;
				$start = $a[3] + $length;
			}
			else {
				$sense = "sense";
				$start = $a[3] + 1;
				$end = $a[3] + $length;
			}
			print $match_fh "query:$a[0] qseq:$seq{$a[0]} hit:$a[2] sense:$sense beg:$start end:$end hseq:$seq{$a[0]}\n";
		}
	}
	$match_fh->close;
	$genomeB->close;

	return %num_hit;
}
#unlink($genome_match);


#################################################################
# retrieve genomic regions base sRNA aligned info		#
# the retrieved sequence is the candidate precursor sequence	#
# Input: the filtered matched small RNA to genome		#
# Output: the precursor sequence for each small RNA		#
#################################################################
my $cmd_retrieve_genomic_regions = "$retrieve_genomic_regions_pl 200 $match $genome_seq $precursor";
system ($cmd_retrieve_genomic_regions) && die "Error in command: $cmd_retrieve_genomic_regions\n";
#unlink($match);

#################################################################
# fold the candidate precursor sequences using RNAfold program	#
# Input: the precursor candidate sequences			#
# Output: the fold info for each precursor candidate sequences	#	
#################################################################
my $cmd_RNAfold = "$rnafold_program -noPS < $precursor > $fold";
system($cmd_RNAfold) && die "Error in command: $cmd_RNAfold\n";

#################################################################
# evaluate the miRNA candidates base one the fold info		#
# Input: RNAfold results					#
# Output: hairpin and mature file				#
#################################################################
my $cmd_evaluate_miRNA = "$evaluate_miRNA_candidates_pl $fold $eval";
system($cmd_evaluate_miRNA) && die "Error in command: $cmd_evaluate_miRNA\n";

#################################################################
# parse the haipin file 					#
# 1. generate hairpin fasta file (miRNA_hairpin_seq)		#
# 2. generate hairpin sql for db (miRNA_hairpin_sql)		#
#################################################################

my %miR_ID;	# uniq miR_ID in hairpin sequence;
		# key: ID of input sequences; value: new ID 

%miR_ID = parse_hairpin($hairpin, $hairpin_seq, $hairpin_sql);

sub parse_hairpin
{
	my ($hairpin, $hairpin_seq, $hairpin_sql) = @_;

	my %miR_ID;

	# the id of hairpin sequence  -- >smallRNAID_Chr1_13252396+       Chr1    13252204        13252419
	open(HAIRPIN, $hairpin) || die "can't open hairpin $hairpin $!\n";
	open(HAIRPIN_SQL, ">$hairpin_sql") || die "can't open hairpin sql $hairpin_sql $!\n";
	open(HAIRPIN_SEQ, ">$hairpin_seq") || die "can't open hairpin seq $hairpin_seq $!\n";

	my ($i, $mID) = (0, 0);

	my ($hairpin_ID, $miR_ID, $start, $strand);

	my %hairpin_sql;

	while(<HAIRPIN>) 
	{
		chomp;
		if (/>/) {
			s/>//;
			# assign new ID to hairpin
			$i++;
			my $x = 6-length($i);
			my $zero = "0"x$x;
		 	$hairpin_ID = "H$zero$i";
			my @a = split "\t";
			$strand = substr($a[0], -1, 1);
			my @b = split("_", $a[0]);
			chop($b[@b -1]);
			if ($strand eq "+") { $start = $b[@b - 1] - $a[2] + 1; }
			else { $start = $a[2] - $b[@b - 1] + 1; }

			# assign new ID for each uniq miR ID in hairpin sequences
			if (!$miR_ID{$b[0]}) {
				$mID++;
				my $x = 5-length($mID);
				my $zero = "0"x$x;
				$miR_ID = "M$zero$mID";
				$miR_ID{$b[0]} = $miR_ID;
			}

			$hairpin_sql{$hairpin_ID} = "$hairpin_ID\t$miR_ID{$b[0]}\t$a[1]\t$a[2]\t$a[3]\t$start\t$strand\t";
		}
		else {
			$hairpin_sql{$hairpin_ID} .= $_;
			print HAIRPIN_SEQ ">$hairpin_ID\n$_\n";
		}
	}
	close(HAIRPIN);
	close(HAIRPIN_SEQ);

	#unlink("output/$hairpin");

	#################################################################
	# generate the fold PS file for each hairpin sequence		#
	# then put the score to hairpin sql				#
	#################################################################

	mkpath("$fd/hairpin", 0, 0777);
	system "$rnafold_program < $hairpin_seq > $hairpin_fold";

	system "mv *.ps $fd/hairpin/";

	open(HAIRPIN_FOLD, "$hairpin_fold") || die "can't open hairpin fold : $hairpin_fold $!\n";
	while(<HAIRPIN_FOLD>) {
		chomp;
		if (/>/) {
			s/>//;
			if ($hairpin_sql{$_}) { print HAIRPIN_SQL $hairpin_sql{$_}; }
		}
		if (/\(/) {
			s/\(/#/g;
			s/\)//g;
			my @a = split "#";
			print HAIRPIN_SQL "\t$a[@a-1]\n";
		}
	}
	close(HAIRPIN_SQL);

	return %miR_ID;
}


=head
unlink("output/$output");
unlink("output/$mature");
unlink("output/$fold");
unlink("output/$hairpin_seq");
unlink("output/$genome_match");
unlink("output/$hairpin_fold");
unlink("output/$uniq_hairpin");
unlink("output/$precursor");
unlink("output/$match");
=cut

my ($miR_seq, $anno) = generate_miRNA(\%miR_ID, $miRNA_temp, $miRNA_seq);
my %miR_seq = %$miR_seq;

generate_sRNA_sql($sRNA_sql, $anno);

=head1 generate miRNA sequence file (miRNA_seq)

 the miRNA_seq was the uniq small RNA sequence that identified 
 with hairpin in its precursor sequences

=cut
sub generate_miRNA
{
	my ($miR_ID, $miRNA_temp, $miRNA_seq) = @_;
	my %miR_ID = %$miR_ID;
	my %anno;
	my %miR_seq;

	open(miR_temp, ">$miRNA_temp") || die "can't open $miRNA_temp $!\n";
	open(miR_Seq,  ">$miRNA_seq")  || die "can't open $miRNA_seq $!\n";
	foreach my $key (sort keys %miR_ID) {
		print miR_temp "$miR_ID{$key}\t$key\t$seq{$key}\n";
		print miR_Seq ">$miR_ID{$key}\n$seq{$key}\n";
		$miR_seq{$miR_ID{$key}} = $seq{$key};
		$anno{$key} = "miRNA";
	}
	close(miR_temp);
	close(miR_Seq);

	return (\%miR_seq, \%anno);
}

=head1 generate sql for all input small RNA

 Format of sRNA_sql
 ID	    Sequence		       Length  Anno		
 Test00003  AAAGGCCGAAAAAAATATGCCCGG   24      miRNA		

=cut
sub generate_sRNA_sql
{
	my ($sRNA_sql, $anno) = @_; 
	my %anno = %$anno;
	open(sRNA_SQL, ">$sRNA_sql") || die "can't open small RNA seq $sRNA_sql $!\n";
	foreach my $key (sort keys %seq) {
		if (!$anno{$key}) { $anno{$key} = "NA"; }
		print sRNA_SQL "$key\t$seq{$key}\t", length($seq{$key}), "\t$anno{$key}\n";
	}
	close(sRNA_SQL);
}

#################################################################
# compare with mirBase using miRNA_seq				#
#################################################################
print "Compare to miRNA database (mature) ...\n";
my $cmd_mirBase = "$bowtie_program -v 2 -a --best --strata -f --norc $mirBase $miRNA_seq $miRNA_match"; 
print $cmd_mirBase."\n";
system($cmd_mirBase) && die "Error in command: $cmd_mirBase\n";
print "Done\n";

# mirBase sequence to hash
my $miRNA_sequence;
my %miRNA_seq;
open(MIRNA, $mirBase) || die "can't open mirBase file $mirBase $!\n";
while(<MIRNA>) { $miRNA_sequence .= $_; }
close(MIRNA);

my @miRNA_seq = split(">", $miRNA_sequence);

for (my $i = 0; $i < @miRNA_seq; $i++) {
	$miRNA_seq[$i] =~ s/\n/#/;
	my @b = split("#", $miRNA_seq[$i]);
	chomp($b[1]);
	$miRNA_seq{$b[0]} = $b[1];
}

=head1 

Input : output/miRNA_match
Output: output/miRNA_conserved_sql

=cut
my ($query_sequence, $hit_sequence, $start_pos, $end_pos, $start_position, $end_position, $to_end, $start_space, $end_space, $start, $end, $input_seq, $count_mismatch);

open(MIRNA_MATCH, $miRNA_match)  || die "can't open miRNA match to mirBase file $miRNA_match $!\n";
open(MIR_conserve, ">$miRNA_conserved_sql")  || die "can't open $miRNA_conserved_sql $!\n";
while(<MIRNA_MATCH>) 
{
	my @a = split "\t";
	my $align;
	if ($a[1] eq "+") {
		$query_sequence = $miR_seq{$a[0]};
		$hit_sequence = $miRNA_seq{$a[2]};
		$start_pos = $a[3];
		$end_pos = length($query_sequence) - 1 + $a[3];
		$start_position = $start_pos + 1;
		$end_position = $end_pos + 1;
		$to_end = length($hit_sequence) - $end_pos - 1;
		$start_space = "x" x $start_pos;
		$end_space = "x" x $to_end;
		$start = "-" x $start_pos;
		$end = "-" x $to_end;
		$input_seq = $start . $query_sequence . $end;
		for (my $k = $start_pos; $k <= $end_pos; $k++) {
			if (substr($input_seq, $k, 1) eq substr($hit_sequence, $k, 1)) { $align .= "|"; }
			else { $align .= "x"; }
		}
		$align = $start_space . $align . $end_space;
		$count_mismatch = $align =~ tr/x/x/;
		print MIR_conserve "$a[0]\t$a[2]\t$input_seq\t$align\t$hit_sequence\n";
	}
}
close(MIRNA_MATCH);
close(MIR_conserve);
#unlink("$miRNA_match");

#################################################################
# compare hairpin sequence with haipin in mirBase               #
#################################################################

#print "Compare to miRNA database (hairpin) ....\n";
#my $cmd_hairpin "$bowtie_program -v 2 -a --best --strata -f --norc $mirBase_hairpin $haipin_seq $haipin_match";
#print $cmd_hairpin."\n";
# system($cmd_hairpin) && die "Error in command: $cmd_hairpin\n";
#print "Done\n";

=head1 miRNA_sql.pl

 function: parse miRNA prediction result to sql info

=cut
if ($miRNA_sql_switch)
{
	#usage: perl miRNA_sql.pl miRNA_conserved_sql miRNA_temp miRNA_mapping miRNA_nomatch output
	my $cmd_miRNA_sql = "$miRNA_sql_pl $miRNA_conserved_sql $miRNA_temp $miRNA_mapping $miRNA_nomatch $miRNA_sql_output";
	print $cmd_miRNA_sql."\n";
	system($cmd_miRNA_sql) && die "Error in command $cmd_miRNA_sql\n";
}

=head1  miRNA_target_pred.pl

Function: predict miRNA target sequence from mRNA/unigene sequence

Usage: perl miRNA_target_pred.pl miRNA_seq_fasta blast_formatted_db source_db output

* miRNA_seq_fasta -- miRNA sequence file identified by miRNA_prediction.pl
* blast_formatted_db -- mRNA DB formatted by blast
* source_db -- mRNA sequence file
* output -- output result (filtered) 

=cut
if ($miRNA_target_switch)
{
	my $cmd_miRNA_target = "$miRNA_target_pred_pl $miRNA_seq $mRNA_seq_index $mRNA_seq $miRNA_target_output";
	print $cmd_miRNA_target."\n";
	system($cmd_miRNA_target) && die "Error in command $cmd_miRNA_target\n";
}

=head1 miRNA_star


=cut
if ($miRNA_star_switch)
{
	my $cmd_miRNA_star = "$miR_star_pl -a $miRNA_sql_output -b $hairpin_sql -d $sRNA_seq -e $sRNA_expr -f $miRNA_star_output";
	print $cmd_miRNA_star."\n";
	system($cmd_miRNA_star) && die "Error in command $cmd_miRNA_star\n";
}

($second, $minute, $hour) = localtime();
$theTime = "$hour:$minute:$second";
print "Finished at $theTime\n";


#!/usr/bin/env perl


# ==========================================================================
# smRNAtarget_Pred.pl: small RNA target prediction tool.
#
# This program returns potential small RNAs for an input 
# gene as well as potential sites for an input small RNA.  
#
# Usage: smRNAtarget_Pred.pl miRNA_seq_fasta blast_formatted_db source_db output
#
# ==========================================================================


###########################################################
# Configuration section

# The path of BLAST tool
$BLAST_PATH = "/usr/local/bin";

# The path of database files 
# 	- Save BLAST formatted database files of mRNAs and small RNAs 
# 	- In case of small RNAs, please make format database files 
#         with their reverse complement

$DATA_PATH = "./";

# The path for temorary file
$TEMP_PATH = "./";

# Default cutoff of target score ( score range: 0 < score < 10 )
$sc_cutoff = 3; 

# Cutoff value of alignment by dynamic programming
$max_score_cutoff = 8; 

# Type of query sequence (smRNA: 0, mRNA: 1)
$seq_type = 0; 

# Direction of mRNA sequence (forward, reverse, both)
$choice_direction = "forward";

###########################################################
my @sel_list_disc = ();
my @sel_list_pos_mRNA_st = ();
my @sel_list_pos_sRNA_st = ();
my @sel_list_pos_sRNA_end = ();
my @sel_list_direction = ();
my @description_list = ();

$blastdb_smRNA = "";	#  small RNA DB
$blastdb_mRNA = "";	#  mRNA DB / unigene DB 

$mRNA_fas = "";		#  mRNA sequence file / unigene sequence file
$smRNA_fas = "";	#  small RNA sequence file

# Temorary file for blast input 
$tempb = $TEMP_PATH."/"."blast.tmp";


=head1 insert_seq

 put fasta sequence to hash. 
 key: seqID ; value: sequence

=cut
sub insert_seq ($) {
	my $fn = $_[0];
	open(INFILE, $fn) or die("Failed to read in $fn");

	my $seq_hash = {};
	my $seq = "";
	my $description;

	$count = 1;
	while($line = <INFILE>) {
		if ($line =~ m/^\>/) {
			if ($count > 1) { 
				$seq_hash{$description} = $seq;	
				$seq = "";	
			}
			$description = split_discription($line);
			chomp $description;
		}
		else {
			chomp  $line;
			$seq .=  $line;
		}
		$count++;
	}
	close(INFILE);
	
	$seq_hash{$description} = $seq;
	return %seq_hash;
}

=head1 query_seq_from_file

 put fasta sequence to hash, 
 key: seqID ; value: sequence 

=cut
sub query_seq_from_file ($) {
	my $fn = $_[0];
	open(INFILE, $fn) or die("Failed to read in $fn");

	my $qseq_hash = {};
	my $seq = "";
	my $description;

	$count = 1;
	while($line = <INFILE>) {
		if ($line =~ m/^\>/) {
			if ($count > 1) { 
				$qseq_hash{$description} = $seq;	
				push @description_list, $description;
				$seq = "";	
			}
			$description = split_discription($line);
			chomp $description;
		}
		else {
			chomp  $line;
			$seq .=  $line;
		}
		$count++;
	}
	close(INFILE);

	chomp($description);
	$seq =~ tr/atcg/ATCG/;
	
	$qseq_hash{$description} = $seq;
	push @description_list, $description;
	return (%qseq_hash);
}

=head1 reverse_complement

 function: get the reverse complement sequence

=cut
sub reverse_complement {
	my $seq = shift;
	$seq =~ s/U/T/g;
	$seq =~ s/u/T/g;
	$seq = reverse($seq);
	$seq =~ tr/atcgATCG/tagcTAGC/;
	return $seq;
}

=head1 lsa (Local Sequence Alignment)

 perform local sequence alignment base on sw-algorithm

=cut
sub lsa {
	my ($adisc, $bdisc, $mRseq, $smRseq, $pos, $direction) = @_;


	my $match = 1;
	my $wobble = 0.5; 
	my $mismatch = -1;
	my $indel = -2;

	my $first_i = 0;
	my $first_j = 0;
	my $last_i = 0;
	my $last_j = 0;

	my $offset = 3;
	my $mRNA_st = $sel_list_pos_mRNA_st[$pos]+0;
	my $sRNA_st = $sel_list_pos_sRNA_st[$pos]+0;
	my $sRNA_end = $sel_list_pos_sRNA_end[$pos]+0;
	my $start = $mRNA_st - ($sRNA_st + $offset);
	if ($direction eq "-") {
		$start = $mRNA_st - ((length($smRseq)-$sRNA_end) + $offset);
	}
	if ($start < 0) { $start = 0; }
	my $seqlen = length($smRseq)+$offset;
	my $end = $start+$seqlen;
	if ($end >= length($mRseq)) { $seqlen = length($mRseq)-$start+1; }

	my $a = substr($mRseq, $start, $seqlen);
	my $b = $smRseq;

	my $S = [];   # An array of scores
	my $R = [];   # An array of backtracking arrows
	my $n = length($a);
	my $m = length($b);

	# We need to work in letters, not in strings.  This is a simple way
	# to turn a string of letters into an array of letters.
	my @a = split // => $a;
	my @b = split // => $b;

	# These are "constants" which indicate a direction in the backtracking array.
	my $UP_AND_LEFT = "\\";
	my $UP          = "|";
	my $LEFT        = "-";
	my $NEITHER     = "o";

	my $max_score = 0;
	my $sc = 0.0;
	my $mc = 0;
	my $wc = 0;
	my $ic = 0;

	my $answer = 0;
	my $bind = "";
	
	# Initialization
	my $score = 0;
	for(my $j = 0; $j <= $m; $j++) { 
		$S->[0][$j] = 0;
		$R->[0][$j] = $NEITHER;
	}

	# This is the main dynamic programming loop that computes the score.
	for(my $i = 1; $i <= $n; $i++) { 
		$S->[$i][0] = 0;
		$R->[$i][0] = $NEITHER;
		for(my $j = 1; $j <= $m; $j++) { 
			$S->[$i][$j] = 0;
			$R->[$i][$j] = $NEITHER;

			if ( $a[$i-1] eq $b[$j-1] ) { $score = $match; }
			elsif ( is_wobble($a[$i-1], $b[$j-1], $direction) eq 1 ) { $score = $wobble; }
			else { $score = $mismatch; }

			if (($S->[$i-1][$j-1]+$score) > $S->[$i][$j]) {
				$S->[$i][$j] = $S->[$i-1][$j-1]+$score;
				$R->[$i][$j] =$UP_AND_LEFT;
			}
			if (($S->[$i][$j-1]+$indel) > $S->[$i][$j]) {
				$S->[$i][$j] = $S->[$i][$j-1]+$indel;
				$R->[$i][$j] =$LEFT;
			}
			if (($S->[$i-1][$j]+$indel) > $S->[$i][$j]) {
				$S->[$i][$j] = $S->[$i-1][$j]+$indel;
				$R->[$i][$j] =$UP;
			}
			if ($S->[$i][$j] > $max_score) {
				$max_score = $S->[$i][$j];
				$last_i = $i;
				$last_j = $j;
			}
		}
	}

	my $seqa = "";
	my $seqb = "";

	if ( $max_score > $max_score_cutoff ) {

		$i = $last_i; 
		$j = $last_j;

		# Trace the backtracking matrix.
		while( $i > 0 and $j > 0 ) {

			if( ($R->[$i][$j] eq $UP_AND_LEFT) and (($a[$i-1] eq $b[$j-1]) or (is_wobble($a[$i-1], $b[$j-1], $direction ) eq 1 )) ) {

				$seqa = $a[$i-1].$seqa;
				$seqb = $b[$j-1].$seqb;
				$i--; $j--;
			}
			elsif( $R->[$i][$j] eq $LEFT ) {
				$seqa = "-".$seqa;
				$seqb = $b[$j-1].$seqb;
				$j--;
			}
			elsif( $R->[$i][$j] eq $UP ) {
				$seqa = $a[$i-1].$seqa;
				$seqb = "-".$seqb;
				$i--;
			}
			else {
				$seqa = $a[$i-1].$seqa;
				$seqb = $b[$j-1].$seqb;
				$i--; $j--;
			}
		}

		$first_i = $i;
		$first_j = $j;

		# Insert postfixes of seqb 
		for(my $j = 1; $j <= length($smRseq)-$last_j; $j++) { 
			$seqb = $seqb.$b[$last_j+$j-1];
		}

		my $lp = get_last_pos($first_i, $seqa);

		# Insert prefixes 
		my $k = 1;
		my $ind = 0;
		for(my $j = $first_j-1; $j >= 0; $j--) { 
				if (($first_i-$k) >= 0) {
					$seqa = $a[$first_i-$k-1].$seqa;
				}	
				else {
					$seqa = "-".$seqa;
					$ind = 1;
				}
				$seqb = $b[$j].$seqb;
				$k++;
		}
		if ($ind eq 1) { $first_i = 0; }
		else { $first_i = $first_i-$k+1; }
		$first_j = $first_j-$k-1;

		# Delete postfix indels
		for(my $i = length($seqa)-1; $i >= 0; $i--) { 
			if ( substr($seqa,$i,1) eq "-" ) { 
				substr($seqa,$i,1) = "";
		 	}
			else { last; }
		}
		chomp($seqa);

		# Insert postfixes
		$k = 1;
		for(my $i = length($seqa); $i < length($seqb); $i++) { 
				if (($lp+$k) < length($a) )    {
					$seqa = $seqa.$a[$lp+$k];
				}
				else {
					$seqa = $seqa."-";
				}
				$k++;
		}

		# Scoring of binding pair
		my $ca = "";
		my $cb = "";
		for(my $i = 0; $i < length($seqb); $i++) { 
			$ca = substr($seqa,$i,1);
			$cb = substr($seqb,$i,1);
			if ( is_match($ca, $cb) eq 1 ) { $bind = $bind."|"; }
			elsif ( is_wobble($ca, $cb, $direction) eq 1 ) { $wc++; $sc=$sc+0.5; $bind = $bind."o"; }
			elsif ( is_indel($ca, $cb) eq 1 ) { $ic++; $sc=$sc+2.0; $bind = $bind."b"}
			else  { $sc=$sc+1.0;  $mc++; $bind = $bind."b"; }
		}
		if ( $sc <= $sc_cutoff ) { $answer = 1; }
	}

	@seqinfo = ($answer, $sc, $mc, $start+$first_i+1, $first_j+1, $seqa, $seqb, $bind, $adisc, $bdisc, $direction, $wc, $ic);
	return @seqinfo;
}

=head1 get_last_pos

 function: get the absolute end position of sequence
 Input: start position, sequence
 Output: end position

=cut
sub get_last_pos (@) {
        my ($start, $seq) = @_;

        my $barc = 0;
	for(my $i = 0; $i < length($seq); $i++) { 
		if (substr($seq,$i,1) eq "-") { $barc++; }
	}
        my $last_pos = $start+length($seq)-$barc-1;
        return $last_pos;
}

=head1 is_match

 function: check if the two base could be match
 Input: baseA and baseB
 Output: 1 (match) or 0 (not match)

=cut
sub is_match (@) {
	my ($ca, $cb) = @_;
	my $answer = 0;
	if ( $ca eq $cb ) { $answer = 1; }
	return $answer;
}

=head1 is_wobble

 function: check if the two base could be wobble
 Input: baseA and baseB
 Output: 1 (wobble) or 0 (not wobble)

=cut
sub is_wobble (@) {
	my ($ca, $cb, $strand) = @_;
	my $answer = 0;
	if ( $strand eq $pos_direction ) {
		if ( (($ca eq "T") and ($cb eq "C")) or (($ca eq "G") and ($cb eq "A" )) ) {
			$answer = 1;
		}
	}
	else {
		if ( (($ca eq "C") and ($cb eq "T" )) or (($ca eq "A") and ($cb eq "G")) ) {
			$answer = 1;
		}
	}
	return $answer;
}

=head1 is_indel

 function: check if the two base could be indel
 Input: baseA and baseB
 Output: 1 (indel) or 0 (not indel)

=cut
sub is_indel (@) {
	my ($ca, $cb) = @_;
	my $answer = 0;
	if ( (($ca eq "-") and ($cb ne "-")) or (($ca ne "-") and ($cb eq "-" )) ) {
		$answer = 1;
	}
	return $answer;
}

=head1 is_mismatch

 function: check if the two base could be mismatch
 Input: baseA and baseB
 Output: 1 (mismatch) or 0 (not mismatch)

=cut
sub is_mismatch (@) {
	my ($ca, $cb) = @_;
	my $answer = 0;
	if ( $ca ne $cb ) { $answer = 1; }
	return $answer;
}

=head1 change_sequence

 function: ??????

=cut
sub change_sequences (@) {
	my ($sad, $sbd, $sa, $sb, $is_mRNA) = @_;
	return 	($sad, $sbd, $sa, $sb);
}

=head1 split_discription

 get the sequence ID  

=cut
sub split_discription ($) {
	my $sd = $_[0];
	my @list = split(" ", $sd);
	my $title = $list[0];
	return $title;
}

sub blast_input_save (@) {
	my ($disc, $qseq, $is_mRNA) = @_;

	my $answer = 1;

        if ( $is_mRNA eq 0 ) {
                $qseq = reverse_complement($qseq);
        }

	chomp($disc);   				
	chomp($qseq);   				
	my $result_save = join "\n", $disc, $qseq;

	open OUTFILE, ">$tempb";
	print OUTFILE "$result_save";
	close OUTFILE;

	return $answer;
}

=head1 blast_run

 run blastp, blast  

=cut
sub blast_run (@) {
	my ($blastdb, $queryfile, $direction) = @_;	

	my $blastdbname = $DATA_PATH."/".$blastdb;
	
	my $blast_out = `$BLAST_PATH/blastall -p blastn -d $blastdbname -i $queryfile -W 7 -q -1 -S $direction -e 100 -m 8`;

	#my $blast_out = `$BLAST_PATH/blastall -p blastn -d $blastdbname -i $queryfile -W 7 -q -1 -S $direction -e 200 -m 8`;
	
	return $blast_out;
}

sub blast (@) {
	my ($query_disc, $query_seq, $queryf, $is_mRNA, $search_direction) = @_;
	my $i = 0;
	my $db = "";
	my $forward = 1;
	my $reverse = 2;
	$pos_direction = "+"; 
	$neg_direction = "-"; 

	if ( $is_mRNA eq 1 ) { 
		$db = $blastdb_smRNA;
		$rv = blast_input_save($query_disc, $query_seq, $is_mRNA);

		if ( ($search_direction eq "forward") or ($search_direction eq "both") ) { 
			$result = blast_run($db, $queryf, $forward);

			my @list = split(/\n/, $result);
			foreach my $line (@list) {
				@entry = split(/\t/,$line);
				$name = ">".$entry[1];
				push @sel_list_disc, $name;
				push @sel_list_pos_mRNA_st, $entry[6];
				push @sel_list_pos_sRNA_st, $entry[8];
				push @sel_list_pos_sRNA_end, $entry[9];
				push @sel_list_direction, $pos_direction;
				$i++;
			}
		}
		if ( ($search_direction eq "reverse") or ($search_direction eq "both") ) { 

			$result = blast_run($db, $queryf, $reverse);

			@list = split(/\n/, $result);
			foreach my $line (@list) {
				@entry = split(/\t/,$line);
				$name = ">".$entry[1];
				push @sel_list_disc, $name;
				push @sel_list_pos_mRNA_st, $entry[6];
				push @f, $entry[9];
				push @sel_list_pos_sRNA_end, $entry[8];
				push @sel_list_direction, $neg_direction;
				$i++;
			}
		}
	}
	else { 
		$db = $blastdb_mRNA;
		$rv = blast_input_save($query_disc, $query_seq, $is_mRNA);
		if ( ($search_direction eq "forward") or ($search_direction eq "both") ) { 

			$result = blast_run($db, $queryf, $forward);

			my @list = split(/\n/, $result);
			foreach my $line (@list) {
				@entry = split(/\t/,$line);
				$name = ">".$entry[1];
				push @sel_list_disc, $name;
				push @sel_list_pos_mRNA_st, $entry[8];
				push @sel_list_pos_sRNA_st, $entry[6];
				push @sel_list_pos_sRNA_end, $entry[7];
				push @sel_list_direction, $pos_direction;
				$i++;
			}		
		}

		if ( ($search_direction eq "reverse") or ($search_direction eq "both") ) { 

			$result = blast_run($db, $queryf, $reverse);

			@list = split(/\n/, $result);
			foreach my $line (@list) {
				@entry = split(/\t/,$line);
				$name = ">".$entry[1];
				push @sel_list_disc, $name;
				push @sel_list_pos_mRNA_st, $entry[9];
				push @sel_list_pos_sRNA_st, $entry[6];
				push @sel_list_pos_sRNA_end, $entry[7];
				push @sel_list_direction, $neg_direction;
				$i++;
			}
		}
	}

	return $i;
}

sub print_tab {
	my $schash, %resulthash =  @_;

	my $summary = "";
	my $i = 1;
	foreach $value (sort {$schash{$a} <=> $schash{$b}} keys %schash) {
		my ($istarget, $sc, $mc, $sa, $sb, $seqa, $seqb, $bind, $adisc, $bdisc, $direction, $wc, $ic) = @{$resulthash{$value}};
		my  $lp = get_last_pos($sa, $seqa);

		$bdisc =~ s/^\>//;
		$adisc =~ s/^\>//;
		my  $sRNA_ID = $bdisc;
		my  $unigene_ID = $adisc;
		
		if ( $direction eq "-") {
			$seqa = reverse_complement($seqa);
			$seqb = reverse($seqb);
        	$seqa =~ tr/tT/uU/;
        	$seqb =~ tr/tT/uU/;
			$bind = reverse($bind);
			my $tmp_pos = $sa;
			$sa =$lp;
			$lp = $tmp_pos;
		}
		else {
			$seqa =~ tr/tT/uU/;
			$seqb =~ tr/atcgATCG/uagcUAGC/;
		}
		my  $target_start = $sa;
		my  $align_string = $seqb;
		my  $align_score = $sc;
		my  $align_query = $seqa;
		my  $align_hit = $bind;
		my  $seq_strand = $direction;
		
		$summary .= "$sRNA_ID\t$unigene_ID\t$target_start\t$lp\t$align_score\t$align_query\t$align_string\t$align_hit\t$seq_strand\n";		

		$i++;
	}
	return $summary;
}

=head1 filter_output_result

=cut
sub filter_output_result
{
	my ($result_tab, $max_mismatch_cutoff, $not_mismatch_position) = @_;

	my @result_tabs = split(/\n/, $result_tab);
	my $result_tab_filtered = "";
	my @position = split(/\t/, $not_mismatch_position);

	foreach my $line (@result_tabs)
	{
		my @a = split(/\t/, $line);
		my $align = reverse $a[7];

		my $mismatch = 0;
		for (my $i = 0; $i <= 8; $i++) 
		{
			if (substr($align, $i, 1) eq "b") { $mismatch++; }
		}

		my $mismatch_meet_position = 0;
		foreach my $pos (@position)
		{
			if ( substr($align,$pos,1) eq "b" ) { $mismatch_meet_position++; }
		}

		if ($mismatch < $max_mismatch_cutoff && $mismatch_meet_position == 0)
		{
			$result_tab_filtered.= $line."\n";
		}
	}

	return $result_tab_filtered;
}

=head1 run_target_score


=cut
sub run_target_score
{

	$USAGE = qq'
Usage: perl $0 miRNA_seq_fasta blast_formatted_db source_db output

* miRNA_seq_fasta -- miRNA sequence file identified by miRNA_prediction.pl
* blast_formatted_db -- mRNA DB formatted by blast 
* source_db -- mRNA sequence file
* output -- output result (filtered)
';
	
	if ($#ARGV < 3) { # one argument input
		print $USAGE;
    		exit;
	}
	
	my $out_file = "./".$ARGV[3];
	open OUTFILE, ">$out_file" || die("Cannot Open File");;
	close OUTFILE;
		
	my $source_file = "";
	my $source_file_name = $ARGV[2];
	my $db_name = $ARGV[1];
	
	if ($seq_type eq 1) 
	{ 
		$source_file = $DATA_PATH."/".$source_file_name; 
		$blastdb_smRNA = $db_name;  #  Tomato smRNA
		
	}
	else 
	{ 
		$source_file = $DATA_PATH."/".$source_file_name;
		$blastdb_mRNA = $db_name;  #  Tomato unigene file
	}

	my %seqsa = insert_seq($source_file);
	
	my $direction = "+";
	my $search_direction = "both";
	
	if ($choice_direction eq "forward") { $search_direction = "forward"; }
	elsif ($choice_direction eq "reverse") { $search_direction = "reverse"; }	

	my %seq_querys = query_seq_from_file($ARGV[0]);
	
	my $seq_num = 1;
	
	my $resulthash = {};
	my $schash = {};
	my @index_list = ();

	foreach my $description (@description_list) {
		
		my $query_n = $description;
		$query_n =~ s/^\>//;
		print $seq_num."-th input query: ".$query_n;


		my $result_tab = "";
		
		my @seqb = ($description, $seq_querys{$description});
	
		if ( (length($seqb[1]) <= 30) and ($seq_type eq 1) ) { print "\nPlease check the query sequence type!";  exit; }
		if ( (length($seqb[1]) > 30) and ($seq_type eq 0) ) { print "\nPlease check the query sequence type!";  exit; }

		my $smRNA = "";
		my $smRNArc = "";
		my $mRNA = "";
		my $mRNAdisc = "";
		my $smRNAdisc = "";
		
		@sel_list_disc=();
		@sel_list_pos_mRNA_st=();
		@sel_list_pos_sRNA_st=();
		@sel_list_pos_sRNA_end=();
		@sel_list_direction=();
		
		my $count = blast($seqb[0], $seqb[1], $tempb, $seq_type, $search_direction);
		
		my $i = 0;

		foreach my $disc (@sel_list_disc) {

			$direction = $sel_list_direction[$i];

			$mRNAdisc = split_discription($mRNAdisc);
			if ( $seq_type eq 1 ) {
				($mRNAdisc, $smRNAdisc, $mRNA, $smRNA) = change_sequences($seqb[0], $disc, $seqb[1], $seqsa{$disc}, $seq_type);
			}
			else {
				($mRNAdisc, $smRNAdisc, $mRNA, $smRNA) = change_sequences($disc, $seqb[0], $seqsa{$disc}, $seqb[1], $seq_type);
			}

			if ( $direction eq "+" ) { $smRNArc = reverse_complement($smRNA); }
			else { $smRNArc = $smRNA; }

			chomp($smRNArc);
			chomp($mRNA);

			my @sinfo = lsa($mRNAdisc, $smRNAdisc, $mRNA, $smRNArc, $i, $direction);

			if ( ($sinfo[0] eq 1) ) {
				my $start_pos=sprintf("%d", $sinfo[3]);
				my $index = join "_", $mRNAdisc, $smRNAdisc, $start_pos;

				$resulthash{$index} = \@sinfo;
				$schash{$index} = $sinfo[1];
				push @index_list, $index;
			}
			$i++;
		}

		$result_tab = print_tab(%schash, %resulthash);

		#################################################
		# filter the result tab info			#
		#################################################
		if ($result_tab) 
		{ 
			$result_tab = filter_output_result($result_tab, 2, "9\t10"); 
			open OUTFILE, ">>$out_file" || die("Cannot Open File");
			print OUTFILE $result_tab;
			close OUTFILE;
		}
		
		$seq_num = $seq_num + 1;
		print "...Done\n";
		
		foreach my $index (@index_list) {
			delete $resulthash{$index}; 
			delete $schash{$index};
		}
		@index_list = ();
	}
}

run_target_score;

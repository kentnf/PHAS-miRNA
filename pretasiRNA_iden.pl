#!/usr/bin/perl

=head1 

 phasing_iden.pl -- identify pahsing genes using small RNA datasets

 Yi Zheng

 11/27/2013
 1. connect the window, than compute the phasing score
 2. generate the phased_tasiRNA

 10/14/2013
 1. correct the phasing score
 2. add cutoff for phasing score
 3. remove repeat calculationof phasing score for saving time
 4. add more parameters for identification of phasing regions
 5. fix the bug for best phasing score
 6. remove some ouput files

 08/12/2013
 init

=cut

use strict;
use warnings;
use PerlIO::gzip;
use IO::File;
use Bio::SeqIO;
use Getopt::Long;

my $usage = qq'
perk $0 [options]
	-i input_sRNAmapping_file (BAM or SAM)  
	-a max multihits	[default 6]
	-x max sRNA length	[default 24]
	-m min sRNA length	[default 21]
	-c cycle size		[default 21]
	-w window cycle		[default 9]
	-s step cycle           [default 3; this is step for window connection]
	-n min uniq reads	[default 10]
	-k min register reads	[default 3]
	-r reads ratio		[default 0.5]
	-p cutoff p-value 	[default 0.001]
	-e cutoff phasing score [default NA]

Output file name: x_connectWindow.txt, 
		  x_phasingscore.bedgraph, 
		  x_phased.fa
		  x is the prefix of input file name

';

my ($help, $input_data, $max_hits, $max_len, $min_len, 
    $cycle_size, $window_cycle, $shift_cycle, 
    $min_uniq, $min_register, $reads_ratio, 
    $cutoff_pvalue, $cutoff_pscore);

GetOptions(
	"h"	=> \$help,
	"i=s"	=> \$input_data,
	"a=i"	=> \$max_hits,
	"x=i"	=> \$max_len,
	"m=i"	=> \$min_len,
	"c=i"	=> \$cycle_size,
	"w=i"	=> \$window_cycle,
	"s=i"	=> \$shift_cycle, 	# steps
	"n=i"	=> \$min_uniq,
	"k=i"	=> \$min_register,
	"r=s"	=> \$reads_ratio,
	"p=s"	=> \$cutoff_pvalue,
	"e=s"	=> \$cutoff_pscore
);

die $usage if $help;
die $usage unless $input_data;

#################################################################
# setting of window and cycle					#
#################################################################
$max_hits ||= 6;
$min_len ||= 21;
$max_len ||= 24;

$cycle_size ||= 21;
$window_cycle ||= 9;
$shift_cycle ||= 3;

$min_uniq ||= 10;
$min_register ||= 3;
$reads_ratio ||= 0.5;

unless ( defined $cutoff_pvalue ) { $cutoff_pvalue ||= 0.001; }

my $window_size = $cycle_size * $window_cycle;
my $shift_size = $cycle_size * $shift_cycle;

#################################################################
# check file name of input dataset 				#
#################################################################
my ($sam_file, $output_prefix);
if ($input_data =~ m/\.bam$/) 
{
	my $bam_file = $input_data;
	$output_prefix = $input_data; 	$output_prefix =~ s/\.bam$//;
	$sam_file = $bam_file; 		$sam_file =~ s/\.bam$/\.sam/;
	bam2sam($bam_file, $sam_file);
} 
elsif ($input_data =~ m/\.sam$/)
{
	$sam_file = $input_data;
	$output_prefix = $input_data;	$output_prefix =~ s/\.sam$//;
}
else
{
	die "Error, file name of input data: $input_data\n";
}


#################################################################
# main:                                                         #
#################################################################

#################################################################
# load small RNA mapping data form sam file			#
#								#
# SAM info were loaded into three hash				#
# $sRNA_map_in_cycle; $sRNA_map_out_cycle, $sRNA_map_all	#
# key1: ref/mRNA ID , position , strand				#
# key2: sRNA ID							#
# value: sequence						#
#								#
# sequence id and length info: $reference_id			#
# key: sequence id						#
# value: length							#
#################################################################
my ($sRNA_map_in_cycle, $sRNA_map_out_cycle, $sRNA_map_all, $reference_id) = load_mapping($sam_file, $cycle_size, $min_len, $max_len, $max_hits);
#print scalar(keys(%$reference_id))."\n"; die;

#################################################################
# identified and connect the window according to below rules	#
# 1. uniq sRNA in window   > $min_uniq  			#
# 2. register in window	   > $min_register			#
# 3. 21nt sRNA / uniq sRNA > $reads_ratio			#
# 4. connected the identified window with overlaps		#
# 5. find if there is any register on window edge to correct	#
#    the window start and end					#
#################################################################
my $connect_window = identify_connect_window(); 
if ( scalar(keys(%$connect_window)) == 0 ) { exit; }

# foreach my $wd (sort keys %$connect_window) { print $wd."\n"; }

#################################################################
# filter the raw windows with P value / Phasing scores		#
# 1. scan connected window with 189bp slide window and 1bp step #
# 2. find the best one according to pvalue			#
# 3. find if there is any register on window edge to correct    #
#    the window start and end					#
# 4. calculate the phasing score for each window		#
# 5. samve the result to file: connect_window, phasingscore.bed #
#################################################################

# create hash for phasing score for all select position, using the globe hash to shave time
# key: refid \t position
# value: pscore
# if the phascore was defined, do not calculate it more
# set output files
my $all_pscore;

my $filtered_window = filter_correct_window();

my $output_conWin = $output_prefix."_connectWindow.txt";
my $outwin = IO::File->new(">".$output_conWin) || die "Can not open connect window $output_conWin $!\n";
print $outwin "#RefId\tstart\tend\ttotal\t".$cycle_size."nt\tregister\tBshift\tp-value\tpScore\n";
print $outwin $filtered_window;
$outwin->close;

my $output_bed = $output_prefix."_phasingscore.bedgraph";
my $outbed = IO::File->new(">".$output_bed) || die "Can not open bedgraph file $output_bed $!\n";
chomp($filtered_window);
my @r = split(/\n/, $filtered_window);
foreach my $region ( @r )
{
        my ($ref_id, $start, $end) = split(/\t/,$region);

        for(my $pos=$start; $pos<=$end; $pos++)
        {
                if ( defined $$all_pscore{$ref_id."\t".$pos} )
                {
                        my $phasing_score = $$all_pscore{$ref_id."\t".$pos};
                        print $outbed $ref_id."\t".$pos."\t".$pos."\t".$phasing_score."\n";
                }
                else
                {
                        print "Error, don not have this phascore before for position $ref_id:$pos\n";
                }
        }
}
$outbed->close;

#################################################################
# generate register reads for each window			#
#################################################################
my $phased_read = generate_phased_sRNA(@r);

my $output_phased_read = $output_prefix."_phased.fa";
my $outpha = IO::File->new(">".$output_phased_read) || die "Can not open phased read file $output_phased_read $!\n";
print $outpha $phased_read;
$outpha->close;

#################################################################
# kentnf: subroutine						#
#################################################################
=head1 bam2sam

 convert bam file to sam file

=cut
sub bam2sam
{
        my ($bam_file, $sam_file) = @_;
        die "Error in suffix of bam file $bam_file\n" unless $bam_file =~ m/\.bam/;
        die "Error in suffix of sam file $sam_file\n" unless $sam_file =~ m/\.sam/;

        my $cmd = "samtools view -h -o $sam_file $bam_file";
        print $cmd."\n";
        system($cmd) && die "Error in command $cmd\n";
}

=head1 load_mapping

 load mapping info from small RNA alignment SAM file

=cut
sub load_mapping
{
	my ($sam_file, $cycle_size, $min_len, $max_len, $max_hits) = @_;

	my (%sRNA_map_in_cycle, %sRNA_map_out_cycle, %sRNA_map_all, %reference_id, %seq_length, %sRNA_count);

	# get the number of hits for each aligned sRNA
	my $fh0 = IO::File->new($sam_file) || die "Can not open input SAM file $sam_file $!\n";
	while(<$fh0>) {
		if($_ =~ m/^@/) { next; }
		my @a = split(/\t/, $_);
		if (defined $sRNA_count{$a[0]}) { $sRNA_count{$a[0]}++; }
		else { $sRNA_count{$a[0]} = 1; }
	}
	$fh0->close;

	# generate hash for sRNA alignemnt
	my $fh = IO::File->new($sam_file) || die "Can not open input SAM file $sam_file $!\n";
	while(<$fh>)
	{
		chomp;
		if($_ =~ m/^@/) 
		{ 
			if ($_ =~ m/^\@SQ/) {
				my @m = split(/\t/, $_);
				my $id = $m[1]; $id =~ s/SN://;
				my $len = $m[2]; $len =~ s/LN://;
				$seq_length{$id} = $len;
				next;
			} else {
				next;
			}
		}

		my @a = split(/\t/, $_);
		my ($sid, $str, $rid, $start, $seq) = ($a[0], $a[1], $a[2], $a[3], $a[9]);

		if      ( $str == 0  ) { $str = 1;   }
		elsif   ( $str == 16 ) { $str = -1;  }
		else	{ next; } 				# filter out unmapped reads
		if ( $sRNA_count{$sid} > $max_hits ) { next; }  # filter out reads with max hits
		if ( length($seq) > $max_len || length($seq) < $min_len ) { next; } # filter out reads with other length
		unless ( defined $seq_length{$rid} ) 		# check the reference length
		{ die "Error, do not have seq length for reference sequnce: $rid\n"; }
		
		if ( length($seq) eq $cycle_size )
		{
			$sRNA_map_in_cycle{$rid.",".$start.",".$str}{$sid} = $seq;
		}
		else
		{
			$sRNA_map_out_cycle{$rid.",".$start.",".$str}{$sid} = $seq;
		}
		$sRNA_map_all{$rid.",".$start.",".$str}{$sid} = $seq;
		$reference_id{$rid} = $seq_length{$rid};
	}
	$fh->close;

	return (\%sRNA_map_in_cycle, \%sRNA_map_out_cycle, \%sRNA_map_all, \%reference_id);
}

=head1 identify_connect_window

=cut
sub identify_connect_window
{
	# identify windows
	my $window_report = "";
	foreach my $ref_id (sort keys %$reference_id)
	{
		my $ref_len = $$reference_id{$ref_id};

		my ($window_start, $window_end, $real_start, $real_end,
		$num_out_cycle_read, $num_in_cycle_read, $num_uniq_read,
		$num_register, $base_shift, $pvalue, $pscore, $report,
 		$window_aligned_reads);

		my %window_mapped;
		for(my $i=1; $i<=$ref_len; $i=$i+$shift_size)
		{
			$window_start = $i;
               		$window_end = $i+$window_size-1;
                	if ($window_end > $ref_len) { $window_end = $ref_len; }

			$num_out_cycle_read = 0;        # read length is equal to cycle length
 			$num_in_cycle_read = 0;         # read length is not equal to cycle length
 			$num_uniq_read = 0;             # all uniq read in window = $num_out_cycle_read + $num_in_cycle_read
			$num_register = 0;              # number of cycle reads fell into register
                
			# checking the report stat of window
                	($num_in_cycle_read, $num_register, $base_shift, $num_uniq_read, $real_start, $real_end) = count_reads($ref_id, $window_start, $window_end, $cycle_size, 0);

			if ($num_uniq_read >= $min_uniq && $num_in_cycle_read/$num_uniq_read > $reads_ratio && $num_register >= $min_register )
                	{
				#my $pvalue = cal_pvalue($num_in_cycle_read, $num_register, $cycle_size, $window_cycle);
				my $pvalue = 1;
                        	$window_report.="$ref_id\t$window_start\t$window_end\t$pvalue\n";
                        	#print "$ref_id\t$window_start\t$window_end\t$num_uniq_read\t$num_in_cycle_read\t$num_register\t$pvalue\n";
                	}
		}
	}

	#################################################################
	# connect report windows if they have overlap                   #
	# input: $window_report                                         #
	# 	it just include reference id, start, end, pvalue	#
	# output: %connect_window                                       #
	# key: referenceID \t start \t end; value = min pvalue		#
	#################################################################
	my %connect_window = ();
	if ($window_report) { %connect_window = connect_window($window_report); }

	return (\%connect_window);
}

=head1 filter window
 
 filter the winow with p value/phasing score

=cut
sub filter_correct_window
{
	my $filter_window = "";
	foreach my $con_window (sort keys %$connect_window)
	{
		my ($ref_id, $con_window_start, $con_window_end) = split(/\t/, $con_window);
		my $ref_length = $$reference_id{$ref_id};

		# scan the connected window with 189bp window, and using one base cycle step
		# find the registers have same base shift with window start, then compute pvalue
		# find the best one from the window according to pvalue
		my $best_pvalue = 1;
		my ($best_report, $best_start, $best_end);

		for(my $i=$con_window_start; $i<=$con_window_end; $i=$i+1)
		{
			my $window_start = $i;
                        my $window_end = $i+$window_size-1;
			my $last_status = 0;
                        if ($window_end >= $con_window_end) { $last_status = 1; }
			if ($window_end >= $ref_length) {
				$window_start = $ref_length-$window_size+1;
				$window_end = $ref_length;
				$last_status = 1;
			}

			# get the stat information for each window
			my $mode = 1; # only identify registers according to the window_start
			my ($num_in_cycle_read, $num_register, $base_shift, $num_uniq_read, $real_start, $real_end) = count_reads($ref_id, $window_start, $window_end, $cycle_size, $mode);

			# cumpute pvalue for each window with suitable stat and generate the best report according to the pvalue
			if ( $num_register >= $min_register  )
			{
				my $pvalue = cal_pvalue($num_in_cycle_read, $num_register, $cycle_size, $window_cycle);
				my $report = "$ref_id, $window_start, $window_end, $real_start, $real_end, $num_uniq_read, $num_in_cycle_read, $num_register, $base_shift, $pvalue";

				if ( $pvalue < $best_pvalue ) {
					$best_pvalue = $pvalue; $best_report = $report; $best_start = $window_start; $best_end = $window_end-$cycle_size+1;
				}
			}
			if ($last_status == 1) { last; }
		}

		# correct the window of best real start and end,
		# scan the correct window with 189bp window
		# get the best one as report results according to pvalue

		if ($best_pvalue < $cutoff_pvalue) {
			my ($correct_start, $correct_end) = correct_window($ref_id, $best_start, $best_end, $cycle_size);
			$correct_end = $correct_end + $cycle_size - 1;

			my ($report_start, $report_end);
			if ( $correct_start  < $best_start  ) { $report_start = $correct_start; }
			else { $report_start = $best_start; }
			
			if ( $correct_end  > $best_end  ) { $report_end = $correct_end; }
			else { $report_end = $best_end; }

			my @a = split(/, /, $best_report);

			# get best phasing score
			my $best_pscore;
			($best_pscore, $all_pscore) = cal_pscore($ref_id, $report_start, $report_end, $all_pscore);
			$filter_window.="$a[0]\t$report_start\t$report_end\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$a[9]\t$best_pscore\n";
		}
	}

	return $filter_window;
}

=head1  generate_phased_sRNA

=cut
sub generate_phased_sRNA
{
	my @window = @_;

	my $phased_sRNA = "";

	foreach my $w (@window)
	{
		my @a = split(/\t/, $w);
		my ($ref_id, $start, $end) = ($a[0], $a[1], $a[2]);
		for(my $i=$start; $i<=$end; $i=$i+$cycle_size)
		{
			my $pos = $i;
			my $anti_pos = $i-2;
			my $key_sense     = $ref_id.",".$pos.",1";
			my $key_antisense = $ref_id.",".$anti_pos.",-1";

			if (defined $$sRNA_map_in_cycle{$key_sense} )
			{
				foreach my $ssid ( sort keys %{$$sRNA_map_in_cycle{$key_sense}} )
				{
					$phased_sRNA.= ">".$ssid."\n";
					$phased_sRNA.= $$sRNA_map_in_cycle{$key_sense}{$ssid}."\n";
				}
			}
		}
	}

	return $phased_sRNA;
}

=head1 correct_window

=cut
sub correct_window
{
	my ($ref_id, $real_start, $real_end, $cycle_size) = @_;

	my ($correct_start, $correct_end) = ($real_start, $real_end);
	while(1)
	{
		if ($correct_start <= $cycle_size) { last; }

		my $upstream = 0;

		$correct_start = $correct_start - $cycle_size;
		my $up1step       = $correct_start;
		my $anti_up1step  = $correct_start - 2;
		my $key_sense1    = $ref_id.",".$up1step.",1";
		my $key_antisense1= $ref_id.",".$anti_up1step.",-1";
		if (defined $$sRNA_map_in_cycle{$key_sense1} || $$sRNA_map_in_cycle{$key_antisense1}) { $upstream = 1; }
		
		if ($upstream == 0)
		{
			$correct_start = $correct_start - $cycle_size;
			my $up2step	  = $correct_start;
			my $anti_up2step  = $correct_start - 2;
			my $key_sense2    = $ref_id.",".$up2step.",1";
			my $key_antisense2= $ref_id.",".$anti_up2step.",-1";
			if (defined $$sRNA_map_in_cycle{$key_sense2} || $$sRNA_map_in_cycle{$key_antisense2}) { $upstream = 1; }
		}

		if ($upstream == 0) { 
			$correct_start = $correct_start + $cycle_size + $cycle_size;
			last; 
		}
	}

	while(1)
	{
		my $downstream = 0;
		
		$correct_end = $correct_end + $cycle_size;
		my $down1step		= $correct_end;
		my $anti_down1step	= $correct_end - 2;
		my $key_sense1		= $ref_id.",".$down1step.",1";
		my $key_antisense1	= $ref_id.",".$anti_down1step.",-1";

		if (defined $$sRNA_map_in_cycle{$key_sense1} || $$sRNA_map_in_cycle{$key_antisense1}) { $downstream = 1; }
		
		if ($downstream == 0)
		{
			$correct_end = $correct_end + $cycle_size;
			my $down2step           = $correct_end;
			my $anti_down2step      = $correct_end - 2;
			my $key_sense2          = $ref_id.",".$down2step.",1";
			my $key_antisense2      = $ref_id.",".$anti_down2step.",-1";

			if (defined $$sRNA_map_in_cycle{$key_sense2} || $$sRNA_map_in_cycle{$key_antisense2}) { $downstream = 1; }
		}

		if ($downstream == 0) { 
			$correct_end = $correct_end - $cycle_size - $cycle_size;
			last;
		}
	}

	#print "#$real_start, $real_end, $correct_start, $correct_end\n";
	return ($correct_start, $correct_end);
}

=head1 count_reads

 get n k m for window connection, P value, and P score;

=cut
sub count_reads
{
	my ($ref_id, $window_start, $window_end, $cycle_size, $mode) = @_;

	my ($num_out_cycle_read, $num_in_cycle_read, $num_uniq_read, $num_register, $base_shift, $window_aligned_reads, $real_start, $real_end);

	$num_out_cycle_read = 0;        # read length is equal to cycle length
	$num_in_cycle_read = 0;         # read length is not equal to cycle length
	$num_uniq_read = 0;             # all uniq read in window = $num_out_cycle_read + $num_in_cycle_read

	$num_register = 0;              # number of cycle reads fell into register
	my @potential_register_position;# position info of cycle reads for calculating register

	# checking the report stat of window
	for(my $j=$window_start; $j<=$window_end; $j++)
	{
		my $aj = $j - 2;
		my $key_sense     = $ref_id.",".$j.",1";
		my $key_antisense = $ref_id.",".$aj.",-1";
		my $sam_antisense = $ref_id.",".$j.",-1";
		my $exp_num = 0;

		if (defined $$sRNA_map_in_cycle{$key_sense} ) {
			foreach my $ssid ( sort keys %{$$sRNA_map_in_cycle{$key_sense}} )
			{
				$window_aligned_reads.=$key_sense."\t".$ssid."\t".$$sRNA_map_in_cycle{$key_sense}{$ssid}."\n";
				my @sid_exp = split(/-/, $ssid);
				$exp_num = $exp_num + $sid_exp[scalar(@sid_exp)-1];
			}
			#print $window_aligned_reads; die;

			$num_in_cycle_read++;
			push(@potential_register_position, $j);
		}

		if (defined $$sRNA_map_in_cycle{$key_antisense} ) {
			foreach my $ssid ( sort keys %{$$sRNA_map_in_cycle{$key_antisense}} )
			{
				$window_aligned_reads.=$key_antisense."\t".$ssid."\t".$$sRNA_map_in_cycle{$key_antisense}{$ssid}."\n";
			}
			#print $window_aligned_reads; die;

			$num_in_cycle_read++;
			push(@potential_register_position, $j);
		}

		if (defined $$sRNA_map_in_cycle{$sam_antisense} ) {
			foreach my $ssid ( sort keys %{$$sRNA_map_in_cycle{$sam_antisense}} )
			{
				my @sid_exp = split(/-/, $ssid);
				$exp_num = $exp_num + $sid_exp[scalar(@sid_exp)-1];
			}
		}

		if (defined $$sRNA_map_out_cycle{$key_sense} ) {
			foreach my $ssid ( sort keys %{$$sRNA_map_out_cycle{$key_sense}} )
			{
				$window_aligned_reads.=$key_sense."\t".$ssid."\t".$$sRNA_map_out_cycle{$key_sense}{$ssid}."\n";
				my @sid_exp = split(/-/, $ssid);
				$exp_num = $exp_num + $sid_exp[scalar(@sid_exp)-1];
			}
                        #print $window_aligned_reads; die;

			my $num_uniq = scalar(keys($$sRNA_map_out_cycle{$key_sense}));
			$num_out_cycle_read = $num_out_cycle_read + $num_uniq;
		}

		if (defined $$sRNA_map_out_cycle{$key_antisense} ) {
			foreach my $ssid ( sort keys %{$$sRNA_map_out_cycle{$key_antisense}} )
 			{
				$window_aligned_reads.=$key_antisense."\t".$ssid."\t".$$sRNA_map_out_cycle{$key_antisense}{$ssid}."\n";
			}
                       	#print $window_aligned_reads; die;

			my $num_uniq = scalar(keys($$sRNA_map_out_cycle{$key_antisense}));
			$num_out_cycle_read = $num_out_cycle_read + $num_uniq;
		}

		if (defined $$sRNA_map_out_cycle{$sam_antisense} ) {
			foreach my $ssid ( sort keys %{$$sRNA_map_out_cycle{$sam_antisense}} )
			{
				my @sid_exp = split(/-/, $ssid);
				$exp_num = $exp_num + $sid_exp[scalar(@sid_exp)-1];
			}
		}

		#if ( $exp_num>0 ) { print $exp "$ref_id\t$j\t$j\t$exp_num\n"; }
	}

	# get register and base-shift information
	my $init_base_shift = "NA";
	if ($mode > 0) { $init_base_shift = $window_start % $cycle_size; }

	if (scalar @potential_register_position > 0 ) {
		($num_register, $base_shift, $real_start, $real_end) = count_register(\@potential_register_position, $cycle_size, $init_base_shift);
	} else {
		($num_register, $base_shift, $real_start, $real_end) = (0, "NA", $window_start, $window_end); 
	}
	$num_uniq_read = $num_out_cycle_read + $num_in_cycle_read;

	return ($num_in_cycle_read, $num_register, $base_shift, $num_uniq_read, $real_start, $real_end);
}

=head1 count_register

=cut
sub count_register
{
	my ($position, $cycle_size, $init_base_shift) = @_;

	# stat the position of cycle-size (21nt) sRNA to hash
	my %uu; # key: base_shift, value: number of register
	foreach my $pos (@$position)
	{
		my $u = $pos % $cycle_size;
		if (defined $uu{$u}) { $uu{$u}++; }
		else { $uu{$u} = 1; }
	}

	# compute the number of register, base shift;
	my $max_reg = 0; my $base_shift = "NA";	

	if ($init_base_shift eq "NA")
	{
		foreach my $k (sort keys %uu) {
			if ($uu{$k} > $max_reg) { $max_reg = $uu{$k}; $base_shift = $k; }
		}
	}
	else
	{
		if ( defined $uu{$init_base_shift} ) {
			$max_reg = $uu{$init_base_shift};
		} else {
			$max_reg = 0;
		}
		$base_shift = $init_base_shift;
	}

	# get the real start and real end according to base shift
	my ($real_start, $real_end);
	foreach my $pos (@$position)
	{
		if ( $pos % $cycle_size == $base_shift )
		{
			if (defined $real_start && defined $real_end) {
				if ($pos < $real_start) { $real_start = $pos; }
				if ($pos > $real_end) { $real_end = $pos; }
			} else {
				$real_start = $pos; $real_end = $pos;
			}
		}
	}

	$real_start = "NA" unless defined $real_start;
	$real_end   = "NA" unless defined $real_end;

	return ($max_reg, $base_shift, $real_start, $real_end);
}

=head1 

=cut
sub connect_window
{
	my $report_window = shift;
	# "ref_id\t window_start\t window_end\t num_uniq_read\t num_in_cycle_read\t num_register\n";

	my %region;

	chomp($report_window);
	my @report_window = split(/\n/, $report_window);

	my $line1 = shift @report_window;

	my @n = split(/\t/, $line1);
	my ($pre_ref_id, $pre_start, $pre_end, $pre_pvalue) = ($n[0], $n[1], $n[2], $n[3]);

	foreach my $line (@report_window)
	{
		my @a = split(/\t/, $line);
		my ($ref_id, $start, $end, $pvalue) = ($a[0], $a[1], $a[2], $a[3]);

		if ($ref_id eq $pre_ref_id )
		{
			if ( $start <= ($pre_end + 1))
			{
				$pre_end = $end;
				if ($pvalue < $pre_pvalue) { $pre_pvalue = $pvalue; }
			}
			else
			{
				$region{$pre_ref_id."\t".$pre_start."\t".$pre_end} = $pre_pvalue;
				$pre_start = $start;
				$pre_end = $end;
				$pre_ref_id = $ref_id;
				$pre_pvalue = $pvalue;
			}
		}
		else
		{
			$region{$pre_ref_id."\t".$pre_start."\t".$pre_end} = $pre_pvalue;
			$pre_start = $start;
			$pre_end = $end;
			$pre_ref_id = $ref_id;
			$pre_pvalue = $pvalue;
		}
	}

	$region{$pre_ref_id."\t".$pre_start."\t".$pre_end} = $pre_pvalue;

	return %region;
}

=head cal_pvalue

=cut
sub cal_pvalue
{
	my ($n, $k, $cycle_size, $window_cycle) = @_;

	my ($pvalue, $pr, $c, $rr, $rw, $m, $m20, $m21);
	$m = $window_cycle * 2;
	$m20 = ($window_cycle * 2) * ($cycle_size - 1);
	$m21 = ($window_cycle * 2) *  $cycle_size;
	$pvalue = 0;

	for ( my $w=$k; $w<=$m; $w++ )
	{
		$c=1; $rr=1; $rw=1;

		for (my $j=0; $j<=$w-1; $j++) { $c=$c*($n-$j)/($j+1); }

                for (my $x=0; $x<=$w-1; $x++) { $rr=$rr*($m-$x)/($m21-$x); }

                for (my $y=0; $y<=$n-$w-1; $y++) { $rw=$rw*($m20-$y)/($m21-$w-$y); }

                $pr = $c*$rr*$rw;
                $pvalue = $pvalue + $pr;
	}
	return $pvalue;
}

=head cal_pscore

 return the best phasing score for a window or region

=cut
sub cal_pscore
{
	my ($ref_id, $start_region, $end_region, $all_pscore) = @_;	
	
	my $best_pscore = 0;

	for(my $pos=$start_region; $pos<=$end_region; $pos++)
	{
		my $phasing_score;

		if ( defined $$all_pscore{$ref_id."\t".$pos} )
		{
			$phasing_score = $$all_pscore{$ref_id."\t".$pos};
		}
		else
		{
			my %nhash;              # hash for count register
			my $cycle = 0;          # init cycle order

			my $start = $pos - ($cycle_size * ($window_cycle - 5));
	                my $end = $pos + ( $cycle_size * ($window_cycle - 4) - 1 );

			my ($knum, $pnum, $unum) = (0,0,0);
		
	                for(my $i=$start; $i<=$end; $i=$i+$cycle_size)
        	        {
                	        $cycle++;
	                        my $key_sense = $ref_id.",".$i.","."1";
        	                if (defined $$sRNA_map_in_cycle{$key_sense}) {
                	                $nhash{$cycle} = 1;
                        	        foreach my $sid (sort keys $$sRNA_map_in_cycle{$key_sense})
                                	{
	                                        my @sid_exp = split(/-/, $sid);
        	                                my ($psid, $exp_num) = ($sid_exp[0], $sid_exp[1]);
                	                        die "Error in sid $sid\n" unless scalar(@sid_exp) == 2;
                        	                die "Error in exp num $sid\n" if $exp_num < 1;
                                	        $pnum = $pnum + $exp_num;
	                                }
        	                }

	                        for (my $o=$i+1; $o<$i+$cycle_size; $o++)
        	                {
                	                my $key = $ref_id.",".$o.","."1";
                        	        if (defined $$sRNA_map_in_cycle{$key})
	                                {
        	                                foreach my $sid (sort keys $$sRNA_map_in_cycle{$key})
                	                        {
                        	                        my @sid_exp = split(/-/, $sid);
                                	                my ($psid, $exp_num) = ($sid_exp[0], $sid_exp[1]);
                                        	        die "Error in sid $sid\n" unless scalar(@sid_exp) == 2;
                                                	die "Error in exp num $sid\n" if $exp_num < 1;
	                                                $unum = $unum + $exp_num;
        	                                }
                	                }
                        	}
                	}

	                $cycle = 0;
        	        for (my $j=$start-2; $j<=$end-2; $j=$j+$cycle_size)
                	{
                        	$cycle++;
	                        my $key_antisense = $ref_id.",".$j.","."-1";
        	                if (defined $$sRNA_map_in_cycle{$key_antisense}) {
                	                $nhash{$cycle} = 1;
                        	        foreach my $sid (sort keys $$sRNA_map_in_cycle{$key_antisense})
                                	{
                                        	my @sid_exp = split(/-/, $sid);
	                                        my ($psid, $exp_num) = ($sid_exp[0], $sid_exp[1]);
        	                                die "Error in sid $sid\n" unless scalar(@sid_exp) == 2;
                	                        die "Error in exp num $sid\n" if $exp_num < 1;
                        	                $pnum = $pnum + $exp_num;
                                	}
                        	}

	                        for(my $k=$j+1; $k<$j+$cycle_size; $k++)
        	                {
                	                my $key = $ref_id.",".$k.","."-1";
                        	        if (defined $$sRNA_map_in_cycle{$key})
                                	{
                                        	foreach my $sid (sort keys $$sRNA_map_in_cycle{$key})
                                        	{
                                                	my @sid_exp = split(/-/, $sid);
	                                                my ($psid, $exp_num) = ($sid_exp[0], $sid_exp[1]);
        	                                        die "Error in sid $sid\n" unless scalar(@sid_exp) == 2;
                	                                die "Error in exp num $sid\n" if $exp_num < 1;
                        	                        $unum = $unum + $exp_num;
                                	        }
                                	}
                        	}
                	}

			foreach my $nk (sort keys %nhash) { if (defined $nhash{$nk}) { $knum++; } }

			# compute $phasing score
			my $nn1 = ($pnum / ( $unum + 1)) * 10 + 1;
	                my $nn2;
        	        if ($knum > 3) {
                	        $nn2 = $nn1;
	                        for(my $ik = 0; $ik < $knum-3; $ik++) { $nn2 = $nn2 * $nn1; }
        	        }
	                elsif ($knum == 3 ) {
        	                $nn2 = $nn1 * 1;
	                }
       		        else {
                	        $nn2 = 1;
                	}

	                $phasing_score = log($nn2);
	                $phasing_score = sprintf("%.2f", $phasing_score);	
			$$all_pscore{$ref_id."\t".$pos} = $phasing_score;
		}

		if ($phasing_score > $best_pscore) { $best_pscore = $phasing_score; }
	}
	
	return ($best_pscore, $all_pscore);
}

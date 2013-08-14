#!/usr/bin/perl

=head1 

 phasing_iden.pl -- identify pahsing genes using small RNA datasets

 Yi Zheng

 08/12/2012

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
	-r reference_fasta
	-c cycle size		[default 21]
	-w window cycle		[default 9]
	-s shift cycle 		[default 3; this is step for window]
	-p cutoff p-value 	[default 0.001]

Output file name: x_pvalue.txt, 
		  x_phasingscore.bedgraph, 
		  x is the prefix of input file name

';

my ($help, $input_data, $ref_fa, $min_len, $max_len, $cycle_size, $window_cycle, $shift_cycle, $cutoff_pvalue);

GetOptions(
	"h"	=> \$help,
	"i=s"	=> \$input_data,
	"r=s"	=> \$ref_fa,
	#"min_len=i"	=> \$min_len,
	#"max_len=i"	=> \$max_len,
	"c=i"	=> \$cycle_size,
	"w=i"	=> \$window_cycle,
	"s=i"	=> \$shift_cycle, 	# steps
	"p=s"	=> \$cutoff_pvalue
);

die $usage if $help;
die $usage unless $input_data;
die $usage unless $ref_fa;

#################################################################
# setting of window and cycle					#
#################################################################
$min_len = 21;
$max_len = 24;

$cycle_size ||= 21;
$window_cycle ||= 9;
$shift_cycle ||= 3;

$cutoff_pvalue ||= 0.001;

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
# load small RNA mapping data form sam file			#
#								#
# SAM info were loaded into three hash				#
# $sRNA_map_in_cycle; $sRNA_map_out_cycle, $sRNA_map_all	#
# key1: ref/mRNA ID , position , strand				#
# key2: sRNA ID							#
# value: sequence						#
#################################################################
my ($sRNA_map_in_cycle, $sRNA_map_out_cycle, $sRNA_map_all, $reference_id) = load_mapping($sam_file, $cycle_size, $min_len, $max_len);
#print scalar(keys(%$reference_id))."\n"; die;

#################################################################
# main: 							#
#################################################################
my $output_report = $output_prefix.".report.gz";
my $output_pvalue = $output_prefix."_pvalue.txt";
my $output_4image = $output_prefix."_4image.txt";

open(OPT, ">:gzip", $output_report) || die "Can not open output report file $output_report $!\n";
my $pva = IO::File->new(">".$output_pvalue) || die "Can not open output pvalue file $output_pvalue $!\n";
my $img = IO::File->new(">".$output_4image) || die "Can not open output pvalue file $output_4image $!\n";

print OPT  "#RefId\tstart\tend\tuniq mapped\tin cycle\tregister\tBshift\tp-value\n";
print $pva "#RefId\tstart\tend\tuniq mapped\tin cycle\tregister\tBshift\tp-value\n";
print $img "#cycle size: $cycle_size\n#window cycle: $window_cycle\n#shift cycle: $shift_cycle\n";
print $img "#RefId\tstart\tend\tuniq mapped\tin cycle\tregister\tBshift\tp-value\n";

my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$ref_fa);
while(my $inseq = $in->next_seq)
{
	my $ref_id = $inseq->id;
	unless ( defined $$reference_id{$ref_id}) { next; }
	my $ref_len = $inseq->length;
	#print "$ref_id $ref_len\n";

	my ($window_start, $window_end, 
	    $num_out_cycle_read, $num_in_cycle_read, $num_uniq_read, 
	    $num_register, $base_shift, $pvalue, $report, 
	    $window_aligned_reads);

	my %window_mapped;
	for(my $i=1; $i<=$ref_len; $i=$i+$shift_size)
	{
		$window_start = $i;
		$window_end = $i+$window_size-1;
		if ($window_end > $ref_len) { $window_end = $ref_len; }
		#print $window_start, "\t", $window_end,"\n";

		$num_out_cycle_read = 0;	# read length is equal to cycle length
		$num_in_cycle_read = 0;		# read length is not equal to cycle length
		$num_uniq_read = 0;		# all uniq read in window = $num_out_cycle_read + $num_in_cycle_read

		$num_register = 0;		# number of cycle reads fell into register
		my @potential_register_position;# position info of cycle reads for calculating register
		
		$report = 0;
		$window_aligned_reads = "";
		#%window_map = (); %window_map_cycle = ();
		
		# checking the report stat of window
		for(my $j=$window_start; $j<=$window_end; $j++)
		{
			my $aj = $j - 2;
			my $key_sense     = $ref_id.",".$j.",1";
			my $key_antisense = $ref_id.",".$aj.",-1";

			if (defined $$sRNA_map_in_cycle{$key_sense} ) {
				foreach my $ssid ( sort keys %{$$sRNA_map_in_cycle{$key_sense}} )
				{
					$window_aligned_reads.=$key_sense."\t".$ssid."\t".$$sRNA_map_in_cycle{$key_sense}{$ssid}."\n";
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
				
			if (defined $$sRNA_map_out_cycle{$key_sense} ) {
				foreach my $ssid ( sort keys %{$$sRNA_map_out_cycle{$key_sense}} )
				{
					$window_aligned_reads.=$key_sense."\t".$ssid."\t".$$sRNA_map_out_cycle{$key_sense}{$ssid}."\n";
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
		}

		$num_uniq_read = $num_out_cycle_read + $num_in_cycle_read;
		($num_register, $base_shift) = count_register(\@potential_register_position, $cycle_size);
		
		$pvalue = "NA"; $report = 0;
		if ($num_uniq_read >= 10 && $num_in_cycle_read/$num_uniq_read > 0.5 && $num_register >= 3)
		{
			$pvalue = cal_pvalue($num_in_cycle_read, $num_register, $cycle_size, $window_cycle);
			#my ($min_pvalue) = phasing_pvalue(\%window_map_cycle, \%sRNA_map, $cycle_size, $window_start, $window_end);
			if ($pvalue < $cutoff_pvalue) { $report = 1; }
		}

		print OPT "$ref_id\t$window_start\t$window_end\t$num_uniq_read\t$num_in_cycle_read\t$num_register\t$base_shift\t$pvalue\n";

		if ($report == 1)
		{
			my $window_report = "$ref_id\t$window_start\t$window_end\t$num_uniq_read\t$num_in_cycle_read\t$num_register\t$base_shift\t$pvalue";
			print $pva $window_report."\n";
			print $img $window_report."\n".$window_aligned_reads;
		}
	}
}

close(OPT);
$pva->close;
$img->close;

# load info from file $outout_pvalue
my $report_pvalue;
my $rpva = IO::File->new($output_pvalue) || die "Can not open output pvalue file $output_pvalue $!\n";
while(<$rpva>)
{
	chomp;
	if ($_ =~ m/^#/) { next; }
	$report_pvalue.=$_."\n";
}
$rpva->close;
unless($report_pvalue) { exit; }

#################################################################
# connect report pvalue region					#
# input: $report_pvalue						#
# output: %report_region					#
# key: referenceID \t start \t end				#
#################################################################
my %report_region = connect_window($report_pvalue);

#################################################################
# call phasing score for report region				#
#################################################################
my $output_bed = $output_prefix."_phasingscore.bedgraph";
my $outbed = IO::File->new(">".$output_bed) || die "Can not open bedgraph file $output_bed $!\n";

foreach my $region (sort keys %report_region)
{
	my ($ref_id, $start_region, $end_region) = split(/\t/,$region);
	#print $region."\n";

	for(my $pos=$start_region; $pos<=$end_region; $pos++)
	{
		my %nhash;		# hash for count register
		my %ahash;
		my $k = 0;		# number of register in cycle-size
		my $cycle = 0;		# init cycle order
		my $start = $pos;	# start of window & end of window
		my $end = $pos + ( $cycle_size * $window_cycle - 1 ); 

		my $pnum;
		my %phasing_num; # key: cycle, value: total num of mapped sRNA reads in each phased cycle

		for(my $i=$start; $i<=$end; $i=$i+$cycle_size)
		{
			$cycle++;
			my $key_sense = $ref_id.",".$i.","."1";
			if (defined $$sRNA_map_in_cycle{$key_sense}) { $nhash{$cycle} = 1; }
		
			$pnum = 0;
			for (my $o=$i; $o<$i+$cycle_size; $o++)
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
						$pnum = $pnum + $exp_num;
					}
				}
			}

			if ( defined $phasing_num{$cycle} ) { $phasing_num{$cycle} = $phasing_num{$cycle} + $pnum; }
			else { $phasing_num{$cycle} = $pnum; }
		}

		$cycle = 0;
		for (my $j=$start-2; $j<=$end-2; $j=$j+$cycle_size)
		{
			$cycle++;
			my $key_antisense = $ref_id.",".$j.","."-1";
			if (defined $$sRNA_map_in_cycle{$key_antisense}) { $ahash{$cycle} = 1; }

			$pnum = 0;
			for(my $k=$j; $k<$j+$cycle_size; $k++)
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
						$pnum = $pnum + $exp_num;
					}
				}
			}

			if ( defined $phasing_num{$cycle} )
			{
				$phasing_num{$cycle} = $phasing_num{$cycle} + $pnum;
			}
			else
			{
				$phasing_num{$cycle} = $pnum;
			}
		}

		foreach my $nk (sort keys %nhash) { if (defined $nhash{$nk}) { $k++; } }
		foreach my $ak (sort keys %ahash) { if (defined $ahash{$ak}) { $k++; } }
		my $phasing_score = phasing_score(\%phasing_num, $k);

		print $outbed $ref_id."\t".$pos."\t".$pos."\t".$phasing_score."\n";

		if ($k > 3)
		{
			#print  $ref_id,"\t",$pos,"\t",$pos,"\t",$phasing_score,"\t$k\n";

			#foreach my $n (sort {$a<=>$b} keys %phasing_num) { 
			#	print $n."\t".$phasing_num{$n}."\n"; 
			#}
		}
	}
}
$outbed->close;

#################################################################
# convert the bedgraph to bigwig file				#
#################################################################


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
	my ($sam_file, $cycle_size, $min_len, $max_len) = @_;

	my (%sRNA_map_in_cycle, %sRNA_map_out_cycle, %sRNA_map_all, %reference_id);

	my $fh = IO::File->new($sam_file) || die "Can not open input SAM file $sam_file $!\n";
	while(<$fh>)
	{
		chomp;
		if($_ =~ m/^@/) { next; }

		my @a = split(/\t/, $_);
		my ($sid, $str, $rid, $start, $seq) = ($a[0], $a[1], $a[2], $a[3], $a[9]);

		if      ($str == 0 ) { $str = 1; }
		elsif   ($str == 16) { $str = -1;}
		else	{ print "Error at strand info in line:\n$_\n"; next; }

		if ( length($seq) > $max_len || length($seq) < $min_len ) { next; }
		
		if ( length($seq) eq $cycle_size )
		{
			$sRNA_map_in_cycle{$rid.",".$start.",".$str}{$sid} = $seq;
		}
		else
		{
			$sRNA_map_out_cycle{$rid.",".$start.",".$str}{$sid} = $seq;
		}
		$sRNA_map_all{$rid.",".$start.",".$str}{$sid} = $seq;
		$reference_id{$rid} = 1;
	}
	$fh->close;

	return (\%sRNA_map_in_cycle, \%sRNA_map_out_cycle, \%sRNA_map_all, \%reference_id);
}

=head1 count_register

=cut
sub count_register
{
	my ($position, $cycle_size) = @_;

	my %uu;
	foreach my $pos (@$position)
	{
		my $u = $pos % $cycle_size;
		if (defined $uu{$u}) { $uu{$u}++; }
		else { $uu{$u} = 1; }
	}

	my $max_reg = 0; my $base_shift = "NA";
	foreach my $k (sort keys %uu)
	{
		if ($uu{$k} > $max_reg) { $max_reg = $uu{$k}; $base_shift = $k; }
	}

	return ($max_reg, $base_shift);
}

=head1 

=cut
sub connect_window
{
	my $report_pvalue = shift;
	# "$ref_id\t$window_start\t$window_end\t$num_uniq_read\t$num_in_cycle_read\t$num_register\t$base_shift\t$pvalue\n";

	my %region;

	my @report_pvalue = split(/\n/, $report_pvalue);

	my $line1 = shift @report_pvalue;
	my @n = split(/\t/, $line1);
	my ($pre_ref_id, $pre_start, $pre_end) = ($n[0], $n[1], $n[2]);

	foreach my $line (@report_pvalue)
	{
		my @a = split(/\t/, $line);
		my ($ref_id, $start, $end) = ($a[0], $a[1], $a[2]);

		if ($ref_id eq $pre_ref_id)
		{
			if ( $start <= ($pre_end + 1))
			{
				$pre_end = $end;
			}
			else
			{
				$region{$pre_ref_id."\t".$pre_start."\t".$pre_end} = 1;
				$pre_start = $start;
				$pre_end = $end;
				$pre_ref_id = $ref_id;
			}
		}
		else
		{
			$region{$pre_ref_id."\t".$pre_start."\t".$pre_end} = 1;
			$pre_start = $start;
			$pre_end = $end;
			$pre_ref_id = $ref_id;
		}
	}
	$region{$pre_ref_id."\t".$pre_start."\t".$pre_end} = 1;

	return %region;
}

=head1 phasing_pvalue

=cut
=head
sub phasing_pvalue
{
	my ($window_map_cycle, $sRNA_map_cycle, $cycle_size, $window_size, $window_start, $window_end);


	my ($n, $k, $w_start, $w_end, $nkey, $key, $p_value, $min_p_value);

	foreach my $skey (sort keys %$window_map_cycle)
	{
		$n=-1;
		$k=-1;

		($ref_id, $position, $strand)=split(/,/, $skey);

		my %sRNA_mapping = ();	# key: sRNA ID, value:
		my %sRNA_label_N = ();	# key: sRNA ID, value: N 
		my %sRNA_label_K = ();  # key: sRNA ID, value: K

		# get n and k number for smallRNA
		if ( $strand == 1 ) # small RNAs on the Watson strand
		{
			$w_start = $position;				# small start
			$w_end   = $position + $window_size - 1;	# small end

			for (my $i=$w_start; $i<=$w_end; $i++) # calculate n on the sense strand
			{
				$nkey = $ref_id.",".$i.","."1";
				if ( defined $$sRNA_map_cycle{$nkey} ) 
				{ 
					$n=$n+1; 
					my @sid = get_sRNA($$sRNA_map_cycle{$nkey});
					foreach my $id (@sid) 
					{ 
						$sRNA_label_N{$id} = 1;
						$key = $id."#".$i."#1";
						$sRNA_mapping{$key} = 1;
					}
				}
			}

			for (my $i=$w_start; $i<=$w_end; $i=$i+$cycle_size) # calculate k on the sense strand
			{
				$nkey = $ref_id.",".$i.","."1";
				if ( defined $$sRNA_map_cycle{$nkey} ) 
				{ 
					$k=$k+1;
					my @sid = get_sRNA($$sRNA_map_cycle{$nkey});
					foreach my $id (@sid) 
					{ 
						$sRNA_label_K{$id} = 1;
						$key = $id."#".$i."#1";
						$sRNA_mapping{$key} = 1;
					}
				}
			}

			for (my $j=$w_start-2; $j<=$w_end-2; $j++) # calculate n on the antisense strand
			{
				$nkey = $ref_id.",".$j.","."-1";
				if ( defined $$sRNA_map_cycle{$nkey} ) 
				{ 
					$n=$n+1; 
					my @sid = get_sRNA($$sRNA_map_cycle{$nkey});
					foreach my $id (@sid) 
					{ 
						$sRNA_label_N{$id} = 1;
						$key = $id."#".$j."#-1";
						$sRNA_mapping{$key} = 1;
					}
				}
			}

			for ($j=$w_start+18; $j<=$w_end-2; $j=$j+$cycle_size) # calculate k on the antisense strand
			{
				$nkey = $ref_id.",".$j.","."-1";
				if ( defined $$sRNA_map_cycle{$nkey} ) 
				{
					$k=$k+1;
					my @sid = get_sRNA($$sRNA_map_cycle{$nkey});
					foreach my $id (@sid) 
					{ 
						$sRNA_label_K{$id} = 1;
						$key = $id."#".$j."#-1";
						$sRNA_mapping{$key} = 1;
					}
				}
			}
		}

		# ? do not have small RNA on crick strand ? using mRNA as ref
		elsif ($str==-1) # small RNAs on the Crick strand
		{
			$w_start = $position;
			$w_end   = $position - $window_size + 1;

			for ($i=$w_start; $i<=$w_end; $i++) # calculate n on the sense strand
			{
				$nkey = $ref_id.",".$i.","."-1";
				if ( defined $$sRNA_map_cycle{$nkey} ) 
				{ 
					$n=$n+1;
					my @sid = get_sRNA($$sRNA_map_cycle{$nkey});
					foreach my $id (@sid) 
					{ 
						$sRNA_label_N{$id} = 1;
						$key = $id."#".$i."#-1";
						$sRNA_mapping{$key} = 1;
					}	
				}
			}

			for ($i=$w_start+20; $i<=$w_end; $i=$i+21) # calculate k on the sense strand
			{
				$nkey = $ref_id.",".$i.","."-1";
				if ( defined $$sRNA_map_cycle{$nkey}) 
				{ 
					$k=$k+1;
					my @sid = get_sRNA($$sRNA_map_cycle{$nkey});
					foreach my $id (@sid) 
					{ 
						$sRNA_label_K{$id} = 1;
						$key = $id."#".$i."#-1";
						$sRNA_mapping{$key} = 1;
					}
				}
			}

			for ($j=$w_start+2; $j<=$w_end+2; $j++) # calculate n on the antisense strand
			{
				$nkey = $ref_id.",".$j.","."1";
				if ( defined $$sRNA_map_cycle{$nkey} ) 
				{ 
					$n=$n+1;
					my @sid = get_sRNA($$sRNA_map_cycle{$nkey});
					foreach my $id (@sid) 
					{ 
						$sRNA_label_N{$id} = 1;
						$key = $id."#".$j."#1";
						$sRNA_mapping{$key} = 1;
					}
				}
			}

			for ($j=$w_start+2; $j<=$w_end+2; $j=$j+21) # calculate k on the antisense strand
			{
				$nkey = $ref_id.",".$j.","."1";
				if ( defined $$sRNA_map_cycle{$nkey} ) 
				{ 
					$k=$k+1;
					my @sid = get_sRNA($$sRNA_map_cycle{$nkey});
					foreach my $id (@sid) 
					{ 
						$sRNA_label_K{$id} = 1;
						$key = $id."#".$j."#1";
						$sRNA_mapping{$key} = 1;
					}
				}
			}	
		}

		# calculate p-value from n and k
		$p_value = cal_pvalue($n, $k);
		
		# output all small RNA clusters with their p values 
		# mRNA   Start   Strand  sRNA    sRNASeq  N  K  P
		my $line = $$window_map_cycle{$skey};

		$output_p_value.="$line\t$n\t$k\t$p_value\n";


		# output mRNA amd small RNA for view.
		# mRNA   Start   Strand  sRNA    sRNASeq
		@ref = split(/\t/, $line);
		$ref_sRNA = $ref[3];

		$out_view = "";

		foreach my $sRNA (sort keys %sRNA_mapping)
		{
	    		($sid, $mapping_cor, $mapping_str) = split(/#/, $sRNA);
	    		$len = $sq{$sid}; 
	    		$end = $mapping_cor+$len-1;
		
	    		if ($sid eq $ref_sRNA) 
	    		{ 
				$label = "R";
				$out_view = "$chr\t$ref_sRNA\t$cor\t$str\t$sid\t$mapping_cor\t$end\t$mapping_str\t$label\n".$out_view;
				#$out_view = "$sid\t$chr\t$cor\t$mapping_cor\t$end\t$mapping_str\t$label\n".$out_view;
	    		}
	    		else 
	    		{
				if ($sRNA_label_N{$sid} && $sRNA_label_K{$sid}) { $label = "K"; }
				elsif ( $sRNA_label_K{$sid}) { $label = "K"; die "Error, a label with K must be labeled with N\n$ref_sRNA\t$sid\n"; }
				elsif ( $sRNA_label_N{$sid}) { $label = "N"; }
				else  { die "Error, $sRNA do not have n or k label for $ref_sRNA\n"; }
				#$out_view.="$sid\t$chr\t$cor\t$mapping_cor\t$end\t$mapping_str\t$label\n";
				$out_view.="$chr\t$ref_sRNA\t$cor\t$str\t$sid\t$mapping_cor\t$end\t$mapping_str\t$label\n";
	    		}
		}

		#select and output small RNA clusters with p<0.001
		if ( $p_value < 0.001 ) { 
			$output_sig_p.="$line\t$n\t$k\t$p\n"; 
			$output_view.= $out_view;
		}
	}
}
=cut
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

=head1 phasing_score

=cut
sub phasing_score
{
	my ($phase_num, $k) = @_;

	my $num_cycle = scalar(keys(%$phase_num)); # number of cycles in window

	my $sum_phase_num = 0; 	# total number of 21nt reads in all cycles phase
	my @p;			# total number of 21nt reads in each cycle phase 
	foreach my $cycle (sort {$a<=>$b} keys %$phase_num )
	{
		my $num = $$phase_num{$cycle};
		push(@p, $num);
		$sum_phase_num = $sum_phase_num + $num;
	}

	# method 1
	my $psum = 0;
	my $usum = 0;

	my @u;
	for(my $i=0; $i<$num_cycle; $i++)
        {
                $psum = $psum + $p[$i];
                my $u = $sum_phase_num - $p[$i];
                push(@u, $u);
                $usum = $usum + $u;

		if ($k > 3) { print "$p[$i]\t$psum\t$u\t$usum\t$sum_phase_num\n"; }
        }

        $usum++;

        my $sum = (10 * $psum/$usum) + 1;
        my $phasing_score = 0;	# set the default value of phasing score
	my $result = 1;		# var for multiply with k

	if ($k >= 3)
        {
                for(my $i=0; $i<$k-2; $i++) { $result = $result * $sum; }
                $phasing_score = log($result);
                $phasing_score = sprintf("%.2f", $phasing_score);
        }


        my $ns = scalar(@p);
        my $nk = scalar(keys(%$phase_num));


	if ($k > 5)
	{
		print "A:$ns B:$nk\tN:$k\tPsum:$psum\tUsum:$usum\tSUM:$sum\tResult:$result\tPscore:$phasing_score\n";
		print "Pi: ",join(" ", @p),"\n";
		print "Ui: ",join(" ", @u),"\n";
	}

	return $phasing_score;
}

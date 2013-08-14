#!/usr/bin/perl

=head1 

 ta-siRNA_drawing1.pl -- drawing ta-siRNA images for reference sequence and sRNA using GD

 Yi Zheng

 08/12/2012

=cut

use strict;
use warnings;
use GD;
use Bio::SeqIO;
use IO::File;

my $usage = qq'
Perl $0  4image.txt  reference_sequence.fa 

* 4image is the output of phasing_iden.pl

';

my $mapping_4image = shift || die $usage;
my $ref_seq = shift || die $usage;

my $flanking_seq_len = 100;

#################################################################
# load parameters                                               #
# * cycle size, window cycle, and shift cycle are automatically #
# * loaded from mapping_4image file                             #
#################################################################
# mapping_4image file to hash                                   #
# key: report_pvalue info                                       #
# value: sRNA align info                                        #
################################################################
my ($report_mapping, $cycle_size, $window_cycle, $step_cycle) = mapping_to_hash($mapping_4image);

#################################################################
# ref fasta seq to hash						#
# key: ref_id, start, end,					#
# value: seq, start_shift, $end_shift				#
#################################################################
my %ref_seq_hash;
my $in = Bio::SeqIO->new(-format=>'fasta', -file=>$ref_seq);
while(my $inseq = $in->next_seq)
{
	my $seq_id = $inseq->id;
	my $seq_len = $inseq->length;

	foreach my $report (sort keys %$report_mapping)
	{
		my ($ref_id, $start, $end, $uniq_mapped, $in_cycle, $register, $base_shift, $pvalue) = split(/\t/, $report);
		if ($ref_id eq $seq_id)
		{
			my $ref_seq_len = $end - $start + 1;
			my ($start_shift, $end_shift) = ($flanking_seq_len, $flanking_seq_len);
			if ($start < $flanking_seq_len) { $start_shift = $start - 1; }
			if ($end > $seq_len) { die "Error in report $report\n"; }
			if ($end + $flanking_seq_len > $seq_len ) { $end_shift = $seq_len - $end; }
			my $real_start = $start - $start_shift;
			my $real_end = $end + $end_shift;
			my $real_length = $real_end - $real_start + 1;
			my $sub_seq = substr($inseq->seq, $real_start-1, $real_length);

			my $key = $ref_id."\t".$start."\t".$end;
			my $val = $sub_seq."\t".$start_shift."\t".$end_shift;
			$ref_seq_hash{$key} = $val;

			#print $key."\n".$val."\n"; 
		}
	}
}

#################################################################
# report info for check						#
#################################################################
print "cycle_size: $cycle_size\nwindow_cycle: $window_cycle\nstep_cycle: $step_cycle\n";
print "Significant windows: ".scalar(keys(%$report_mapping))."\n";
print "Reference sequences: ".scalar(keys(%ref_seq_hash))."\n";

#################################################################
# Set vars and Constants					#
#################################################################
use constant IMAGE_SIZE_X => 750;
use constant BG_COLOR  => (255, 255, 255);
use constant XMARGIN  => 45;
use constant YMARGIN  => 50;
use constant MULPIXEL => 6;  #length of each character (Sequence Text)
#use constant MULPIXEL => 1;  # Without Sequence Text (with SEQLINEHIGHT > 0)
use constant SEQHIGHT => 40;
use constant SEQLINEHIGHT => 0; # Draw line (>0)
use constant SEQLINEMARGIN => 4;
use constant SRNAHIGHT => 1;
use constant SRNAYMARGIN => 10;
use constant MAXSRNALEN => 50;

#################################################################
# Parse the report_mapping hash					#
#################################################################
foreach my $report (sort keys %$report_mapping)
{
	my $info = $$report_mapping{$report};
	chomp($report);
	chomp($info);

	# clean the hash
	my %hithash = ();

	# new values from report info
	my ($ref_id, $start, $end, $uniq_mapped, $in_cycle, $register, $base_shift, $pvalue) = split(/\t/, $report);
	# get reference sequence in window size with flanking sequence 
	die "Error, can not get reference seq for:\n$report\n" unless defined $ref_seq_hash{$ref_id."\t".$start."\t".$end};
	my ($ref_seq, $start_shift, $end_shift)  = split(/\t/, $ref_seq_hash{$ref_id."\t".$start."\t".$end});

	# put the sRNA mapping info to hash
	my @lines = split(/\n/, $info);
	foreach my $line (@lines)
	{
		my @n = split(/\t/, $line);
		my @m = split(/,/,  $n[0]);
		my ($rid, $sRNA_start, $sRNA_strand, $sRNA_id, $sRNA_seq) = ($m[0], $m[1], $m[2], $n[1], $n[2]);
		my $sRNA_len = length($sRNA_seq);
		
		# get label info for c, r, or o
		# c : reads in cycle
		# r : reads in register
		# o : reads out of cycle
		my $cro = get_cro($sRNA_start, $sRNA_strand, $start, $base_shift, $cycle_size, $sRNA_len);
		die "Error, can not assign cro to info: $line\n$report\n" if $cro eq "NA";

		my $map_cor = $sRNA_start - $start + 1 + $start_shift; 

		my $key = join "\t", $rid, $map_cor;
		my $val = join "\t", $sRNA_id, $sRNA_len, $sRNA_strand, $cro;

		if (exists $hithash{$key}) { $hithash{$key}.=$val."\n"; }
		else { $hithash{$key} = $val."\n"; };
		#print "$key\n$val\n";
	}
	build_image($ref_id, $start, $end, $start_shift, $end_shift, $ref_seq, \%hithash);
}

#################################################################
# kentnf: subroutine						#
#################################################################

=head1 mapping_to_hash

 convert mapping info to hash, and get the cycle_size, window_size, shift_cycle form mapping_4image file

=cut
sub mapping_to_hash
{
        my $mapping_4image = shift;

        my %report_mapping;
        my ($cycle_size, $window_cycle, $shift_cycle, $key, $value);

        my $fh = IO::File->new($mapping_4image) || die "Error, Can not open mapping 4image file $mapping_4image $!\n";
        while(<$fh>)
        {
                chomp;
                if ($_ =~ m/^#/)
                {
                        my @a = split(/:\s+/, $_, 2);
                        if ($a[0] =~ "cycle size")      { $cycle_size = $a[1]; }
                        if ($a[0] =~ "window cycle")    { $window_cycle = $a[1]; }
                        if ($a[0] =~ "shift cycle")     { $shift_cycle = $a[1]; }
                }
                else
                {
                        my @a = split(/\t/, $_);

                        if (scalar(@a) == 8)
                        {
                                if ($key && $value) { $report_mapping{$key} = $value; }
                                $key = "";
                                $value = "";
                                $key = $_;
                        }
                        elsif (scalar(@a) == 3)
                        {
                                $value.=$_."\n";
                        }
                        else
                        {
                                die "Error in line: $_\n";
                        }
                }

                if (eof)
                {
                        if ($key && $value) { $report_mapping{$key} = $value; }
                }
        }
        $fh->close;

        die "Error in cycle size\n" unless $cycle_size;
        die "Error in window cycle\n" unless $window_cycle;
        die "Error in shift cycle\n" unless $shift_cycle;

        return(\%report_mapping, $cycle_size, $window_cycle, $shift_cycle);
}

=head1 get_cro


=cut
sub get_cro
{
	my ($sRNA_start, $sRNA_strand, $start, $base_shift, $cycle_size, $sRNA_len) = @_;

	my $cpos; # cpos : compute position, it is pos - 2 when strand is minus
	if ($sRNA_strand eq '-1') { $cpos = $sRNA_start + 2; }
	else { $cpos = $sRNA_start; }

	my $cro = "NA";

	if ($sRNA_len == $cycle_size)
	{
		if ( (($cpos - ($start + $base_shift - 1)) % $cycle_size ) == 0 ) { $cro = "r"; }
		else { $cro = "c"; }
	}
	else
	{
		$cro = "o";
	}

	return $cro;
}

=head1 build_image

 build image for mapped sRNA with reference sequences 

=cut
sub build_image 
{ 
	#my ($ref_id, $ref_seq,$sid) = @_;
	my ($ref_id, $start, $end, $start_shift, $end_shift, $ref_seq, $hithash) = @_;

	# defined the reference sequence length in window with flanking seq
	my $ref_seq_len = length($ref_seq);

	# define output file names
	my $output_png_file  = $ref_id."-".$start."-".$end.".png";

	# split sequence to array
       	my @seq = split // => $ref_seq;

	my $gap = 1;
	
	my $max_level = 500;

        my $E = [];
        for(my $i = 0; $i <= $max_level; $i++) { $E->[$i] = -100; }

	my $LS = []; my $LE = [];
        for(my $i = 0; $i <= $ref_seq_len; $i++) {  $LS->[$i] = 0; $LE->[$i] = 0; }

	my $max_depth = 1;

	# convert hit hash to lev hash
	my %levhash = ();
	for(my $i = 0; $i < $ref_seq_len; $i++) 
	{
		my $key = $ref_id."\t".$i;

		my $depth = 0;

		if (defined $$hithash{$key} ) 
		{
			my $sRNA_list = $$hithash{$key}; 
			my @sRNAs = split(/\n/, $sRNA_list);

			foreach my $S (@sRNAs) 
			{
				my ($sRNA_id, $sRNA_len, $sRNA_strand, $cro) = split(/\t/, $S);

        			for(my $j = 0; $j <= $max_level; $j++)
				{
					my $e_pos = $E->[$j]+$gap;

        				if ($LE->[$i+$sRNA_len] < $j) { # $i: 1 mRNA base
        					$LE->[$i+$sRNA_len] = $j;
					}

					if ($i > $e_pos ) {
						$E->[$j] = $i+$sRNA_len; 
        					my $key_lev = $ref_id."\t".$i."\t".$j;
					    	$levhash{$key_lev} = $S;

						if($max_depth < $j) {$max_depth = $j;}

						$depth = $j;
					}

					last if ($i > $e_pos ); 
				}

			}
		}
        	$LS->[$i] = $depth;
	}

	# get max depth from pre code
	#print "$max_depth, $ref_seq_len, $ref_seq, $ref_id\n";

	display_sRNA_seq($max_depth, $ref_seq_len, $ref_seq, $ref_id, \%levhash, $LS, $LE, $output_png_file, $start_shift, $end_shift); 

	######  For debug......
	#print "\nLS:\n";
        #for(my $i = 0; $i <= $seq_len; $i++) { 
	#	if ($LS->[$i] > 0) {
	#		print "$i:$LS->[$i] ";
	#	}
	#}
	#print "\nLE:\n";
        #for(my $i = 0; $i <= $seq_len; $i++) { 
	#	if ($LE->[$i] > 0) {
	#		print "$i:$LE->[$i] ";
	#	}
	#}

	# Get depth information
        #for(my $i = 0; $i <= $seq_len; $i++) { 
        #        my $max_d = $LS->[$i];
	#	print "$i: $LS->[$i]\n";

       	# 	for(my $j = 0; $j <= $max_d; $j++) {
	#		$key = $description."_".$i."_".$j;
	#		if ($levhash{$key}=~ m/\S+/) {
	#			print "key: $key\n";
	#		}
	#	}
	#}


       	#for(my $j = 0; $j <= $max_depth; $j++) {
        #	print "Level $j";
        #	for(my $i = 0; $i <= $seq_len; $i++) { 
	#		$key = $description."_".$i."_".$j;
	#		if ($levhash{$key}=~ m/\S+/) {
	#			print " ($key, $levhash{$key}) ";
	#		}
	#	}
        #	print "\n";
	#}
}

=head1 dispaly_sRNA_seq

 dispaly sRNA sequence 

=cut
sub display_sRNA_seq 
{ 
	my ($max_depth, $ref_seq_len, $ref_seq, $ref_id, $levhash, $LS, $LE, $output_png_file, $start_shift, $end_shift) = @_;
	#print "$max_depth, $ref_seq_len, $ref_seq, $ref_id, $levhash, $LS, $LE, $ouput_png_file\n";

	my @imgmaps = ();
	my $img_pos;  
	my $seq_start;
        my $seq_end;
        my $subseq_len;
        my $subseq;
        my $y_pos_end;
        my $last_seq_end;

	my $str_seq_end;
	my $str_len;

	my $y_pos = YMARGIN;  				# 50 pix of ymargin
	my $max_x = $ref_seq_len * MULPIXEL;		# length of sequence pix length, 6pix/base
	my $line_x = IMAGE_SIZE_X - (XMARGIN * 2);	# 700 - 20 x 2 = 660, 
	my $seq_line_num = int($max_x / $line_x);	# get number of lines 
	my $max_line_num = $seq_line_num;		# max line num = seq line number
	my $last_line_x = int($max_x % $line_x);	# the num of pix in last line

	$max_depth = $max_depth + 1;

	if ($last_line_x > 0) { $max_line_num = $seq_line_num + 1; }

	my $LL = [];	# L: length, x axis
	my $YY = [];	# Y: height, y axis

        for(my $i = 0; $i < $max_line_num; $i++) {	# a lot of 0 to array, init pixel of Y axis
		$YY->[$i] = 0; 
	}
	for(my $i = 0; $i < $ref_seq_len; $i++) {		# a lot of 0 to array, init pixel of X axis
		$LL->[$i] = 0; 
	}

	######## Get the depth and y postions in each line
        my $max_d;
        my $line_max_d;

	$subseq_len = int($line_x/MULPIXEL);		# 660/6 = 110
	$YY->[0] = YMARGIN;

        for(my $i = 0; $i < $seq_line_num; $i++) {	# seq_line_num = int(1052/110)+1 = 10
		$seq_start = $subseq_len*$i;
		$seq_end = $seq_start+$subseq_len;
		$line_max_d = 0;
        	for(my $j = $seq_start; $j <= $seq_end; $j++) { 
        		$max_d = $LS->[$j];
			if ($max_d > $line_max_d) { $line_max_d = $max_d; }
        		$max_d = $LE->[$j];
			if ($max_d > $line_max_d) { $line_max_d = $max_d; }
		}
        	for(my $j = $seq_start; $j <= $seq_end; $j++) { 
			$LL->[$j] = $i;
		}
		$YY->[$i+1] = $YY->[$i]+SEQLINEHIGHT+SEQLINEMARGIN;  
		$YY->[$i+1] = $YY->[$i+1]+((SRNAHIGHT+SRNAYMARGIN)*$line_max_d);  
		$YY->[$i+1] = $YY->[$i+1]+SEQHIGHT;  
		$y_pos_end = $YY->[$i+1];	
	}

	if ($last_line_x > 0) {
		$seq_start = $subseq_len*$seq_line_num;
		$seq_end = $ref_seq_len;
		$line_max_d = 0;
        	for(my $j = $seq_start; $j <= $seq_end; $j++) { 
        		$max_d = $LS->[$j];
			if ($max_d > $line_max_d) { $line_max_d = $max_d; }
        		$max_d = $LE->[$j];
			if ($max_d > $line_max_d) { $line_max_d = $max_d; }
		}
        	for(my $j = $seq_start; $j <= $seq_end; $j++) { 
			$LL->[$j] = $seq_line_num;
		}
		$YY->[$seq_line_num+1] = $YY->[$seq_line_num]+SEQLINEHIGHT+SEQLINEMARGIN;  
		$YY->[$seq_line_num+1] = $YY->[$seq_line_num+1]+((SRNAHIGHT+SRNAYMARGIN)*$line_max_d);  
		$YY->[$seq_line_num+1] = $YY->[$seq_line_num+1]+SEQHIGHT;  
		$y_pos_end = $YY->[$seq_line_num+1];	
	}

	#########################################################
	# draw images using GD					#
	#########################################################

	my $image_size_y = $y_pos_end;

	#print "\n\n sRNA viewer\n";
	#print "\nLine Num: $seq_line_num\n";
	#print "\nLast Line X: $last_line_x\n";
	#print "\nimage_size_y: $image_size_y\n";

	my $im = new GD::Image(IMAGE_SIZE_X, $image_size_y);

	my ($white, $black, $red, $blue, $green, $gray, $orange);
	$white = $im->colorAllocate(255,255,255);
	$black = $im->colorAllocate(0,0,0);
	$red = $im->colorAllocate(255,0,0);
	$blue = $im->colorAllocate(0,0,255);
	$green = $im->colorAllocate(0,255,0);
	$gray = $im->colorAllocate(128,128,128);
	$orange = $im->colorAllocate(255,165,0);

	$im->transparent($white);
	$im->interlaced('true');

	######## Draw referece sequence 
	$subseq_len = int($line_x/MULPIXEL);

	my $start_i = int($start_shift/$subseq_len);
	my $start_x = $start_shift % $subseq_len;

	my $end_i = int( ($ref_seq_len - $end_shift) / $subseq_len);
	my $end_x = ($ref_seq_len - 1 - $end_shift) % $subseq_len;
	#print "$subseq_len $start_shift $start_i $start_x $end_shift $end_i $end_x\n"; #die;

        for(my $i = 0; $i < $seq_line_num; $i++) { 

		if (SEQLINEHIGHT > 0) { 
			$im->rectangle(XMARGIN, $YY->[$i], XMARGIN+$line_x, $YY->[$i]+SEQLINEHIGHT, $gray);
			$im->fill(XMARGIN+1, $YY->[$i]+1, $gray); 
		}
		
		$seq_start = int($line_x/MULPIXEL)*$i;
		$subseq = substr $ref_seq, $seq_start, $subseq_len;

		$str_seq_end = sprintf("%d", ($seq_start+$subseq_len) );
		$str_len = length($str_seq_end);

		$im->string(gdSmallFont, XMARGIN, $YY->[$i]-25, $seq_start+1, $gray);
		$im->string(gdSmallFont, XMARGIN+$line_x-(MULPIXEL*$str_len), $YY->[$i]-25, $str_seq_end, $gray);
		if (MULPIXEL eq 6) {
			$im->string(gdSmallFont, XMARGIN, $YY->[$i]-12, $subseq, $gray);
		}

		# draw start and end point
		if ($i eq $start_i) {
			$im->string(gdSmallFont, XMARGIN+MULPIXEL*$start_x, $YY->[$i]-25, "S", $black);
		}

		if ($i eq $end_i) {
			$im->string(gdSmallFont, XMARGIN+MULPIXEL*$end_x, $YY->[$i]-25, "E", $black);
		}
	}

	if ($last_line_x > 0) {
		if (SEQLINEHIGHT > 0) { 
			$im->rectangle(XMARGIN, $YY->[$seq_line_num], XMARGIN+$last_line_x, $YY->[$seq_line_num]+SEQLINEHIGHT, $gray);
			$im->fill(XMARGIN+1, $YY->[$seq_line_num]+1, $gray); 
		}

		$seq_start = int($line_x/MULPIXEL)*$seq_line_num;
		$subseq_len = length($ref_seq)-$seq_start;
		$subseq = substr $ref_seq, $seq_start, $subseq_len;

		$last_seq_end = $subseq_len*MULPIXEL;
		$str_seq_end = sprintf("%d", ($seq_start+$subseq_len) );
		$str_len = length($str_seq_end);

		$im->string(gdSmallFont, XMARGIN, $YY->[$seq_line_num]-25, $seq_start+1, $gray);
		$im->string(gdSmallFont, XMARGIN+$last_seq_end-(MULPIXEL*$str_len), $YY->[$seq_line_num]-25, $str_seq_end, $gray);
		if (MULPIXEL eq 6) {
			$im->string(gdSmallFont, XMARGIN, $YY->[$seq_line_num]-12, $subseq, $gray);
		}
	}

	if ($end_i > $seq_line_num) {
		$im->string(gdSmallFont, XMARGIN+MULPIXEL*$end_x, $YY->[$seq_line_num]-25, "E", $black);
	}

	######## Draw mapping sRNAs ########
	
	my $cro_exp;	# lables, c:in cycle, r:register, o: out cycle, exp: number of sRNA reads
	my $ys_pos;

	foreach my $key (sort keys %$levhash )
	{
		my ($rid, $seq_pos, $level) = split(/\t/, $key); $seq_pos = $seq_pos - 1;
		my ($sRNA_id, $sRNA_len, $sRNA_strand, $cro) = split(/\t/, $$levhash{$key});

		# set x and y of small RNA
		my $cx_pos = $seq_pos*MULPIXEL;					# $seq_pos = start of mRNA
		my $xs_pos = int( $cx_pos % $line_x);				# get start of X, line_x = 660 pixel , qiu yu shu %
		my $xe_pos = $xs_pos+($sRNA_len*MULPIXEL);			# get end of X,
		$ys_pos = $YY->[$LL->[$seq_pos]]+((SRNAHIGHT+SRNAYMARGIN)*$level)+SRNAYMARGIN;  # get Y of small RNA

		# compute the expression of sRNA_num count
		my $sRNA_num; if ($sRNA_id =~ m/-(\d+)$/) { $sRNA_num = $1; }
		# set color of small RNA
		my $sR_color;
		my $color_bright;
		if ($sRNA_num > 255) {  $color_bright = 255; } 
		else { $color_bright = $sRNA_num; }

		if ($cro eq "o") 
		{
			if ($sRNA_strand eq "+" || $sRNA_strand eq "1") { $sR_color = $im->colorAllocate(230,200,200); }
			else { $sR_color = $im->colorAllocate(200,235,200); }
		}
		elsif ($cro eq "r") 
		{
			if ($sRNA_strand eq "+" || $sRNA_strand eq "1") { $sR_color = $im->colorAllocate(200,0,0); } 
			else { $sR_color = $im->colorAllocate(0,0,200); }
		}
		else
		{
			if ($sRNA_strand eq "+" || $sRNA_strand eq "1") { $sR_color = $im->colorAllocate(250,165,15); }
			else { $sR_color = $im->colorAllocate(15,255,50); }
		}

	
		$cro_exp = "$cro ($sRNA_num)";

		# sRNA mapped to the mRNA range
		# im->rectangle will draw small RNA
		# img_pos will stroe small RNA position to image as map link
		# can add N or K before the small RNA
		if ( $ref_seq_len > ($seq_pos+$sRNA_len) ) 
		{
			if ($LL->[$seq_pos] eq $LL->[$seq_pos+$sRNA_len] ) 
			{
				$im->rectangle(XMARGIN+$xs_pos, $ys_pos, XMARGIN+$xe_pos, $ys_pos+SRNAHIGHT, $sR_color);
				$im->string(gdSmallFont, XMARGIN+$xs_pos-30, $ys_pos-5, $cro_exp, $black);
				$img_pos = join(",", $sRNA_id, XMARGIN+$xs_pos, $ys_pos, XMARGIN+$xe_pos, $ys_pos+SRNAHIGHT);
				push @imgmaps, $img_pos;
			}
			else 
			{
				$im->rectangle(XMARGIN+$xs_pos, $ys_pos, XMARGIN+$line_x, $ys_pos+SRNAHIGHT, $sR_color);
				$im->string(gdSmallFont, XMARGIN+$xs_pos-30, $ys_pos-5, $cro_exp, $black);
				$img_pos = join(",", $sRNA_id, XMARGIN+$xs_pos, $ys_pos, XMARGIN+$line_x, $ys_pos+SRNAHIGHT);
				push @imgmaps, $img_pos;

				$xe_pos = $xe_pos - $line_x;
				$ys_pos = $YY->[$LL->[$seq_pos+$sRNA_len]]+((SRNAHIGHT+SRNAYMARGIN)*$level)+SRNAYMARGIN;  

				$im->rectangle(XMARGIN, $ys_pos, XMARGIN+$xe_pos, $ys_pos+SRNAHIGHT, $sR_color);

				$img_pos = join(",", $sRNA_id, XMARGIN, $ys_pos, XMARGIN+$xe_pos, $ys_pos+SRNAHIGHT);
				push @imgmaps, $img_pos;
			}
		}
		else 
		{
				$xe_pos = $xs_pos+(($ref_seq_len-$seq_pos)*MULPIXEL);
				$im->rectangle(XMARGIN+$xs_pos, $ys_pos, XMARGIN+$xe_pos, $ys_pos+SRNAHIGHT, $sR_color);
				$im->string(gdSmallFont, XMARGIN+$xs_pos-30, $ys_pos-5, $cro_exp, $black);
				$img_pos = join(",", $sRNA_id, XMARGIN+$xs_pos, $ys_pos, XMARGIN+$xe_pos, $ys_pos+SRNAHIGHT);
				push @imgmaps, $img_pos;
		}	
	}

	# output image with png format 
	my $pic = IO::File->new(">$output_png_file") || die "Can not open png file $output_png_file $!\n";
	binmode $pic;
	print $pic $im->png;
	$pic->close;

	# return imgmap for link in map
	#return (\@imgmaps);
}

=head
sub html_out {
	print qq'
        <br><h3 align=center>siRNA Viewer of <a href=gene.cgi?ID=$gene_ID>$gene_ID</a></h3>
        <table border=1 cellpadding=0 cellspacing=0 align=center bordercolor=#999999 style=border-collapse:collapse>
        <tr><td>
        <img src=/TFGD/tmp/$png_file usemap="#sRNA" border="0">
	<map name="sRNA">
	';

	foreach my $L (@imgmaps) {
		my ($sRNA_ID, $x, $y, $xx, $yy) = split(",",$L);
		print qq'<area shape="rect" coords="$x,$y,$xx,$yy" href="$sRNA_URL$sRNA_ID">
		';
	}

	print qq'
	</map>	
        </td></tr>
	</table>
	';
}
=cut



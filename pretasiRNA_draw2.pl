#!/usr/bin/perl

=head1

 draw exprssion images for identified ta-siRNA base on sRNA mapping info

 Yi Zheng

 08/12/2013

=cut

use strict;
use warnings;
use IO::File;

my $usage = qq'
perl $0 4image.txt output_prefix

';

my $mapping_4image = shift || die $usage;
my $output_prefix = shift || die $usage;

#################################################################
# load parameters						#
# * cycle size, window cycle, and shift cycle are automatically	#
# * loaded from mapping_4image file				#
#################################################################
# mapping_4image file to hash					#
# key: report_pvalue info					#
# value: sRNA align info					#
#################################################################
my ($report_mapping, $cycle_size, $window_cycle, $step_cycle) = mapping_to_hash($mapping_4image);
#print "$cycle_size, $window_cycle, $step_cycle\n";
#print scalar(keys(%$report_mapping))."\n";

#################################################################
# main								#
#################################################################
foreach my $report (sort keys %$report_mapping)
{
	my $mapping_info = $$report_mapping{$report};

	# get the plot info form mapping info
	# outout pdf: ref_id-start-end.pdf 
	# plot info: has number of lines equal to window size
	# position \t plus_register \t plus_in_cycle \t minus_register \t minus_incycle \n
	my ($output_pdf, $plot_info) = parse_line( $report, $mapping_info, $cycle_size, $output_prefix);

	# plot the image for only in cycle and register sRNA
	plot_it($output_pdf, $plot_info);

	# convert pdf file to png file
	my $output_png = $output_pdf;
	$output_png =~ s/\.pdf$/\.png/;
	system("convert $output_pdf $output_png");
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
			if ($a[0] =~ "cycle size")	{ $cycle_size = $a[1]; }
			if ($a[0] =~ "window cycle")	{ $window_cycle = $a[1]; }
			if ($a[0] =~ "shift cycle")	{ $shift_cycle = $a[1]; }
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

=head1 parse_line

 parse lines of mapping info for one window, output file name and position info for plot

=cut
sub parse_line
{
	my ($report_info, $mapping_info, $cycle_size, $output_prefix) = @_;
	chomp($report_info);
	chomp($mapping_info);
	my @report = split(/\t/, $report_info);
	die "Error in $report_info\n" unless scalar(@report) == 8;
	my ($ref_id, $start, $end, $uniq_mapped, $in_cycle, $register_num, $base_shift, $pvalue) = @report;
	my @line = split(/\n/, $mapping_info);

	my $output_pdf = $output_prefix."-".$ref_id."-".$start."-".$end.".pdf";
	my %position_plus_in_cycle;
	my %position_plus_register;
	my %position_minus_in_cycle;
	my %position_minus_register;

	my ($pos, $cpos, $id, $exp, $strand, $register, $seq, $length);

	foreach my $line (@line)
	{
		#ref_id, position, strand	sRNA_id-count	seq
		#Solyc01g005780.1.1,948,1	SLY1A1168334-1	AGTGTTAGACTTGAGCAATAA
		my @n = split(/\t/, $line);
		my @m = split(/,/, $n[0]);
		
		$pos = $m[1];
		$strand = $m[2];
		$id = $n[1];
		my @p = split(/-/, $id);
		$exp = $p[scalar(@p)-1];
		if ($exp > 1) { $exp = log($exp)/log(10); }
		else { $exp = 0.1; }
		$seq = $n[2];
		$length = length($n[2]);

		# only output the reads with length equal cycle size (21nt or 24nt)
		if ($length ne $cycle_size) { next; }

		# cpos : compute position, it is pos - 2 when strand is minus
		if ($strand eq '-1') { $cpos = $pos + 2; }
		else { $cpos = $pos;}

		# check if the read is on register
		$register = 0;
		if ( (($cpos - ($start + $base_shift - 1)) % $cycle_size ) == 0 ) { $register = 1; }
		
		if ($strand eq '-1')
		{
			if  ($register == 1) { $position_minus_register{$pos} = "-".$exp; } 
			else                 { $position_minus_in_cycle{$pos} = "-".$exp; }
		}
		elsif ($strand eq '1')
		{
			if   ($register == 1) { $position_plus_register{$pos} = $exp; }
			else                  { $position_plus_in_cycle{$pos} = $exp; }
		}
	}

	my ($position_info, 
		$mc,  # minus_in_cycle number
		$mr,  # minus_register_number
		$pc,  # plus_in_cycle_number
		$pr   # plus_register_number
		);

	for(my $i=$start-2; $i<=$end; $i++)
	{
		if (defined $position_minus_in_cycle{$i} ) { $mc = $position_minus_in_cycle{$i}; } else { $mc = 0; }
		if (defined $position_minus_register{$i} ) { $mr = $position_minus_register{$i}; } else { $mr = 0; }
		if (defined $position_plus_in_cycle{$i} )  { $pc = $position_plus_in_cycle{$i};  } else { $pc = 0; }
		if (defined $position_plus_register{$i} )  { $pr = $position_plus_register{$i};  } else { $pr = 0; }

		$position_info.= "$i\t$pr\t$pc\t$mr\t$mc\n";
	}
	#print $position_info; die;

	return ($output_pdf, $position_info);
}


=head1 plot_it

 plot images using R

=cut
sub plot_it
{
	my ($pdf, $table) = @_;
	
	open(DD, ">temp.table") || die "Can not open temp.table $!\n";
	print DD $table;
	close(DD);

	#if ($strand eq "1") { $segment = "segments(0, 0, 21, 0, col= 'red4')";} 
	#else { $segment = "segments(210, 0, 231, 0, col= 'red4')";  }

	my $name = $pdf; $name =~ s/\.pdf//;

my $R_LD =<< "END";
a<-read.table("temp.table")
x<-a[,1]
plusk<-a[,2]
plusn<-a[,3]
minusk<-a[,4]
minusn<-a[,5]
pdf("$pdf",width=12,height=6)
y_range<-range(plusk, plusn, minusk, minusn)
plot(minusn, col="green", type="h", ylim=y_range, main="$name", xlab="Position", ylab="Expression LOG10(Count)")
lines(minusk, col="blue", type="h")
lines(plusn, col="orange", type="h")
lines(plusk, col="red", type="h")
invisible(dev.off())
END

	open R,"|/usr/bin/R --vanilla --slave" or die $!;
	print R $R_LD;
	close R;

	unlink("temp.table");
}


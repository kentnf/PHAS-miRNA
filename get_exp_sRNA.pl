#!/usr/bin/perl

=head1

 get expression of small RNA uniq reads;

=cut
use strict;
use warnings;
use IO::File;

my $usage  = qq'
usage: perl $0 list_uniq_read output_prefix 

* the input of this script is the uniq read of small RNA (cleaned)
  it mut have the expression of reads separated by \'-\' in ID.
 
>IP31L_0000001-23814  	# 23814 is the number of read TTCCACAGCTTTCTTGAACTG
TTCCACAGCTTTCTTGAACTG

* this script will use a lot of memory depend on the number of input files

* the 1st output is raw expression ( prefix_sRNA_expr )
* the 2nd output is library size (prefix_sRNA_libsize )
* the 3rd output is normalized exp ( prefix_sRNA_expTPM )
* the 4th output is sequence with > 5 TPM 
';

my $list_uniq_read = shift || die $usage;
my $output = shift || die $usage;

my $count_cutoff = 10;
my $norm_cutoff = 5;

my $output1 = $output."_sRNA_expr";
my $output2 = $output."_sRNA_libsize";
my $output3 = $output."_sRNA_expTPM";
my $output4 = $output."_sRNA_seq";

# put list of uniq small RNA reads to array
my @list;
my $fh = IO::File->new($list_uniq_read) || die "Can not open list file $list_uniq_read $!\n";
while(<$fh>)
{
	chomp;
	push(@list, $_);
}
$fh->close;

# main expression and libsize value to hash
my %uniq_read;
my %libsize;

foreach my $file (@list)
{
	my $total_read = 0;
	my $fu = IO::File->new($file) || die "Can not open uniq read file $file $!\n";
	while(<$fu>)
	{
		chomp;
		my $id = $_;
		my @a = split(/-/, $id);
		my $exp = $a[1];
		my $seq = <$fu>; chomp($seq);
		$uniq_read{$seq}{$file} = $exp;
		$total_read = $total_read + $exp;
	}
	$fu->close;
	$libsize{$file} = $total_read;
}

# output the libsize 
my $out2 = IO::File->new(">$output2") || die "Can not open output libsize $output2 $!\n";
foreach my $k (sort keys %libsize) { print $out2 $k."\t".$libsize{$k}."\n"; }
$out2->close;

# output the expression 
my $out1 = IO::File->new(">$output1") || die "Can not open output expression $output1 $!\n";
my $out3 = IO::File->new(">$output3") || die "Can not open output expression TPM $output3 $!\n";
my $out4 = IO::File->new(">$output4") || die "Can not open output small RNA $output4 $!\n";

print $out1 "UniqRead";
print $out3 "UniqRead";
foreach my $file (@list)
{
	print $out1 "\t".$file;
	print $out3 "\t".$file;
}
print $out1 "\n";
print $out3 "\n";

my $seq_order = 0;

my $length = length(scalar(keys(%uniq_read)));

foreach my $seq (sort keys %uniq_read)
{
	my $line = $seq;
	my $line_tpm = $seq;
	my $is_express = 0;
	foreach my $file (@list)
	{
		if ( defined $uniq_read{$seq}{$file} )
		{
			$line.="\t$uniq_read{$seq}{$file}";

			my $tpm = $uniq_read{$seq}{$file} / $libsize{$file} * 1000000;
			$tpm = sprintf("%.2f", $tpm);
			$line_tpm.="\t$tpm";

			$is_express = 1;
			#if ($uniq_read{$seq}{$file} >= 50 ) { $is_express = 1; }
			#if ($tpm > $norm_cutoff){ $is_express = 1; }
		}
		else
		{
			$line.="\t0";
			$line_tpm.="\t0";
		}
	}

	if ($is_express) { 
		print $out1 $line."\n";
		print $out3 $line_tpm."\n"; 
	}

	$seq_order++;
	my $zlen = $length - length($seq_order);
	my $z = "0"x$zlen;
	my $seq_id = "sRNA".$z.$seq_order;

	print $out4 $seq_id."\n".$seq."\n";
}

$out1->close;
$out3->close;
$out4->close;

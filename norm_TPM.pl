#!/usr/bin/perl

=head1

 the raw expression of small RNA is normalized using TPM method

=cut

use strict;
use warnings;
use IO::File;

my $usage = qq'
perl $0 sRNA_expr libsize[option] > output

* this script will normalize the raw expression to TPM
* if the libsize do not provided, this script will sum
  all raw count of sample as libsize
';

my $raw_count = shift || die $usage;
my $libsize = shift;

#################################################################
# load libsize to hash						#
#################################################################
my %lib;

#################################################################
# library size from file or expression file			#
#################################################################
if ($libsize) 
{
	my $fh = IO::File->new($libsize) || die "Can not open libsize file $libsize $!\n";
	while(<$fh>)
	{
        	chomp;
		my @b = split(/\t/, $_);
        	$lib{$b[0]} = $b[1];
	}
	$fh->close;
}
else
{
	my $fh = IO::File->new($raw_count) || die "Can not open raw count file $raw_count $!\n";
	my $title = <$fh>; chomp($title);
	my @title = split(/\t/, $title);
	while(<$fh>)
	{
		chomp;
		my @a = split(/\t/, $_);
		for(my $i=1; $i<@a; $i++)
		{
			my $tt = $title[$i];
			if (defined $lib{$tt})
			{
				$lib{$tt} = $lib{$tt} + $a[$i];
			}
			else
			{
				$lib{$tt} = $a[$i];
			}
		}
	}
	$fh->close;
}

#################################################################
# parse the expression again to save memory			#
#################################################################

my $in = IO::File->new($raw_count) || die "Can not open raw count file $raw_count $!\n";
my $title = <$in>; chomp($title);
my @title = split(/\t/, $title);

# output The results title
print $title."\n";

while(<$in>)
{
	chomp;
	my @a = split(/\t/, $_);

	print $a[0];

	for(my $i=1; $i<@a; $i++)
	{
		my $tt = $title[$i];	
		my $tpm;
		if ( $a[$i] > 0 )
		{
			$tpm = $a[$i]/$lib{$tt} * 1000000;
			$tpm = sprintf("%.2f", $tpm);
		}
		else
		{
			$tpm = 0;
		}

		print "\t".$tpm;
	}
	print "\n";
}
$in->close;

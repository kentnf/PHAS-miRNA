#!/usr/bin/perl

=head1

 function: combine list of sRNA_expr into one file

=cut

use strict;
use warnings;
use IO::File;

my $usage = qq'
perl $0 list > output

* the list file include sRNA_expr file per line
* this script will use large memory
';

my $list = shift || die $usage;

#################################################################
# checking the list file					#
#################################################################
my @list;
my $fh = IO::File->new($list) || die "Can not open list file $list $!\n";
while(<$fh>)
{
	chomp;
	my $file = $_;
	if (-s $file) { push(@list, $file); }
	else { die "Error, can not find sRNA_expr file $file $!\n"; }
}
$fh->close;

#################################################################
# expression info to hash					#
#################################################################
my %title;
my %uniq_expr;

foreach my $file (@list)
{
	my $fh = IO::File->new($file) || die "Can not open sRNA_expr file $file $!\n";
	my $title = <$fh>; chomp($title);
	my @title = split(/\t/, $title);

	for(my $i=1; $i<@title; $i++) {
		my $title_name = $title[$i];
		unless (defined $title{$title_name}) { $title{$title_name} = 1; }
	}

	while(<$fh>)
	{
		chomp;
		my @a = split(/\t/, $_);
		my $seq = $a[0];
		for(my $i=1; $i<@a; $i++)
		{
			my $t = $title[$i];
			if ($a[$i] > 0 ) { $uniq_expr{$seq}{$t} = $a[$i]; }
		}
	}
	$fh->close;
}


#################################################################
# output the combined results					#
#################################################################

# print title
print "sRNA";
foreach my $t (sort keys %title) { print "\t".$t; }
print "\n";

# print expression table
foreach my $seq (sort keys %uniq_expr)
{
	my $line = $seq;
	my $print_line = 0;

	foreach my $t (sort keys %title)
	{
		if (defined $uniq_expr{$seq}{$t})
		{
			$line.="\t".$uniq_expr{$seq}{$t};
		}
		else
		{
			$line.="\t0";
		}
	}	
	print $line."\n";
}

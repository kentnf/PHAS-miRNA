#!/usr/bin/perl

#===============================================================#
# This script processes miRNA_temp output file from		#
# the analysis pipeline and generate miRNA mysql table 		#
#===============================================================#

my $usage = qq'
usage: perl miRNA_sql.pl miRNA_conserved_sql miRNA_temp miRNA_mapping miRNA_nomatch output

';


if (@ARGV < 5) {
	print $usage;
	exit(0);
}

$input1 = $ARGV[0];
$input2 = $ARGV[1];
$input3 = $ARGV[2];
$nomatch= $ARGV[3];
$output = $ARGV[4];


# load miRNA mapping info to hash
open(IN3, $input3) || die "can't open miRNA_temp file $ARGV[2] $!\n";
while(<IN3>) {
	chomp;
	@a = split "\t";
	$miRNA_ID{$a[0]} = $a[1];
}
close(IN3);

# using mapping info to create nomatch information
# load the miRNA_conserved_sql to map hash
open(IN1, $input1) || die "can't open miRNA_conserved_sql $ARGV[0] $!\n";
open(NO, ">$nomatch");
while(<IN1>) {
	chomp;;
	@a = split "\t";
	$a[1] =~ s/\-/#/;
	@b = split("#", $a[1]);
	if (!$miRNA_ID{$b[1]}) { print NO $_, "\n"; }
	if (!$map{$a[0]}) { $map{$a[0]} = $miRNA_ID{$b[1]}; }
	else {
		if ($map{$a[0]} !~ /$miRNA_ID{$b[1]}/) { $map{$a[0]} .= "/".$miRNA_ID{$b[1]}; }
	}
}
close(IN1);
close(NO);

# output information
open(IN2, $input2) || die "can't open $ARGV[1] $!\n";
open(OUT, ">$output") || die "can't open $ARGV[4] $!\n";

while(<IN2>) {
	chomp;
	@a = split "\t";
	if (!$map{$a[0]}) { $map{$a[0]} = "NA";	}
	print OUT "$a[0]\t$a[1]\t$map{$a[0]}\n";
}
close(IN2);
close(OUT);


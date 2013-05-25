#!/usr/bin/perl

=head1

predict miRNA star sequence

=cut

use FindBin;
use IO::File;
use Getopt::Long;

my $usage = qq'
usage: ex_miR_star_bowtie.pl [-a miRNA_sql_output -b miRNA_hairpin_sql | -c hairpin ] -d sRNA_seq -e sRNA_expr -f output_file
 
	-a miRNA_sql		miRNA_sql_output file	
	-b miRNA_hairpin_sql 	miRNA_hairpin_sql file
	-c hairpin 		prase miRNA_sql_output and miRNA_hairpin_sql to generate hairpin file
	-d sRNA_seq 		input small RNA sequence [required]
	-e sRNA_expr 		expression of small RNA sequence [required]
	-f output_file		output file [required]
	-h 			print help info

';

my ($miRNA_sql, $hairpin_sql, $hairpin, $sRNA_seq, $sRNA_expr, $output);

GetOptions(
        "h|?|help"		=> \$help,
        "a|miRNA-sql=s"		=> \$miRNA_sql,
        "b|hairpin-sql=s"	=> \$hairpin_sql,
        "c|hairpin=s"		=> \$hairpin,
        "d|small-RNA=s"		=> \$sRNA_seq,
        "e|sRNA-expr=s"		=> \$sRNA_expr,
        "f|output=s"		=> \$output
);

die $usage if $help;
die $usage unless $sRNA_seq;
die $usage unless $sRNA_expr;
die $usage unless $output;
if ($miRNA_sql && $hairpin_sql) {} 
elsif ($hairpin) { }
else { die $usage; }

#################################################################
# Configuration section                                         #
#################################################################

$BOWTIE_PATH = ${FindBin::RealBin};
$TEMP_PATH = ".";
$gap = 2;    # distance from end position of hairpin

#################################################################


my $star_candidate = "Star_candidate";

my $temp_fas = $TEMP_PATH."/"."bowtie.tmp.fas.";
my $temp_db = $TEMP_PATH."/"."bowtie.tmp.db.";
unless ($hairpin) { $hairpin = $TEMP_PATH."/"."HAIRPINC"; }


my ($second, $minute, $hour);
($second, $minute, $hour) = localtime(); 
$fid = $second.$minute.$hour;
$temp_fas .= $fid; 
$temp_db .= $fid;

#################################################################
# constructure temp hairpin sequence file 			#
#################################################################
if ($miRNA_sql, $hairpin_sql)
{
	my %m; # key: miRNA new ID; value: miRNA old ID from small RNA sequence;
	my $in1 = IO::File->new($miRNA_sql) || die "Can not open miRNA sql file $miRNA_sql $!\n";
	while(<$in1>)
	{
		my @a = split(/\t/, $_);
		$m{$a[0]} = $a[1];
	}
	$in1->close;

	# H000001 Test01029       Chr1    13252204        13252419        193     +       Hairpin_Sequence        -191.60
	my $out = IO::File->new(">".$temp_fas) || die "Can not open temp fasta file $temp_fas $!\n";	
	my $out2 = IO::File->new(">".$hairpin) || die "Can not open hairpin file $hairpin $!\n";
	my $in2 = IO::File->new($hairpin_sql) || die "Can not open miRNA hairpin sql file $hairpin_sql $!\n";
	while(<$in2>)
	{
        	my @b = split(/\t/, $_);
        	$b[1] = $m{$b[1]};
		print $out2 join("\t", @b);
		my $len_seq = length($b[7]);
                print $out ">$b[0]_$b[1]_$b[5]_$len_seq\n";
                my $seq = $b[7];
                $seq =~ tr/uU/tT/;
                print $out "$seq\n";
	}
	$in2->close;
	$out->close;
	$out2->close;

}
elsif ($hairpin)
{
	# H000001 Test01029       Chr1    13252204        13252419        193     +       Hairpin_Sequence        -191.60
	my $out = IO::File->new(">".$temp_fas) || die "Can not open temp fasta file $temp_fas $!\n";
	my $in = IO::File->new($hairpin) || die "Can not open converted hairpin file $hairpin $!\n";
	while(<$in>)
	{
		chomp;
		my $list = $_;
		my @a = split("\t", $list);
		my $len_seq = length($a[7]);
		print $out ">$a[0]_$a[1]_$a[5]_$len_seq\n";
		my $seq = $a[7];
  		$seq =~ tr/uU/tT/;
		print $out "$seq\n";
	}
	$in->close;
	$out->close;
}
else
{
	print $usage;
	exit(0);
}

#################################################################
# compare small RNA with hairpin fasta sequence			#
#################################################################
my $bowtie_db_out = `$BOWTIE_PATH/bowtie-build $temp_fas $temp_db`;
my $bowtie_run_out = `$BOWTIE_PATH/bowtie -v 0 -a -f $temp_db $sRNA_seq`;

#################################################################
# convert bowtie run output to star candidate 			#
#################################################################
open OUTFILE, ">$star_candidate" || die "Cannot Open start candidate file $!\n";

my @list = split(/\n/, $bowtie_run_out);
foreach my $line (@list) {
      chomp $line;
      my ($sRNA, $strand, $miR, $loc, $seq, $sc, $mismatch) = split(/\t/, $line);
      if ($strand eq "+") {
  		$sRNA_s = $loc+1;
  		$sRNA_e = $loc+length($seq);

  		#$miR = $L[0];
  		#$sRNA_b = $L[1];

  		#$sRNA_s = $L[2];
  		#$sRNA_e = $L[3];

  		@miR_L = split("_", $miR);
  		$miR_id = $miR_L[0];
  		$sRNA_a = $miR_L[1];
  		$miR_s = $miR_L[2];
  		$miR_l = $miR_L[3];

  		$miR_e = $miR_l - $gap;

        	if (int($miR_s) <=  $gap) {      
                	if ($sRNA_e >= $miR_e) {
                        	print OUTFILE "$miR_id\t$sRNA_a\t$sRNA\t$miR_s\t$miR_l\t$sRNA_s\t$sRNA_e\n";
                	}
        	}       
        	else {
                	if ($sRNA_s <=  $gap) {
                        	print OUTFILE "$miR_id\t$sRNA_a\t$sRNA\t$miR_s\t$miR_l\t$sRNA_s\t$sRNA_e\n";
                	}
        	}
  	}
}
close OUTFILE;

system "rm $TEMP_PATH/bowtie.tmp.*";

#################################################################
# Ratio between miR mature read and miR star read		#
# Minmum read abundence						#
# Difference between miR mature length and miR star length 	#
#################################################################

$ratio_cutoff = 2.0;
$min_freq_a_cutoff = 2;
$min_freq_b_cutoff = 1;
$len_diff_cutoff = 3;     

#################################################################

my %miRNA_seq = {};
my %sRNA_seq = {};
my %sRNA_freq = {};

#################################################################
# $hairpin sequence to hash					#
#################################################################
print "Loading miR hairpin....\n";
%miRNA_seq = load_hairpin($hairpin);

=head1 load_hairpin

 function: load miR hairpin sequence to hash

=cut
sub load_hairpin
{
	my $hairpin = shift;

	my %miRNA_seq;

	my $fh = IO::File->new($hairpin) || die "Can not open hairpin sequence file $hairpin $!\n";
	while(<$fh>)
	{
		chomp;
		my $line = $_;
		my @a = split("\t", $line);
		my $m_id = $a[0];
		my $seq = $a[7];
		$miRNA_seq{$m_id} = $seq;
	}
	$fh->close;

	return %miRNA_seq;
}

#################################################################
# small RNA sequence to hash					#
#################################################################
print "Loading sRNAs....\n";
%sRNA_seq = load_sRNA($sRNA_seq);

=head1 load_sRNA

 function: load small RNA sequence to hash

=cut
sub load_sRNA
{
	my $sRNA_seq = shift;

	my $s_id = "";
	my %sRNA_seq;

	my $fh = IO::File->new($sRNA_seq) || die "Can not open small RNA sequence $sRNA_seq $!\n";
	while(<$fh>)
	{
		chomp;
		$line = $_;
		
		if ($line =~ m/^\>/) 
		{
			$s_id = $line;
     			$s_id =~ s/^\>//;
  		}
  		else
		{
     			$sRNA_seq{$s_id} = $line;
  		}
	}
	$fh->close;

	return %sRNA_seq;
}

#################################################################
# load small RNA expression number to hash			#
#################################################################
print "Loading read numbers....\n";
%sRNA_freq = load_read_number($sRNA_expr);

=head1 load_read_number

 function: load read expression number to hash

=cut

sub load_read_number
{
	my $sRNA_expr = shift;
	
	my %sRNA_freq;

	my $fh = IO::File->new($sRNA_expr) || die "Can not open small RNA expr file $sRNA_expr $!\n";
	while(<$fh>)
	{
		chomp;
  		my $line = $_;
  		my @a = split("\t", $line);
  		my $s_id = $a[0];
  		shift @a;
  		$sRNA_freq{$s_id} = join("\t", @a);
	}
	$fh->close;

	return %sRNA_freq;
}

#################################################################
# checking miRNA start						#
#################################################################
print "Checking miRNA star....\n";

open(IN, $star_candidate) || die "can't open $star_candidate $!\n";
open(OUT, ">$output") || die "can't open $output $!\n";

my @exp = split(/\t/,$sRNA_freq{'sRNA'});
my $exp_title = '';
foreach my $sample (@exp) { $exp_title.="\t$sample:A\t$sample:B\t$sample:ratio"; }

print OUT "miR_id\tsRNA_a\tsRNA_b\tmiR_s\tmiR_l\tsRNA_s\tsRNA_e$exp_title\n";

while(<IN>)
{
  	chomp;
  	$line = $_;
	my $exp_line = "";
  	@a = split("\t", $line);
  	$miR_id = $a[0];
  	$sRNA_id_a = $a[1];
  	$sRNA_id_b = $a[2];

	# check the hash
	unless(defined $miRNA_seq{$miR_id}) { die "Error in miRNA_seq"; }
	unless(defined $sRNA_seq{$sRNA_id_a}) { die "Error in sR_seq_a"; }
	unless(defined $sRNA_seq{$sRNA_id_b}) { die "Error in sR_seq_b"; }
	$miR_seq = $miRNA_seq{$miR_id};
	$sR_seq_a = $sRNA_seq{$sRNA_id_a};
	$sR_seq_b = $sRNA_seq{$sRNA_id_b};
	#"$sR_seq_a\n$sR_seq_b\n";
	
	unless(defined $sRNA_freq{$sR_seq_a}) { die "Error in sRNA_freq\n>$sRNA_id_a\n$sR_seq_a\n"; }
	unless(defined $sRNA_freq{$sR_seq_b}) { die "Error in sRNA_freq\n>$sRNA_id_b\n$sR_seq_b\n"; }
	#print "$sRNA_id_a\n$sRNA_freq{$sR_seq_a}\n$sRNA_id_b\n$sRNA_freq{$sR_seq_b}\n"; die;

	@freq_a = split("\t", $sRNA_freq{$sR_seq_a});
	@freq_b = split("\t", $sRNA_freq{$sR_seq_b});

  	my $ratio_avail = 0;
	my $ratio;
	# code for question
	for ($i = 0; $i < @freq_a; $i++)
	{
		if ($freq_b[$i] > 0 && $freq_a[$i] / $freq_b[$i] < $ratio_cutoff)
		{
			#$ratio_avail = 0;
			#last;
		}

    		if ($freq_b[$i] > 0 && $freq_a[$i] / $freq_b[$i] >= $ratio_cutoff)
		{
			$ratio_avail = 1;
		}

		if ($freq_b[$i] > 0 && $freq_a[$i] > 0) 
		{
			$ratio = $freq_a[$i] / $freq_b[$i];
			$ratio = sprintf("%2.f", $ratio);
		}
		else
		{
			if ($freq_a[$i] > 0)
			{
				$ratio = "max";
			}
			else 
			{
				$ratio = 0;
			}
		}
		$exp_line.="\t".$freq_a[$i]."\t".$freq_b[$i]."\t".$ratio;
	}

  	$len_a = length($sR_seq_a);
  	$len_b = length($sR_seq_b);
  	$len_diff = abs($len_a-$len_b);

	if ( $ratio_avail == 1 && $len_diff <= $len_diff_cutoff){
    		print OUT $line.$exp_line."\n";		
		#print OUT join("\t", @freq_a), "\t";
		#print OUT join("\t", @freq_b), "\n";
  	}
}
close(IN);
close(OUT);


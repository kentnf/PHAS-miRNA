#!/usr/bin/perl -w

# script that calls miRcheck to evaluate miRNA candidates
# usage: evaluate_miRNA_candidates.pl sample_candidates_g_folded sample_candidates_g_folded_mirs
# input: takes output from RNAfold as input.  However
# output: 

use FindBin;
use lib ${FindBin::RealBin};
use miRcheck;

$hitfile =$ARGV[0]; # output from RNAfold, in which header lines contain sequence of putative miRNA and coordinate of putative miRNA in hairpin
$output = $ARGV[1];


#parameters for miRNA checking
$MAX_UNPAIR = 4; #max # unpaired bases in putative mir
$MAX_STAR_UNPAIR = 5;  #max # unpaired bases in putative mir*
$MAX_SIZEDIFFERENCE = 3; # max size difference between mir and mir*
$MAX_MIR_GAP = 2; # longest acceptable run of unpaired bases in mir
$MAX_STAR_GAP = 3; # longest acceptable run of unpaired bases in mir*
$MIN_FBACK_SIZE = 60; # shortest acceptable length of folback including mir and mir*
$MAX_MIR_AS_BULGE = 2; # maximum total # assymetrically unpaired bases in mir
$MIN_UNPAIR = 1; # minimum number of unpair bases in acceptable extended mirs/mir*s
$BP_EXTENSION = 4; # number of nt to extend mirs and mir*s

while (@ARGV){
    $thisarg = shift @ARGV;
    if ($thisarg eq "-unpair" ) {$MAX_UNPAIR=shift @ARGV;} 
    if ($thisarg eq "-star_unpair" ) {$MAX_STAR_UNPAIR=shift @ARGV;}
    if ($thisarg eq "-size_diff" ) {$MAX_SIZEDIFFERENCE=shift @ARGV;}
    if ($thisarg eq "-mir_bulge") {$MAX_MIR_GAP=shift @ARGV;}
    if ($thisarg eq "-star_bulge") {$MAX_STAR_GAP=shift @ARGV;}
    if ($thisarg eq "-fback_min") {$MIN_FBACK_SIZE=shift @ARGV;}
    if ($thisarg eq "-ass") {$MAX_MIR_AS_BULGE=shift @ARGV;}
    if ($thisarg eq "-min_unpair") {$MIN_UNPAIR=shift @ARGV;}
    if ($thisarg eq "-bp_ext") {$BP_EXTENSION=shift @ARGV;}
}
open (O,">$output");
open (M,">$output\_mature");
open (H,">$output\_hairpin");

open(F,$hitfile);
while(<F>){
    if (/>(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)/){
	# assumes input file contains an ID, a potential miRNA sequence ,and the 1st coordinate of potential miRNA within a genomic sequence
	$ID = $1;
	$chr = $2;
	$genomic_start = $3;
	$genomic_stop = $4;
	$miR_seq = $5;
	$miR_start = $6;
	$miR_stop = $miR_start + length($miR_seq)-1;
	$fold = '';
	$genomic_seq = '';
    }
    elsif (/(\S+) \(\S+\)/){
	$fold = $1;
#	print "$ID\t$miR_start\t$miR_stop\t$miR_seq\n$genomic_seq\n";
	#judge quality of foldback
	$fback = '';
	($fback,$fback_start,$fback_stop)  = miR_check($fold,$miR_start,$miR_stop,
						       "-unpair",$MAX_UNPAIR,
						       "-star_unpair",$MAX_STAR_UNPAIR,
						       "-size_diff",$MAX_SIZEDIFFERENCE,
						       "-mir_bulge",$MAX_MIR_GAP,
						       "-star_bulge",$MAX_STAR_GAP,
						       "-fback_min",$MIN_FBACK_SIZE,
						       "-ass",$MAX_MIR_AS_BULGE,
						       "-min_unpair",$MIN_UNPAIR,
						       "-bp_ext",$BP_EXTENSION
						       );
	print O  "$ID $fback $fback_start $fback_stop\n";
	$hairpin_seq = get_subseq($genomic_seq,$fback_start+1,$fback_stop+1);
	if ($genomic_stop > $genomic_start){
	    $fback_x = $genomic_start + $fback_start;
	    $fback_y = $genomic_start + $fback_stop;
	    $mir_x = $genomic_start + $miR_start;
	    $mir_y = $genomic_start + $miR_stop;
	}
	else{
	    $fback_x = $genomic_start - ($fback_start);
	    $fback_y = $genomic_start - ($fback_stop);
	    $mir_x = $genomic_start - $miR_start;
            $mir_y = $genomic_start - $miR_stop;
	}
	if ($fback =~ /prime/){
	    print H ">$ID\t$chr\t$fback_x\t$fback_y\n";
	    print H "$hairpin_seq\n";
	    
	    print M ">$ID\t$chr\t$mir_x\t$mir_y\n";
	    print M "$miR_seq\n";
	}	
    }
    elsif (/(\S+)/){$genomic_seq = $1;}
}


	

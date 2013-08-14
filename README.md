
=================
Quality Control for small RNA
=================

clip adapter sequence, remove low quality and N reads, as well as remove rRNA using bwa

	$(Shan's method)

remove rRNA and tRNA using bowtie

	convert the fastq reads to fasta format
	$fq2fa -i list
	
	the suffix of fasta file must be 'fasta'
	list all the fasta file to a list file, then remove suffix name
	
	remove the rRNA or other contamination using bowtie
	$rRNA_rm.pl -i list -s SS -r your_rRNA_index_file -v 0 -p 24

check small RNA length distribution

	$smallRNA_len_dist.pl input_seq

Remove redundancy and report count information for each unique sRNA.

	$get_uniq_read.pl input_seq

=================
Identify new miRNA from small RNA
=================

Combine all cleaned small RNA sequences. Update the count information.
Normalize the expression to TPM. Only keep those with TPM >= 5 in at least one sample

	$get_exp_sRNA.pl list_uniq_read[output of step 3] output_prefix

* combine the raw expression of small RNA.

	$combine_sRNA_expr.pl list > output

* normalize the raw expression to TPM.

	$norm_TPM.pl sRNA_expr libsize[option] > output

Run miRNA identification pipeline script. This pipeline script will:
Process to generate miRNA MySQL table file.
Identify miRNA targets.
Finding miRNAStar sequences.

	$miRNA_prediction.pl -s sRNA_test.fa -g input/cucumber_genome -b \
                        -t -u input/cucumber_mRNA_v2.fa -d input/cucumber_mRNA_v2.fa \
                        -r -e sRNA_exp_TPM

	* now the old step 5-8 are combined except cut -f1,3,6,7 mir_star_output > miRNA_star,
   	because I do not see the real star result

+-----------------------------------------------+
+ the identification could be run step by step  +
+-----------------------------------------------+

1) Run miRNA identification pipeline.

	$miRNA_prediction.pl -s sRNA_test.fa -g input/cucumber_genome

2) Process to generate miRNA MySQL table file.

   	$./script/miRNA_sql.pl miRNA_conserved_sql miRNA_temp mirBase/miRNA_mapping miRNA_nomatch miRNA_sql_output

3) Identify miRNA targets.

	$./script/$miRNA_target_pred.pl miRNA_seq mRNA_seq_index[formatdb] mRNA_seq miRNA_target_output

4) Finding miRNAStar sequences.

	$./script/$miR_star.pl -a miRNA_sql_output -b hairpin_sql -d sRNA_seq -e sRNA_expr -f miRNA_star_output
   	or
   	$./script/$miR_star.pl -c hairpin -d sRNA_seq -e sRNA_expr -f miRNA_star_output

* hairpin is the file after pasing miRNA_sql_output and hairpin_sql using script hairpin_process.pl

usage: perl miRNA_prediction.pl  [options]

  -s|small-RNA          small RNA sequence (required)
  -g|genome             genome sequnece & genome index created by bowtie-build (required)
  -m|multi-hit          maximum of multiple hits for one small RNAs (default = 20)
  -b|sql                generate miRNA sql file (default = off)
  -t|target-predict     predict target from mRNA (default = off)
  -u|mRNA               mRNA sequence (required if -t|target-predict is on)
  -d|mRNA-index         mRNA index created by formatdb
  -r|miRNA-star         predict miRNA star (default = off)
  -e|sRNA-expr          expression of small RNA (required if -r|miRNA-star is on)
  -o|output-prefix      the prefix of output file (default = output)
  -h|?|help             print help info

* the input small RNA is : cleaned, uniq, remove rRNA and tRNA
* genome is the genome fasta file, as well as index build by bowtie-build
	
=================
View the expression of small RNA with fixed/all different length
=================

For uniq cleaned reads alignment with genome/reference, add reads alignment info base on expression
For cleaned reads (not uniq) alignment with genome/reference, skip this step.

	$filter_SAM.pl -i sample.sam -e -o sample_exp.sam

For view the expression of small RNA with fixed length (set 21nt as example), keep the 21nt alignment.

	$filter_SAM.pl -i sample.sam -l 21 -o sample_21.sam

Sort the output sam file, than convert alignment files to bigwig format.
File Reference size include ID and length of reference sequence.
Constant factor is for normalization of raw count.

	$genomeCoverageBed -bg -ibam sample_exp.bam -g reference_size -scale constant_factor > sample_exp.bedgraph
	$wigToBigWig sample_exp.bedgraph reference_size sample_exp.bw
		
Or we can convert list of bam files to bw format using script.

	$sam2bw.pl -i 

* All the bw files, bam files (indexed) could be viewed in IGV.

=================
Identify pre-ta-siRNA using small RNA
=================

align uniq cleaned small RNA reads to genome/reference

	$rRNA_rm.pl -i list -r reference -t 2 -v 0 -k 7

remove the unmapped reads and multi=hit reads (>6)

	$filter_SAM.pl -i input.sam -n
	$filter_SAM.pl -i input.sam -m 6

Convert sam to bam file, sort bam file, save the sorted bam file.
Then identify pre-ta-siRNA using small RNA

	$pretasiRNA_iden.pl -i input.sam -r reference

The output of identification include four files:
1) input_pvalue.txt : 
2) input_phasingscore.bedgraph : 
3) input.report.gz : 
4) input_4image.txt :_

Draw images for identified pre-ta-siRNA with reference and mapped sRNA.

	$pretasiRNA_draw1.pl input_4image.txt reference

Draw images for identified pre-ta-siRNA with expression of in cycle and register reads.

	$pretasiRNA_draw2.pl input_4image.txt

=================
To be continue
=================

=================
Old Steps to Identify new miRNA from small RNA
=================
5. Run miRNA identification pipeline (all genome and sRNA sequence files are 
   in the "input" folder)
   perl miRNA_prediction.pl sRNA_seq genome_seq output
   Note: in miRBase, change all "U" to "T" before bowtie-build

   * done by JG

6. Process to generate miRNA MySQL table file
   perl miRNA_sql.pl miRNA_conserved_sql miRNA_temp miRNA_mapping output

7. Identify miRNA targets (plus strand only - $choice_direction = "forward" 
   in the program)
   formatdb the cDNA database

   perl miRNA_target_pred.pl miRNA_seq_fasta formatdb_name cDNA_db_name output

   filter out some targets that do not meet the requirement (this can be integrated 
   into the above script)

   perl miRNA_target_filter.pl miRNA_target >output

8. Finding miRNA* sequences. (codes for this part may need update)

   perl hairpin_process.pl miRNA_hairpin >HAIRPIN

   perl miR_star.pl HAIRPIN[could combine] sRNA_seq[input small RNA] sRNA_expr[expression] mir_star[output]

   cut -f1,3,6,7 mir_star >miRNA_star

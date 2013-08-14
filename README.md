
=================
Identify pre-ta-siRNA using small RNA
=================





1. Remove adaptor sequences. (Shan's method)

2. check small RNA length distribution
   
   $smallRNA_len_dist.pl input_seq

3. Remove redundancy and report count information for each unique sRNA. 

   $get_uniq_read.pl input_seq

4. Combine all cleaned small RNA sequences. Update the count information. 
   Normalize the expression to TPM. Only keep those with TPM >= 5 in at least one sample

   $get_exp_sRNA.pl list_uniq_read[output of step 3] output_prefix

   * combine the raw expression of small RNA.
   $combine_sRNA_expr.pl list > output

   * normalize the raw expression to TPM.
   $norm_TPM.pl sRNA_expr libsize[option] > output

5. Filter the rRNA and tRNA

   using bowtie to compare the output of step 4 with tRNA and rRNA database.
   Keep the un-mapped reads.

6. Run miRNA identification pipeline.
   Process to generate miRNA MySQL table file.
   Identify miRNA targets.
   Finding miRNA* sequences.
   ------ using miRNA_prediction.pl script

   $miRNA_prediction.pl -s sRNA_test.fa -g input/cucumber_genome -b \
                        -t -u input/cucumber_mRNA_v2.fa -d input/cucumber_mRNA_v2.fa \
                        -r -e sRNA_exp_TPM

   * now the old step 5-8 are combined except cut -f1,3,6,7 mir_star_output > miRNA_star,
   because I do not see the real star result

   +------------------------------------------------------------------------------------+
   + the step 6 could be run step by step						+
   +------------------------------------------------------------------------------------+

   6.1 Run miRNA identification pipeline.
   $miRNA_prediction.pl -s sRNA_test.fa -g input/cucumber_genome
   
   6.2 Process to generate miRNA MySQL table file.
   $./script/miRNA_sql.pl miRNA_conserved_sql miRNA_temp mirBase/miRNA_mapping miRNA_nomatch miRNA_sql_output

   6.3 Identify miRNA targets.
   $./script/$miRNA_target_pred.pl miRNA_seq mRNA_seq_index[formatdb] mRNA_seq miRNA_target_output
 
   6.4 Finding miRNA* sequences.   
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

    
===================================================================================
Old Step
===================================================================================
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
===================================================================================

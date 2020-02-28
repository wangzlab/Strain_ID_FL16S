# Strain_ID_FL16S
A simplified, automated pipeline for strain identification using Pacbio full-length 16S data.

1. Ensure blastn and Rscript is in your system path.

2. Make sure perl module Getopt::Long is install.

   $ cpanm Getopt::Long
   
3. Download 16S_DB.fa, decompress the file, and index the BLASTn database:

   $ unzip 16S_DB.fa.zip
   
   $ makeblastdb -dbtype nucl -in 16S_DB.fa
   
   The 16S_DB.fa contains unique 16S sequences from all 14,062 complete bacterial genomes in Genbank (accessed May 9, 2019)
   
4. Download the 16S_ratio.txt and Strain_ID_FL16S.pl files and place them in the same directory.
   
   The 16S_ratio.txt contains the intragenomic 16S copy number ratio for each complete bacterial genomes with SSU IDs matching the 16S_DB.fa
   
5. Format input files (example files provided): 

   profile.txt contains the abundances (absolute count) for each ASV in each sample

   seqs.fa contains the full-length 16S sequences for each ASV
   
6. Run the script:
   running "perl Strain_ID_FL16S.pl" should print the help message.

   $ perl Strain_ID_FL16S.pl -profile profile.txt -seqs seqs.fa > output.txt

7. Optional parameters can be set include: (1) -threads: number of threads for BLASTn search (default 1), (2) -coef: Pearson correlation coefficient cutoff (default 0.7), (3) -deviation: allowed deviations between observed and genuine 16S copy number ratios (default 0.3)

8. The output.txt file is a table in which each line represents a pair of ASVs with 1) identical match to the 16S alleles of same bacterial genome, 2) showed correlation patterns across all samples and 3) in integral copy number ratio with the intragenomic 16S copy number ratio for the bacterial genome (±0.3).

 
   



# Strain_ID_FL16S
Strain_ID_FL16S is a simplified, automated pipeline for strain identification using Pacbio full-length 16S data. The script generates strain-level (or genome-level) amplicon sequence variant (ASV) bins by leveraging intra-genomic variations of 16S rRNA gene for bacterial genomes. For this purpose, the 16S sequences and their copy number ratios for 14,062 complete bacterial genomes were curated.

The detailed usage for the script is as follows:

1. Ensure blastn and Rscript is in your system path.

2. Make sure perl module Getopt::Long is install.

   $ cpanm Getopt::Long
   
3. Download 16S_DB.fa, decompress the file, and index the BLASTn database:

   $ tar -zxvf 16S_DB.tar.gz
   
   $ makeblastdb -dbtype nucl -in 16S_DB.fa
   
   The 16S_DB.fa contains unique 16S sequences from all 14,062 complete bacterial genomes in Genbank (accessed May 9, 2019) and will be updated periodically.
   
4. Download the 16S_ratio.txt and Strain_ID_FL16S.pl files and place them in the same directory.
   
   The 16S_ratio.txt contains the intragenomic 16S copy number ratio for each complete bacterial genomes with SSU IDs matching the 16S_DB.fa, along with species and strain-level taxonomy of each genome.
   
5. Format input files (example files provided): 

   profile.txt contains the abundances (absolute count) for each ASV in each sample, obtained from QIIME2/dada2 (with or without header info #Constructed from biom file, with or without taxonomy in the last column).

   seqs.fa contains the full-length 16S sequences for each ASV in FASTA format.
   
6. Run the script:
   running "perl Strain_ID_FL16S.pl" should print the help message.

   $ perl Strain_ID_FL16S.pl -profile profile.txt -seqs seqs.fa > output.txt

7. Optional parameters include: (1) -threads: number of threads for BLASTn search (default 1), (2) -coef: Pearson correlation coefficient cutoff (default 0.7), (3) -deviation: allowed deviations between observed and genuine 16S copy number ratios (default 0.3).

8. The output.txt file is a table in which each line represents a pair of ASVs with 1) identical match to the 16S alleles of same bacterial genome, 2) showed correlation patterns across all samples and 3) in integral copy number ratio with the intragenomic 16S copy number ratio for the bacterial genome (±0.3).

9. Note: If there are more than one pair of ASVs corresponding to the same genome (meaning there are more than two ASVs recovered from the genome), one can manually aggregate them into a single copy number ratio (i.e. SSU1:SSU2=2:1, SSU4:SSU1=3:1 => SSU4:SSU1:SSU2=6:2:1) and report the results.

10. The genomic abundances of each strain in each sample can be calculated as the average abundances of ASVs normalized by their copy numbers in the genome.
   



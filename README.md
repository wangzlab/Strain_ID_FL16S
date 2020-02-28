# Strain_ID_FL16S: an automated pipeline for strain identification using Pacbio full-length 16S data

1. Ensure blastn and Rscript is in your system path.
2. Download 16S_DB.fa, decompress the file, and index the BLASTN database:
   
   unzip 16S_DB.fa.zip
   makeblastdb -dbtype nucl -in 16S_DB.fa
   
   ### The 16S_DB.fa contains unique 16S sequences from all 14,062 complete bacterial genomes in Genbank (accessed May 9, 2019) ###
3. Download the 16S_ratio.txt and Strain_ID_FL16S.pl files and place them in the same directory.
   ### The 16S_ratio.txt contains the intragenomic 16S copy number ratio for each complete bacterial genomes with SSU IDs matching to the 16S_DB.fa ### 
4. Format your input files: 
   profile.txt which contains the abundances (absolute count) for each ASV in each sample
   seqs.fa which contains the full-length 16S sequences for each ASV
   
5. Run the script:
   running "perl Strain_ID_FL16S.pl" should print the help message:
   
   Description: Given the ASV abundance profile and sequences, this script will generate strain-level bins when possible.
   Usage: Strain_ID_FL16S.pl <options>
        Options:
                -(Required) profile the ASV abundance profile from i.e. QIIME2/dada2
                -(Required) seqs the ASV sequences in Fasta format
                -(Optional) coef Pearson correlation coefficient (default 0.7)
                -(Optional) threads number of threads for BLASTn (default 1)
                -(Optional) deviation max deviations between ASV and genuine 16S ratio (default 0.3)
   
   Make sure blastn and Rscript is in your system path
   Contact:     Zhang Wang (wangz@m.scnu.edu.cn)

   perl Strain_ID_FL16S.pl -profile profile.txt -seqs seqs.fa > output.txt
   
 6. The output file is a table in which each line represent a pair of ASVs with 1) identical match to the same bacterial genome, 2) showed correlation patterns across all samples (Pearson's R>0.7) and 3) in integral copy number ratio with the intragenomic 16S copy number ratio for the bacterial genome (±0.3).
 
 
   



#! bin/perl -w
use Getopt::Long;
use strict;

my $usage = qq~
Description: Given the ASV abundance profile and sequences, this
        script will generate strain-level bins when possible.
Usage: $0 <options>
        Options:
                -(Required) profile the ASV abundance profile from i.e. QIIME2/dada2
                -(Required) seqs the ASV sequences in Fasta format
                -(Optional) coef Pearson correlation coefficient (default 0.7)
		-(Optional) threads number of threads for BLASTn (default 1)
		-(Optional) deviation max deviations between ASV and genuine 16S ratio (default 0.3)

	Make sure blastn and Rscript is in your system path
	Contact:     Zhang Wang (wangz\@m.scnu.edu.cn)
~;

my ($profile,$seqs,$coef,$threads,$deviation)=undef;
$coef=0.7 unless $coef;
$threads=1 unless $threads;
$deviation=0.3 unless $deviation;

GetOptions(     'profile=s'=>\$profile,
		'seqs=s'=>\$seqs,
		'coefficient=s'=>\$coef,
		'threads=s'=>\$threads,
		'deviation=s'=>\$deviation)  || die $usage;
$profile || die $usage;
$seqs || die $usage;

my (%abundance, %nonsingle, %organism, %species, %strain, %copy, %genome, %revgenome)=();
print STDERR "Step 1: read in ASV abundance table\n";
open (IN, $profile);
my $header=<IN>;
if ($header =~ /Constructed from biom file/) {
	$header=<IN>;
}
chop $header;
my @headers=split("\t",$header);
while (<IN>) {
	chop;
	my @a=split("\t",$_);
	if (/k__Bacteria|k__Archaea/) {
		pop @a; ## filter taxonomy column
	}
	my $count=0;
	for my $i (1..$#a) {
		$abundance{$headers[$i]}{$a[0]}=$a[$i];
		if ($a[$i]>0) {
			$count ++;
		}
	}
	if ($count>=2) {
		$nonsingle{$a[0]}=1; # flag singletons
	}
}

print STDERR "Step 2: read in intragenomic 16S ratio\n";
open (IN, "16S_ratio.txt");
my $dump=<IN>;
while (<IN>) {
	chop;
	my @a=split("\t",$_);
	my @b=split(":",$a[1]);
	my @c=split(":",$a[2]);
	if (scalar @b>1) {
		for my $i (0..$#b) {
			for my $j (0..$#b) {
				$copy{$a[0]}{$b[$i]."_".$b[$j]}=$c[$i]/$c[$j];
			}
		}
	}
	else {
		$copy{$a[0]}{$a[1]}=$a[2];
	}
	$organism{$a[0]}=$a[3];
	$species{$a[0]}=$a[4];
	$strain{$a[0]}=$a[5];
}
print STDERR "Step 3: Perform BLASTn against database\n";
system ("/software/ncbi-blast-2.8.1+/bin/blastn -query $seqs -db 16S_DB.fa -outfmt 6 -out $$.btab -perc_identity 100 -qcov_hsp_perc 100 -num_threads $threads");

print STDERR "Step 4: Obtain initial taxonomic bins\n";
open (IN, "$$.btab");
while (<IN>) {
	chop;
	my @a=split("\t",$_);
	next if (! exists $nonsingle{$a[0]}); # remove singletons
	my ($gn, $variant)=($a[1]=~ /(GCA_\d+)_(SSU_\d+)/);
	$genome{$gn}{$a[0]}=$variant;
	$revgenome{$a[0]}{$gn}=$variant;
}
open (OUT1, ">$$.R"); # print temporary R code for cor.test
print OUT1 "data<-read.table(\"$$.txt\",row.names=1)\n";
print OUT1 "out<-data.frame(coef=1,slope=1)\n";
print OUT1 "out[1,1]<-unname(cor.test(data\$V2,data\$V3)\$estimate)\n";
print OUT1 "out[1,2]<-lm(data\$V2~data\$V3)\$coefficient[\"data\$V3\"]\n";
print OUT1 "print (out)\n";

print STDERR "Step 5: Obtain refined bins\n";
print "Genome_ID\tSpecies\tStrain\tASV1\tASV2\tASV1:ASV2\tSSU1\tSSU2\tSSU1:SSU2\n"; # print header
for my $key (keys %genome) {
	if (scalar keys %{$genome{$key}}>1) {
		my %select=();
		for my $key2 (keys %{$genome{$key}}) {
			for my $key3 (keys %{$genome{$key}}) {
				next if (exists $select{$key3}{$key2} or $key2 eq $key3);
				open (OUT, ">$$.txt"); # output abundance to temp file for cor test #
				for my $key4 (keys %abundance) {
					print OUT $key4."\t".$abundance{$key4}{$key2}."\t".$abundance{$key4}{$key3}."\n";
				}
				$select{$key2}{$key3}=1;
				system("Rscript $$.R > $$.result"); # correlation test
				open (IN, "$$.result");
				my $dump=<IN>;
				my $data=<IN>;
				my @data=split(" ",$data);
				next if ($data[1]<=$coef); # filter correlations by coefficient
				my $ssu1=$genome{$key}{$key2};
				my $ssu2=$genome{$key}{$key3};
				if ($data[2]<1) { # get correlation slope
					my $slope=1/$data[2]; #invert the slope if it is less than 1.0
					if (abs($slope-$copy{$key}{$ssu2."_".$ssu1})<$deviation*$copy{$key}{$ssu2."_".$ssu1}) {
						print $key."\t".$species{$key}."\t".$strain{$key}."\t".$key3."\t".$key2."\t".$slope."\t".$ssu2."\t".$ssu1."\t".$copy{$key}{$ssu2."_".$ssu1}."\n";
					}
				}
				else {
					my $slope=$data[2];
					if (abs($slope-$copy{$key}{$ssu1."_".$ssu2})<$deviation*$copy{$key}{$ssu2."_".$ssu1}) {
						print $key."\t".$species{$key}."\t".$strain{$key}."\t".$key2."\t".$key3."\t".$slope."\t".$ssu1."\t".$ssu2."\t".$copy{$key}{$ssu1."_".$ssu2}."\n";
					}

				}
			}
		}
	}
	else { # also print if perfect match between single ASV and genome with only one copy of 16S
		my @asv=keys %{$genome{$key}};
		if (scalar keys %{$copy{$key}} == 1 and scalar keys %{$revgenome{$asv[0]}}==1) { ## remove spurious cases with multiple identical strain-level hits
			print $key."\t".$species{$key}."\t".$strain{$key}."\t".$asv[0]."\t".$asv[0]."\t1\t".$genome{$key}{$asv[0]}."\t".$genome{$key}{$asv[0]}."\t".$copy{$key}{$genome{$key}{$asv[0]}}."\n";
		}
	}
}
system ("rm $$.*"); #clean up
print STDERR "Step 6: Done\n";

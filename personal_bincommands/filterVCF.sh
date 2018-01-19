#!/bin/bash -l
#SBATCH -J filterVCF
#SBATCH -t 12:00:00
#SBATCH -A snic001-12-255
#SBATCH -p core

# This script filters vcf files for SNPs with a minimum of reads in a minimum of individuals. As if 28.11.2013 it also eliminates all-het sites (sites with all individuals being heterozygotes).
# inputs are
# 1) number of individuals in the file (noInd), 
# 2) minimal number of individuals (indFilt) for which
# 3) a minimal number of reads (readFilt) is required,
# 4) the input vcf file (infile)
# run example: sbatch filterVCF 20 10 7 gt_I.vcf

# RUN: sbatch /gulo/glob/reto/scripts/bash/filterVCF.ssh 20 10 7 /proj/b2010010/repos/variation/200genomes/vcf/population_genotypes/gt_E.vcf

# It outputs a vcf that is automatically named according to input file name, and filtering criteria in the working folder

noInd=$1
indFilt=$2
readFilt=$3
infile=$4

b=$(basename $infile .vcf)
out=$b.R$readFilt.I$indFilt.vcf

#set population coordinates in vcf file
min=9
max=$(( $min+$noInd ))


head -50000 $infile | grep '^#' > $out
cat $infile | grep -v '^#' | perl -ane 'BEGIN{$READFILT=shift; $NOIND=shift; $INDFILT=shift; $MIN=shift; $MAX=shift;}  $mc=0; $ct=0; for ($i=$MIN;$i<$MAX;$i++) {@s=split(/\:/,$F[$i]); if($s[0] eq "\.\/\." || $s[1]< $READFILT ) {$mc++} }  $ct=$NOIND-$mc;  print if ($ct>=$INDFILT);'  $readFilt $noInd $indFilt $min $max | perl -ane 'BEGIN{$MIN=shift; $MAX=shift; $NOIND=shift;} $ch=0; for($i=$MIN;$i<$MAX;$i++){ @s=split(/\:/,$F[$i]); @g=split(/\//,$s[0]); if($g[0] ne $g[1]){$ch++};} print if ($ch < $NOIND)' $min $max $noInd >> $out
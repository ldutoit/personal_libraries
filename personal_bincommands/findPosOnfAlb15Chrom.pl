#!/usr/bin/perl

my $usage = "
# # # # # #
# findPosOnfAlb15Chrom.pl
# Version 1.3, Sept 2013 (original written in Dec 2012)
# Author: Linnéa Smeds 
# NEW v1.3: Added options to include unassigned scaffolds, and choose input
# format gtf and bed (set the position columns automatically)
# NEW v1.2: BUGFIXED! Before 12 Sept 2013, all positions on reverse scaffolds
# was one base too small (pos 1000 should be 1001, etc).    
# NEW v1.1 (Feb 2013): THIS VERSION IS UPDATED TO INCLUDE THE CORRECT SCAFFOLD N00498 
# INSTEAD OF THE WRONG N00546.
# ===========================================================================
# Takes a list with positions and finds the location on the concatenated 
# chromosome, given a certain gap size between chromosomes (and what scaffold
# that should be included - all or only the anchored or linked ones).
# 
# USAGE:
# perl findPosOnfAlb15Chrom.pl -in=file [-posCol=\"2,3\" -scafCol=1, -gap=5000 
#	-level=all|strict -out=newfile -h]

# All input flags: [default in brackets]
-in=file\t Any kind of table with scaffolds and positions (in version fAlb15!)
-posCol=\"n,n\"\t Column or columns with positions to translate. separated by
\t\t a comma \",\". [2]
-scafCol=n\t The column that contains the scaffold names. [first column]
-gap=n\t\t How many Ns that are counted between the scaffolds [5000]
-level=all|strict\t What scaffolds to use for the chromosomes, strict means
\t\t only scaffolds where we are absolutely sure or position/direction, all
\t\t means all possible scaffolds.
-un\t\t If set, regions in unassigned scaffolds will be printed to output,
\t\t unchanged. [not used]
-format=gtf|bed\t When used, scaffold and position columns are chosen automatically
\t\t according to the format standard, and for gtf the strand (col 7) is
\t\t reversed if the scaffold has a negative orientation. [not used]
-out=file\t Name of output file. [add suffix \"chrompos\" to infile]
-h\t\t Print this help message.
#============================================================================\n";

use strict;
use warnings;
use Getopt::Long;

my ($in,$posCol,$scafCol,$level,$unassigned,$format,$out,$help);
my $gapsize = 5000;	#Default value (updated this 20130227 because gap=0 was
					#interpreted as "undefined" when testing unless($gapsize)
GetOptions(
  	"in=s" => \$in,
   	"posCol=s" => \$posCol,
  	"scafCol=s" => \$scafCol,
	"gap=i" => \$gapsize,
	"level=s" => \$level,
	"un" => \$unassigned,
	"format=s" => \$format,
	"out=s" => \$out,
	"h" => \$help);

my $strictfile = "/proj/b2010010/private/Linkage/fAlb15_chromosomes_allWithLinksToCertain_FLIPPED_20130221.txt";
my $allfile = "/proj/b2010010/private/Linkage/fAlb15_chromosomes_withIntermediates_FLIPPED_20130221.txt";
my $choice;
#-------------------------------------------------------------------------------
#Checking input

if($help) {
	die $usage . "\n";
}
unless($in) {
	die "ERROR: No infile was given! \n" . $usage ."\n";
}
if($posCol) {
	if($format) {
		if($format eq "gtf"  && $posCol ne "4,5") {
			die "ERROR: Parameter conflict: -gtf is not compatible with -posCol=$posCol.\n\tLeave -posCol out and it will be set to \"4,5\"\n";
		}
		if($format eq "bed" && $posCol ne "2,3") {
			die "ERROR: Parameter conflict: -bed is not compatible with -posCol=$posCol.\n\tLeave -posCol out and it will be set to \"2,3\"\n";
		}
	}
}
else {
	$posCol = 2;
	
}	
if($scafCol) {
	if($format) {
		if($format eq "gtf" && $scafCol != 1) {
			die "ERROR: Parameter conflict: -gtf is not compatible with -scafCol=$scafCol.\n\tLeave -scafCol out and it will be set to 1\n";
		}
		if($format eq "bed" && $scafCol != 1) {
			die "ERROR: Parameter conflict: -bed is not compatible with -scafCol=$scafCol.\n\tLeave -scafCol out and it will be set to 1\n";
		}
	}
	$scafCol=$scafCol-1;
	
}
else {
	$scafCol = 0;
}

if($level) {
	if($level eq "strict"){
		$choice = $strictfile;
	}
	elsif($level eq "all"){
		$choice = $allfile;
	}
	else {
		die "ERROR: Flag -level must be either strict or all! \n".$usage."\n";
	}
}else {
	$choice = $strictfile;
}

if($unassigned){
	$unassigned="yes";
	#print "you've chosen to include unassigned\n";
}
else {
	$unassigned="no";
}
if($format) {
	#print "you've set format to $format!\n";
	if($format eq "gtf") {
		$posCol = "4,5";
	}
	elsif($format eq "bed") {
		$posCol = "2,3";
	}
	else {
		die "ERROR: Unknown format $format. Known formats are gtf and bed.\n".$usage."\n";
	}
}
else {
	$format="no";
}		 

unless($out) {
	$out= $in.".chrompos";
}
unless(-e $in) {
	die "ERROR: File $in doesn't exist!\n";
}
# ----------------------------------------------------------------------------------------
# 


# Make hash of the scaffolds
my %scaffs = ();
open(IN, $choice);
my $prev = "";
my $tmpstart=0;
while(<IN>) {
	my ($chr, $scaf, $len, $dir, $type, $color, $link) = split(/\s+/, $_);
	
	$scaffs{$scaf}{'chr'}=$chr;
	$scaffs{$scaf}{'len'}=$len;
	$scaffs{$scaf}{'dir'}=$dir;
	
	if($prev ne $chr) {
		$tmpstart = 0;
	}
	$scaffs{$scaf}{'offset'}=$tmpstart;

	$prev=$chr;
	$tmpstart+=$len;
	$tmpstart+=$gapsize;
#	print "adding gapsize $gapsize\n";
}


# Open outfile 
open(OUT, ">$out");


# Go through the infile, change the name and position columns
# and print the line with the new scaffold name and positions 
open(IN, $in);
while(<IN>) {
	chomp($_);
	my @tab = split(/\t+/, $_);

	my @col = split(/,/, $posCol);
	my $scaff = $tab[$scafCol];
	my $printFlag = "on";

	if(defined $scaffs{$scaff}) {

		my $chrom = $scaffs{$scaff}{'chr'};

		foreach my $c (@col) {
			my $c=$c-1;	#index is real col_no-1

			# CALCULATE THE NEW POSITION BASED ON THE OFFSET AND THE DIRECTION!
			my $newpos;
			if($scaffs{$scaff}{'dir'} eq "+") {
				$newpos = $scaffs{$scaff}{'offset'}+$tab[$c];		
			}
			elsif($scaffs{$scaff}{'dir'} eq "-") {
				$newpos = $scaffs{$scaff}{'offset'}+$scaffs{$scaff}{'len'}-$tab[$c]+1;		#BUGFIXED 20130912 (added +1)	
			}
			else {
				print "WARNING: $scaff has no assigned direction, assuming +!\n";
				$newpos = $scaffs{$scaff}{'offset'}+$tab[$c];		
			}	
			$tab[$c]=$newpos;

		}
		$tab[$scafCol]=$scaffs{$scaff}{'chr'};
		# Switch start and stop on negative scaffolds for gtf files
		if($format eq "gtf" && $scaffs{$scaff}{'dir'} eq "-") {
			#print "DEBUG: format is gtf and scaffold orientation is negative!!\n";
			my $tempstart=$tab[4];	
			$tab[4]=$tab[3];
			$tab[3]=$tempstart;
			if($tab[6] eq "+") {
				$tab[6]="-";
			}
			elsif($tab[6] eq "-") {
				$tab[6]="+";
			}
		}
		print OUT join("\t", @tab) . "\n";
				
	}
	else {
		if($unassigned eq "yes") {
			print OUT join("\t", @tab) . "\n";	
		}
		else {
			print $scaff." is not anchored to a chromosome, skipping line...\n";	
		}	
	}
}
close(IN);

#!/usr/bin/perl -w
use strict;
use threads;
use Thread::Queue;
use File::Basename;

my $size_of_each_piece = 50000000; #50 Mb per piece
my $cut = "/opt/LTR_FINDER.x86_64-1.0.5/cut.pl"; #the program to cut sequence
my $ltr_finder = "/opt/LTR_FINDER.x86_64-1.0.5/ltr_finder"; #the program to ltr_finder
my $seq_file = ""; #specify the sequence file
my $threads = 4; #specify thread number to run ltr_finder
my $finder_para = "-w 2 -C -D 15000 -d 1000 -L 7000 -l 100 -p 20 -M 0.85"; #specify parameters to run ltr_finder. The "-w 2" parameter is required.
my $timeout = 3600; #set maximum time for a thread to run. After $timeout seconds, the child thread is killed.
my $try1 = 0; #0, further split to 50 Kb regions if thread killed by timeout. 1, no further split.
my $next = 0;
my $harvest_format = 0; #indicate if the print out format is in LTR_finder (0, default) table output format or in LTRharvest (1) screen output format

my $usage = "
~ ~ ~ Run LTR_finder with multi-threading ~ ~ ~

Author: Shujun Ou (shujun.ou.1\@gmail.com)
Date: 09/19/2018

Usage: perl run_ltrfinder_in_pieces.pl -seq [file] -size [int] -threads [int]
Options:	-seq	[file]	Specify the sequence file
		-size	[int]	Specify the size you want to split the seq.
				Please make it large enough to avoid spliting too much of the seq. Default 50000000 (bp)
		-time	[int]	Specify the maximum time to run a subregion (a thread).
				This helps to skip simple repeat regions that take a substantial of time to run. Default: 3600 (seconds).
				Suggestion: 300 for -size 1000000. Increase -time when -size increased.
		-harvest_out	Output LTRharvest screen output format if specified. Default: output LTR_finder table format
		-next		Only summarize the results for previous jobs without rerunning LTR_finder.
		-threads	[int]	Indicate how many CPU/threads you want to run LTR_finder
		-cut	[file]	The path to the program cut.pl to split the sequence
		-finder	[file]	The path to the program LTR_finder
";

my $script_path = dirname(__FILE__);

my $i=0;
foreach (@ARGV){
	$seq_file = $ARGV[$i+1] if /^-seq$/i;
	$cut = $ARGV[$i+1] if /^-cut$/i;
	$ltr_finder = $ARGV[$i+1] if /^-finder$/i;
	$size_of_each_piece = $ARGV[$i+1] if /^-size$/i;
	$timeout = $ARGV[$i+1] if /^-time$/i;
	$try1 = 1 if /^-try1$/i;
	$threads = $ARGV[$i+1] if /^-threads|-t$/i;
	$next = 1 if /^-next$/i;
	$harvest_format = 1 if /^-harvest_out$/i;
	$i++;
	}

open SEQ, "<$seq_file" or die $usage;
my %seq; #a hash to store seq name and length info
my %order; #store seq order in genome
$i=0;
$/ = "\n>";
while (<SEQ>){
	s/>//g;
	my ($id, $seq) = (split /\n/, $_, 2);
	$seq =~ s/\s+//g;
	my $len = length $seq;
	$seq{$id} = $len;
	$order{$id} = $i;
	$i++;
	}
close SEQ;
$/="\n";

goto Next if $next == 1; #run the next step if all LTR_finder processes are done

`mkdir $seq_file.finder` unless -d "$seq_file.finder";
chdir "$seq_file.finder";
`perl $cut ../$seq_file -s -l $size_of_each_piece`;

##multi-threading using queue, create a worker module for parallel computation
my $process_q=Thread::Queue->new();
sub worker {
	while (my $seq = $process_q -> dequeue()){
		chomp $seq;
		$seq =~ s/\s+//;
		print localtime() ." CPU".threads -> self() -> tid().": running on $seq\n";
		`timeout $timeout $ltr_finder $finder_para $seq > $seq.finder.scn`; #set the max running time for a thread as $timeout seconds
		if ($? ne 0 and $try1 ne 1){
			my $in=`perl $script_path/run_ltrfinder_in_pieces.pl -seq $seq -size 50000 -time 10 -try1 -threads 1 -cut $cut -finder $ltr_finder`;
			`mv $seq.finder.combine.scn $seq.finder.scn`;
			`rm -rf ${seq}.finder`;
			}
		}
	}

#insert seq names into the worker queue
open List, "<../$seq_file.list" or die $!;
$process_q -> enqueue (<List>);
$process_q -> end(); #stop adding items to the queue
close List;

#work and finish
for (1..$threads){
	threads -> create ( \&worker );
	}
foreach my $thr (threads -> list()){
	$thr -> join();
	}

chdir "../";

Next:
#combine split ltr_finder results
open Out, ">$seq_file.finder.combine.scn" or die $!;

#print out headers
if ($harvest_format == 0){
	print Out "index	SeqID	Location	LTR len	Inserted element len	TSR	PBS	PPT	RT	IN (core)	IN (c-term)	RH	Strand	Score	Sharpness	Similarity\n";
	} else {
	print Out "#perl run_ltrfinder_in_pieces.pl -seq $seq_file -size $size_of_each_piece -time $timeout -harvest_out -threads $threads -cut $cut -finder $ltr_finder
# LTR_FINDER args=$finder_para
# predictions are reported in the following way
# s(ret) e(ret) l(ret) s(lLTR) e(lLTR) l(lLTR) s(rLTR) e(rLTR) l(rLTR) sim(LTRs) seq-nr chr
# where:
# s = starting position
# e = ending position
# l = length
# ret = LTR-retrotransposon
# lLTR = left LTR
# rLTR = right LTR
# sim = similarity
# seq-nr = sequence order\n";
	}

open List, "<$seq_file.list" or die $!;
foreach my $seq (<List>){
	$seq =~ s/\s+//;
	my ($base, $order) = ($1, $2) if $seq =~ /(.*)_sub([0-9]+)$/;
	my $coord_adj = ($order - 1) * $size_of_each_piece;
	next unless -e "$seq_file.finder/$seq.finder.scn";
	open Scn, "<$seq_file.finder/$seq.finder.scn" or die $!;
	while (<Scn>){
		next unless /^\[/;
		s/^\[\s+[0-9]+\]/\[NA\]/;
		my ($index, $id, $loc, $ltr_len, $ele_len, $TSR, $PBS, $PPT, $RT, $IN_core, $IN_cterm, $RH, $strand, $score, $sharpness, $sim) = (split);
		my @coord = ($loc, $PBS, $PPT, $RT, $IN_core, $IN_cterm, $RH);
		my $i = -1;
		foreach (@coord) {
			$i++;
			next if /^N-N$/;
			my ($start, $end) = ($1+$coord_adj, $2+$coord_adj) if /([0-9]+)\-([0-9]+)/;
			$coord[$i] = "$start-$end";
			}
		if ($harvest_format == 0){
			print Out "[NA]\t$base\t$coord[0]\t$ltr_len\t$ele_len\t$TSR\t$coord[1]\t$coord[2]\t$coord[3]\t$coord[4]\t$coord[5]\t$coord[6]\t$strand\t$score\t$sharpness\t$sim\n";
			} else {
			my ($start, $end) = ($1, $2) if $coord[0] =~ /([0-9]+)\-([0-9]+)/;
			my ($lltr, $rltr) = ($1, $2) if $ltr_len=~/([0-9]+),([0-9]+)/; 
			my ($lltr_e, $rltr_s) = ($start+$lltr-1, $end-$rltr+1);
			$sim*=100;
			print Out "$start $end $ele_len $start $lltr_e $lltr $rltr_s $end $rltr $sim $order{$base} $base\n";
			}
		}
	}
close List;
print localtime() ." Job finished! Check out $seq_file.finder.combine.scn\n";

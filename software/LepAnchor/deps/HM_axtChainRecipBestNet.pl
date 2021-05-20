#!/usr/bin/perl  

## Version 2.50

## Introduction.
## To create chain (all.chain.gz), reciprocal-best net/chain (all.rbest.Net/all.rbest.chain.gz)
## from raw blastz file (raw.axt.gz) files.
## By doing so, create a suitable data set for the further tba procedures. 
## Each time only one pair of species can be processed.

## Recommand log file: species1.species2.result/_axtChainRecipBestNet.log

## Assumptions.
## 1)The current/home directory is the working directory.
## 2)size files.
## 3)raw.axt.gz, all.chain.gz, all.rbest.net.gz
## 4)This script is in its very preliminary form, so users of this script take their own risk. 
## 5)the script does not check for argument syntax, the users need to ensure correct syntax. 

## Dependences.
## axtChain, chainAntiRepeat, chainMergeSort, chainStitchId, chainNet, NetSyntenic, chainSwap ...

## Update
## V1.00. 
## creation of  this script.
## V1.50.
## Add --zeroMinSpaceNet function.
## The code has been rewriten to suit HaploMerger.
## Program name is changed to HM_axtChainRecipBestNet.pl.
## V2.50.
## fix the scoring issue after rbest chainNeting.
## bug fixing #3: 
## fix the chainSort command in HM_axtChainRecipBestNet.pl to avoid writing /dev/stdout error in some systems.

use strict; 
use warnings;
use threads;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 1 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit; }

	This script is to perform chainNet and create reciprocal-best chain/net.

  Required excutables http://genome.ucsc.edu/admin/jksrc.zip,
including axtChain, chainAntiRepeat, chainMergeSort, chainStitchId, chainNet, 
NetSyntenic, chainSwap, netChainSubset, chainPreNet, chainSort, and netFilter.

   --help or no arguments,   show this message. Undefined inputs are ignored.

Main procedure: (!!!case sensitive)
   --axtChain - produce chain from raw blast alignemnt (axt files)
   --rbestNet - the step next to --axtChain, produce reciprol-best net

Options:   
   --Species -   mandatory, to provide species names which are wanted to be 
                 processed
                 NOTE THAT the species name should be given 
                 in ORDER, because the species name comes first will be used
                 as the target species.
                 Note: each run only allow two species.
                 Example: --Species human chimp
                 
   --linearGap=<medium|loose> - mandatory, set the -linearGap for axtChain
   --minScore=N - set the -minScore for axtChain (default=1000)
   
   --minScore2=N - set -minScore for chainNet (def=2000)
   --minSpace=N  - set -minSpace for chainNet (def=20)
   
   --axtSuffix=<name> - set suffix of raw axt files (def=axt, means *.axt),
                    
   --chainFile=<file_name> - specify the chain.gz file name for --rbestNet 
                 in order to create the 'all.rbest.net.gz' file
                 Do not include the path! The default is 'all.chain.gz'.
   --netFile=<net_name.net.gz> - specify the net file name for --rbestNet in order
                 to create 'net_name.rbest.net.gz' instead of 'all.rbest.net.gz'
                 Do not include the path!
                 This option requires that --chainFile have been specified.
                 A chain file named 'net_name.chain.gz' will be first
                 generated using the specified net and chain files, and 
                 then the 'net_name.rbest.net.gz' and 'net_name.rbest.chain.gz'
                 file will be generated subsequently.
   --qrbest -    to produce rbest net and chain files with the query
                 sequence as the reference (*.qrbest.net/chain.gz)
                 Note that this option can not be run alone!
   --tbest -     save tbest.net/chain.gz files for further use
   --zeroMinSpaceNet - create net with -minSpace=0 for further use
   
   --threads=N - number of CPUs used for axtChain (default=1)                                             
   --Force -     allow over-writing existing files (default=no)
   --Delete -    delete the temporary files (default=no)
   
   Required files:
   species1/2.sizes,
   species1/2.seq/*.nib,
   ~result/raw.axt/*.axt.*.
   
   Output files includes:
   ~result/all.chain.gz - results of axtChain,
   ~result/all.rbest.chain.gz - reciprocal-best chain,
   ~result/all.rbest.net.gz - rbest net,
   ~result/zeroMinSpace.rbest.net.gz - net with -minSpace=0.
   ~result/all.tbest.chain.gz - target-best chain,
   ~result/all.tbest.net.gz - target-best net, 
   ~result/all.qrbest.chain.gz - query-rbest chain,
   ~result/all.qrbest.net.gz - query-rbest net.    

Notes:
1)The current/home directory is the working directory.
2)The raw axt file name must be "*.axt.*" and stored in
  direcotry "species1.species2.result/raw.axt/".
3)--axtChain and --rbestNet can be run separately or together,
  but remember their running order.
4)The results are put into ./species1.species2.result.  
	 
USAGES

print "\n\n========== Start at "; system("date");
print "$0 $Arg_list\n\n";
my $timer = time();


#set the output_field_separator to " " and autoflushing
$,=' ';
$|=1;

# set filter score
my $filter_ali= 100;
if($Arg_list =~ m/--filterAli=(\d+)/){
	$filter_ali=$1;
}

# set $axt_suffix
my $axt_suffix="axt";
if($Arg_list =~ m/--axtSuffix=([\.\w]+)/){
	$axt_suffix=$1;
}

# set the over-write flag, to over-write any existing files
my $Force=0;
if ($Arg_list =~ m/--Force/) {
	$Force = 1;
	print "Set to OVER-WRITING mode!\n";
}

# Check and store species names
my @Species;
unless ($Arg_list =~ m/--Species\s+([^-]*)/) {die "No --Species argument or species_names found! die!\n" }
unless (@Species = $1 =~ m/(\w+)/g)  {die "No species names found! Die!\n" };
unless (scalar(@Species) == 2) { die "only two species are allowed! Die! \n"; } 
print "Species included: ", @Species, "\n";

print "Thread number for axtChain is set to ... ";
my $thread_count = 1;
if ($Arg_list =~ m/--threads=(\d+)/){
	$thread_count = $1 > 0 ? $1 : 1;
}
print "$thread_count !\n";

sub axtChain;
sub rbestNet;
sub netFilter;

axtChain if ($Arg_list =~ m/--axtChain/);
rbestNet if ($Arg_list =~ m/--rbestNet/);
netFilter if ($Arg_list =~ m/--netFilter/);

print "\n\n========== Time used = ", time()-$timer, " seconds or ", (time()-$timer)/3600, " hours.\n\n";

######################### subroutine defination #####################################

sub axtChain_workhorse($$$$$){ #target, query, axtfile, thread_id

	my ($target, $query, $axt_file, $option_string, $thread_id) = @_;
	
	my $cmd;
	$cmd  = " axtChain $option_string $axt_file $target.seq $query.seq /dev/stdout ";
	$cmd .= " | chainAntiRepeat $target.seq $query.seq /dev/stdin $axt_file.chain ";
	print "Thread_ID $thread_id: $cmd \n";
	system($cmd);
	
	return($thread_id);
		
}

sub axtChain {
	
	print "=== axtChain start ===\n";
	
	# checking sizes files
	print "Checking missing sizes files ...\n";
	my $missing_sizes =0;
	unless (-f "$Species[0].sizes") { $missing_sizes =1; print "$Species[0].sizes\n"; }
	unless (-f "$Species[1].sizes") { $missing_sizes =1; print "$Species[1].sizes\n"; }
	die "File missing! Die! \n" if ($missing_sizes == 1);
	
	# checking nib files
	my $missing_nib = 0;
	my $sizeFH;
	foreach my $temp (@Species) {
		print "checking nib files in $temp.seq/ ...\n";
		open($sizeFH, "<$temp.sizes") or die "Cant not open directory $temp.sizes! die!\n";
		while (<$sizeFH>) {
			next unless m/([-\.\w]+)\t/;
			unless (-f "$temp.seq/$1.nib") { $missing_nib=1; print "$1.nib\n"; }
		}
		close $sizeFH;
	}
	die "Some nib files are missing! Die!\n" if ($missing_nib == 1);	
	
	# checking axt files
	#my $missing_axt = 0;
	#print "checking axt files in $Species[0].$Species[1].result/raw.axt/ ...\n";
	#open($sizeFH, "<$Species[0].sizes") or die "Cant not open size file $Species[0].sizes! die!\n";
	#while (<$sizeFH>) {
	#	next unless m/([-\.\w]+)\t/;
	#	unless (-f "$Species[0].$Species[1].result/raw.axt/$1.axt") { $missing_axt=1; print "$1.axt is not found !\n"; }
	#}
	#close $sizeFH;
	#print "These axt files are missing! Going on anyway!\n" if ($missing_axt == 1);	

	# checking old all.chain.gz files
	print "checking chain file $Species[0].$Species[1].result/all.chain.gz ...\n";
	if (-f "$Species[0].$Species[1].result/all.chain.gz" ) {
		$Force ? print "File existed, will be over-written!\n" : die "File existed! die!\n";
	}
	
	# creating all.chain.gz files
	my $path = "$Species[0].$Species[1].result";
	my $option_string;
	my $min_score=1000;
	die "--linearGap is mandatory! die!\n" unless (($option_string) = $Arg_list =~ m/-(-linearGap=medium|-linearGap=loose)/);
	if ($Arg_list =~ m/--minScore=(\d+)/) { $min_score=$1; $option_string .= ' '.'-minScore='.$min_score; }
	$option_string .= " -verbose=0 ";
	if(-f "$Species[0].$Species[1].q"){
		$option_string .= " -scoreScheme=$Species[0].$Species[1].q ";
	}elsif(-f "scoreMatrix.q"){
		$option_string .= " -scoreScheme=scoreMatrix.q ";
	}
	
	my @axt_files = glob "$path/raw.axt/*.$axt_suffix";
		
	my $thread_id = 1; # 0=no available threads, >0 refer to id of an available thread
	my @thread_stat = (0) x ($thread_count+1); #0 - free, >=1 = occupied or in used
	my @workhorses;
	my $rv; 	

	foreach my $temp (@axt_files){
				
		$workhorses[$thread_id]	= threads->create('axtChain_workhorse', $Species[0], $Species[1], $temp, $option_string, $thread_id);
		$thread_stat[$thread_id]=1;	
		$thread_id=0;		
		
		while($thread_id<1){
			for (my $i=1;$i<=$thread_count;$i++){
				if ($thread_stat[$i]==0){
					$thread_id=$i;
					last;
				}
			}
			if ($thread_id<1){
				for (my $i=1;$i<=$thread_count;$i++){
					if ($workhorses[$i]->is_joinable()){
						$rv = $workhorses[$i]->join();
						print "WARNING: an error occurred to thread ID $i !!!!!!!\n" if $rv ne $i;
						$thread_stat[$i]=0;							
					}
				}
				sleep(1);
			}
		}
		
	}
		
	for (my $i=1;$i<=$thread_count;$i++){
		if ($thread_stat[$i]>0){
			$rv = $workhorses[$i]->join();
			print "WARNING: an error occurred to thread ID $i !!!!!!!\n" if $rv ne $i;							
		}
	}
	
	my $cmd = "chainMergeSort $path/raw.axt/*.$axt_suffix.chain | gzip -c >$path/all.chain.gz";
	print $cmd,"\n";
	system($cmd);
	
	print "=== axtChain done! ===\n\n";
	if ($Arg_list =~ m/--Delete/){
		unlink glob "$path/raw.axt/*.$axt_suffix.chain";
	}
}

####

sub rbestNet {
	
	print "=== rbestNet start ===\n";
	
	my ($min_score,$min_space)=(2000,20);
	if($Arg_list =~ m/--minScore2=(\d+)/){
		$min_score=$1;
		print "--minScore2 is set to $min_score.\n";
	}
	if($Arg_list =~ m/--minSpace=(\d+)/){
		$min_space=$1;
		print "--minSpace is set to $min_space.\n";
	}
		
	# checking sizes files
	print "Checking missing sizes files ...\n";
	my $missing_sizes = 0;
	unless (-f "$Species[0].sizes") { $missing_sizes =1; print "$Species[0].sizes\n"; }
	unless (-f "$Species[1].sizes") { $missing_sizes =1; print "$Species[1].sizes\n"; }
	die "File missing! Die! \n" if ($missing_sizes == 1);
	
	# Store all.chain.gz name and check its existence
	print "checking chain.gz files ...\n";
	my $in_chain_file;
	if ($Arg_list =~ m/--chainFile=([-\.\w]+\.chain\.gz)/) {
		$in_chain_file =  "$Species[0].$Species[1].result/$1";
	}else{
		$in_chain_file =  "$Species[0].$Species[1].result/all.chain.gz";
	} 
	print $in_chain_file, "\n";
	unless (-f $in_chain_file) { die "File not found! Die!\n"; }
	
	# Store net file name and check its existence
	print "checking the net file and if present, produce the corresponding chain file ...\n";
	my $net_name = 'all';
	if ($Arg_list =~ m/--netFile=([-\.\w]+)\.net/) { $net_name = $1; }
	if ($net_name ne 'all') {
		if ($net_name =~ m/rbest/) { die "The input netfile name should not contain the string \'rbest\'! Die!\n"; }
		my $net_file =  "$Species[0].$Species[1].result/$net_name.net";
		print $net_file, "\n";
		unless (-f $net_file) { die "File not found! Die!\n"; }
		
		if ( -f "$Species[0].$Species[1].result/$net_name.chain.gz"){
			$Force ? print "$net_name.chain.gz is already existed, will be over-written!\n" : die "$net_name.chain.gz is already existed! die!\n";
		}
		
		my $cmd;
		$cmd  = "gunzip -c $in_chain_file ";
		$cmd .= "| netChainSubset $net_file /dev/stdin /dev/stdout 2>>$Species[0].$Species[1].result/_netChainSubset.log ";
		$cmd .= "| chainStitchId /dev/stdin /dev/stdout 2>>$Species[0].$Species[1].result/_ChainStitchId.log ";
		$cmd .= "| gzip -c >$Species[0].$Species[1].result/$net_name.chain.gz";
		print $cmd,"\n";
		system($cmd);
		$in_chain_file = "$Species[0].$Species[1].result/$net_name.chain.gz";
	}					
		
	
	# checking old rbest.net and rbest.chain.gz files
	print "checking old rbest.net/chain.gz and qrbest.net/chain.gz file ...\n";
	if (-f "$Species[0].$Species[1].result/$net_name.rbest.net.gz" 
				or -f "$Species[0].$Species[1].result/$net_name.rbest.chain.gz"
				or -f "$Species[0].$Species[1].result/$net_name.qrbest.net.gz"
				or -f "$Species[0].$Species[1].result/$net_name.qrbest.chain.gz"
				or -f "$Species[0].$Species[1].result/$net_name.tbest.net.gz" ) {
		$Force ? print "File existed, will be over-written!\n" : die "File existed! Die!\n";
	}

	####################### create the all.rbes.net file  ###########################
	my ($target,$query) = @Species[0,1];
	my $path = "$target.$query.result";
	my $cmd;
	
	print "creating $target.$query.$net_name.tbest.net and chain ...\n";
	$cmd  = "gunzip -c $in_chain_file ";
	$cmd .= "| chainPreNet /dev/stdin $target.sizes $query.sizes /dev/stdout 2>>$path/_ChainPreNet.log ";
	$cmd .= "| chainNet -minSpace=$min_space -minScore=$min_score /dev/stdin $target.sizes $query.sizes /dev/stdout /dev/null 2>>$path/_ChainNet.log ";
	$cmd .= "| netSyntenic /dev/stdin $path/$target.$query.$net_name.tbest.net 2>>$path/_netSyntenic.log";
	print $cmd,"\n";
	system($cmd);
	$cmd  = "gunzip -c $in_chain_file ";
	$cmd .= "| netChainSubset $path/$target.$query.$net_name.tbest.net /dev/stdin /dev/stdout 2>>$path/_netChainSubset.log ";
	$cmd .= "| chainStitchId /dev/stdin $path/$target.$query.$net_name.tbest.chain 2>>$path/_ChainStitchId.log";
	print $cmd,"\n";
	system($cmd);
	
	print "creating all.rbest.net.gz ... \n";
	# Swap target-referenced tbest.chain to query-referenced tbest.chain
	$cmd  = "chainSwap $path/$target.$query.$net_name.tbest.chain /dev/stdout 2>>$path/_ChainSwap.log";
	$cmd .= "| chainSort /dev/stdin $path/$query.$target.$net_name.tbest.chain 2>>$path/_ChainSort.log"; #bug fixing #3
	print $cmd,"\n";
	system($cmd);
	# Net the query-referenced tbest.chain to get the reciprocal-best net
	$cmd  = "chainPreNet $path/$query.$target.$net_name.tbest.chain $query.sizes $target.sizes /dev/stdout 2>>$path/_ChainPreNet.log ";
	$cmd .= "| chainNet -minSpace=$min_space -minScore=$min_score /dev/stdin $query.sizes $target.sizes /dev/stdout /dev/null 2>>$path/_ChainNet.log ";
	$cmd .= "| netSyntenic /dev/stdin $path/$query.$target.$net_name.rbest.net 2>>$path/_netSyntenic.log";
	print $cmd,"\n";
	system($cmd);
	# Extract query-referenced reciprocal-best chain
	$cmd  = "netChainSubset $path/$query.$target.$net_name.rbest.net $path/$query.$target.$net_name.tbest.chain /dev/stdout 2>>$path/_netChainSubset.log ";
	$cmd .= "| chainStitchId /dev/stdin $path/$query.$target.$net_name.rbest.chain 2>>$path/_ChainStitchId.log";
	print $cmd,"\n";
	system($cmd);
	# Swap query-refenced rbest.chain to target-referenced rbest.chain
	$cmd  = "chainSwap $path/$query.$target.$net_name.rbest.chain /dev/stdout 2>>$path/_ChainSwap.log "; 
	$cmd .= "| chainSort /dev/stdin $path/$target.$query.$net_name.rbest.chain 2>>$path/_ChainSort.log"; #bug fixing #3
	print $cmd,"\n";
	system($cmd);
	# Net the target-referenced rbest.chain to get the target-referenced reciprocal-best (or target-query-best) net
	$cmd  = "chainPreNet $path/$target.$query.$net_name.rbest.chain $target.sizes $query.sizes /dev/stdout 2>>$path/_ChainPreNet.log ";
	$cmd .= "| chainNet -minSpace=$min_space -minScore=$min_score /dev/stdin $target.sizes $query.sizes /dev/stdout /dev/null 2>>$path/_ChainNet.log ";
	$cmd .= "| netSyntenic /dev/stdin $path/$target.$query.$net_name.rbest.net 2>>$path/_netSyntenic.log ";
	print $cmd,"\n";
	system($cmd);
	
	# net to axt
	$cmd  = "netToAxt -maxGap=99999999 $path/$target.$query.$net_name.rbest.net $path/$target.$query.$net_name.rbest.chain ";
	$cmd .= "$target.seq $query.seq $path/$target.$query.$net_name.rbest.axt"; ###########fixing??
	print $cmd,"\n";
	system($cmd);
	# axt to chain
	my $options;
	die "--linearGap is mandatory! die!\n" unless (($options) = $Arg_list =~ m/-(-linearGap=medium|-linearGap=loose)/);
	$options .= ' -minScore=0 ';
	$options .= " -verbose=0 ";
	if(-f "$Species[0].$Species[1].q"){
		$options .= " -scoreScheme=$Species[0].$Species[1].q ";
	}elsif(-f "scoreMatrix.q"){
		$options .= " -scoreScheme=scoreMatrix.q ";
	}	
	$cmd  = "axtChain $options $path/$target.$query.$net_name.rbest.axt ";
	$cmd .= "$target.seq $query.seq  /dev/stdout | gzip -c > $path/$net_name.rbest.chain.gz";
	print $cmd,"\n";
	system($cmd);	
	# chainNet
	$cmd  = "gunzip -c $path/$net_name.rbest.chain.gz ";
	$cmd .= "| chainPreNet /dev/stdin $target.sizes $query.sizes /dev/stdout 2>>$path/_ChainPreNet.log ";
	$cmd .= "| chainNet -minSpace=$min_space -minScore=$min_score /dev/stdin $target.sizes $query.sizes /dev/stdout /dev/null 2>>$path/_ChainNet.log ";
	$cmd .= "| netSyntenic /dev/stdin $path/$net_name.rbest.net 2>>$path/_netSyntenic.log";
	print $cmd,"\n";
	system($cmd);
	$cmd  = "gzip -f $path/$net_name.rbest.net";
	system($cmd);
	print "reciprocal-best $target-reference chain file $path/$net_name.rbest.chain.gz created!\n";
	print "reciprocal-best $target-reference net file $path/$net_name.rbest.net created!\n";	

	if ($Arg_list =~ m/--qrbest/){

		$cmd  = "gunzip -c $path/$net_name.rbest.chain.gz ";
		$cmd .= "| chainSwap /dev/stdin /dev/stdout 2>>$path/_ChainSwap.log "; 
		$cmd .= "| chainSort /dev/stdin  $path/$query.$target.$net_name.qrbest.chain 2>>$path/_ChainSort.log "; #bug fixing #3
		print $cmd,"\n";
		system($cmd);

		$cmd  = "chainPreNet $path/$query.$target.$net_name.qrbest.chain $query.sizes $target.sizes /dev/stdout 2>>$path/_ChainPreNet.log ";
		$cmd .= "| chainNet -minSpace=$min_space -minScore=$min_score /dev/stdin $query.sizes $target.sizes /dev/stdout /dev/null 2>>$path/_ChainNet.log ";
		$cmd .= "| netSyntenic /dev/stdin $path/$net_name.qrbest.net 2>>$path/_netSyntenic.log";
		print $cmd,"\n";
		system($cmd);		
		
		$cmd  = "netChainSubset $path/$net_name.qrbest.net $path/$query.$target.$net_name.qrbest.chain /dev/stdout 2>>$path/_netChainSubset.log ";
		$cmd .= "| chainStitchId /dev/stdin /dev/stdout 2>>$path/_ChainStitchId.log ";
		$cmd .= "| gzip -c > $path/$net_name.qrbest.chain.gz ";
		print $cmd,"\n";
		system($cmd);
		$cmd  = "gzip -f $path/$net_name.qrbest.net ";
		system($cmd);					
		print "reciprocal-best $query-reference chain file $path/$net_name.qrbest.chain.gz created!\n";
		print "reciprocal-best $query-reference net file $path/$net_name.qrbest.net created!\n";
	}

	if ($Arg_list =~ m/--tbest/) {
		print "Save $net_name.tbest.net and $net_name.tbest.chain.gz ...\n";
		system("gzip -f -c $path/$target.$query.$net_name.tbest.net >$path/$net_name.tbest.net.gz");
		system("gzip -f -c $path/$target.$query.$net_name.tbest.chain > $path/$net_name.tbest.chain.gz");
	}
	
	if ($Arg_list =~ m/--zeroMinSpaceNet/) {
		print "create zeroMinSpace.rbest.net.gz file ...\n";
		$cmd  = "gunzip -c $path/$net_name.rbest.chain.gz ";
		$cmd .= "| chainPreNet /dev/stdin $target.sizes $query.sizes /dev/stdout 2>>$path/_ChainPreNet.log ";
		$cmd .= "| chainNet -minSpace=0 -minScore=$min_score /dev/stdin $target.sizes $query.sizes /dev/stdout /dev/null 2>>$path/_ChainNet.log ";
		$cmd .= "| netSyntenic /dev/stdin /dev/stdout 2>>$path/_netSyntenic.log | gzip -f -c >$path/zeroMinSpace.rbest.net.gz";
		print $cmd,"\n";
		system($cmd);
		print "Finished creating zeroMinSpace.rbest.net.gz file.\n";
	}	
	
	if ($Arg_list =~ m/--Delete/) {
		print "Deleting imtermediate net and chain files ... \n";
		unlink "$path/$target.$query.$net_name.tbest.net";
		unlink "$path/$target.$query.$net_name.tbest.chain";
		unlink "$path/$query.$target.$net_name.tbest.chain";
		unlink "$path/$query.$target.$net_name.rbest.net";
		unlink "$path/$query.$target.$net_name.rbest.chain";
		unlink "$path/$target.$query.$net_name.rbest.chain";
		unlink "$path/$target.$query.$net_name.rbest.net";
		unlink "$path/$target.$query.$net_name.rbest.axt";
		unlink "$path/$query.$target.$net_name.qrbest.chain";
	}

	print "=== rbestNet done! ===\n\n";
}

sub netFilter{
	
	#########################################
	#### create filtered net (and chain.gz files)
	print "create filtered.rbest.net.gz file ...\n";
	unless (-f "$Species[0].$Species[1].result/all.rbest.net.gz" and -f "$Species[0].$Species[1].result/all.rbest.chain.gz")  {die "$Species[0].$Species[1].result/all.rbest.net(chain.gz) files not found!\n"; }
	#
	# set filter score
	my $out_netfile = "filtered.rbest.net.gz";
	if($Arg_list =~ m/--outputNet=([-\.\w]+)/){
		$out_netfile=$1;
	}
	#
	print "Set minScore to $filter_ali for net filtering ...\n";
	#
	my $cmd = "gunzip -c $Species[0].$Species[1].result/all.rbest.net.gz | netFilter -minAli=$filter_ali /dev/stdin | gzip -c >$Species[0].$Species[1].result/$out_netfile";
	print $cmd,"...\n";
	system($cmd);
}		

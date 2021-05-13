#!/usr/bin/perl   

## Version 1.60

## Introduction.
## This script runs lasz between each species pair and produce raw.axt files.
## By doing so, create a suitable data set for the further axtChainNet and tba procedures. 

## Recommand log file for HM_all_lastz.pl             : ./_all_lastz.log
## Recommand log file for normal lastz (buildin)   : ./species1.species2.result/_lastz.log

## Assumptions.
## 1)The current/home directory is the working directory.
## 2)nib file directory (species.seq/*.*) and ./species.fa.gz is present (required by lastz all to all).
## 2)Executable programs required for lastz can be found through the $PATH.
## 3)All pertinent files, such as nib files, ctl file are in correct format, in proper place, naming and other condition.
## 5)HM_all_lastz.pl requires control files.
##   'species1.species2.all_lastz.ctl' is the control file for all_lastz.pl for each comparison pair.
##   if such file is not found. all_lastz.ctl is the default control file for all_lastz.pl.
## 6)HM_all_lastz.pl may require scoreMatrix files, but not a must.
##   'species1.species2.q' is the score file for each comparison pair.
## 4)This script is in its very preliminary form, so users of this script take their own risk. 
## 5)the script does not check for argument syntax, the users need to ensure correct syntax. 

## Dependences.
## lastz, gzip
## faFilter -v -name=sequence_name input.fa output.fa ; from ucsc jksrc.zip

## Update
## V1.00. 
## creation of  this script.
## V1.10. 
## inferOnly, set --step=20 and --notransition.
## add an option: --identity=<range..range>.
## add an option: --noself to supress self comparison in real all_lastz processes, which means
## chromosome with the same name in target and query is skipped. 
## V1.20
## all_lastz.pl was rewritten to HM_all_lastz.pl (HM for HaploMerger) to handle and SPEED UP
## thousands of target sequences versus thousands of query sequences.
## (all_lastz.pl runs much slower than this script).
## V1.50
## The code is rewritten to acommadate the use of multiple CPU.
## The code is rewritten to suit the use of HaploMerger.
## The name also has been changed to HM_all_lastz_mThreads.pl. 
## V1.60
## add the option to close soft masking

use strict; 
use warnings;
use threads;

my $Arg_list = join " ", @ARGV;
if (@ARGV < 1 || $Arg_list =~ m/--help|-\?/) { print <<'USAGES' and exit }

	 This script is to run lastz on all species pairs, the result is put into
result directory (specie1.species2.result/raw.axt/*.axt).

   This new script can perform target*query parallel jobs.

   Require the excutable lastz from http://www.bx.psu.edu/~rsharris/lastz/.
   Requrie excutables from  http://genome.ucsc.edu/admin/jksrc.zip,
including faSplit.

Options: (!!!case sensitive)
   --help or no arguments - show this message
   --Species -   mandatory, provide species names which are wanted to be
                 processed
                 NOTE THAT the species name should be given
                 in ORDER, because the species name comes first will be used
                 as the target species.
                 Example: --Species human chimp rhesus
                 
   --notrivial - produce axt files (*.axt.self.notrivial) with
                 trivial self-self alignments removed (def=no, radius=5000 bp)
   --radius=N - used with --notrivial, to remove trivial self-alignment
                along the diagonal (def=5000 bp)              
   
   --noself - produce axt files (*.axt.noself) containing no self-self 
              alignments (default=no)
   --deleteSelfAlignment - delete those axt files (*.axt.self.notrivial, 
                 *.axt.selfOnly.notrivial) containing 
                 self alignments (def=no, works only with --noself) 
   
   --identity=<min>[..<max>] - specify identity range (default=no range)                             

   --targetSize=N - split target fasta file by size N (default=50000000)
   --querySize=N - split query fasta file by size N (default=1600000000)

   --threads=N - specified N threads/CPUs for use (default=1)   
   --Force -     over-write existing files (default=no)
   --Delete -    delete temp *.fa files (default=no)
   
   --unmask - for coverting lowercases to uppercases; def=not
   
   Required files:
   species.fa.gz,
   species.sizes,
   species.seq/*.nib,
   scoreMatrix.q or species1.species2.q,
   all_lastz.ctl or species1.species2.all_lastz.ctl.
   
   Results are put into ~result/raw.axt/*.axt.*:
   *.axt - original axt,
   *.axt.self.notrivial - axt without trivial self-alignment,
   *.axt.noself - axt without self-alignment,
   *.axt.selfOnly.notrivial - axt with only nontrivial self-alignment.   

Notes:
1)This script is suitable to deal with thousands of target sequences
	against thousands of query sequences (it is done partially by collecting all query 
	sequences in one tba-fasta file[species.fa.gz]).
2)The current directory is the working directory.
3)nib file directory (species.seq), and ./species.fa.gz are required by lastz.
4)Executable programs required for lastz can be found through the $PATH.
5)All pertinent files, such as nib files, ctl files, size files, .q files,
  are in correct format, proper place, naming and other condition.
  This evironment can be created by using script named "initiation.pl".
6)HM_all_lastz.pl requires control files.
  'species1.species2.all_lastz.ctl' is the control file for HM_all_lastz.pl for 
  each comparison pair. if such file is not found, all_lastz.ctl is the default 
  control file for HM_all_lastz.pl.
7)HM_all_lastz.pl may require scoreMatrix files, but not a must.
  'species1.species2.q' is the score file for each comparison pair.
8)The results are put into "./species1.species2.result/raw.axt/".  
	 
USAGES

print "\n\n========== Start at "; system("date");
print "$0 $Arg_list\n\n";
my $timer = time();

#set the output_field_separator to " "
$,=' ';

# set the over-write flag, to over-write any existing files
my $Force=0;
if ($Arg_list =~ m/--Force/) {
	$Force = 1;
	print "Set to OVER-WRITING mode!\n";
}

# set the deletion flag
my $Delete=0;
if ($Arg_list =~ m/--Delete/) {
	$Delete = 1;
	print "Set to Delete mode!\n";
}

# Store species names and to check if all prerequesite data are in proper conditions.
my @Species;
unless ($Arg_list =~ m/--Species\s+([^-]*)/) {die "No --Species argument or species_names found!\n"; }
unless (@Species = $1 =~ m/(\w+)/g)  {die "No species names found!\n"; }
print "Species included: ", @Species, "\n";

print "checking missing gz-compressed multiple-fasta files ...\n";
my $missing_gzFa = 0;
foreach my $temp (@Species) {
	unless (-f "$temp.fa.gz") { $missing_gzFa = 1; print "$temp.fa.gz\n"; }
}

if ( $missing_gzFa > 0 ) {
	die "Some files/directories are missing! Die!\n"
}else{
	print "Every required files seem ok, going on ...\n\n";
}

print "Check existed raw.axt directory ...\n";
my $has_axt = 0;
for (my $i=0;$i < scalar(@Species)-1;$i++) {
	for (my $j=$i+1;$j < scalar(@Species);$j++) {
		my $temp = "$Species[$i].$Species[$j].result/raw.axt";
		print "checking if $temp is existing ...\n";
		if (-d $temp)   { $has_axt=1; print "$temp is already existed !\n"; }
	}
}
die   "existing raw.axt directory found! Die!\n" if ($has_axt>0 and $Force == 0);
print "existing raw.axt directory found! \nForce to " if ( $has_axt>0 and $Force == 1);

print "cleaning the evironment anyway ...\n";
for (my $i=0;$i < scalar(@Species)-1;$i++) {
	for (my $j=$i+1;$j < scalar(@Species);$j++) {
		my $temp = "$Species[$i].$Species[$j].result";
		unless (-d $temp) {	mkdir $temp or die "Can not create directory $temp !\n"; }
		if (-d "$temp/raw.axt") { 
			unlink glob "$temp/raw.axt/*.axt";
			unlink glob "$temp/raw.axt/*.chain";
			unlink glob "$temp/raw.axt/*.log";  
		}else{
			mkdir "$temp/raw.axt" or die "Can not create directory $temp/raw.axt !\n"; 
		}
	}
}

print "Thread number is set to ... ";
my $thread_count = 1;
if ($Arg_list =~ m/--threads=(\d+)/){
	$thread_count = $1 > 0 ? $1 : 1;
}
print "$thread_count !\n";

print "Target file size is ... ";
my $target_size =  50000000;
if($Arg_list =~ m/--targetSize=(\d+)/){$target_size = $1; }
print "$target_size !\n";

print "Query file size is ... ";
my $query_size = 1600000000;
if($Arg_list =~ m/--querySize=(\d+)/){ $query_size = $1; }
print "$query_size !\n";

print "Is_noself status is ... ";
my $is_noself = 0;
$is_noself = 1 if $Arg_list =~ m/--noself/;
print "$is_noself !\n";

print "--notrivial status is ... ";
my $notrivial=0;
if($Arg_list =~ m/--notrivial/){	$notrivial = 1;	}
print "$notrivial !\n";

print "--radius for notrivial function is ... ";
my $radius=5000;
if($Arg_list =~ m/--radius=(\d+)/){	$radius = $1;	}
print "$radius !\n";

print "--unmask for covert lowcase letters to upcase letters ... ";
my $unmask='';
if($Arg_list =~ m/--unmask/){	$unmask = "unmask";	}
print "$unmask !\n";

print "go to lastz all to all step ...\n\n";

########################## main body of this all_lastz #################################

#the subroutine to read the all_lastz.ctl and scoreMatrix.q file
sub read_ctl_q($$);

#the workhorse for lastz (multithreads)
sub lastz_workhorse($$$$);

############################################################
##doing lastz all to all

my @target_fas;
my @query_fas;
my ($target, $query);

for (my $i=0;$i < scalar(@Species)-1;$i++) {
	for (my $j=$i+1;$j < scalar(@Species);$j++) {
		
		($target, $query) = @Species[$i, $j];
		
		#### split target fasta files into HM_all_lastz_target_N
		#### split query fasta files into HM_all_lastz_query_N
		my $rr=int(rand(99999))+1;
		system("faSplit about $target.fa.gz $target_size HM_all_lastz_target_$rr"."_");
		system("faSplit about $query.fa.gz  $query_size  HM_all_lastz_query_$rr"."_");
		@target_fas = glob("HM_all_lastz_target_$rr"."_*.fa");
		@query_fas  = glob("HM_all_lastz_query_$rr"."_*.fa"); 
		print "Finished splitting the target fasta file ...\n";
		print "$_\n" foreach (@target_fas); 
		print "Finished splitting the query fasta file ...\n";
		print "$_\n" foreach (@query_fas);
		
		####
		####		
		my $option_string = read_ctl_q($target, $query);
		
		my $thread_id = 0; # 0=no available threads, >0 refer to id of an available thread
		my @thread_stat = (0) x ($thread_count+1); #0 - free, >=1 = occupied or in used
		my @workhorses;
		my $rv; 	
		
		print "lastz $target to $query ... \n";	
		foreach my $temp (@target_fas) {
		foreach my $temp1 (@query_fas) {	
			while ($thread_id<1){ #get an availabe thread
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
					sleep(5);
				}
			}
			
			$workhorses[$thread_id]	= threads->create('lastz_workhorse', $temp, $temp1, $thread_id, $option_string);
			$thread_stat[$thread_id]=1;	
			$thread_id=0;						
		}
		}
		
		print "\n";
		
		#join all threads
		for (my $i=1;$i<=$thread_count;$i++){
			if ($thread_stat[$i]>0){
				$rv = $workhorses[$i]->join();
				print "WARNING: an error occurred to thread ID $i !!!!!!!\n" if $rv ne $i;							
			}
		}
		
		print "Finished lastz all $target to all $query !\n";
		
		#merge axt file by target_fa name
		foreach my $temp (@target_fas) {
		foreach my $temp1 (@query_fas) {
			system("cat $target.$query.result/raw.axt/$temp.$temp1.axt >>$target.$query.result/raw.axt/$temp.axt");
			unlink "$target.$query.result/raw.axt/$temp.$temp1.axt";				
		}
		}		
		
		#### notrivial filter
		if($notrivial > 0 or $is_noself >0){
			foreach my $tt (@target_fas){
				print "To delete self-self trivial alignment from $target.$query.result/raw.axt/$tt.axt ...\n";
				open(my $inFH,"<$target.$query.result/raw.axt/$tt.axt") or die "Can't read $target.$query.result/raw.axt/$tt.axt.\n";
				open(my $outFH,">$target.$query.result/raw.axt/$tt.axt.self.notrivial") or die "Can't open $target.$query.result/raw.axt/$tt.axt.self.notrivial.\n";
				while(<$inFH>){
					if(m/^#|^\s/){ print $outFH $_; next; }
					m/^\d+\s+([-\.\w]+)\s+(\d+)\s+(\d+)\s+([-\.\w]+)\s+(\d+)\s+(\d+)\s+([+-])/;
					if($1 eq $4 and (($2>$5-$radius and $2<$5+$radius) or ($3>$6-$radius and $3<$6+$radius))){
						$_=<$inFH>; $_=<$inFH>; 
					}else{
						print $outFH $_;
						$_=<$inFH>; print $outFH $_; $_=<$inFH>; print $outFH $_; 
					} 
				}
				close $inFH; close $outFH;
			}
		}				
		
		if($is_noself > 0){
			my $delete_self = $Arg_list =~ m/--deleteSelfAlignment/ ? 1 : 0;
			foreach my $tt (@target_fas){
				print "To delete self-self alignment from $target.$query.result/raw.axt/$tt.axt.self.notrivial ...\n";
				#rename("$target.$query.result/raw.axt/$tt.axt","$target.$query.result/raw.axt/$tt.self.axt") or die "Can't rename axt files!\n";
				open(my $inFH,"<$target.$query.result/raw.axt/$tt.axt.self.notrivial") or die "Can't read $target.$query.result/raw.axt/$tt.axt.self.notrivial.\n";
				open(my $outFH,">$target.$query.result/raw.axt/$tt.axt.noself") or die "Can't open $target.$query.result/raw.axt/$tt.axt.noself.\n";
				open(my $s_outFH,">$target.$query.result/raw.axt/$tt.axt.selfOnly.notrivial") or die "Can't open $target.$query.result/raw.axt/$tt.axt.selfOnly.notrivial\n";
				while(<$inFH>){
					if(m/^#|^\s/){ print $outFH $_; print $s_outFH $_; next; }
					m/^\d+\s+([-\.\w]+)\s+\d+\s+\d+\s+([-\.\w]+)\s+/;
					if($1 eq $2){
						print $s_outFH $_;
						$_=<$inFH>; print $s_outFH $_; $_=<$inFH>; print $s_outFH $_;
					}else{
						print $outFH $_;
						$_=<$inFH>; print $outFH $_; $_=<$inFH>; print $outFH $_; 
					} 
				}
				close $inFH; close $outFH; close $s_outFH;
				unlink "$target.$query.result/raw.axt/$tt.axt.self.notrivial" if $delete_self==1;
				unlink "$target.$query.result/raw.axt/$tt.axt.selfOnly.notrivial" if $delete_self==1;
			}
		}	
		unlink @target_fas if $Delete==1;
		unlink @query_fas if $Delete==1;
	}
}

print "\n\n========== Time used = ", time()-$timer, " seconds or ", (time()-$timer)/3600, " hours.\n\n";


########################## subroutines of this all_lastz #################################

#the subroutine to read the all_lastz.ctl and scoreMatrix.q file
sub read_ctl_q($$) {
	my $target = shift;
	my $query  = shift;
	my $ctlFH;
	my $option_string=' ';
	if (-f "$target.$query.all_lastz.ctl") {
		open($ctlFH, "<$target.$query.all_lastz.ctl") or die "Can not open $target.$query.all_lastz.ctl!\n";
	}elsif (-f "all_lastz.ctl") {
		open($ctlFH, "<all_lastz.ctl") or die "Can not open all_lastz.ctl!\n";
	}else{
		die "No control file (*.all_lastz.ctl) found! Die!\n";
	}
	while (<$ctlFH>) {
		s/#.*//;
		if (m/(--[^\s]+)/) { $option_string .= $1.' '; }
	}
	if (-f "$target.$query.q") { 
		$option_string .= "--scores=$target.$query.q "; 
	}elsif( -f "scoreMatrix.q"){
		$option_string .= "--scores=scoreMatrix.q ";
	}
	if ($Arg_list =~ m/(--identity=)(\d+\.\.\d+|\d+)/){ $option_string .= $1 . $2 ." "; }
	return $option_string;
}

#the workhorse for lastz (multithreads)
sub lastz_workhorse($$$){ #parameters in this order: $target, $query, $chrT, $thread_id, $cat_string, $option_string

	my ($target_fa, $query_fa, $thread_id, $option_string) = @_;
	
	my $cmd;
	## unmasking or not
	if($unmask eq "unmask"){
	  $cmd  = "lastz $target_fa"."[multiple,$unmask] ";
	  $cmd .= "$query_fa"."[$unmask] $option_string ";
	}else{
	  $cmd  = "lastz $target_fa"."[multiple] ";
	  $cmd .= "$query_fa $option_string ";	
	}
	$cmd .= "1>> $target.$query.result/raw.axt/$target_fa.$query_fa.axt 2>> $target.$query.result/raw.axt/$target_fa.$query_fa.log";
	print "Thread_ID $thread_id : $cmd \n";
	system("$cmd");
	
	return($thread_id);
}

#! /usr/bin/env bash

cat <<EOF
### This file is MY default parameter setting file for all_lastz.pl!

##(!!!actually this section is never useful)if want to infer score matrix, uncomment the following two line 
#--inferonly				 	# when infering score, lastz_D is used instead of lastz, gap penalty should not be inferred in present. 
#--infscores[=<output_file>], 	# inferred score output to this file. 

##dynamic masking
# equivalent to the BLASTZ_M, default is 0, when required, set to 50 normally, some set to 254
--masking=254					
# --notrivial is required for haploMerger, should not changed!
--notrivial

##scoring parameters for tuning
#--scores=<scoring_file>       	# should not specify when using scoreMatrix inference procedure, i.e. --infer
#--inner=2000					# equivalent to the BLASTZ_H, normally H=2000
--hspthresh=3000				# equivalent to the BLASTZ_K, 3000 by default
--ydrop=3400					# equivalent to the BLASTZ_Y, O+300E by default
--gappedthresh=3000				# equivalent to the BLASTZ_L, L=K by default

##################################################
##scoring parameters usually don't need change

# equivalent to BLASTZ_O and BLASTZ_E, 400 and 30, by default
--gap=400,30					
#--nogapped

#Offset between the starting positions of successive target words considered for potential seeds. 
#(But this does not apply to the query words, which always use a step size of 1.) 
# equivalent to the BLASTZ_Z, z=1 by default
# set to 20 can reduce the memory consumtion by a factor of 3
--step=20

# by default --seed=12of19	, or --seed=14of22					
--seed=12of19					

# default setting=--transition, alternative --notransition
# set to --notranstion can shorten the time by a factor of 10
--notransition					

################################################## 
##non scoring parameters usually don't need change
--format=axt					# lav (default), axt, maf, and so on
--markend						# Just before normal completion, write "# lastz end-of-file" to the output file.
--verbosity=1					# logfile details
EOF
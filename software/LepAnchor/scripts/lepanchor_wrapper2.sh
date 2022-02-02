#!/bin/bash
##########################################################
#
#  Lep-Anchor wrapper, version 2
#
#  usage: lepanchor_wrapper2.sh -f ref.fasta -e extra_cut_site_file -n num_chr -c chain_file -p paf_file -r prox_file -m map_file1 -m map_file2, ... 
#
#  output: REF_LA.fa.gz anchored reference genome
#          REF_LA.agp   agp file describing the genome
#          REF_LA_scaffolds.fa.gz anchored reference genome in scaffolds
#          REF_LA_scaffolds.agp anchored reference genome in scaffolds agp file
#          marey*.png   Marey maps for visual verification
#
#          chr*.agp     agp files for each chromosome
#          scaffolds_chr*.agp agp files for each chromosome in scaffolds (each block of linked contigs as a scaffold)
#
#
#  Pasi Rastas, (c) 2021, pasi.rastas@gmail.com
#
##########################################################

#if [ "$#" -ne 3 ]; then
#    echo "At least three input parameters must be provided"
#    exit 1
#fi


#parse parameters

function print_usage()
{
echo "##########################################################"
echo "#"
echo "#  Lep-Anchor wrapper2"
echo "#"
echo "#  usage: lepanchor_wrapper2.sh -T thread_per_run -t threads -f ref.fasta -n num_chr -e extra_cut_site_file -c chain_file -r prox_file -p paf_file -m map_file1 -m map_file2, ... "
echo "#"
echo "#  download Lep-Anchor by"
echo "#  wget https://sourceforge.net/projects/lep-anchor/files/binary%2Bcode.zip/download -O la.zip;unzip la.zip"
echo "#"
echo "#  Pasi Rastas, (c) 2021, pasi.rastas@gmail.com"
echo "#"
echo "##########################################################"
}

while getopts ":e:n:c:p:m:f:t:r:T:" OPTION; do
	case ${OPTION} in
	T)
	THREADS2=$OPTARG;;
	t)
	THREADS=$OPTARG;;
	n)
	CHR=$OPTARG;;
	c)
	CHAIN=$OPTARG;;
	p)
	PAF=$OPTARG;;
        e)
	CUT=$OPTARG;;
	r)
	PROX=$OPTARG;;
	m)
	MAP="$MAP $OPTARG";;
	f)
	REF=$OPTARG;;
	*)
        echo "Incorrect options provided"
        exit 1;;
    esac
done

if [[ $MAP =~ ^$ ]];then
	print_usage
	echo "Please provide at least one map"
	exit 1
fi

if [[ $REF =~ ^$ ]];then
    if [[ ! -e "contigs.length" ]];then
	print_usage
	echo "Please provide either reference fasta or contigs.length file"
	exit 1
    fi
fi

if [[ $CHR =~ ^[0-9]+$ ]];then
	echo "Number of chromosomes = $CHR"
else
	print_usage
	echo "Please provide parameter n (number of chromosomes)"
	exit 1
fi

if [[ $CHAIN =~ ^$ ]];then
	echo "No chain provided"
	CHAIN="/dev/null"
fi

if [[ $PROX =~ ^$ ]];then
        echo "No proximity data provided"
        PROX="/dev/null"
fi


if [[ $PAF =~ ^$ ]];then
	echo "No paf provided" 
	PAF="/dev/null"
fi

#number of threads
if [[ $THREADS =~ ^$ ]];then
	THREADS=8
fi

if [[ $THREADS2 =~ ^$ ]];then
	THREADS2=1
fi

#get Lep-Anchor
#wget https://sourceforge.net/projects/lep-anchor/files/binary%2Bcode.zip/download -O la.zip
#unzip la.zip

#where Lep-Anchor binaries are located
LABIN=bin/

#java runtime located here
JAVA=java

#parallel 
if ! command -v "parallel"
then
	echo "command parallel not found, using only one thread"
	PARALLEL=bash
else
	PARALLEL="parallel --jobs $THREADS"
fi

if [[ ! $REF =~ ^$ ]];then
	echo "calculating contigs.length file..."
	echo "gunzip -fc $REF|awk -f contigLength.awk >contigs.length"|bash
fi


echo "finding full haplotypes..."
echo "gunzip -fc $CHAIN|awk -f findFullHaplotypes.awk >fullHaplotypes50.txt"|bash
wc -l fullHaplotypes50.txt


echo "running liftoverHaplotypes for all input maps..."
for i in $MAP
do
echo "gunzip -fc $CHAIN|$JAVA -cp $LABIN LiftoverHaplotypes map=$i haplotypes=fullHaplotypes50.txt chain=- >$i.liftover"
done|$PARALLEL

#store lift overed maps to $MAPL
for i in $MAP
do
	MAPL="$MAPL $i.liftover"
done

#make input for CleanMap
for i in $MAPL
do
cat $i
done|sort -V -k 1,1 -k 2,2n >map_all_sorted.liftover


#CleanMap
echo "running CleanMap..."
$JAVA -cp $LABIN CleanMap map=map_all_sorted.liftover >map_all.clean

#Map2Bed
echo "running Map2Bed..."
$JAVA -cp $LABIN Map2Bed map=map_all.clean contigLength=contigs.length >map.bed

if [[ ! $CUT =~ ^$ ]];then
echo "adding extra cuts to map.bed..."
awk '{print $1"\t"$2"\t"$3}' $CUT|sort -V -k 1,1 -k 2,2n|awk -f cutBed.awk map.bed - >map.bed.tmp
awk -f cutBed_fixN.awk contigs.length map.bed.tmp >map.bed
fi


#find contigs not put into chromosomes
cut -f 1 contigs.length|grep -v -w -F -f <(cut -f 2 fullHaplotypes50.txt; cut -f 1 map.bed) >not_used.txt
#cut -f 1 contigs.length|grep -v -w -F -f <(cut -f 1 map.bed) >not_used.txt


grep -w -F -f not_used.txt contigs.length|awk -vn=$CHR '{s=$1"\t1\t"$2"\t?\t"; for (i=1;i<=n;++i) print s i}' >chr0.bed
cat map.bed chr0.bed >map_extra.bed

#PlaceAndOrientContigs
echo "running PlaceAndOrientContigs (first iteration)..."
for i in $(seq $CHR)
do
echo "gunzip -fc $CHAIN|$JAVA -cp $LABIN PlaceAndOrientContigs numThreads=$THREADS2 bed=map_extra.bed chromosome=$i map=$MAPL chain=- paf=$PAF proximity=$PROX keepEmptyIntervals=0 >chr$i.la 2>chr$i.la.err"
done|$PARALLEL

#propagate4
echo "running propagate4..."
awk -f propagate4.awk pass=1 chr*[0-9].la pass=2 chr*[0-9].la.err|awk -f pickbed.awk - map_extra.bed >map_extra_prop.bed

#PlaceAndOrientContigs
echo "running PlaceAndOrientContigs (second iteration)..."
for i in $(seq $CHR)
do
echo "gunzip -fc $CHAIN|$JAVA -cp $LABIN PlaceAndOrientContigs numThreads=$THREADS2 $(awk -f pickorientation.awk chr$i.la) bed=map_extra_prop.bed chromosome=$i map=$MAPL chain=- paf=$PAF proximity=$PROX keepEmptyIntervals=1 >ichr$i.la 2>ichr$i.la.err"
done|$PARALLEL

#propagate
echo "running propagate..."

awk -f propagate.awk ichr*.la >tmp1.la
awk -f propagate.awk tmp1.la >tmp2.la
i=2

while ! cmp -s "tmp$i.la" "tmp$(( $i-1 )).la" ;do
	awk -f propagate.awk tmp$i.la >tmp$[$i+1].la
	i=$[$i+1]
done

#create prop*.la, take only contig put to uniquely to a chromosome
#use propagate2 to include possible bridge contigs as well...
awk -f propagate2.awk tmp$i.la|awk '(/^[^#]/ && NF>=8){++d[$1"\t"($7+0)"\t"($8+0)]; data[++line]=$0}END{for (i=1; i<=line; ++i) {$0=data[i];if (d[$1"\t"($7+0)"\t"($8+0)] == 1) {fn="prop"$5".la";print $0>fn}}}'

#create a new bed by combining prop[1-9]*.la and map_extra_prop.bed
awk '{print $1"\t"($7+0)"\t"($8+0)"\t?\t"$5}' prop[1-9]*.la|awk -f pickbed.awk - map_extra_prop.bed >map_extra_prop2.bed

#third iteration could be just improveOrder

#PlaceAndOrientContigs
echo "running PlaceAndOrientContigs (third iteration)..."
for i in $(seq $CHR)
do
echo "gunzip -fc $CHAIN|$JAVA -cp $LABIN PlaceAndOrientContigs numThreads=$THREADS2 $(awk -f pickorientation.awk chr$i.la) bed=map_extra_prop2.bed chromosome=$i map=$MAPL chain=- paf=$PAF proximity=$PROX keepEmptyIntervals=1 evaluateAnchoring=prop$i.la improveAnchoring=1 numRuns=1 >iichr$i.la 2>iichr$i.la.err"
done|$PARALLEL

#pruning contig blocks without map support
for i in $(seq $CHR)
do
        awk -f prune.awk iichr$i.la >iichr${i}_pruned.la
done 2>pruned.la

#remove overlap(s)
awk -f removeOverlaps.awk map_extra_prop2.bed iichr*_pruned.la >iiall.la

#construct agp files
for i in $(seq $CHR)
do
awk -vn=$i '($5==n)' iiall.la|awk -vprefix="LG" -vlg=$i -f makeagp_full2.awk - >chr$i.agp
awk -vn=$i '($5==n)' iiall.la|awk -vprefix="LG" -vlg=$i -f makeagp2.awk - >scaffolds_chr$i.agp
done

#find contigs not used
cut -f 1 contigs.length|grep -v -w -F -f <(cut -f 2 fullHaplotypes50.txt;awk '($5!="U"){print $6}' chr*.agp) >not_used_final.txt

grep -F -w -f not_used_final.txt contigs.length|awk '{print $1,1,$2,1,"W",$1,1,$2,"+"}' >not_used.agp

cat chr*.agp not_used.agp >REF_LA.agp
#one could use scaffolds_chr*.agp as well instead of chr*.agp
cat scaffolds_chr*.agp not_used.agp >REF_LA_scaffolds.agp

#make final fasta
if [[ ! $REF =~ ^$ ]];then
	echo "constructing final fasta (REF_LA.fa.gz)..."
	echo "gunzip -fc $REF|awk -f makefasta.awk - REF_LA.agp|gzip >REF_LA.fa.gz"|bash

	echo "constructing final fasta (REF_LA_scaffolds.fa.gz)..."
	echo "gunzip -fc $REF|awk -f makefasta.awk - REF_LA_scaffolds.agp|gzip >REF_LA_scaffolds.fa.gz"|bash
fi

#construct Marey map
echo "constructing Marey maps... (marey*.png)"
j=1
for m in $MAP
do
	for c in $(seq $CHR)
	do
#	awk -vn=$c '($3==n)' $m.liftover|awk -f liftover.awk chr$c.agp -|awk -vm=$j '(/LG/ && NF>=4){if (NF==4) $5=$4;print $1"\t"$2"\t"$3"\t"m"\t"$4"\t"$5}'
	awk -vn=$c '($3==n)' $m.liftover|awk -f liftover.awk chr$c.agp -|awk -vm=$j '(/LG/ && NF>=4){if (NF==4) {print $1"\t"$2"\t"$3"\t"m"\t"$4"\t"$4} else for (i=4;i<NF;i+=2) print $1"\t"$2"\t"$3"\t"m"\t"$(i)"\t"$(i+1)}'
	done
	j=$[$j + 1]
done|gzip >marey.data.gz

Rscript plot_marey.R

#TODO: lepanchor_wrapper_step2.sh


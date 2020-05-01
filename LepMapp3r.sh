#! /usr/bin/env bash
## This bash wrapper is will act like a liason to run LepMap3's various modules. Strongly recommended to run this in a screen environment. 

parentcall(){
    printf "->"
    printf "\033[01;33m"
    printf "\nParentCall"
    printf "\033[0m"
    printf 'What is the name of the filtered vcf file?  '
    read -r VCFFILE
    if [ ! -f $VCFFILE ]; then
        echo "Error: $VCFFILE not found in working directory. Please check spelling." && exit 1
    fi
    java $LM3PATH ParentCall2 data=${PEDIGREE} vcfFile=${VCFFILE} removeNonInformative=1 | gzip >data.call.gz

    echo -e "\n-- Initiating Filtering2 module --"
    printf '\nSpecify your data tolerance (0.0001 to 0.01):  '
    read -r DTOLERANCE
    zcat data.call.gz | java $LM3PATH Filtering2 data=- dataTolerance=$DTOLERANCE | gzip >data_f.call.gz
}

export -f parentcall

separatechromosomes(){
    printf "\nParentCall->"
    printf "\033[01;33m" 
    printf " SeparateChromosomes\n" 
    printf "\033[0m"
    echo -e "\nChromosome separation can be iterated over a range of LOD score limits to find the best map"
    printf 'What LOD limit do you want to start with?  '
    read -r LODSTART
    printf 'What LOD limit do you want to end with?  '
    read -r LODEND 
    echo -e "\nThis may take a while depending on your data and the range of LOD scores you're exploring."
    printf 'How many CPUs would you like to use per iteration (max=%s)?  ' "$(nproc)"
    read -r NBPROCS
    mkdir -p maps.splitchrom
    for i in $(seq $LODSTART $LODEND)
        do
        printf "\033[01;33m" 
        printf "Running SeparateChromosomes 2 with LOD limit=%s\n" "$i" 
        printf "\033[0m"
        zcat data_f.call.gz | java $LM3PATH SeparateChromosomes2 data=- lodLimit=$i distortionLod=1 numThreads=$NBPROCS > maps.splitchrom/map.$i
        # this exhausted pipe summarizes the maps, removes leading whitespaces, and sorts by LG
        sed '1,1d' maps.splitchrom/map.$i | sort | uniq -c | sed 's/^[[:space:]]*//' | sort -k2n > maps.splitchrom/.map.$i.summary.txt
        # prepend column names
        sed  -i "1i map.$i LG" maps.splitchrom/.map.$i.summary.txt
    done
    
    # initialize max number of LG's in summary file
    cut -f2 -d " " maps.splitchrom/.map.$LODEND.summary.txt > maps.splitchrom/maps.summary.txt
    
    # merge summaries into one file
    for summfile in $(ls ./maps.splitchrom/.map.*.summary.txt | sort -V )
        do
        # append the marker numbers onto the summary file
        cut -f1 -d " " $summfile | paste maps.splitchrom/maps.summary.txt -  > maps.splitchrom/summtemp && mv maps.splitchrom/summtemp maps.splitchrom/maps.summary.txt && rm $summfile
    done
    # replace all spaces with tabs
    sed -i 's/ /\t/g' maps.splitchrom/maps.summary.txt

    echo -e "\nExamine the maps produced ("maps.splitchrom/maps.summary.txt") and decide on the best map before proceeding"
    echo "if using a screen environment, press ctrl+a then d to detach the screen and when ready return to it with the command "screen -r" "
    echo "if not using a screen, press ctrl+z to pause followed by the "bg" command to put in the background. When ready, enter the "fg" command to bring LepMapp3r into the foreground again"
    echo "Press Enter to proceed"
    read PROCEED
}

export -f separatechromosomes

joinsingles(){
    printf  "\n ParentCall->SeparateChromosomes->"
    printf "\033[01;33m" 
    printf " JoinSingles\n" 
    printf "\033[0m"
    echo -e "\n------------List of map files in maps.splitchrom/------------"
    find ./maps.splitchrom -name "map.*" -printf "%f\n"
    echo -e "\n------------------------------------------------------\n"
    printf 'Which map would you like to use (just the filename, exclude directory)?  '
    read -r BESTMAP
    if [ ! -f maps.splitchrom/$BESTMAP ]; then
        echo "Error: maps.splitchrom/$BESTMAP not found in working directory. Please check spelling." && exit 1
    fi
    printf 'LOD limit (4 is common)?  '
    read -r LODLIMIT
    printf 'LOD difference cutoff (2 is common)?  '
    read -r LODDIFF
    printf 'How many CPUs would you like to use per iteration (max=%s)?  ' "$(nproc)"
    read -r PROCS
    DATAMAP=$(echo "$BESTMAP.master")
    zcat data_f.call.gz | java $LM3PATH JoinSingles2All map=maps.splitchrom/$BESTMAP data=- lodLimit=$LODLIMIT lodDifference=$LODDIFF iterate=1 numThreads=$PROCS > $DATAMAP
    cut -f 1 $DATAMAP | sort | uniq -c | sort -k2n
    echo "Your filtered map can be found in ./$BESTMAP.master"
}

export -f joinsingles

ordermarkers(){
    printf "\nParentCall->SeparateChromosomes->JoinSingles->"
    printf "\033[01;33m" 
    printf " OrderMarkers\n" 
    printf "\033[0m"
    echo -e "\nThis step will order the markers on linkage groups 1:N"
    printf 'What distance-calculating method would you like to use (default, morgan, or kosambi)?  '
    read -r DISTTYPE
    case $DISTTYPE in
        kosambi)
            echo "using the Kosambi distance function"
            DISTMETHOD="useKosambi=1"
        ;;
        morgan)
            DISTMETHOD="useMorgan=1"
            echo "using Morgan linear distance"
        ;;
        *)
            echo "using the default LepMap3 distance calculation method"
            DISTMETHOD=""
    esac
    printf 'You want to order the markers on linkage groups 1 to... ?  '
    read -r NUMCHROM
    printf 'How many iterations per linkage group (>30 preferred)? '
    read -r NUMITER
    printf 'This step runs in parallel. The number of chromosomes being run at once must be <= number CPUs '
    printf '\nSubtract 1.5x from how many linkage groups you would like to order in parallel at once. i.e. if you input "10", it will use 15. (max CPUs=%s)?  ' "$(nproc)"
    read -r NUMJOBS
    THREADS=$(echo $(nproc) $NUMCHROM 1.5 | awk '{ printf("%.0f\n", $1/($2/$3)) }')
    NUMLOCI=$(tail -n +2 ./$DATAMAP | wc -l)
    NUMINDS=$(( $(head -n 1 $PEDIGREE | awk '{print NF}') - 2 ))
    SCALEVAL=$(echo $NUMINDS $NUMLOCI | awk '{ print $1/$2 }')
    printf 'Data scaling will be set to %s based on %s markers and %s individuals\n' "$SCALEVAL" "$NUMLOCI" "$NUMINDS"
    echo -e "Depending on  the number of input CPUs iterations, this may take a WHILE. \nFeel free to detach the screen."
    mkdir -p ordermarkers/bestlikelihoods
    mkdir -p ordermarkers/logs
    DATAMAP="$(ls map.*.master)"
    export DATAMAP
    export THREADS
    export DISTMETHOD
    export SCALEVAL
    export NUMITER
    runordermarkers() {
        zcat data_f.call.gz | java $LM3PATH OrderMarkers2 map=$DATAMAP data=- numThreads=$THREADS $DISTMETHOD scale=$SCALEVAL chromosome=$2 &> ordermarkers/logs/ordered.$2.$1.log
        grep -A 100000 \*\*\*\ LG\ \= ordermarkers/logs/ordered.$2.$1.log > ordermarkers/ordered.$2.$1.txt
    }
    export -f runordermarkers
    seq 1 $NUMCHROM | parallel --jobs $NUMJOBS runordermarkers {1} {2} ::: $(seq 1 $NUMITER) ::: $(seq 1 $NUMCHROM)

    # get likelihoods
    for i in $(seq 1 $NUMCHROM)
        do
        for j in $(seq 1 $NUMITER)
            do
            LG="ordered.$i"
            ITERUN="$j"
            LIKELIHOOD=$(head -1 ordermarkers/ordered.$i.$j.txt | tail -1 | cut -c 27-) 
            echo -e "$LG\t$ITERUN\t$LIKELIHOOD" >> ordermarkers/likelihoods.txt
        done
    done
    sort ordermarkers/likelihoods.txt -k1,1 -k3,3nr > ordermarkers/likelihoods.sorted.txt
}

export -f ordermarkers

trimming(){    
    # pull out best maps for each linkage group
    TOTALMAPS=$(cat ordermarkers/likelihoods.txt | wc -l)
    NUMITER=$(ls ordermarkers/ordered.1.*.txt | wc -l) 
    echo "Best ordered maps:"
    for i in $(seq 1 $NUMITER $TOTALMAPS); 
        do
        LIKELYMAP=$(sed -n ${i}p ordermarkers/likelihoods.sorted.txt | cut -f1,2 | awk '{print $0, $1 "." $NF}' | cut -d ' ' -f2)
        echo "$LIKELYMAP.txt"
        cp ordermarkers/$LIKELYMAP.txt ordermarkers/bestlikelihoods
    done
    printf "\nParentCall->SeparateChromosomes->JoinSingles->OrderMarkers->"
    printf "\033[01;33m" 
    printf " Trimming Ends\n" 
    printf "\033[0m"

    # run trimming script
    printf '\nAt what map distance cutoff would you like to trim end-clusters of markers? (10 is conservative)? '
    read -r TRIMTHRESH
    mkdir -p ./ordermarkers/best.trimmed
    Rscript /bin/LepMapp3rQA.r $(pwd)/ordermarkers/bestlikelihoods ordered $TRIMTHRESH > Trimming.log
    mv ./ordermarkers/bestlikelihoods/trimmed* ./ordermarkers/bestlikelihoods/bad.markers.txt ./ordermarkers/bestlikelihoods/trimming.plots.pdf ./ordermarkers/best.trimmed
}    

export -f trimming
    
# reorder the trimmed markers
reordering(){
    mkdir -p reordermarkers/bestlikelihoods reordermarkers/best.trimmed
    mkdir -p reordermarkers/logs
    printf "\nParentCall->SeparateChromosomes->JoinSingles->OrderMarkers->Trimming Ends->"
    printf "\033[01;33m" 
    printf " Reordering\n" 
    printf "\033[0m"
    NUMITER=$(ls ordermarkers/ordered.1.*.txt | wc -l)
    TOTALFILES=$(ls ordermarkers/ordered.*.txt | wc -l) 
    NUMCHROM=$(echo $NUMITER $TOTALFILES | awk '{ print $2/$1 }')
    THREADS=$(echo $(nproc) $NUMCHROM 1.5 | awk '{ printf("%.0f\n", $1/($2/$3)) }')
    NUMINDS=$(( $(head -n 1 $PEDIGREE | awk '{print NF}') - 2 ))
    NUMLOCI=$(tail -n +2 ./$DATAMAP | wc -l)
    SCALEVAL=$(echo $NUMINDS $NUMLOCI | awk '{ print $1/$2 }')
    printf '\nWhat distance-calculating method would you like to use (default, morgan, or kosambi)?  '
    read -r DISTTYPE
    case $DISTTYPE in
        kosambi)
            echo "using the Kosambi distance function"
            DISTMETHOD="useKosambi=1"
        ;;
        morgan)
            DISTMETHOD="useMorgan=1"
            echo "using Morgan linear distance"
        ;;
        *)
            echo "using the default LepMap3 distance calculation method"
            DISTMETHOD=""
    esac
    printf '\nHow many linkage groups would you like to order in parallel at once (max CPUs=%s)?  ' "$(nproc)"
    read -r NUMJOBS
    echo "This will take about as long as before."
    export DATAMAP
    export THREADS
    export DISTMETHOD
    export SCALEVAL
    export NUMITER
    export NUMCHROM
    export NUMJOBS
    reordermarkers() {
	    LOGFILE=$(echo "reordered.$(basename $1).$2.log" | cut -d'.' -f1,4,7,8)
	    OUTFILE=$(echo "reordered.$(basename $1).$2.txt" | cut -d'.' -f1,4,7,8)
            zcat data_f.call.gz | java $LM3PATH OrderMarkers2 evaluateOrder=$1\
            map=$DATAMAP \
            data=- numThreads=$THREADS $DISTMETHOD \
            scale=$SCALEVAL &> reordermarkers/logs/$LOGFILE
	    grep -A 100000 \*\*\*\ LG\ \= reordermarkers/logs/$LOGFILE > reordermarkers/$OUTFILE
    }
    export -f reordermarkers
    
    parallel --group --jobs $NUMJOBS reordermarkers {2} {1} ::: $(seq 1 $NUMITER) ::: $(find ./ordermarkers/best.trimmed -name "trimmed.*.txt")

    # get reordered linkage group likelihoods
    for i in $(seq 1 $NUMCHROM)
        do
        for j in $(seq 1 $NUMITER)
            do
            LG="reordered.$i"
            ITERUN="$j"
            LIKELIHOOD=$(head -1 reordermarkers/$LG.$j.txt | cut -c 27-) 
            echo -e "$LG\t$ITERUN\t$LIKELIHOOD" >> reordermarkers/reordered.likelihoods.txt
        done
    done

    # sort by linkage group and likelihood
    sort reordermarkers/reordered.likelihoods.txt -k1,1 -k3,3nr > reordermarkers/reordered.likelihoods.sorted.txt

    # pull out best maps for each linkage group
    TOTALMAPS=$(cat reordermarkers/reordered.likelihoods.sorted.txt | wc -l)
    echo "Finding best maps by likelihoods"
    echo "Best reordered maps:"
    for i in $(seq 1 $NUMITER $TOTALMAPS); 
        do
        LIKELYMAP=$(sed -n ${i}p reordermarkers/reordered.likelihoods.sorted.txt | cut -f1,2 | awk '{print $0, $1 "." $NF}' | cut -d ' ' -f2)
        echo "$LIKELYMAP"
        cp reordermarkers/$LIKELYMAP.txt reordermarkers/bestlikelihoods
    done

    # trim markers with 2 zero positions
    ## REMOVED
    #for i in reordermarkers/bestlikelihoods/reordered* ; do
        # remove markers with position 0
    #FILE=$(basename $i)
    #    sed '/0.000\t0.000/d' ./reordermarkers/bestlikelihoods/$FILE > reordermarkers/best.trimmed/$FILE
        # append the trimmed markers to a log file
    #    sed -n '/0.000\t0.000/p' ./reordermarkers/bestlikelihoods/$FILE  | cut -f1 >> ./reordermarkers/best.trimmed/bad.markers2.txt
    #done
   
    ## merge all trimmed markers
    #cat ./ordermarkers/best.trimmed/bad.markers.txt ./reordermarkers/best.trimmed/bad.markers2.txt > Trimmed.Markers.log
    cat ./ordermarkers/best.trimmed/bad.markers.txt > Trimmed.Markers.log

    ## calculate map distance, sex specific map distance, and intervals
    mkdir -p distances distances_sexAveraged intervals
    distancecalc(){
        zcat data_f.call.gz | java $LM3PATH OrderMarkers2 data=- evaluateOrder=$1  improveOrder=0 $DISTMETHOD > distances/$(basename ${1%.txt}).distance
        zcat data_f.call.gz | java $LM3PATH OrderMarkers2 data=- evaluateOrder=$1  improveOrder=0 sexAveraged=1 $DISTMETHOD > distances_sexAveraged/$(basename ${1%.txt}).sexdistance
        zcat data_f.call.gz | java $LM3PATH OrderMarkers2 data=- evaluateOrder=$1  calculateIntervals=intervals/$(basename ${1%.txt}).intervals $DISTMETHOD 
    }
    
    export -f distancecalc
    find ./reordermarkers/bestlikelihoods -maxdepth 1 -name "reordered.*.txt" | parallel --jobs $NUMJOBS distancecalc &> distance.log

    echo -e "\n LepMapp3r is all done! \n"
}

export -f reordering

# Run everything
if [ $# -eq 0 ]; then
    LM3PATH="-cp /bin/LM3"
    echo -e "\n------------List of all files in directory------------"
    ls -p | grep -v /
    echo -e "------------------------------------------------------"
    printf '\nWhat is the name of the pedigree file?  '
    read -r PEDIGREE
    if [ ! -f $PEDIGREE ]; then
        echo "Error: $PEDIGREE not found in working directory. Please check spelling." && exit 1
    fi
    export LM3PATH
    export PEDIGREE
    if [ ! -f "data_f.call.gz" ]; then
        parentcall && separatechromosomes && joinsingles && ordermarkers && trimming && reordering
    else
        if [ ! -d "maps.splitchrom" ]; then
            echo "file "data_f.call.gz" detected, skipping ParentCall..."
            separatechromosomes && joinsingles && ordermarkers && trimming && reordering
        elif test -n "$(find ./ -maxdepth 1 -name 'map.*.master' -print -quit)"; then
            echo -e "\n------------List of all maps in working directory------------"
            ls map.*.master
            echo -e "----------------------------------------------------------------"
            printf 'Which map would you like to use for ordering (map.*.master)? '
                read -r BESTMAP
            if [ ! -f $BESTMAP ]; then
                echo "Error: $DATAMAP not found in working directory. Please check spelling." && exit 1
            else
                DATAMAP=$(echo "$BESTMAP")
            fi
            if [ ! -d "ordermarkers" ]; then
                ordermarkers && trimming && reordering
            elif [ "$(ls -A "ordermarkers/best.trimmed" )" ]; then
                echo "ordered maps in "ordermarkers/best.trimmed" detected, skipping OrderMarkers..."
                reordering
            else
                echo "directory "ordermarkers" detected, skipping OrderMarkers..."
                trimming && reordering
            fi
        else
            echo "directory "maps.splitchrom" detected, skipping SeparateChromosomes..."
            joinsingles && ordermarkers && trimming && reordering
        fi
    fi
    exit 
else
            cat <<EOF

LepMapp3r is a wrapper for LepMap3 (Rasta 2017) intended to link the modules together into a single workflow. To use it correctly, LepMap3 and LepMapperQA.r needs to be installed in /bin/LM3, however you can always modify LepMapp3r to point to where your installation is by editing Lines 173 and 267. To use LepMapp3r, simply run the command with no arguments.

                             [LepMapp3r workflow]
ParentCall -> SeparateChromosomes -> JoinSingles -> OrderMarkers -> Trimming Ends -> Reordering

By default, if you run LepMapp3r without arguments it will start at ParentCall and work through until the end. Ordering markers may take a while, so it is recommended to run LepMapp3r in a screen environment. LepMapp3r creates several folders during operation, and for safety, it is made to identify if these folders exist so as to skip that step and not overwrite data. Please remove or rename these folders as necessary to avoid unintentionally skipping steps: 
maps.splitchrom
ordermarkers/best.trimmed 
ordermarkers/bestlikelihoods
reordermarkers/bestlikelihoods

                                [Disclaimer] 
LepMap3 is a very comprehensive software, and LepMapp3r does not incorporate all the features and nuances within the various modules. Your study is unique, so you are encouraged to fork the GitHub repo and adapt LepMapp3r to your needs! If using LepMapp3r and publishing, cite Pasi Rastas for his work on LepMap3, and if you like using it, give Pavel Dimens a shout out on Twitter @pvdimens =)

                                 [Citation] 
Pasi Rastas, Lep-MAP3: robust linkage mapping even for low-coverage whole genome sequencing data, Bioinformatics, Volume 33, Issue 23, 01 December 2017, Pages 3726–3732, https://doi.org/10.1093/bioinformatics/btx494

EOF
fi

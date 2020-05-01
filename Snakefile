from os import path
import glob
#FILENAMES = [os.path.basename(i) for i in glob.iglob("test/*.R1.fq.gz")]
#FILEBASENAMES = [i.replace('.R1.fq.gz', '') for i in FILENAMES]

vcf_file = os.path.basename(glob.glob("./*.vcf"))

rule all:
    input:
        "data_f.call.gz"


rule parentcall:
    input:
        vcf = expand("{vcffile}", vcffile = vcf_file),
        pedigree = "pedigree.txt"
    output:
        "data.call.gz"
    shell:
        "java -cp ./LM3 ParentCall2 data={input.pedigree} vcfFile={input.vcf} removeNonInformative=1 | gzip > data.call.gz"

rule filtering:
    input:
        "data.call.gz"
    output:
        "data_f.call.gz"
    shell:
        "echo -e -n "\nSpecify your data tolerance (0.0001 to 0.01):  " "
        "read -r "
        "zcat data.call.gz | java -cp ./LM3 Filtering2 data=- dataTolerance=$REPLY | gzip > data_f.call.gz"




#rule separatechromosomes:
#    input:

#    output:

#    shell:




"""
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
"""
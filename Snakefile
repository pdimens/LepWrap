from os import path
import glob
#FILENAMES = [os.path.basename(i) for i in glob.iglob("test/*.R1.fq.gz")]
#FILEBASENAMES = [i.replace('.R1.fq.gz', '') for i in FILENAMES]


rule all:
    input:
        "data.call.gz"


rule parentcall:
    input:
        vcf = "YFT70_maxmiss80.recode.vcf"
        pedigree = "pedigree.txt"
    output:
        "data.call.gz"
    shell:
        "java -cp ./LM3 ParentCall2 data={input.pedigree} vcfFile={input.vcf} removeNonInformative=1 | gzip > data.call.gz"

#rule filtering:
#    input:
#        "data.call.gz"
#    output:
#        "data_f.call.gz"
#    shell:
#        'echo -en "\nSpecify your data tolerance (0.0001 to 0.01):  " '
#        "read -r "
#        "zcat data.call.gz | java -cp ./LM3 Filtering2 data=- dataTolerance=$REPLY | gzip > data_f.call.gz"




#rule separatechromosomes:
#    input:

#    output:

#    shell:
rule parentcall:
    input:
        vcf = {config["vcf"]},
        pedigree = {config["predigree"]}
    output:
        "data.call.gz"
    message:
        """
        Creating Lep-Map3 data file from VCF and pedigree files
        """
    shell:
        "java -cp LM3 ParentCall2 data={input.pedigree} vcfFile={input.vcf} removeNonInformative=1 | gzip > data.call.gz"

rule filtering:
    input:
        "data.call.gz"
    output:
        "data_f.call.gz"
    message:
        """
        Filtering the data
        """
    shell:
        """
        echo -e -n '\nSpecify your data tolerance (0.0001 to 0.01):  '
        read -r
        zcat {input} | java -cp LM3 Filtering2 data=- dataTolerance=$REPLY | gzip > data_f.call.gz
        """
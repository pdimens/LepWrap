rule parent_call:
    input:
        vcf = vcf,
        pedigree = pedigree
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
    params:
        data_tolerance = data_tol
    shell:
        """
        if [ {params} == 0 ]; then
            ln -sr {input} {output}
        else
            zcat {input} | java -cp LM3 Filtering2 data=- dataTolerance={params} | gzip > {output}
        fi
        """

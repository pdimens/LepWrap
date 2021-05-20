rule parent_call:
    input:
        vcf = vcf,
        pedigree = pedigree
    output:
        "1_ParentCall/data.lepmap3.gz"
    message: "Creating Lep-Map3 data file from {input.vcf} and {input.pedigree}"
    shell:
        "java -cp LM3 ParentCall2 data={input.pedigree} vcfFile={input.vcf} removeNonInformative=1 | gzip > {output}"

rule filtering:
    input: "1_ParentCall/data.lepmap3.gz"
    output: "2_Filtering/data.filtered.lepmap3.gz"
    message: "Filtering {input}"
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

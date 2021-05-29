rule parent_call:
    input:
        vcf = vcf,
        pedigree = pedigree
    output:
        "1_ParentCall/data.lepmap3.gz"
    message: "Creating Lep-Map3 data file from {input.vcf} and {input.pedigree}"
    params:
        extra = parentcall_extra
    shell:
        "java -cp software/LepMap3 ParentCall2 data={input.pedigree} vcfFile={input.vcf} {params} | gzip > {output}"

rule filtering:
    input: "1_ParentCall/data.lepmap3.gz"
    output: "2_Filtering/data.filtered.lepmap3.gz"
    message: "Filtering {input}"
    params:
        data_tolerance = data_tol,
        extra = filtering_extra
    shell:
        """
        if [ {params.data_tolerance} == 0 ]; then
            echo "Skipping Filtering2 and creating symlink {output} instead"
            ln -sr {input} {output}
        else
            zcat {input} | java -cp software/LepMap3 Filtering2 data=- dataTolerance={params.data_tolerance} {params.extra} | gzip > {output}
        fi
        """

![logo](.misc/logo.png)

[![alt text](https://img.shields.io/badge/docs-wiki-75ae6c?style=for-the-badge&logo=Read%20The%20Docs)](https://github.com/pdimens/LepMak3r/wiki) 

# LepMak3r

LepMak3r is a reusable pipeline to use the linkage map software [Lep-Map3](https://sourceforge.net/projects/lep-map3/). It is the Snakemake-based successor to [LepMapp3r](https://github.com/pdimens/LepMapp3r). Check out [the wiki](https://github.com/pdimens/LepMak3r/wiki) for detailed installation, usage, and workflow information.



### How to install
You will need a `Snakemake` installation, along with cloning this repository locally.
#### 1. Installing snakemake
Assuming you have `anaconda` or `miniconda` installed:
```bash
conda create -n environment_name -c bioconda snakemake
```
where `environment_name` is whatever name you want to call the the snakemake conda environment.
#### 2. Cloning LepMak3r
Download a zip of this repository using the "Code" button on the top-right and unzip it on your machine or:
```bash
git clone https://github.com/pdimens/LepMak3r.git
```

### How to run
You will need to modify `config.yaml` to suit your needs, then you can simply run the pipeline with the wrapper:
```bash
LepMak3r <number_of_cores>
```
where `<number_of_cores>` is an integer of the maximum number of cores/threads you want the pipeline to use.

### Something to keep in mind

Lep-Map3 is a **very** comprehensive software, and LepMak3r does incorporate all the features and nuances. Your study is unique, so you are encouraged to clone or fork this repo and adapt LepMak3r to your needs! All of the code in LepMak3r is written in human-readable bash or aggressively annotated R, so give it a shot and adapt it to your workflow. If using LepMak3r in a publication, cite **Pasi Rastas** for his work on Lep-Map3 and please include a link to this repository. If you like using it, give me (Pavel) a shout out on Twitter [@pvdimens](https://twitter.com/PVDimens) [![alt text](http://i.imgur.com/wWzX9uB.png)](https://twitter.com/PVDimens)  =)



## Citation

Pasi Rastas, Lep-MAP3: robust linkage mapping even for low-coverage whole genome sequencing data, Bioinformatics, Volume 33, Issue 23, 01 December 2017, Pages 3726–3732,https://doi.org/10.1093/bioinformatics/btx494

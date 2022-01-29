![logo](.misc/logo.png)

_It's Lep-Map3 (and Lep-Anchor), but with snakes 🐍🐍_

[![alt text](https://img.shields.io/badge/docs-wiki-75ae6c?style=for-the-badge&logo=Read%20The%20Docs)](https://github.com/pdimens/LepWrap/wiki) 

# LepWrap

LepWrap is a reusable pipeline to use the linkage map software [Lep-Map3](https://sourceforge.net/projects/lep-map3/) and the genome assembly map-based anchoring and orienting software [Lep-Anchor](https://sourceforge.net/p/lep-anchor/wiki/Home/). It is the Snakemake-based successor to [LepMapp3r](https://github.com/pdimens/LepMapp3r). Check out [the wiki](https://github.com/pdimens/LepWrap/wiki) for detailed installation, usage, and workflow information.



### How to install
You will need a `conda` installation ([Anaconda](https://docs.anaconda.com/anaconda/install/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html), I recommend Miniconda), along with downloading the latest release or cloning this repository locally. Using the latest release is a more stable solution.

#### 1. Cloning LepWrap
Download a zip of this repository using the "Code" button on the top-right and unzip it on your machine or:
```bash
git clone https://github.com/pdimens/LepWrap.git
```

#### 2. Installing other dependencies
Assuming you have `anaconda` or `miniconda` installed:
```bash
cd LepWrap
conda env create -f conda_setup.yml
```
This will create an environment called `lepwrap` that can be activated with:
```bash
conda activate lepwrap
```

### How to run
You will need to modify `config.yml` to suit your needs, then you can simply run the pipeline with the wrapper:
```bash
./LepWrap <number_of_cores> <config.yml>
```
where `<number_of_cores>` is an integer of the maximum number of cores/threads you want the pipeline to use and `<config.yml>` (optional!) is the name of the config file, if it's different than `config.yml`. 

**Examples**
```bash
./LepWrap 15                    # assumes config.yml
./LepWrap 32 nojoinsingles.yml  # specific config file
```
### Something to keep in mind
LepWrap does things a certain way, employing the most common/reasonable way of using Lep-Map3 (and LepAnchor more or less). Current versions are **a lot** more flexible that the predecessors, but might still lack something you need. Your study is unique, and I encourage you to clone/fork this repository and adapt LepWrap to it! All of the code in LepWrap is written in human-readable bash or aggressively annotated R, so give it a shot and adapt it to your workflow. PR's always welcome!


## Citation
If using LepWrap in a publication, cite **Pasi Rastas** for their work on Lep-Map3/Lep-Anchor and please include a link to this repository. If you like using it, please Star the repository and/or give me (Pavel) a shout out on Twitter [@pvdimens](https://twitter.com/PVDimens) [![alt text](http://i.imgur.com/wWzX9uB.png)](https://twitter.com/PVDimens)  =)

> Pasi Rastas, Lep-MAP3: robust linkage mapping even for low-coverage whole genome sequencing data, Bioinformatics, Volume 33, Issue 23, 01 December 2017, Pages 3726–3732,https://doi.org/10.1093/bioinformatics/btx494

> Pasi Rastas, Lep-Anchor: automated construction of linkage map anchored haploid genomes, Bioinformatics, Volume 36, Issue 8, 15 April 2020, Pages 2359–2364, https://doi.org/10.1093/bioinformatics/btz978

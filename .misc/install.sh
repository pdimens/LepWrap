# install LepWrap into conda PATH
mkdir -p $CONDA_PREFIX/bin
# LepWrap executable
cp LepWrap $CONDA_PREFIX/bin/
chmod +x $CONDA_PREFIX/bin/LepWrap
# associated scripts
chmod +x scripts/*
cp scripts/* $CONDA_PREFIX/bin/
# LepMap3 modules and scripts
cp software/LepMap3/*.class $CONDA_PREFIX/bin
cp software/LepMap3/scripts/* $CONDA_PREFIX/bin
# LepAnchor modules and scripts
cp software/LepAnchor/*.class software/LepAnchor/lepanchor_wrapper.sh $CONDA_PREFIX/bin
cp software/LepAnchor/scripts/* $CONDA_PREFIX/bin
cp software/LepAnchor/deps/ucsc_binaries/* $CONDA_PREFIX/bin
cp software/LepAnchor/deps/*.pl software/LepAnchor/deps/Red software/LepAnchor/deps/all_lastz.ctl software/LepAnchor/deps/scoreMatrix.q software/LepAnchor/deps/step* $CONDA_PREFIX/bin
# Snakemake rules
cp rules/LepAnchor/*.smk rules/LepMap3/*.smk $CONDA_PREFIX/bin

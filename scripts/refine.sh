#! /usr/bin/env bash

mkdir -p refine_map

for i in {22..70}; do
    zcat data_f.call.gz | java -cp LM3 SeparateChromosomes2 data=- map=maps.splitchrom/map.26 lg=1 sizeLimit=30 lodLimit=$i distortionLod=1 numThreads=10 > refine_map/map.$i
done

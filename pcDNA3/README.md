Workflow for detecting integration events using isling: https://github.com/aehrc/isling/tree/master

<h3>Instructions:</h3>
This Nextflow workflow comes with a test data set and can be run as follows:

1. Checkout this version of the code in this repository under the master branch: `git checkout https://github.com/mclaugsf/mgc-public.git`
2. Checkout isling: `git checkout https://github.com/aehrc/isling.git` we used git hash `66f983fe7c8f5f2cb26ccf2cd23c6b3603adcb2f` for processing these data.  The parts of isling that get called here are the two python scripts: `${isling_dir}/scripts/find_ints.py` and `${isling_dir}/scripts/filter.py` and we use their Docker container in the workflow `szsctt/isling:latest` but call the code outside of the docker container.
3. Download hsa1 from UCSC: https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/ and create indexes on it.  

#!/bin/bash
MAX_TASKS=500
export TAG=13JAN2023
export SHAKE_SCRATCH=/cluster/scratch/kucerat/proteinshake/$TAG
export SHAKE_STORE=borg:/links/scratch/borg/scratch/Datasets/proteinshake/$TAG
export HTTP_PROXY=proxy.ethz.ch:3128
export HTTPS_PROXY=proxy.ethz.ch:3128
export BATCHSIZE=20000
LOGDIR=./logs/$TAG

mkdir -p $LOGDIR
exec &>> $LOGDIR/submit.txt

datasets=("GeneOntologyDataset" "EnzymeCommissionDataset" "PfamDataset" "ProteinProteinInterfaceDataset" "ProteinLigandInterfaceDataset" "TMAlignDataset" "SCOPDataset" "RCSBDataset" "AlphaFoldDataset")

organisms=("arabidopsis_thaliana" "caenorhabditis_elegans" "candida_albicans" "danio_rerio" "dictyostelium_discoideum" "drosophila_melanogaster" "escherichia_coli" "glycine_max" "homo_sapiens" "methanocaldococcus_jannaschii" "mus_musculus" "oryza_sativa" "rattus_norvegicus" "saccharomyces_cerevisiae" "schizosaccharomyces_pombe" "zea_mays" "swissprot")

#datasets=("GeneOntologyDataset" "EnzymeCommissionDataset" "PfamDataset" "RCSBDataset" "AlphaFoldDataset")
#organisms=("methanocaldococcus_jannaschii")

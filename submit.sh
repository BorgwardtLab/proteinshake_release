#!/bin/bash
MAX_TASKS=1000
export TAG=13JAN2023
export SHAKE_SCRATCH=/cluster/scratch/kucerat/proteinshake/$TAG
export SHAKE_STORE=borg:/links/scratch/borg/scratch/Datasets/proteinshake/$TAG
export HTTP_PROXY=proxy.ethz.ch:3128
export HTTPS_PROXY=proxy.ethz.ch:3128
LOGDIR=./logs/$TAG

mkdir -p $LOGDIR
exec &>> $LOGDIR/submit.txt

datasets=("GeneOntologyDataset" "EnzymeCommissionDataset" "PfamDataset" "ProteinProteinInterfaceDataset" "ProteinLigandInterfaceDataset" "TMAlignDataset" "SCOPDataset" "RCSBDataset" "AlphaFoldDataset")

organisms=("arabidopsis_thaliana" "caenorhabditis_elegans" "candida_albicans" "danio_rerio" "dictyostelium_discoideum" "drosophila_melanogaster" "escherichia_coli" "glycine_max" "homo_sapiens" "methanocaldococcus_jannaschii" "mus_musculus" "oryza_sativa" "rattus_norvegicus" "saccharomyces_cerevisiae" "schizosaccharomyces_pombe" "zea_mays" "swissprot")

for DATASET in "${datasets[@]}"; do
    if [[ $DATASET == "AlphaFoldDataset" ]]; then
        loop=("${organisms[@]}");
    else
        loop=("none")
    fi
    for ORGANISM in "${loop[@]}"; do
        if [[ $DATASET == "AlphaFoldDataset" ]]; then
            NAME=${DATASET}_${ORGANISM};
        else
            NAME=$DATASET
        fi
        echo $NAME
        # download
        sbatch --ntasks 1 --cpus-per-task 20 --mem-per-cpu 3G --time 23:59:00 -o $LOGDIR/${NAME}_download.log -J $NAME --wait --wrap "python -m main --njobs 20 --dataset $DATASET --organism $ORGANISM"
        echo "$DATASET downloaded."
        # parse
        if [[ $DATASET != "AlphaFoldDataset" && $DATASET != "RCSBDataset" ]]; then
            # sequence clustering
            sbatch --ntasks 1 --cpus-per-task 100 --mem-per-cpu 1G --time 3:59:00 -o $LOGDIR/${NAME}_sequence.log -J $NAME --wait --wrap "python -m cluster_sequence --njobs 100 --dataset $DATASET --organism $ORGANISM"
            echo "$DATASET sequence clustering done."
            # structure clustering: prepare
            sbatch --ntasks 1 --cpus-per-task 1 --mem-per-cpu 50G --time 3:59:00 -o $LOGDIR/${NAME}_structure_prepare.log -J $NAME --wait --wrap "python -m cluster_structure --prepare --njobs 1 --dataset $DATASET --organism $ORGANISM"
            # structure clustering: compute
            NUM_TASKS=$(find $SHAKE_SCRATCH/$NAME/jobs -name "*.txt" | wc -l)
            echo "Submitting $NUM_TASKS jobs for structure clustering."
            JOBID=$(sbatch --parsable --ntasks 1 --cpus-per-task 1 --mem-per-cpu 10G --time 3:59:00 -o $LOGDIR/${NAME}_structure_compute.log -J $NAME --array 0-$NUM_TASKS%$MAX_TASKS --wrap "python -m cluster_structure --compute --njobs 1 --dataset $DATASET --organism $ORGANISM")
            # structure clustering: collect
            sbatch --ntasks 1 --cpus-per-task 1 --mem-per-cpu 50G --time 3:59:00 -o $LOGDIR/${NAME}_structure_collect.log -J $NAME --dependency $JOBID --wrap "python -m cluster_structure --collect --njobs 1 --dataset $DATASET --organism $ORGANISM"
        fi
    done
done

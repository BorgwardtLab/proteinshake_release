
CLUSTER_NJOBS=10
DOWNLOAD_NJOBS=10
MEM=10 # 100GB
MEM_PER_CPU=$(($MEM*1024/$CLUSTER_NJOBS))
STRUCTURE_TIME="20:00" #"30-0"
SEQUENCE_TIME="20:00" #"23:00:00"
PARSE_TIME="20:00" #"23:00:00"

export TAG="12JAN2023"
export GLOBAL_SCRATCH="/cluster/scratch/kucerat/proteinshake/$TAG"
export RELEASE_DIR="borg:/links/scratch/borg/scratch/Datasets/proteinshake/$TAG"
LOGDIR="./logs/$TAG"


datasets=("RCSBDataset" "GeneOntologyDataset" "EnzymeCommissionDataset" "PfamDataset" "ProteinProteinInterfaceDataset" "ProteinLigandInterfaceDataset" "TMAlignDataset" "SCOPDataset" "AlphaFoldDataset")
organisms=("arabidopsis_thaliana" "caenorhabditis_elegans" "candida_albicans" "danio_rerio" "dictyostelium_discoideum" "drosophila_melanogaster" "escherichia_coli" "glycine_max" "homo_sapiens" "methanocaldococcus_jannaschii" "mus_musculus" "oryza_sativa" "rattus_norvegicus" "saccharomyces_cerevisiae" "schizosaccharomyces_pombe" "zea_mays" "swissprot")

datasets=("RCSBDataset" "AlphaFoldDataset")
organisms=("methanocaldococcus_jannaschii")

mkdir -p $LOGDIR
mkdir -p $GLOBAL_SCRATCH

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
        # download
        python -m main --njobs $DOWNLOAD_NJOBS --dataset $DATASET --organism $ORGANISM --download
        # parse
        JOBID1=$(sbatch --parsable --ntasks 1 --cpus-per-task $CLUSTER_NJOBS --mem-per-cpu $MEM_PER_CPU --tmp $MEM --time $PARSE_TIME -o $LOGDIR/${NAME}_parse.log -J $NAME --wrap "python -m main --njobs $CLUSTER_NJOBS --dataset $DATASET --organism $ORGANISM")
        # sequence clustering
        JOBID2=$(sbatch --parsable --ntasks 1 --cpus-per-task $CLUSTER_NJOBS --mem-per-cpu $MEM_PER_CPU --tmp $MEM --time $SEQUENCE_TIME -o $LOGDIR/${NAME}_sequence.log -J $NAME --dependency=afterok:$JOBID1 --wrap "python -m cluster_sequence --njobs $CLUSTER_NJOBS --dataset $DATASET --organism $ORGANISM")
        # structure clustering
        JOBID3=$(sbatch --parsable --ntasks 1 --cpus-per-task $CLUSTER_NJOBS --mem-per-cpu $MEM_PER_CPU --tmp $MEM --time $STRUCTURE_TIME -o $LOGDIR/${NAME}_structure.log -J $NAME --dependency=afterok:$JOBID2 --wrap "python -m cluster_structure --njobs $CLUSTER_NJOBS --dataset $DATASET --organism $ORGANISM")
        # cleanup
        JOBID4=$(sbatch --parsable --ntasks 1 --cpus-per-task 1 --mem-per-cpu 5G --tmp 5G --time 3:59:00 -o $LOGDIR/${NAME}_clean.log -J $NAME --dependency=afterok:$JOBID3 --wrap "rm -r $GLOBAL_SCRATCH; rm -r $LOCAL_SCRATCH")
    done
done

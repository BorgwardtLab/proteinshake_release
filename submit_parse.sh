source ./vars.sh

for DATASET in "${datasets[@]}"; do
    if [[ $DATASET == "AlphaFoldDataset" ]]; then
        loop=("${organisms[@]}");
    else
        loop=("none")
    fi
    for ORGANISM in "${loop[@]}"; do
        NAME=${DATASET}_${ORGANISM};
        # download
        echo "Submitting parse:"
        sbatch --ntasks 1 --cpus-per-task 100 --mem-per-cpu 3G --time 23:59:00 -o $LOGDIR/${NAME}_download.log -J $NAME --wrap "python -m main --njobs 100 --dataset $DATASET --organism $ORGANISM"
        echo "$NAME parsed."
    done
done

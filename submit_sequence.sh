source ./vars.sh

for DATASET in "${datasets[@]}"; do
    if [[ $DATASET != "AlphaFoldDataset" && $DATASET != "RCSBDataset" ]]; then
        # sequence clustering
        echo "Submitting sequence clustering:"
        sbatch --ntasks 1 --cpus-per-task 100 --mem-per-cpu 1G --time 3:59:00 -o $LOGDIR/${DATASET}_sequence.log -J $DATASET --wrap "python -m cluster_sequence --njobs 100 --dataset $DATASET"
        echo "$DATASET sequence clustering done."
    fi
done

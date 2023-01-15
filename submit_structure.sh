source ./vars.sh

for DATASET in "${datasets[@]}"; do
        if [[ $DATASET != "AlphaFoldDataset" && $DATASET != "RCSBDataset" ]]; then
            # structure clustering: prepare
            echo "Submitting structure clustering preparation:"
            sbatch --ntasks 1 --cpus-per-task 1 --mem-per-cpu 30G --time 3:59:00 -o $LOGDIR/${DATASET}_structure_prepare.log -J $DATASET --wait --wrap "python -m cluster_structure --prepare --njobs 1 --dataset $DATASET"
            # structure clustering: compute
            echo "Submitting structure clustering computation:"
            NUM_TASKS=$(find $SHAKE_SCRATCH/$DATASET/jobs -name "*.txt" | wc -l)
            echo "Submitting $NUM_TASKS jobs for structure clustering."
            JOBID=$(sbatch --parsable --ntasks 1 --cpus-per-task 1 --mem-per-cpu 10G --time 3:59:00 -o $LOGDIR/${DATASET}_structure_compute.log -J $DATASET --array 0-$NUM_TASKS%$MAX_TASKS --wrap "python -m cluster_structure --compute --njobs 1 --dataset $DATASET")
            # structure clustering: collect
            echo "Submitting structure clustering collection:"
            sbatch --ntasks 1 --cpus-per-task 1 --mem-per-cpu 100G --time 3:59:00 -o $LOGDIR/${DATASET}_structure_collect.log -J $DATASET --dependency $JOBID --wrap "python -m cluster_structure --collect --njobs 1 --dataset $DATASET"
        fi
done

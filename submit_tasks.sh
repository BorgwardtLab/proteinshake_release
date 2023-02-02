source ./vars.sh

for i in "${!tasks[@]}"; do
    TASK="${tasks[i]}"
    DATASET="${task_roots[i]}"
    # download
    sbatch --ntasks 1 --cpus-per-task 2 --mem-per-cpu 10G --time 3:59:00 -o $LOGDIR/${TASK}.log -J $TASK --wrap "python -m main --njobs 2 --dataset $DATASET --task $TASK"
done

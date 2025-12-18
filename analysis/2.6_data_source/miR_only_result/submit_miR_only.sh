#!/bin/bash
#SBATCH --job-name=miR_only
#SBATCH --account=project_2011179
#SBATCH --time=6:00:00
#SBATCH --partition=large
#SBATCH --nodes=1                      # 仅用一个节点
#SBATCH --ntasks=1                     # 每个 array 任务是一个独立任务
#SBATCH --cpus-per-task=16             # 每个任务申请 10 个 CPU 核心
#SBATCH --mem=32G                      # 申请 20GB 内存
#SBATCH --array=1-16              # 运行 100 个子任务 (每个子任务独立运行) 496
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

module load r-env/442


cancer_list=("BLCA" "BRCA" "CESC" "COAD" "ESCA" "KIRC" "KIRP" "LIHC" "LUAD" "LUSC" "PAAD" "PRAD" "READ" "STAD" "THCA" "UCEC")
ncore=16


# 计算任务索引
total_cancer=${#cancer_list[@]}
task_id=$((SLURM_ARRAY_TASK_ID - 1))

cancer=${cancer_list[$task_id]}

# 定义日志文件名
log_dir="/scratch/project_2011179/code/miREA/analysis/2.6_data_source/miR_only_result/slurm_logs"
mkdir -p $log_dir  # 确保 logs 目录存在
log_out="${log_dir}/${cancer}.out"
log_err="${log_dir}/${cancer}.err"

echo "Running enrichment analysis for cancer: $cancer, seed: $seed"
echo "Log output: $log_out"
echo "Error log: $log_err"
seff $SLURM_JOBID

# 运行命令并将标准输出和错误输出分别重定向
srun apptainer_wrapper exec Rscript --no-save /scratch/project_2011179/code/miREA/analysis/2.6_data_source/miR_only_result/get_miR_only_result.R "$cancer" > "$log_out" 2> "$log_err"

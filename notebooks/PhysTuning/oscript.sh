#!/bin/bash

#PBS -N comparison_scripts
#PBS -S /bin/bash
#PBS -q smp
#PBS -l walltime=24:00:00
# email when the job [b]egins and [e]nds, or is [a]borted
#PBS -m bea
#PBS -M sallen@eoas.ubc.ca
#PBS -l nodes=1:ppn=4
#PBS -l pmem=2000mb
# stdout and stderr file paths/names
#PBS -o /scratch/sallen/sallen/202111/stdout
#PBS -e /scratch/sallen/sallen/202111/stderr

WORK_DIR="$HOME/MEOPAR/analysis-susan/notebooks/PhysTuning"

cd ${WORK_DIR}
echo "working dir: $(pwd)"

echo "Starting run at $(date)"

__conda_setup="$('/home/Software/system/Miniconda/3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/Software/system/Miniconda/3/etc/profile.d/conda.sh" ]; then
        . "/home/Software/system/Miniconda/3/etc/profile.d/conda.sh"
    else
        export PATH="/home/Software/system/Miniconda/3/bin:$PATH"
    fi
fi
unset __conda_setup

conda activate py38_sqlAlchemy

echo "Setup Complete"

bash script_runner.sh
echo "Ended run at $(date)"

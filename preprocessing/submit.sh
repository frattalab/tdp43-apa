#!/bin/bash
#Submit to the cluster, give it a unique name
#$ -S /bin/bash

#$ -cwd
#$ -V
#$ -l h_vmem=1.9G,h_rt=20:00:00,tmem=1.9G
#$ -pe smp 2

# join stdout and stderr output
#$ -j y
#$ -R y


if [[ ( $@ == "--help") ||  $@ == "-h" ]]; then
    echo "Usage: source submit.sh CONFIG_PATH (BIND_PATH) (RUN_NAME)"
    echo
    echo "Submit script for UCL cluster"
    echo
    echo "CONFIG_PATH - Required - path to config file for snakemake run (e.g. config/config.coverage.yaml)"
    echo "BIND_PATH - Optional - path to externally mount to singularity containers. Default is /SAN/vyplab"
    echo "RUN_NAME - Optional - argument to name run. Config file for run will be copied to folder containing cluster log files (.submissions/<date><time>/) with run name prefixed"
    echo "-h/--help - print this help message and exit"
    exit 0
fi

# define config path
CONFIG_PATH=$1

# Define bind directory, using default if not provided
if [ "$2" == "" ]; then
    IN_BIND_PATH="/SAN/vyplab"
else
    IN_BIND_PATH=$2
fi

# Define run name, using default if not provided
if [ "$3" == "" ]; then
    RUN_NAME="pas_window_coverages"
else
    RUN_NAME=$3
fi

echo "Constructed config file path: "$CONFIG_PATH
echo "Constructed singularity bind command: --bind "$IN_BIND_PATH

# Create directory for cluster job submission script outputs & copy config file
FOLDER=submissions/$(date +"%Y%m%d%H%M")

mkdir -p ${FOLDER}/${RUN_NAME}
cp $1 ${FOLDER}/${RUN_NAME}/

snakemake \
-p \
-s coverage.smk \
--configfile=$CONFIG_PATH \
--use-singularity \
--singularity-args="--bind $IN_BIND_PATH" \
--jobscript ucl_cluster_qsub.sh \
--cluster-config config/cluster.yaml \
--cluster-sync "qsub -l tmem={cluster.tmem},h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -o $FOLDER {cluster.submission_string}" \
-j 20 \
--rerun-incomplete \
--latency-wait 100 \
--keep-going

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
    echo "Usage: source submit.sh CONFIG_PATH (RUN_NAME)"
    echo
    echo "Submit script for UCL cluster"
    echo
    echo "CONFIG_PATH - Required - path to config file for snakemake run (e.g. config/config.QAPA.yaml)"
    echo "RUN_NAME - Optional - argument to name run. Config file for run will be copied to folder containing cluster log files (.submissions/<date><time>/) with run name prefixed"
    echo "-h/--help - print this help message and exit"
    exit 0
fi

# define config path
CONFIG_PATH=$1



# Define run name, using default if not provided
if [ "$2" == "" ]; then
    RUN_NAME="PEKA"
else
    RUN_NAME=$2
fi

echo "Constructed config file path: "$CONFIG_PATH

# Create directory for cluster job submission script outputs & copy config file
FOLDER=submissions/$(date +"%Y%m%d%H%M")

mkdir -p ${FOLDER}/${RUN_NAME}
cp $1 ${FOLDER}/${RUN_NAME}/

snakemake \
-p \
--configfile=$CONFIG_PATH \
--use-conda \
--jobscript ucl_cluster_qsub.sh \
--cluster-config config/cluster.yaml \
--cluster-sync "qsub -l tmem={cluster.tmem},h_vmem={cluster.h_vmem},h_rt={cluster.h_rt} -o $FOLDER/$RUN_NAME {cluster.submission_string}" \
-j 40 \
--nolock \
--rerun-incomplete \
--latency-wait 100 \
--keep-going

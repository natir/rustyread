#!/bin/bash 

set -e
set -u
set -o pipefail
set -x

function main() {
    if [ "$#" -ne 11 ]; then
	usage
	exit
    fi

    coverage=$1
    thread=$2
    ref=$3
    len_avg=$4
    len_std=$5
    identity=$6
    max=$7
    stdev=$8
    error_model=$9
    qscore_model=${10}
    output=${11}

    pid=$BASHPID
    
    cov_per_thread=$(($coverage / $thread))
    for i in $(seq 1 $thread); do
	run ${cov_per_thread} ${ref} ${len_avg} ${len_std} ${identity} ${max} ${stdev} ${error_model} ${qscore_model} ${pid} $i &
    done

    wait

    cat temp_${pid}_*.fastq > ${output}
    rm temp_${pid}_*.fastq

    cat temp_${pid}_*.log
    rm temp_${pid}_*.log
}

function run() {
    coverage=$1
    ref=$2
    len_avg=$3
    len_std=$4
    identity=$5
    max=$6
    stdev=$7
    error_model=$8
    qscore_model=$9
    output=temp_${10}_${11}.fastq
    log=temp_${10}_${11}.log
    seed=$((42 + ${11}))
    
    badread simulate --reference ${ref} --quantity ${coverage}x --length ${len_avg},${len_std} --identity ${identity},${max},${stdev} --error_model ${error_model} --qscore_model ${qscore_model} --seed ${seed} 2> ${log} > ${output}
}

function usage() {
    echo "multi-thead-badread.sh usage:"
    echo "multi-thead-badread.sh {coverage} {thread} {reference} {length_average} {length_stdev} {identity} {identity_max} {identity_stdev} {error_model} {qscore_model} {output}"
}

main $@


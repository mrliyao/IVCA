#!/bin/bash

# =========== set the test dataset path here ==============
# test_set_dir="./testset/"
test_set_dir="/mnt/c/DatasetFast/"

predictor_path="./vca"
summarizer_path="./sum"
stat_dir="./csvs/"
if [ ! -d "${stat_dir}" ]; then
    mkdir "${stat_dir}"
fi
log_path="./predictor.log"

for file in "${test_set_dir}"/*.{yuv,y4m}; 
do
    seq_name=$(basename -- "${file}")
    seq_name="${seq_name%.*}"

# =========== set the resolution here ==============
    seq_width=3848
    seq_height=2160
    seq_csp=420

    seq_log_path="${stat_dir}/${seq_name}.csv"
    file_temp=$(mktemp)
    # file_temp=./log.log
    "${predictor_path}" --input-res "${seq_width}x${seq_height}" --input-csp "${seq_csp}" --input "${file}" --complexity-csv "${seq_log_path}" > "${file_temp}"
    # "${predictor_path}.exe" --input-res "${seq_width}x${seq_height}" --input-csp "${seq_csp}" --input "${file}" --complexity-csv "${seq_log_path}"
    prediction=$("${summarizer_path}" "${seq_log_path}")
    time=$(cat "${file_temp}" | grep time | awk -F 'time ' '{print $2}')
    echo "${seq_name},${prediction},${time}" >> "${log_path}"
    echo "${seq_name} --- ${prediction} --- ${time}"
done

 
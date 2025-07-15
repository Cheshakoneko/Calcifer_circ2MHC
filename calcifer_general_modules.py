#!/usr/bin/env python

import os
from pathlib import Path
import subprocess

# module for general used functions i.e. file structure and trimming

def file_structure(file_path: str, dataset: str, ending: str) -> Path:
    base_path = Path(file_path)
    dataset_path = base_path / dataset

    if not dataset_path.exists():
        dataset_path.mkdir(parents=True, exist_ok=True)
        print(f"Creating directory {dataset_path} for dataset {dataset}...")

    if ending == "pe":
        fastq_files = [f"{dataset}_1.fastq", f"{dataset}_2.fastq"]
    else:
        fastq_files = [f"{dataset}.fastq"]

    for fq_file in fastq_files:
        orig_fastq = base_path / fq_file
        dest_path = dataset_path / orig_fastq.name
        if orig_fastq.exists() and not dest_path.exists():
            orig_fastq.rename(dest_path)
        elif not orig_fastq.exists():
            print(f"Warning: {orig_fastq.name} does not exist. Skipping...")
        else:
            print(f"Warning: {orig_fastq.name} already exists in {dataset_path}. Skipping...")
    return dataset_path

# single-end read trimming the raw read-data with flexbar #
# flexbar trimming based on a pre trim phred of 20 and a min read len of 50 #
def se_flexbar(cores, dataset, working_dir):
    trim_gz_exist = os.path.isfile(working_dir + 'trim_' + dataset + '.fastq.gz')
    trim_fq_exist = os.path.isfile(working_dir + 'trim_' + dataset + '.fastq')
    if trim_gz_exist:
        trimmed_data = working_dir + 'trim_' + dataset + '.fastq.gz'
        unzip_cmd = 'gunzip ' + trimmed_data
        os.system(unzip_cmd)
        unzip_trimmed_data = trimmed_data[:-3]
        return unzip_trimmed_data
    elif trim_fq_exist:
        unzip_trimmed_data = working_dir + 'trim_' + dataset + '.fastq'
        return unzip_trimmed_data
    else:
        trimmed_data = working_dir + 'trim_' + dataset + '.fastq.gz'
        trimming_cmd = 'flexbar -r ' + working_dir + dataset + '.fastq -t ' + working_dir + 'trim_' + dataset \
                       + ' --zip-output GZ --qtrim-format i1.8 --pre-trim-phred 20 --min-read-length 50 -n ' + cores + ' --output-reads ' + trimmed_data[:-3]
        os.system(trimming_cmd)
        unzip_cmd = 'gunzip ' + trimmed_data
        os.system(unzip_cmd)
        unzip_trimmed_data = trimmed_data[:-3]
        return unzip_trimmed_data


# paired-end read trimming the raw read-data with flexbar #
# trimming based on default values and a min read len of 20 #
def pe_flexbar(cores, dataset, working_dir):
    trim_gz_exist = os.path.isfile(working_dir + 'trim_' + dataset + '_1.fastq.gz')
    trim_fq_exist = os.path.isfile(working_dir + 'trim_' + dataset + '_1.fastq')
    if trim_gz_exist:
        trimmed_data_1 = working_dir + 'trim_' + dataset + '_1.fastq.gz'
        trimmed_data_2 = working_dir + 'trim_' + dataset + '_2.fastq.gz'
        unzip_cmd_1 = 'gunzip ' + trimmed_data_1
        os.system(unzip_cmd_1)
        unzip_trimmed_data_1 = trimmed_data_1[:-3]
        unzip_cmd_2 = 'gunzip ' + trimmed_data_2
        os.system(unzip_cmd_2)
        unzip_trimmed_data_2 = trimmed_data_2[:-3]
        return unzip_trimmed_data_1, unzip_trimmed_data_2
    elif trim_fq_exist:
        unzip_trimmed_data_1 = working_dir + 'trim_' + dataset + '_1.fastq'
        unzip_trimmed_data_2 = working_dir + 'trim_' + dataset + '_2.fastq'
        return unzip_trimmed_data_1, unzip_trimmed_data_2
    else:
        trimmed_data_1 = working_dir + 'trim_' + dataset + '_1.fastq.gz'
        trimmed_data_2 = working_dir + 'trim_' + dataset + '_2.fastq.gz'
        trimming_cmd = 'flexbar -r ' + working_dir + dataset + '_1.fastq -p ' + working_dir + dataset + '_2.fastq '\
                       + '-t ' + working_dir + 'trim_' + dataset \
                       + ' -n ' + cores + ' --zip-output GZ --qtrim-format i1.8 --pre-trim-phred 20 --min-read-length 50 --output-reads '\
                       + trimmed_data_1[:-3] + ' --output-reads2 ' + trimmed_data_2[:-3]
        os.system(trimming_cmd)
        unzip_cmd_1 = 'gunzip ' + trimmed_data_1
        os.system(unzip_cmd_1)
        unzip_trimmed_data_1 = trimmed_data_1[:-3]
        unzip_cmd_2 = 'gunzip ' + trimmed_data_2
        os.system(unzip_cmd_2)
        unzip_trimmed_data_2 = trimmed_data_2[:-3]
        return unzip_trimmed_data_1, unzip_trimmed_data_2


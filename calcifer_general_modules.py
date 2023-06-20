#!/usr/bin/env python

import os


# module for general used functions i.e. file structure and trimming


# creates a directory for each dataset if there is not already one #
def se_file_structure(file_path, dataset):
    dir_exists = os.path.exists(file_path + dataset + '/')
    if dir_exists:
        return file_path + dataset + '/'
    else:
        make_data_dir = 'mkdir ' + file_path + dataset
        os.system(make_data_dir)
        os.rename(file_path + dataset + '.fastq', file_path + dataset + '/' + dataset + '.fastq')
        return file_path + dataset + '/'
     
        
def pe_file_structure(file_path, dataset):
    dir_exists = os.path.exists(file_path + dataset + '/')
    if dir_exists:
        return file_path + dataset + '/'
    else:
        make_data_dir = 'mkdir ' + file_path + dataset
        os.system(make_data_dir)
        os.rename(file_path + dataset + '_1.fastq', file_path + dataset + '/' + dataset + '_1.fastq')
        os.rename(file_path + dataset + '_2.fastq', file_path + dataset + '/' + dataset + '_2.fastq')
        return file_path + dataset + '/'
   
   
# single-end read trimming the raw read-data with flexbar #
# trimming based on default values and a min read len of 20 #
def se_flexbar(dataset, working_dir):
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
                       + ' --zip-output GZ --qtrim-format i1.8 --min-read-length 20 -n 20 --output-reads ' + trimmed_data[:-3]
        os.system(trimming_cmd)
        unzip_cmd = 'gunzip ' + trimmed_data
        os.system(unzip_cmd)
        unzip_trimmed_data = trimmed_data[:-3]
        return unzip_trimmed_data


# paired-end read trimming the raw read-data with flexbar #
# trimming based on default values and a min read len of 20 #
def pe_flexbar(dataset, working_dir):
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
                       + ' -n 20 --zip-output GZ --qtrim-format i1.8 --min-read-length 20 --output-reads '\
                       + trimmed_data_1[:-3] + ' --output-reads2 ' + trimmed_data_2[:-3]
        os.system(trimming_cmd)
        unzip_cmd_1 = 'gunzip ' + trimmed_data_1
        os.system(unzip_cmd_1)
        unzip_trimmed_data_1 = trimmed_data_1[:-3]
        unzip_cmd_2 = 'gunzip ' + trimmed_data_2
        os.system(unzip_cmd_2)
        unzip_trimmed_data_2 = trimmed_data_2[:-3]
        return unzip_trimmed_data_1, unzip_trimmed_data_2


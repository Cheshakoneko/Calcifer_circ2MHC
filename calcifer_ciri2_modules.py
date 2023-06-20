#!/usr/bin/env python

import os


# module with functions for circRNA detection with CIRI2 #
# CIRI2 needs to be installed manually as there is no conda version!!! #


# se bwa mapping for the input for ciri2 #
def se_bwa_mapping(working_dir, unzip_trimmed_data, bwa_index, cores):
    bwa_exist = os.path.isfile(working_dir + 'aln.sam')
    if bwa_exist:
        return
    else:
        bwa_cmd = 'bwa mem -t ' + cores + ' -T 19 ' + bwa_index + ' ' + unzip_trimmed_data + ' > ' + working_dir + \
                  'aln.sam '
        os.system(bwa_cmd)


# pe bwa mapping for the input for ciri2 #
def pe_bwa_mapping(working_dir, unzip_trimmed_data_1, unzip_trimmed_data_2, bwa_index, cores):
    bwa_exist = os.path.isfile(working_dir + 'aln.sam')
    if bwa_exist:
        return
    else:
        bwa_cmd = 'bwa mem -t ' + cores + ' -T 19 ' + bwa_index + ' ' + unzip_trimmed_data_1 + ' ' + \
                  unzip_trimmed_data_2 + ' > ' + working_dir + 'aln.sam'
        os.system(bwa_cmd)


# finding putative circRNAs with ciri2, >= 2 bs-reads is baseline for the tool #
def ciri2_find(ciri_path, working_dir, ref_path, genome_path):
    ciri2_cmd = 'perl ' + ciri_path + ' -I ' + working_dir + 'aln.sam -F ' + genome_path + ' -A ' + ref_path + ' -O ' \
                + working_dir + 'ciri_2_output.txt -G ' + working_dir + 'ciri_logs -T 10'
    os.system(ciri2_cmd)


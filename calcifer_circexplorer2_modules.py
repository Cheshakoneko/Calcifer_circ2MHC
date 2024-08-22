#!/usr/bin/env python

import os


# module with functions for circRNA detection with CE2 #
# CE2 can directly be installed with conda #


# se star-mapping for the chimeric junctions and as input for ce2 #
def se_star_mapping(working_dir, unzip_trimmed_data, star_index, cores):
    chimeric_exist = os.path.isfile(working_dir + 'Chimeric.out.junction')
    if chimeric_exist:
        return
    else:
        star_mapping_cmd = "STAR --runThreadN " + cores + " --runMode alignReads --readFilesIn " + unzip_trimmed_data \
                       + " --genomeDir " + star_index + " --outFileNamePrefix " + working_dir \
                       + " --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate" \
                       + " --alignSJDBoverhangMin 15 --alignSJoverhangMin 15 --chimSegmentMin 15 --chimScoreMin 15" \
                       + "--chimScoreSeparation 10 --chimJunctionOverhangMin 15"
        os.system(star_mapping_cmd)


# pe star-mapping for the chimeric junctions and as input for ce2 #
def pe_star_mapping(working_dir, unzip_trimmed_data_1, unzip_trimmed_data_2, star_index, cores):
    chimeric_exist = os.path.isfile(working_dir + 'Chimeric.out.junction')
    if chimeric_exist:
        return
    else:
        star_mapping_cmd = "STAR --runThreadN " + cores + " --runMode alignReads --readFilesIn " + unzip_trimmed_data_1 + " " \
                       + unzip_trimmed_data_2 + " --genomeDir " + star_index + " --outFileNamePrefix " + working_dir \
                       + " --outFilterMultimapNmax 1 --outFilterMismatchNmax 2 --outSAMtype BAM SortedByCoordinate" \
                       + " --alignSJDBoverhangMin 15 --alignSJoverhangMin 15 --chimSegmentMin 15 --chimScoreMin 15" \
                       + " --chimScoreSeparation 10 --chimJunctionOverhangMin 15"
        os.system(star_mapping_cmd)


# CIRCexplorer2 detection of circRNAs #
# parse all back-spliced-junctions out of the chimeric junctions #
def ce2_parse(working_dir):
    ce2_parse_cmd = "CIRCexplorer2 parse -t STAR " + working_dir + "Chimeric.out.junction -b " + working_dir \
                    + "back_spliced_junction.bed"
    os.system(ce2_parse_cmd)


# get the putative circRNAs from the back-spliced-junctions #
def ce2_annotate(working_dir, ref_file, genome_file):
    ce2_annotate_cmd = "CIRCexplorer2 annotate -r " + ref_file + " -g " + genome_file + " -b " + working_dir \
                       + "back_spliced_junction.bed -o " + working_dir + "circularRNA_full.txt"
    os.system(ce2_annotate_cmd)


# filter ce2 results for at least 2 bs-reads per circRNA #
def ce2_initial_filter(working_dir):
    path_results = working_dir + "circularRNA_full.txt"
    awk_filter_cmd = "awk '$13 > 1' " + path_results + " > " + working_dir + "filtered_circularRNA_full.txt"
    os.system(awk_filter_cmd)


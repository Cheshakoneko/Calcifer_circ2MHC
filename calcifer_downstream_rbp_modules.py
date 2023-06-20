#!/usr/bin/env python

import os
import numpy as np


# module for functions of the rbp binding downstream analysis #


# use fimo on the pseudo circular circRNA sequence #
# filter results for q-val < 0.1 #
def rbp_analysis_circ(working_dir, rbp_db):
    output_dir = working_dir + "all_circs/"
    rbp_file = rbp_db
    fimo_cmd = "fimo -o " + output_dir + "fimo_circ_out/ " + rbp_file + " " + output_dir + "pseudo_circular_seq.fasta"
    fimo_file = output_dir + "fimo_circ_out/fimo.tsv"
    # use fimo to get rbp binding on circ exon seq first
    os.system(fimo_cmd)

    # filter the fimo results for q-val < 0.1
    with open(output_dir + "filtered_fimo_circ_res.txt", "w") as fimo_out:
        with open(fimo_file, "r") as fimo_in:
            for line in fimo_in:
                line_content = line.split()
                if len(line_content) == 10 and line_content[0] != "#":
                    if line_content[0] == "motif_id":
                       fimo_out.write(line)
                    else:
                        if float(line_content[8]) <= 0.1:
                            fimo_out.write(line)

    fimo_res_dict = {}
    with open(output_dir + "filtered_fimo_circ_res.txt", "r") as rbp_in:
        for line in rbp_in:
            line_content = line.split()
            if line_content[2] in fimo_res_dict:
                fimo_res_dict[line_content[2]].append(line_content[1])
            else:
                fimo_res_dict[line_content[2]] = [line_content[1]]

    with open(output_dir + "rbp_analysis_circ_res.tab", "w") as rbp_out:
        for key in fimo_res_dict.keys():
            values, counts = np.unique(fimo_res_dict[key], return_counts=True)
            list_val = list(values)
            list_counts = list(counts)
            rbp_output_string = ""
            for i in range(len(list_val)):
                val = list_val[i]
                count = str(list_counts[i])
                out_count = val + ":" + count + ";"
                rbp_output_string += out_count
            final_rbp_out = key + "\t" + rbp_output_string[:-1] + "\n"
            rbp_out.write(final_rbp_out)


# repeat fimo analysis for sequence around the bsj (and 25 bp into the circRNA on both junction sites #
def rbp_analysis_bsj(working_dir, rbp_db):
    output_dir = working_dir + "all_circs/"
    rbp_file = rbp_db
    fimo_cmd = "fimo -o " + output_dir + "fimo_bsj_out/ " + rbp_file + " " + output_dir + "circ_bsj_seq.fasta"
    fimo_file = output_dir + "fimo_bsj_out/fimo.tsv"
    os.system(fimo_cmd)

    # filter the fimo results for q-val < 0.1
    with open(output_dir + "filtered_fimo_bsj_res.txt", "w") as fimo_out:
        with open(fimo_file, "r") as fimo_in:
            for line in fimo_in:
                line_content = line.split()
                if len(line_content) == 10 and line_content[0] != "#":
                    if line_content[0] == "motif_id":
                       fimo_out.write(line)
                    else:
                        if float(line_content[8]) <= 0.1:
                            fimo_out.write(line)

    fimo_res_dict = {}
    with open(output_dir + "filtered_fimo_bsj_res.txt", "r") as rbp_in:
        for line in rbp_in:
            line_content = line.split()
            if line_content[2] in fimo_res_dict:
                fimo_res_dict[line_content[2]].append(line_content[1])
            else:
                fimo_res_dict[line_content[2]] = [line_content[1]]

    with open(output_dir + "rbp_analysis_bsj_res.tab", "w") as rbp_out:
        for key in fimo_res_dict.keys():
            values, counts = np.unique(fimo_res_dict[key], return_counts=True)
            list_val = list(values)
            list_counts = list(counts)
            rbp_output_string = ""
            for i in range(len(list_val)):
                val = list_val[i]
                count = str(list_counts[i])
                out_count = val + ":" + count + ";"
                rbp_output_string += out_count
            final_rbp_out = key + "\t" + rbp_output_string[:-1] + "\n"
            rbp_out.write(final_rbp_out)


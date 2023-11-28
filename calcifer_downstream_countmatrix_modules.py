#!/usr/bin/env python

import os


# module with function to create count matrixes and coldata for putative DESeq2 analysis #


# generate count matrix for linear and circular RNAs, also generate coldata based on conditions #
# count matrix is also generated if it is a single condition dataset #
# counts can also be used for other analysis beside DGE! #
def deseq2_analysis(working_dir, datasets, conditions, condition_name, read_type, gtf_file, strand):
    col_path = working_dir + "coldata.csv"
    counts_path = working_dir + "count_matrix.csv"
    id_gene_path = working_dir + "all_circs/id_gene_names.tab"
    clr_path = working_dir + "all_percent_circularized.csv"
    dataset_count = 0
    condition_count = 0

    replicate_list = ["replicate"]
    condition_list = ["condition"]

    count_matrix_dict = {}
    gene_count_matrix_dict = {}

    # create header for count_matrix.csv
    count_matrix_header = "identifier"
    for rep in datasets:
        count_matrix_header += "," + rep

    # create output for coldata.csv
    for cons in conditions:
        cons = int(cons)
        for i in range(dataset_count, cons + dataset_count):
            replicate_list.append(datasets[i])
            condition_list.append(condition_name[condition_count])
        condition_count += 1
        dataset_count += cons

    # write out coldata.csv
    with open(working_dir + "coldata.csv", "w") as col_out:
        for i in range(len(replicate_list)):
            if i == len(replicate_list):
                col_out_line = replicate_list[i] + "," + condition_list[i]
            else:
                col_out_line = replicate_list[i] + "," + condition_list[i] + "\n"
            col_out.write(col_out_line)

    # write count data of circ rnas and htseq count of every replicate in dict
    # also write out count data for circRNAs on parental gene level for additional analysis
    rep_count = 0
    for rep in datasets:
        htseq_count_exist = os.path.isfile(working_dir + rep + "/count_" + rep + ".txt")
        # if htseq results exist pass the creation
        if htseq_count_exist:
            if os.stat(working_dir + rep + "/count_" + rep + ".txt").st_size != 0:
                pass
        else:
            indexing_cmd = "samtools index " + working_dir + rep + "/Aligned.sortedByCoord.out.bam"
            os.system(indexing_cmd)
            if read_type == "se":
                htseq_count_cmd = "htseq-count -f bam --stranded=" + strand + " " + working_dir + rep \
                                  + "/Aligned.sortedByCoord.out.bam " + \
                                  gtf_file + " > " + working_dir + rep + "/count_" + rep + ".txt "
                os.system(htseq_count_cmd)
            else:
                htseq_count_cmd = "htseq-count -f bam --stranded=" + strand + " --order=pos " + working_dir + rep \
                                  + "/Aligned.sortedByCoord.out.bam " + gtf_file + " > " \
                                  + working_dir + rep + "/count_" + rep + ".txt "
                os.system(htseq_count_cmd)

        circ_rnas = working_dir + rep + "/filtered_circs.txt"
        htseq_counts = working_dir + rep + "/count_" + rep + ".txt"
        with open(circ_rnas, "r") as circ_count_input:
            for line in circ_count_input:
                line_content = line.split()
                circ_ident = line_content[0]
                circ_count = int(line_content[1])
                gene_circ_ident = line_content[5]
                if circ_ident in count_matrix_dict:
                    count_matrix_dict[circ_ident][rep_count] += circ_count
                else:
                    count_matrix_dict[circ_ident] = [0] * dataset_count
                    count_matrix_dict[circ_ident][rep_count] += circ_count
                if gene_circ_ident in gene_count_matrix_dict:
                    gene_count_matrix_dict[gene_circ_ident][rep_count] += circ_count
                else:
                    gene_count_matrix_dict[gene_circ_ident] = [0] * dataset_count
                    gene_count_matrix_dict[gene_circ_ident][rep_count] += circ_count

        with open(htseq_counts, "r") as htseq_count_input:
            for line in htseq_count_input:
                line_content = line.split()
                gene_ident = line_content[0]
                gene_count = int(line_content[1])
                if gene_ident in count_matrix_dict:
                    count_matrix_dict[gene_ident][rep_count] += gene_count
                    gene_count_matrix_dict[gene_ident][rep_count] += gene_count
                else:
                    # write dict for the full count matrix
                    count_matrix_dict[gene_ident] = [0] * dataset_count
                    count_matrix_dict[gene_ident][rep_count] += gene_count
                    gene_count_matrix_dict[gene_ident] = [0] * dataset_count
                    gene_count_matrix_dict[gene_ident][rep_count] += gene_count

        rep_count += 1

    with open(working_dir + "count_matrix.csv", "w") as count_out:
        count_out.write(count_matrix_header)
        for key in count_matrix_dict.keys():
            full_count_output = "\n" + key
            for count in count_matrix_dict[key]:
                full_count_output += "," + str(count)
            count_out.write(full_count_output)
    
    # matrix with counts per parental gene for circRNAs
    with open(working_dir + "gene_count_matrix.csv", "w") as count_out:
        count_out.write(count_matrix_header)
        for key in gene_count_matrix_dict.keys():
            full_count_output = "\n" + key
            for count in gene_count_matrix_dict[key]:
                full_count_output += "," + str(count)
            count_out.write(full_count_output)


#!/usr/bin/env python

import os


# module with function for the mirna binding analysis #


# miRNA binding analysis with miranda on the pseudo circular sequence #
def mirna_analysis(working_dir, mirna_run):

    output_dir = working_dir + "all_circs/"
    seq_dict = {}
    miranda_cmd = "miranda " + mirna_run + " " + output_dir + "pseudo_circular_seq.fasta -sc 150 -strict -out " + \
                  output_dir + "miranda_circ_res.txt"

    mirna_exist = os.path.isfile(output_dir + "miranda_circ_res.txt")
    if not mirna_exist:
        os.system(miranda_cmd)

    with open(output_dir + "pseudo_circular_seq.fasta", "r") as fasta_in:
        for line in fasta_in:
            seq_key = line.split(">")[1][:-1]
            nextline = next(fasta_in)
            seq_len = len(nextline[:-1])
            seq_dict[seq_key] = [seq_len, 0, {}]

    with open(output_dir + "miranda_circ_res.txt", "r") as miranda_in:
        for line in miranda_in:
            if line[:1] == ">" and line[:2] != ">>":
                search_key = line.split()[1]
                mirna_name = line.split()[0][1:]
                if search_key in seq_dict:
                    seq_dict[search_key][1] += 1
                    if mirna_name in seq_dict[search_key][2]:
                        seq_dict[search_key][2][mirna_name] += 1
                    else:
                        seq_dict[search_key][2][mirna_name] = 1

    with open(output_dir + "analysed_miranda_circ_res.txt", "w") as mir_out:
        for key in seq_dict.keys():
            if seq_dict[key][1] != 0:
                max_mirna = max(seq_dict[key][2], key=seq_dict[key][2].get)
                max_mirna_binding_count = seq_dict[key][2][max_mirna]
                max_mirna_percent_of_all = str(round(float(max_mirna_binding_count) / float(seq_dict[key][1]), 5))
                norm_binding = str(round(float(seq_dict[key][1]) / float(seq_dict[key][0]), 5))
                max_mirna_output = max_mirna + ":" + str(max_mirna_binding_count) + ":" + max_mirna_percent_of_all
                mir_out.write(key + "\t" + str(seq_dict[key][0]) + "\t" + str(seq_dict[key][1]) + "\t" + norm_binding +
                              "\t" + max_mirna_output + "\n")

    mirna_density_list = []
    mirna_dict = {}
    id_gene_dict = {}

    with open(output_dir + "id_gene_names.tab", "r") as ids:
        for line in ids:
            line_content = line.split()
            id_gene_dict[line_content[0]] = line_content[1]

    with open(output_dir + "analysed_miranda_circ_res.txt", "r") as mirna_res:
        for line in mirna_res:
            line_content = line.split()
            mirna_density_list.append(float(line_content[3]))
            if float(line_content[3]) >= 0.1:
                if line_content[0] in id_gene_dict:
                    line_content.append(id_gene_dict[line_content[0]])
                mirna_dict[line_content[0]] = line_content

    mirna_density_list.sort()

    with open(output_dir + "circ_mirna_results.txt", "w") as res_out:
        header = "circRNA\tlength\tbinding_sites\tbs_density\tgene_name"
        res_out.write(header)
        for key in mirna_dict.keys():
            converted_output_list = [str(element) for element in mirna_dict[key]]
            output_string = "\t".join(converted_output_list)
            mirna_results = "\n" + output_string
            res_out.write(mirna_results)


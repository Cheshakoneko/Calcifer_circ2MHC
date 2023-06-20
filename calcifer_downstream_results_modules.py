#!/usr/bin/env python

import os


# modules to collect downstream analyses and create result folder structure #
# also includes a function for the final clean up #

# create data structure for the downstream analysis results #
def data_structure_filtering(working_dir):

    create_out_dir_all = "mkdir " + working_dir + "all_circs/"
    create_out_dir_normoxia = "mkdir " + working_dir + "conditions_circs/"

    dir_all_exists = os.path.exists(working_dir + "all_circs/")
    dir_conditions_exists = os.path.exists(working_dir + "conditions_circs/")

    if not dir_all_exists:
        os.system(create_out_dir_all)

    if not dir_conditions_exists:
        os.system(create_out_dir_normoxia)


# collect all different analyses results and write important end results in the final calcifer output #
def final_output(working_dir, conditions, mirna_run):
    
    # path to all different results
    two_unique_file = working_dir + "all_circs/two_unique_filtered.txt"
    orf_file = working_dir + "all_circs/circ_orfs.tab"
    miranda_file = working_dir + "all_circs/analysed_miranda_circ_res.txt"
    unique_pep_file = working_dir + "all_circs/unique_circ_pep.tab"
    rbp_circ_analysis_file = working_dir + "all_circs/rbp_analysis_circ_res.tab"
    rbp_bsj_analysis_file = working_dir + "all_circs/rbp_analysis_bsj_res.tab"
    path_to_cons = working_dir + "conditions_circs"
    number_conditions = len(conditions)

    result_output = {}
    clr_dict = {}
    
    # read in general circRNA results
    with open(two_unique_file, "r") as circ_in:
        for line in circ_in:
            line_content = line.split()
            circ_id = line_content[0]
            max_unique_bsr = line_content[2]
            circ_origin = line_content[3]
            parental_gene = line_content[7]
            result_output[circ_id] = [parental_gene, circ_origin, max_unique_bsr]
            clr_dict[circ_id] = [0] * number_conditions
    
    # read in ORF results
    orf_dict = {}
    with open(orf_file, "r") as orf_in:
        next(orf_in)
        for line in orf_in:
                line_content = line.split()
                circ_id = line_content[0]
                orf_dict[circ_id] = [line_content[1], line_content[2], line_content[3]]

    if mirna_run != "n":
        miranda_res_dict = {}
        with open(miranda_file, "r") as miranda_in:
            for line in miranda_in:
                line_content = line.split()
                circ_id = line_content[0]
                mirna_bs_density = line_content[3]
                max_mirna_results = line_content[4]
                miranda_res_dict[circ_id] = [mirna_bs_density, max_mirna_results]

    # condition number list for flexible output header
    con_clr_header = ""
    for con in range(number_conditions):
        condition_number = con + 1
        con_clr_header += "\tcon_" + str(condition_number) + "_pc"
        con_circ_file = path_to_cons + "/con_" + str(condition_number) + "_filtered.txt"
        with open(con_circ_file, "r") as con_in:
            for line in con_in:
                line_content = line.split()
                circ_id = line_content[0]
                circ_clr = line_content[8]
                if circ_id in clr_dict:
                    clr_dict[circ_id][con] = float(circ_clr)
    
    # read in unique peptide results
    pep_dict = {}
    with open(unique_pep_file, "r") as pep_in:
        for line in pep_in:
            line_content = line.split()
            if len(line_content) > 1:
                pep_dict[line_content[0]] = line_content[1]
    
    # read in rbp binding results
    rbp_circ_dict = {}
    with open(rbp_circ_analysis_file, "r") as rbp_in:
        for line in rbp_in:
            line_content = line.split()
            rbp_circ_dict[line_content[0]] = line_content[1]

    rbp_bsj_dict = {}
    with open(rbp_bsj_analysis_file, "r") as rbp_in:
        for line in rbp_in:
            line_content = line.split()
            rbp_bsj_dict[line_content[0]] = line_content[1]

    # check and collect results for all circRNA IDs
    for circ_rna in result_output.keys():
        if circ_rna in clr_dict:
            result_output[circ_rna].extend(clr_dict[circ_rna])
    
        if mirna_run != "n":
            if circ_rna in miranda_res_dict:
                result_output[circ_rna].append(miranda_res_dict[circ_rna][0])
                result_output[circ_rna].append(miranda_res_dict[circ_rna][1])
            else:
                result_output[circ_rna].append("no_binding_sites")
                result_output[circ_rna].append("no_max_mirna")

        if circ_rna in rbp_circ_dict:
            result_output[circ_rna].append(rbp_circ_dict[circ_rna])
        else:
            result_output[circ_rna].append("no_rbp_binding")

        if circ_rna in rbp_bsj_dict:
            result_output[circ_rna].append(rbp_bsj_dict[circ_rna])
        else:
            result_output[circ_rna].append("no_rbp_binding")

        if circ_rna in orf_dict:
            result_output[circ_rna].append(orf_dict[circ_rna][0])
            result_output[circ_rna].append(orf_dict[circ_rna][1])
            result_output[circ_rna].append(orf_dict[circ_rna][2])
        else:
            result_output[circ_rna].append(3 * "0:0:NA\t")

        if circ_rna in pep_dict:
            result_output[circ_rna].append(pep_dict[circ_rna] + "\t")
        else:
            result_output[circ_rna].append("non_unique\t")

    # calcifer out write in working dir
    with open(working_dir + "calcifer_output.tab", "w") as calcifer_out:
        file_header = "circID\tparental_gene\ttype\tunique_bsr" + con_clr_header + "\tmirna_binding_site_density\tmost_mirna" \
                      "\trbp_circ_binding\trbp_bsj_binding\tlinear_seq_orf\tpseudo_circular_seq_orf\tmulti_cycle_seq_orf" \
                      "\tunique_region"
        calcifer_out.write(file_header)
        for circ_rna in result_output.keys():
            converted_output_list = [str(element) for element in result_output[circ_rna]]
            output_string = "\t".join(converted_output_list)
            calcifer_results = "\n" + circ_rna + "\t" + output_string
            calcifer_out.write(calcifer_results)
            
 
# remove all not needed files which were created in the analysis #
# move all relevant results into specific folders for better organisation # 
def clean_up(working_dir, mirna_run):
    dir_exists = os.path.exists(working_dir + "all_circs/sub_result_dir/")
    if not dir_exists:
        res_dir = "mkdir " + working_dir + "all_circs/sub_result_dir"
        os.system(res_dir)
    if mirna_run != "n":
        dir_exists = os.path.exists(working_dir + "all_circs/sub_result_dir/mirna_miranda_res_files/")
        if not dir_exists:
            mirna_res_dir = "mkdir " + working_dir + "all_circs/sub_result_dir/mirna_miranda_res_files"
            os.system(mirna_res_dir)
        move_mirna_res_cmd = "mv " + working_dir + "all_circs/*mir* " + working_dir + "all_circs/sub_result_dir" \
                                                                                      "/mirna_miranda_res_files/ "
        os.system(move_mirna_res_cmd)
    
    dir_exists = os.path.exists(working_dir + "all_circs/sub_result_dir/rbp_fimo_res_files/")
    if not dir_exists:
        rbp_res_dir = "mkdir " + working_dir + "all_circs/sub_result_dir/rbp_fimo_res_files"
        os.system(rbp_res_dir)
    else:
        remove_old_res = "rm -r " + working_dir + "all_circs/sub_result_dir/rbp_fimo_res_files/*fimo*"
        os.system(remove_old_res)
    move_rbp_res_cmd_1 = "mv " + working_dir + "all_circs/*fimo* " + working_dir + "all_circs/sub_result_dir" \
                                                                                   "/rbp_fimo_res_files/ "
    move_rbp_res_cmd_2 = "mv " + working_dir + "all_circs/*rbp* " + working_dir + "all_circs/sub_result_dir" \
                                                                                  "/rbp_fimo_res_files/ "
    os.system(move_rbp_res_cmd_1)
    os.system(move_rbp_res_cmd_2)
    
    dir_exists = os.path.exists(working_dir + "all_circs/sub_result_dir/orf_res_files/")
    if not dir_exists:
        orf_res_dir = "mkdir " + working_dir + "all_circs/sub_result_dir/orf_res_files"
        os.system(orf_res_dir)
    move_orf_res_cmd_1 = "mv " + working_dir + "all_circs/*seq* " + working_dir + "all_circs/sub_result_dir" \
                                                                                  "/orf_res_files/ "
    move_orf_res_cmd_2 = "mv " + working_dir + "all_circs/*orf* " + working_dir + "all_circs/sub_result_dir" \
                                                                                  "/orf_res_files/ "
    move_orf_res_cmd_3 = "mv " + working_dir + "all_circs/*prediction* " + working_dir + "all_circs/sub_result_dir" \
                                                                                         "/orf_res_files/ "
    os.system(move_orf_res_cmd_1)
    os.system(move_orf_res_cmd_2)
    os.system(move_orf_res_cmd_3)
    
    dir_exists = os.path.exists(working_dir + "all_circs/sub_result_dir/unique_peptide_res_files/")
    if not dir_exists:
        unique_peptide_res_dir = "mkdir " + working_dir + "all_circs/sub_result_dir/unique_peptide_res_files"
        os.system(unique_peptide_res_dir)
    move_unique_peptide_res_cmd_1= "mv " + working_dir + "all_circs/unique* " + working_dir + "all_circs/sub_result_dir/unique_peptide_res_files/"
    move_unique_peptide_res_cmd_2= "mv " + working_dir + "all_circs/*Unique* " + working_dir + "all_circs/sub_result_dir/unique_peptide_res_files/"
    move_unique_peptide_res_cmd_3= "mv " + working_dir + "all_circs/*.mums " + working_dir + "all_circs/sub_result_dir/unique_peptide_res_files/"
    move_unique_peptide_res_cmd_4= "mv " + working_dir + "all_circs/Query.bed " + working_dir + "all_circs/sub_result_dir/unique_peptide_res_files/"
    os.system(move_unique_peptide_res_cmd_1)
    os.system(move_unique_peptide_res_cmd_2)
    os.system(move_unique_peptide_res_cmd_3)
    os.system(move_unique_peptide_res_cmd_4)


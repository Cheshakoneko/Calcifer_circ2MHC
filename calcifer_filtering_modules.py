#!/usr/bin/env python

from pyensembl import Genome
import os


# module for filtering the raw circRNA results from CE2 and CIRI2


# filtering ce2 and ciri2 results based on chimeric junctions found by star #
# additional compution clr for each circRNA based on bsj and sj read counts #
# uses PyEnsembl to convert ENSG/ENST from CE2 and CIRI2 into gene names #
def chimeric_filtering(working_dir, datasets, genome_fasta, gtf_file):

    path_ce2 = "/filtered_circularRNA_full.txt"
    path_ciri2 = "/ciri_2_output.txt"
    path_chim = "/Chimeric.out.junction"
    list_datasets = datasets

    # create file for clr output for all replicates/datasets
    with open(working_dir + "/all_percent_circularized.csv", "w") as clr_out:
        clr_out.write("circRNA,replicate,clr")

    # create db for getting gene names
    data = Genome(reference_name='ref_annotation', annotation_name='genome_features', gtf_path_or_url=gtf_file)
    # parse GTF and construct database of genomic features
    data.index()

    for dataset in list_datasets:
        ciri2_circ_pos = {}
        ce2_circ_pos = {}
        chim_filter_dict = {}
        filtered_circ = {}
        sj_dict = {}

        # read in sj.out.tab for CLR calculation
        with open(working_dir + dataset + "/SJ.out.tab", "r") as sj_info:
            for line in sj_info:
                sj_pos_1 = line.split()[1]
                sj_pos_2 = line.split()[2]
                sj_read_count = int(line.split()[6])
                if sj_pos_1 in sj_dict:
                    sj_dict[sj_pos_1] += sj_read_count
                else:
                    sj_dict[sj_pos_1] = sj_read_count
                if sj_pos_2 in sj_dict:
                    sj_dict[sj_pos_2] += sj_read_count
                else:
                    sj_dict[sj_pos_2] = sj_read_count

        # read in circexplorer2 results #
        with open(working_dir + dataset + path_ce2, "r") as ce2_input:
            for line in ce2_input:
                line_content = line.split()
                # search linear sj count
                junction_pos_1 = str(line_content[1])
                junction_pos_2 = str(int(line_content[2]) + 1)
                linear_junction_count = 0
                if junction_pos_1 in sj_dict:
                    linear_junction_count += sj_dict[junction_pos_1]
                if junction_pos_2 in sj_dict:
                    linear_junction_count += sj_dict[junction_pos_2]
                ce2_pos_key = str(line_content[0]) + ":" + str(line_content[1]) + "-" + str(
                    line_content[2])
                ce2_circ_pos[ce2_pos_key] = ["exon", line_content[15], ":" + str(line_content[5]),
                                             linear_junction_count]

        # read in ciri2 results #
        with open(working_dir + dataset + path_ciri2, "r") as ciri2_input:
            next(ciri2_input)
            for line in ciri2_input:
                line_content = line.split()
                # search linear sj count
                junction_pos_1 = str(int(line_content[2]) - 1)
                junction_pos_2 = str(int(line_content[3]) + 1)
                linear_junction_count = 0
                if junction_pos_1 in sj_dict:
                    linear_junction_count += sj_dict[junction_pos_1]
                if junction_pos_2 in sj_dict:
                    linear_junction_count += sj_dict[junction_pos_2]
                ciri2_pos_key = str(line_content[1]) + ":" + str(int(line_content[2]) - 1) + "-" + str(
                    line_content[3])
                ciri2_circ_pos[ciri2_pos_key] = [line_content[8], line_content[9][:-1], ":" + str(line_content[10]),
                                                 linear_junction_count]

        # read in chimeric junctions found by STAR for filtering, filter the junctions before #
        # filtering after same chromosome, no compassing junctions, junction lenghts < 100kb #
        # also consider the strand for noting the positions in the directory #
        chim_data = working_dir + dataset + path_chim
        unique_cigars = {}

        with open(chim_data, "r") as chim_in:
            for line in chim_in:
                line_content = line.split()
                junc_len = abs(int(line_content[1]) - int(line_content[4]))

                if line_content[0] == line_content[3] and line_content[2] == line_content[5] and int(
                        line_content[6]) != -1 and junc_len <= 100000:
                    junc_strand = line_content[2]

                    if junc_strand == "+":
                        line_content[1] = str(int(line_content[1]) - 1)
                        chim_key = str(line_content[0]) + ":" + (line_content[4]) + "-" + str(
                            line_content[1])
                    elif junc_strand == "-":
                        line_content[4] = str(int(line_content[4]) - 1)
                        chim_key = str(line_content[0]) + ":" + str(line_content[1]) + "-" + str(
                            line_content[4])

                    # only bs-reads with unique cigar are counted!
                    if chim_key in chim_filter_dict:
                        cigars = [line_content[11], line_content[13]]
                        chim_filter_dict[chim_key][0] += 1
                        if cigars not in unique_cigars[chim_key]:
                            chim_filter_dict[chim_key][1] += 1
                            unique_cigars[chim_key].append(cigars)
                    else:
                        chim_filter_dict[chim_key] = [1, 1]
                        unique_cigars[chim_key] = [[line_content[11], line_content[13]]]

        # filter the chimeric junctions for canonical splice sites  with samtools faidx and genome.fasta #
        gene_fasta = genome_fasta
        chim_splice_site_filtered = {}
        canonical_splice_motifs = ["GT/AG", "GC/AG", "AT/AC", "CT/AC", "CT/GC", "GT/AT"]
        for key in chim_filter_dict.keys():
            chromosome = key.split(":")[0]
            start_pos = key.split(":")[1].split("-")[0]
            end_pos = key.split(":")[1].split("-")[1]
            start_splice_pos = [int(start_pos) - 1, int(start_pos)]
            end_splice_pos = [int(end_pos) + 1, int(end_pos) + 2]
            start_motif_search = '"' + chromosome + ":" + str(start_splice_pos[0]) + "-" + str(
                start_splice_pos[1]) + '"'
            end_motif_search = '"' + chromosome + ":" + str(end_splice_pos[0]) + "-" + str(end_splice_pos[1]) + '"'
            
            if "chr" in chromosome:
                start_mo = os.popen('samtools faidx ' + gene_fasta + ' ' + start_motif_search).read().split()
                end_mo = os.popen('samtools faidx ' + gene_fasta + ' ' + end_motif_search).read().split()
                if len(start_mo) > 1 and len(end_mo) > 1:
                    start_motif = start_mo[1]
                    end_motif = end_mo[1]
                    splice_site_motif = end_motif + "/" + start_motif
                    if splice_site_motif in canonical_splice_motifs:
                        chim_splice_site_filtered[key] = chim_filter_dict[key]

        # filter ciri2 results with remaining chimeric junctions #
        for key in ciri2_circ_pos.keys():
            if key in chim_splice_site_filtered:
                n_key = key + ciri2_circ_pos[key][2]
                # clr computation
                linear_count = float(ciri2_circ_pos[key][3])
                circular_count = float(chim_splice_site_filtered[key][0])
                clr = round(circular_count / (linear_count + circular_count), 10)
                filtered_circ[n_key] = [chim_splice_site_filtered[key][0], chim_splice_site_filtered[key][1],
                                        ciri2_circ_pos[key][0], ciri2_circ_pos[key][1], clr]

        # filter circexplorer2 results with remaining chimeric junctions #
        for key in ce2_circ_pos.keys():
            if key in chim_splice_site_filtered and key not in filtered_circ:
                n_key = key + ce2_circ_pos[key][2]
                # clr computation
                linear_count = float(ce2_circ_pos[key][3])
                circular_count = float(chim_splice_site_filtered[key][0])
                clr = round(circular_count / (linear_count + circular_count), 10)
                filtered_circ[n_key] = [chim_splice_site_filtered[key][0], chim_splice_site_filtered[key][1],
                                        ce2_circ_pos[key][0], ce2_circ_pos[key][1], clr]

        # write out all filtered circRNAs
        with open(working_dir + dataset + "/filtered_circs.txt", "w") as filtered_out:
            for key in filtered_circ.keys():
                # get gene name by ensembl
                id_of_circ = str(filtered_circ[key][3])
                only_id = id_of_circ.split(",")[0]
                id_split = key.split(":")
                chromosome = id_split[0][3:]
                start_pos = int(id_split[1].split("-")[0])
                end_pos = int(id_split[1].split("-")[1])
                gene_start_name = data.gene_names_at_locus(contig=chromosome, position=start_pos)
                gene_end_name = data.gene_names_at_locus(contig=chromosome, position=end_pos)
                gene_names = list(set(gene_start_name + gene_end_name))
                if len(gene_names) == 0:
                    gene_names.append("NA")
                output_name = "/".join(gene_names)
                if len(output_name) == 0:
                    output_name = "NA"
                if output_name[0] == "/":
                    output_name = output_name[1:]
                if output_name[-1] == "/":
                    output_name = output_name[:-1]

                out_string = key + "\t" + str(filtered_circ[key][0]) + "\t" + str(filtered_circ[key][1]) + "\t" \
                                + filtered_circ[key][2] + "\t" + str(filtered_circ[key][3]) + "\t" \
                                + output_name + "\t" + str(filtered_circ[key][4]) + "\n"
                if "chr" in key:
                    filtered_out.write(out_string)

        with open(working_dir + "/all_percent_circularized.csv", "a") as clr_out:
            for key in filtered_circ.keys():
                output = "\n" + key + "," + dataset + "," + str(100 * filtered_circ[key][4])
                if "chr" in key:
                    clr_out.write(output)
                
                
# merging the resulting circRNAs from same condition datasets and all datasets #
def merging_results(working_dir, datasets, conditions, ubsjr_filter):

    all_output = working_dir + "all_circs/"
    conditions_output = working_dir + "conditions_circs/"
    all_datasets = []
    for dataset in datasets:
        all_datasets.append(dataset)
    all_conditions = conditions
    all_filtered = {}

    # merge and write out results for datasets of one condition #
    condition_count = 0
    for con in all_conditions:
        condition_count += 1
        dataset_count = int(con)
        condition_circs = {}
        for i in range(dataset_count):
            actual_dataset = all_datasets.pop(0)
            path_to_dataset = working_dir + actual_dataset
            with open(path_to_dataset + "/filtered_circs.txt", "r") as act_data:
                for line in act_data:
                    line_content = line.split()
                    if int(line_content[1]) > 0:
                        if line_content[0] in all_filtered:
                            all_filtered[line_content[0]][3] += 1
                            all_filtered[line_content[0]][0] += int(line_content[1])
                            all_filtered[line_content[0]][1] += int(line_content[2])
                            if int(line_content[2]) > all_filtered[line_content[0]][4]:
                                all_filtered[line_content[0]][4] = int(line_content[2])
                        else:
                            all_filtered[line_content[0]] = [int(line_content[1]), int(line_content[2]),
                                                             line_content[3], 1, int(line_content[2]), line_content[4],
                                                             line_content[5]]

                    if int(line_content[1]) > 0:
                        if line_content[0] in condition_circs:
                            condition_circs[line_content[0]][3] += 1
                            condition_circs[line_content[0]][0] += int(line_content[1])
                            condition_circs[line_content[0]][1] += int(line_content[2])
                            condition_circs[line_content[0]][7] += float(line_content[6])
                            if int(line_content[2]) > all_filtered[line_content[0]][4]:
                                all_filtered[line_content[0]][4] = int(line_content[2])
                        else:
                            condition_circs[line_content[0]] = [int(line_content[1]), int(line_content[2]),
                                                                line_content[3], 1, int(line_content[2]),
                                                                line_content[4], line_content[5],
                                                                float(line_content[6])]

        with open(conditions_output + "con_" + str(condition_count) + "_filtered.txt", "w") as output:
            for key in condition_circs.keys():
                mean_clr = round((condition_circs[key][7] / dataset_count), 5)
                out_string = key + "\t" + str(condition_circs[key][0]) + "\t" + str(condition_circs[key][1]) + "\t" \
                             + str(condition_circs[key][2]) + "\t" + str(condition_circs[key][3]) + "\t" \
                             + str(condition_circs[key][4]) + "\t" + str(condition_circs[key][5]) + "\t" \
                             + str(condition_circs[key][6]) + "\t" + str(mean_clr) + "\n"
                output.write(out_string)

    # write out circID with gene_name for downstream analysis (DESeq2/GO)
    with open(all_output + "id_gene_names.tab", "w") as id_and_gene_names:
        first_line = "ID\tgene_name"
        id_and_gene_names.write(first_line)

    # merge and write out results for all datasets without regarding conditions #
    with open(all_output + "all_filtered.txt", "w") as out:
        for key in all_filtered.keys():

            # write id and gene_name in tab file for downstream analysis
            with open(all_output + "id_gene_names.tab", "a") as id_and_gene_names:
                add_line = "\n" + key + "\t" + str(all_filtered[key][6])
                id_and_gene_names.write(add_line)

            out_string = key + "\t" + str(all_filtered[key][0]) + "\t" + str(all_filtered[key][1]) + "\t" \
                         + str(all_filtered[key][2]) + "\t" + str(all_filtered[key][3]) + "\t" \
                         + str(all_filtered[key][4]) + "\t" + str(all_filtered[key][5]) + "\t" \
                         + str(all_filtered[key][6]) + "\n"
            out.write(out_string)

    # filter all circRNAs for >= 2 unique bs-reads (CIGAR-based, default) or if possible found in >= 2 datasets #
    with open(all_output + "two_unique_filtered.txt", "w") as out:
        for key in all_filtered.keys():
            if all_filtered[key][4] >= ubsjr_filter:
                out_string = key + "\t" + str(all_filtered[key][0]) + "\t" + str(all_filtered[key][1]) + "\t" \
                         + str(all_filtered[key][2]) + "\t" + str(all_filtered[key][3]) + "\t" \
                         + str(all_filtered[key][4]) + "\t" + str(all_filtered[key][5]) + "\t" \
                         + str(all_filtered[key][6]) + "\n"
                out.write(out_string)



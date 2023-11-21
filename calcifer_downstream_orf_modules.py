#!/usr/bin/env python

import os
from Bio import SeqIO
import pybedtools


# translate the rna sequence into peptide sequence to detect ORFs
def translate(seq):
    protein_code_dict = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                         'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R',
                         'AGG': 'R', 'CTA': 'L', 'CTC': 'L', 'CTG': 'L',
                         'CTT': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P', 'CAC': 'H', 'CAT': 'H', 'CAA': 'Q',
                         'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'GTA': 'V', 'GTC': 'V', 'GTG': 'V',
                         'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'GAC': 'D', 'GAT': 'D', 'GAA': 'E',
                         'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S',
                         'TCT': 'S', 'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L', 'TAC': 'Y', 'TAT': 'Y', 'TAA': '*',
                         'TAG': '*', 'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'}
    protein = ""
    if len(seq) % 3 == 0:
        # Ns in the sequence are rare but can happen, then the peptide sequence can't be retrieved
        if "N" in seq:
            protein = '"NA"'
        else:
            for i in range(0, len(seq), 3):
                codon = seq[i:i + 3]
                protein += protein_code_dict[codon]
        return protein


# fast check if an amino acid is M or * for ORF detection
def find_amino_acid(pep_seq, amino_acid):
    return [i for i, ltr in enumerate(pep_seq) if ltr == amino_acid]


# execute ORF finder function on all three sequence types
def orf_detection(working_dir, min_aa):
    translation_cycles = ["linear_seq", "pseudo_circular_seq", "multi_cycle_seq"]
    for cycle in translation_cycles:
        out_file = working_dir + "all_circs/" + str(cycle) + "_complete_orf.pep"
        circ_seq = working_dir + "all_circs/" + str(cycle) + ".fasta"
        orf_finder(circ_seq, out_file, min_aa)


# detect all putative ORfs on all three reading frames and generate output
def orf_finder(circ_seq, out_file, min_aa):
    with open(circ_seq, "r") as seq_in:
        with open(out_file, "w") as orf_out:
            line_count = 0
            for line in seq_in:
                header = line.replace("\n", "")
                forward_sequence = next(seq_in).replace("\n", "")

                one_frame_forward_seq = forward_sequence
                remaining_bases = len(one_frame_forward_seq) % 3
                if remaining_bases != 0:
                    one_frame_forward_seq = one_frame_forward_seq[:-remaining_bases]
                one_frame_forward_pep = translate(one_frame_forward_seq)

                two_frame_forward_seq = forward_sequence[1:]
                remaining_bases = len(two_frame_forward_seq) % 3
                if remaining_bases != 0:
                    two_frame_forward_seq = two_frame_forward_seq[:-remaining_bases]
                two_frame_forward_pep = translate(two_frame_forward_seq)

                three_frame_forward_seq = forward_sequence[2:]
                remaining_bases = len(three_frame_forward_seq) % 3
                if remaining_bases != 0:
                    three_frame_forward_seq = three_frame_forward_seq[:-remaining_bases]
                three_frame_forward_pep = translate(three_frame_forward_seq)

                pep_sequences = [one_frame_forward_pep, two_frame_forward_pep, three_frame_forward_pep]
                counter = 0
                orf_counter = 1
                for pep_seq in pep_sequences:

                    if counter == 1 or counter == 4:
                        add_base = 1
                    elif counter == 2 or counter == 5:
                        add_base = 2
                    else:
                        add_base = 0

                    complete_orf_pos_list = []
                    complete_orf_pep_list = []
                    start_pos_list = find_amino_acid(pep_seq, "M")
                    end_pos_list = find_amino_acid(pep_seq, "*")
                    if len(start_pos_list) > 0:
                        for start_pos in start_pos_list:
                            orf_pos = 0
                            start_pos = start_pos + 1
                            for end_pos in end_pos_list:
                                if orf_pos == 0:
                                    if end_pos > start_pos and (end_pos - start_pos) < min_aa:
                                        orf_pos = 1
                                    elif (end_pos - start_pos) >= min_aa:
                                        orf_pos = str(start_pos) + ":" + str(end_pos)
                                        nuc_start_pos = (start_pos * 3) - 2 + add_base
                                        nuc_end_pos = (end_pos * 3) + 3 + add_base
                                        complete_orf_pos_list.append(str(nuc_start_pos) + ":" + str(nuc_end_pos))
                                        complete_orf_pep_list.append(pep_seq[start_pos - 1:end_pos + 1])

                    if len(complete_orf_pos_list) != 0:
                        for orf_count in range(0, len(complete_orf_pos_list)):
                            if line_count != 0:
                                orf_pos = complete_orf_pos_list[orf_count]
                                orf_pep = "\n" + complete_orf_pep_list[orf_count]
                                orf_header = "\n>lcl|ORF" + str(orf_counter) + "_" + header[1:] + ":" + orf_pos + \
                                             " complete ORF "
                                orf_counter += 1
                                line_count += 1
                            else:
                                orf_pos = complete_orf_pos_list[orf_count]
                                orf_pep = "\n" + complete_orf_pep_list[orf_count]
                                orf_header = ">lcl|ORF" + str(orf_counter) + "_" + header[1:] + ":" + orf_pos + \
                                             " complete ORF"
                                orf_counter += 1
                                line_count += 1

                            orf_out.write(orf_header)
                            orf_out.write(orf_pep)

                    counter += 1


# merge ORF results for the three different sequence types
def longest_orf_filtering(working_dir):
    translation_cycles = ["linear_seq", "pseudo_circular_seq", "multi_cycle_seq"]
    cycle_count = 0
    circ_orf_dict = {}
    unique_orf_dict = {}
    seq_type_dict = {}
    for cycle in translation_cycles:
        in_file = working_dir + "all_circs/" + str(cycle) + "_complete_orf.pep"
        with open(in_file, "r") as orf_in:
            pep_seq = ""
            for line in orf_in:
                if line[0] == ">":
                    if len(pep_seq) > 0 and "partial" not in header_line:
                        # check if orf for this circRNA is unique or redundant (multi cycle sequence short orfs)
                        if circ_id in unique_orf_dict:
                            if pep_seq not in unique_orf_dict[circ_id]:
                                unique_orf_dict[circ_id].append(pep_seq)
                                # save if orf is from linear, pseudo circular or multi cycle circ rna sequence
                                seq_type_dict[circ_id + pep_seq] = [cycle, orf_position]
                                if circ_id in circ_orf_dict:
                                    circ_orf_dict[circ_id][cycle_count][0] += 1
                                    if len(pep_seq) > circ_orf_dict[circ_id][cycle_count][1]:
                                        circ_orf_dict[circ_id][cycle_count][1] = len(pep_seq)
                                        circ_orf_dict[circ_id][cycle_count][2] = orf_number
                                        circ_orf_dict[circ_id][cycle_count][3] = pep_seq
                                else:
                                    circ_orf_dict[circ_id] = [[0, 0, "NA", ""], [0, 0, "NA", ""], [0, 0, "NA", ""]]
                                    circ_orf_dict[circ_id][cycle_count][0] += 1
                                    circ_orf_dict[circ_id][cycle_count][1] = len(pep_seq)
                                    circ_orf_dict[circ_id][cycle_count][2] = orf_number
                                    circ_orf_dict[circ_id][cycle_count][3] = pep_seq
                        else:
                            unique_orf_dict[circ_id] = [pep_seq]
                            seq_type_dict[circ_id + pep_seq] = [cycle, orf_position]
                            if circ_id in circ_orf_dict:
                                circ_orf_dict[circ_id][cycle_count][0] += 1
                                if len(pep_seq) > circ_orf_dict[circ_id][cycle_count][1]:
                                    circ_orf_dict[circ_id][cycle_count][1] = len(pep_seq)
                                    circ_orf_dict[circ_id][cycle_count][2] = orf_number
                                    circ_orf_dict[circ_id][cycle_count][3] = pep_seq
                            else:
                                circ_orf_dict[circ_id] = [[0, 0, "NA", ""], [0, 0, "NA", ""], [0, 0, "NA", ""]]
                                circ_orf_dict[circ_id][cycle_count][0] += 1
                                circ_orf_dict[circ_id][cycle_count][1] = len(pep_seq)
                                circ_orf_dict[circ_id][cycle_count][2] = orf_number
                                circ_orf_dict[circ_id][cycle_count][3] = pep_seq
                    pep_seq = ""
                    new_header = (line.split("_")[1]).split(":")
                    header_line = line
                    circ_id = new_header[0] + ":" + new_header[1] + ":" + new_header[2]
                    orf_number = line.split("_")[0][5:]
                    orf_position = new_header[3] + "-" + new_header[4].split()[0]
                else:
                    pep_seq += line[:-1]
        cycle_count += 1
    with open(working_dir + "all_circs/circ_orfs.tab", "w") as orf_res_out:
        orf_res_out.write("circID\tlinear_orf\tpseudo_circular_orf\tmulti_cycle_orf")
        for key in circ_orf_dict.keys():
            out_line = "\n" + key + "\t" + str(circ_orf_dict[key][0][0]) + ":" + str(
                circ_orf_dict[key][0][1]) + "\t" + str(circ_orf_dict[key][1][0]) + ":" + str(
                circ_orf_dict[key][1][1]) + "\t" + str(circ_orf_dict[key][2][0]) + ":" + str(circ_orf_dict[key][2][1])
            orf_res_out.write(out_line)

    with open(working_dir + "all_circs/orf_seq.pep", "w") as pep_out:
        line_count = 0
        for key in unique_orf_dict.keys():
            counter = 1
            for orf_seq in unique_orf_dict[key]:
                if line_count == 0:
                    header_line = ">" + key + ".ORF" + str(counter) + "." + seq_type_dict[key + orf_seq][0] + "." + \
                                  seq_type_dict[key + orf_seq][1]
                else:
                    header_line = "\n>" + key + ".ORF" + str(counter) + "." + seq_type_dict[key + orf_seq][0] + "." + \
                                  seq_type_dict[key + orf_seq][1]
                pep_out.write(header_line + "\n" + orf_seq)
                counter += 1
                line_count += 1


# check the sequence in front of each ORF for IRES-like or DRACH motifs
def ires_m6a_prediction(working_dir):
    circ_seq_file = working_dir + "all_circs/multi_cycle_seq.fasta"
    circ_pep_file = working_dir + "all_circs/orf_seq.pep"

    circ_seq_dict = {}
    circ_orf_dict = {}

    # create dict with all multi cycle nucleotide sequences for each circRNA
    with open(circ_seq_file, "r") as seq_in:
        for line in seq_in:
            circ_id = line[:-1]
            seq = next(seq_in)[:-1]
            circ_seq_dict[circ_id] = seq

    # create dict with all found complete ORFs from linear, pseudo circular and multi cycle seq
    with open(circ_pep_file, "r") as orf_in:
        for line in orf_in:
            if line[0] == ">":
                line_content = line.split(".")
                circ_id = line_content[0]
                orf_pos = line_content[3][:-1]
                circ_orf_dict[line[:-1]] = [circ_id, orf_pos]

    front_orf_seq_dict = {}

    # get sequence of x nucleotides in front of the start codon
    # If the ORF is near the BSJ the remaining nt are taken from the end of the sequence
    for key in circ_orf_dict.keys():
        # pseudo circular sequence contains linear sequence with additional 25 nt from the end at the start
        # these 25 nt are substracted so the multi cycle sequence can also be used for it
        if "pseudo" not in key:

            end_pos_seq = int(circ_orf_dict[key][1].split("-")[0]) - 1
            start_pos_seq = end_pos_seq - 10
            circ_seq = circ_seq_dict[circ_orf_dict[key][0]]
            # depending on the start and end position the string slicing is different
            # as we always need x nt in front of a given ORF
            # so the sequence can also be ranging over the bsj
            if start_pos_seq < 0:
                front_orf_seq_1 = circ_seq[start_pos_seq:]
                front_orf_seq_2 = circ_seq[:end_pos_seq]
                front_orf_seq = front_orf_seq_1 + front_orf_seq_2
            else:
                front_orf_seq = circ_seq[start_pos_seq:end_pos_seq]
        else:
            end_pos_seq = int(circ_orf_dict[key][1].split("-")[0]) - 26
            start_pos_seq = end_pos_seq - 10
            circ_seq = circ_seq_dict[circ_orf_dict[key][0]]
            if start_pos_seq < 0 and end_pos_seq >= 0:
                front_orf_seq_1 = circ_seq[start_pos_seq:]
                front_orf_seq_2 = circ_seq[:end_pos_seq]
                front_orf_seq = front_orf_seq_1 + front_orf_seq_2
            elif start_pos_seq < 0 and end_pos_seq < 0:
                front_orf_seq = circ_seq[start_pos_seq:end_pos_seq]
            else:
                front_orf_seq = circ_seq[start_pos_seq:end_pos_seq]

        front_orf_seq_dict[key] = front_orf_seq

    # ires_like and drach motifs which should be searched in front of all circRNA ORFs
    ires_like_motif_list = ['AATATA', 'AAAATA', 'AAGATA', 'AACATA', 'AACATT', 'AAACAT', 'ACATAA', 'GAGATA', 'GGAGAT',
                            'TGACAT', 'GACATA', 'AAATAT', 'ATATAT', 'ACATAT', 'AGATAT', 'AATATC', 'ATATCT', 'TAATAT',
                            'TAATCT', 'AAAAAT', 'AAAATT', 'AAATAC', 'AAATTC', 'AAATCC', 'AATAAA', 'TCAAGC', 'ATCAAG',
                            'AATCAA', 'ATAAAG', 'ATAAAT', 'AAAAAA', 'CAAAAA', 'ACAAAA', 'TAAAAA', 'ATAAAA', 'ATAAAC',
                            'CGAAAC', 'AATACA', 'ATACAA', 'AAACAA', 'TAAACA', 'TATACA', 'ATATAA', 'CATATA', 'ATATAG',
                            'ATATAC', 'TATATA', 'TATATT', 'TATTTT', 'TATATG', 'TTATAT', 'TATAAA', 'ATTTAA', 'TTTAAA',
                            'AAATTA', 'AATTAT', 'TAATTA', 'AATTAA', 'AATTTA', 'AATTCA', 'AACTGA', 'ATATTA', 'TATTAA',
                            'ATTAAT', 'ATTATT', 'TAGATT', 'AGATTA', 'ATTAGG', 'CATTAG', 'ATTCGA', 'AATAGA', 'AAATAA',
                            'ATAAGA', 'AAAAGA', 'AAAGAC', 'TAAGAC', 'AGAAGA', 'GAAGAA', 'AAGAAG', 'AAGAAT', 'TAAGAA',
                            'AATAAG', 'ATAAGT', 'AATAAT', 'AATATT', 'AATACT', 'TATACT', 'ATACTG', 'ATACTA', 'ATAATA',
                            'TATAAT', 'TTATAA', 'ATTATA', 'TAATAA', 'TAAATA', 'TTAATA', 'TGAATA']

    drach_motif_list = ['GGACA', 'GGACC', 'GGACT', 'GAACA', 'GAACC', 'GAACT', 'AGACA', 'AGACC', 'AGACT', 'AAACA',
                        'AAACC', 'AAACT', 'TGACA', 'TGACC', 'TGACT', 'TAACA', 'TAACC', 'TAACT']

    pentamer_dict = {}
    hexamer_dict = {}

    # get all pentamer and hexamer sequences and their position included in any given circRNA ORF
    # scan these for ires_like and drach motifs
    # later also include identifier for each occurancy of a given k-mer
    # so ires_like and drach motifs can be assigned to specific circRNA ORFs
    for key in front_orf_seq_dict.keys():
        fasta_seq = front_orf_seq_dict[key]
        key_content = key.split(">")[1].split(".")
        output_key = key_content[0] + "." + key_content[1]

        for i in range(0, len(fasta_seq) - 4):
            pentamer_seq = fasta_seq[i:i + 5]
            if pentamer_seq in pentamer_dict:
                pentamer_dict[pentamer_seq].append(output_key + "|" + str(i) + ":" + str(i + 5))
            else:
                pentamer_dict[pentamer_seq] = [output_key + "|" + str(i) + ":" + str(i + 5)]

        for i in range(0, len(fasta_seq) - 5):
            hexamer_seq = fasta_seq[i:i + 6]
            if hexamer_seq in hexamer_dict:
                hexamer_dict[hexamer_seq].append(output_key + "|" + str(i) + ":" + str(i + 6))
            else:
                hexamer_dict[hexamer_seq] = [output_key + "|" + str(i) + ":" + str(i + 6)]

    ires_m6a_dict = {}

    with open(working_dir + "all_circs/ires_like_binding_prediction.tsv", "w") as ires_out:
        ires_out.write("circID\tORF_number\tires_like_binding_position [10 nt in front of start codon]")
        for ires in ires_like_motif_list:
            if ires in hexamer_dict:
                for pos in hexamer_dict[ires]:
                    circ_id = pos.split(".")[0]
                    orf_number = pos.split(".")[1].split("|")[0]
                    circ_orf_id = circ_id + "." + orf_number
                    ires_pos = pos.split("|")[1]
                    ires_out.write("\n" + circ_id + "\t" + orf_number + "\t" + ires_pos)
                    if circ_orf_id in ires_m6a_dict:
                        ires_m6a_dict[circ_orf_id][0] += 1
                    else:
                        ires_m6a_dict[circ_orf_id] = [1, 0]

    with open(working_dir + "all_circs/drach_site_prediction.tsv", "w") as drach_out:
        drach_out.write("circID\tORF_number\tdrach_site_position [10 nt in front of start codon]")
        for drach in drach_motif_list:
            if drach in pentamer_dict:
                for pos in pentamer_dict[drach]:
                    circ_id = pos.split(".")[0]
                    orf_number = pos.split(".")[1].split("|")[0]
                    circ_orf_id = circ_id + "." + orf_number
                    drach_pos = pos.split("|")[1]
                    drach_out.write("\n" + circ_id + "\t" + orf_number + "\t" + drach_pos)
                    if circ_orf_id in ires_m6a_dict:
                        ires_m6a_dict[circ_orf_id][1] += 1
                    else:
                        ires_m6a_dict[circ_orf_id] = [0, 1]

    return ires_m6a_dict


# detect unique peptide subsequences on putative translated circRNA peptides
# (original algorithm implemented by Justin Murtagh from schulz lab)
def unique_peptides_analysis(working_dir, pep_ref, ires_m6a_dict):
    file_path = working_dir + "all_circs/"
    db_file = pep_ref
    query_file = file_path + "orf_seq.pep"
    matchlength = 10
    mmumer_res = file_path + "Maxmatch.mums"
    non_un = file_path + "Non_Unique_Regions.bed"
    quer = file_path + "Query.bed"
    unique = file_path + "Unique_Regions.bed"
    out_name = file_path + "Unique_Regions_merged.bed"
    mmumer_cmd = "mummer -maxmatch -l " + str(matchlength) + " " + db_file + " " + query_file + " > " + mmumer_res

    os.system(mmumer_cmd)

    infile = open(mmumer_res, "r")
    lines = infile.readlines()

    with open(non_un, "w") as out1:
        name = ""
        for i in lines:
            if i[0] == ">":
                name = i[2:-1]
            else:
                columns = i.split()
                start = int(columns[2]) - 1
                end = start + int(columns[3])
                out1.write(name + "\t" + str(start) + "\t" + str(end) + "\n")

    reference = SeqIO.parse(query_file, "fasta")
    with open(quer, "w") as out2:
        for line in reference:
            out2.write(line.id + "\t" + "0" + "\t" + str(len(str(line.seq))) + "\n")

    a = pybedtools.BedTool(non_un)
    b = pybedtools.BedTool(quer)
    b.subtract(a).saveas(unique)

    infile = open(unique, "r")
    lines = infile.readlines()
    gap = int(matchlength)
    previous = 3 * [""]
    with open(out_name, "w") as out:
        for i in lines:
            columns = i.split()
            if columns[0] == previous[0]:
                if int(columns[1]) - int(previous[2]) <= gap:
                    previous[2] = columns[2]
                else:
                    if int(previous[2]) - int(previous[1]) >= gap:
                        out.write(previous[0] + "\t" + previous[1] + "\t" + previous[2] + "\n")
                    previous = columns
            else:
                if previous[0] != "":
                    if int(previous[2]) - int(previous[1]) >= gap:
                        out.write(previous[0] + "\t" + previous[1] + "\t" + previous[2] + "\n")
                previous = columns
        if int(previous[2]) - int(previous[1]) >= gap:
            out.write(previous[0] + "\t" + previous[1] + "\t" + previous[2] + "\n")

    unique_regions = {}

    with open(out_name, "r") as unique_in:
        for line in unique_in:
            line_content = line.split()
            circ_id = line_content[0].split(".")[0]
            orf_number = line_content[0].split(".")[1]

            # add ires_like and m6a site prediction results to the unique peptides
            circ_orf_id = circ_id + "." + orf_number
            if circ_orf_id in ires_m6a_dict:
                ires_m6a = ires_m6a_dict[circ_orf_id]
                ires_m6a_res = "ireslike_" + str(ires_m6a[0]) + ":m6a_" + str(ires_m6a[1])
            else:
                ires_m6a_res = "none"

            starting_pos = line_content[1]
            end_pos = line_content[2]
            unique_region = orf_number + ":" + starting_pos + "-" + end_pos + ":" + ires_m6a_res
            # only unique peptides with at least one ires_like or drach motif 10 nt in front
            # of the start codon are considered
            if circ_id in unique_regions and ires_m6a_res != "none":
                unique_regions[circ_id] += "|" + unique_region
            else:
                if ires_m6a_res != "none":
                    unique_regions[circ_id] = unique_region

    with open(file_path + "unique_circ_pep.tab", "w") as unique_out:
        for key in unique_regions.keys():
            circ_key = key
            unique_region = unique_regions[key]
            unique_out.write(circ_key + "\t" + unique_region + "\n")

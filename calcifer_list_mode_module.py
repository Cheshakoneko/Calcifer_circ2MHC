#!/usr/bin/env python

import os
import random
from pyfaidx import Fasta
from pyfaidx import Faidx
from pyensembl import Genome
import numpy as np
from Bio import SeqIO
import pybedtools


# get parental gene name for each circRNA identifier
def parental_gene_name(working_dir, crna_list_file, gtf_file):
    data = Genome(reference_name='ref_annotation', annotation_name='genome_features', gtf_path_or_url=gtf_file)
    # parse GTF and construct database of genomic features
    data.index()
    with open(working_dir + "circrna_name_list.tsv", "w") as crna_out:
        header_line = "parental_gene\tcircRNA"
        crna_out.write(header_line)
        with open(working_dir + crna_list_file, "r") as crna_in:
            for line in crna_in:
                line_content = line.split(":")
                chromosome = line_content[0][3:]
                start_pos = int(line_content[1].split("-")[0])
                end_pos = int(line_content[1].split("-")[1])
                strand = line_content[2].split("\n")[0]
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
                crna_out.write("\n" + output_name + "\t" + "chr" + str(chromosome) + ":" + str(start_pos) + "-"
                               + str(end_pos) + ":" + strand)


# get dict with the annotation for all exons, cds and 3'utrs
def mirna_annotation(gtf_file):
    cds_annotation = {}
    three_utr_annotation = {}
    exon_annotation = {}
    exon_endings = {}
    # transform generic identifier to chromosome names
    with open(gtf_file, "r") as anno:
        for line in anno:
            if line[0] != "#":
                line_content = line.split()

                if line_content[0] == "NC_000023.11":
                    line_content[0] = "chrX"
                elif line_content[0] == "NC_000024.10":
                    line_content[0] = "chrY"
                elif line_content[0] == "NC_000022.11":
                    line_content[0] = "chr22"
                elif line_content[0] == "NC_000021.9":
                    line_content[0] = "chr21"
                elif line_content[0] == "NC_000020.11":
                    line_content[0] = "chr20"
                elif line_content[0] == "NC_000019.10":
                    line_content[0] = "chr19"
                elif line_content[0] == "NC_000018.10":
                    line_content[0] = "chr18"
                elif line_content[0] == "NC_000017.11":
                    line_content[0] = "chr17"
                elif line_content[0] == "NC_000016.10":
                    line_content[0] = "chr16"
                elif line_content[0] == "NC_000015.10":
                    line_content[0] = "chr15"
                elif line_content[0] == "NC_000014.9":
                    line_content[0] = "chr14"
                elif line_content[0] == "NC_000013.11":
                    line_content[0] = "chr13"
                elif line_content[0] == "NC_000012.12":
                    line_content[0] = "chr12"
                elif line_content[0] == "NC_000011.10":
                    line_content[0] = "chr11"
                elif line_content[0] == "NC_000010.11":
                    line_content[0] = "chr10"
                elif line_content[0] == "NC_000009.12":
                    line_content[0] = "chr9"
                elif line_content[0] == "NC_000008.11":
                    line_content[0] = "chr8"
                elif line_content[0] == "NC_000007.14":
                    line_content[0] = "chr7"
                elif line_content[0] == "NC_000006.12":
                    line_content[0] = "chr6"
                elif line_content[0] == "NC_000005.10":
                    line_content[0] = "chr5"
                elif line_content[0] == "NC_000004.12":
                    line_content[0] = "chr4"
                elif line_content[0] == "NC_000003.12":
                    line_content[0] = "chr3"
                elif line_content[0] == "NC_000002.12":
                    line_content[0] = "chr2"
                elif line_content[0] == "NC_000001.11":
                    line_content[0] = "chr1"
                elif "chr" not in line_content[0]:
                    line_content[0] = "chr" + line_content[0]

                information_type = line_content[2]
                if int(line_content[3]) > int(line_content[4]):
                    feature_start = line_content[4]
                    feature_end = line_content[3]
                else:
                    feature_start = line_content[3]
                    feature_end = line_content[4]
                key_pos = line_content[0] + ":" + str(feature_start) + "-" + str(feature_end)
                # get all annotated cds
                if information_type == "CDS":
                    cds_annotation[key_pos] = '"' + key_pos + '"'
                # get all annotated 3'utr
                elif information_type == "three_prime_utr":
                    three_utr_annotation[key_pos] = '"' + key_pos + '"'

                # get position of all exons
                elif information_type == "exon":
                    if line_content[0] + ":" + str(feature_start) in exon_annotation:
                        if line_content[0] + ":" + str(feature_end) not in exon_annotation[line_content[0] + ":" +
                                                                                           str(feature_start)]:
                            exon_annotation[line_content[0] + ":" + str(feature_start)].append(line_content[0] + ":" +
                                                                                               str(feature_end))
                    else:
                        exon_annotation[line_content[0] + ":" + str(feature_start)] = [line_content[0] + ":" +
                                                                                       str(feature_end)]
                    exon_endings[line_content[0] + ":" + str(feature_end)] = str(feature_start)

    return cds_annotation, three_utr_annotation, exon_annotation, exon_endings


# get the circRNA sequences
def circ_exon_seq(working_dir, gene_fasta, exon_anno, exon_endings):
    genes = Faidx(gene_fasta)
    circ_file = working_dir + "circrna_name_list.tsv"
    circ_rnas = {}
    # keys are start or end position, every start/end can have multiple circRNAs (isoforms) #
    # read in circRNAs with identifier as key #
    with open(circ_file, "r") as circs_in:
        next(circs_in)
        with open(working_dir + "circ_bsj_seq.fasta", "w") as bsj_seq:
            for line in circs_in:
                line_content = line.split()
                all_pos = line_content[0].split(":")
                chromosome = all_pos[0]
                strand = line_content[0].split(":")[2]
                positions = all_pos[1].split("-")
                positions[0] = str(int(positions[0]) + 1)
                circ_rnas[line_content[0]] = [positions[0], positions[1], line_content[1:]]
                # get sequence around bsj +/- 250bp
                header = "\n>" + line_content[0] + "\n"
                up_bsj = chromosome + ":" + str(int(positions[0]) - 250) + "-" + str(int(positions[0]) + 24)
                down_bsj = chromosome + ":" + str(int(positions[1]) - 24) + "-" + str(int(positions[1]) + 250)
                if strand == "+":
                    fasta_seq_up = str(genes.fetch(chromosome, (int(positions[0]) - 250), (int(positions[0]) + 24)))
                    full_fasta_seq_up = ""
                    full_fasta_seq_up += ''.join(fasta_seq_up)
                    fasta_seq_down = str(genes.fetch(chromosome, (int(positions[1]) - 24), (int(positions[1]) + 250)))
                    full_fasta_seq_down = ""
                    full_fasta_seq_down += ''.join(fasta_seq_down)
                # -i option does not work maybe implement something else!
                elif strand == "-":
                    fasta_seq_up = str(genes.fetch(chromosome, (int(positions[0]) - 250), (int(positions[0]) + 24)).reverse.complement)
                    full_fasta_seq_up = ""
                    full_fasta_seq_up += ''.join(fasta_seq_up)
                    fasta_seq_down = str(genes.fetch(chromosome, (int(positions[1]) - 24), (int(positions[1]) + 250)).reverse.complement)
                    full_fasta_seq_down = ""
                    full_fasta_seq_down += ''.join(fasta_seq_down)
                bsj_seq.write(header)
                bsj_seq.write(full_fasta_seq_up)
                bsj_seq.write(header)
                bsj_seq.write(full_fasta_seq_down)

    circ_rna_seq_pos = {}
    # set nonsense value for iteration failures
    end_point = -1
    start_exon = -1
    pos = -1
    for key in circ_rnas.keys():
        circ_rna_starts = []
        circ_rna_exons = {}
        start = int(circ_rnas[key][0])
        end = int(circ_rnas[key][1])
        chrom = key.split(":")[0]
        for pos in range(start, end):
            # check for partial sequence of exon at circRNA start
            if len(circ_rna_exons) == 0:
                ident_pos_end = chrom + ":" + str(pos)
                if ident_pos_end in exon_endings:
                    real_exon_start = exon_endings[ident_pos_end]
                    if int(real_exon_start) < pos:
                        start_exon = start
                        end_point = pos
                        circ_rna_starts.append(int(start_exon))
                        circ_rna_exons[start_exon] = [int(end_point), chrom + ":" + str(start_exon) + "-" +
                                                      str(end_point)]
            end_point = -1
            ident_pos = chrom + ":" + str(pos)
            # check for full exons in circRNA
            if ident_pos in exon_anno:
                start_exon = pos
                end_point = 0
                exon_ending_list = exon_anno[ident_pos]
                # only take largest version of a given exon (i.e. largest end point inside a circ)
                for i in exon_ending_list:
                    end_pos = int(i.split(":")[1])
                    if end >= end_pos > end_point:
                        end_point = end_pos
                if end_point > pos:
                    circ_rna_starts.append(int(pos))
                    circ_rna_exons[pos] = [int(end_point), chrom + ":" + str(pos) + "-" + str(end_point)]
        # add partial exon sequence at the end of circRNA
        if end_point == 0:
            circ_rna_starts.append(int(start_exon))
            circ_rna_exons[start_exon] = [int(pos), chrom + ":" + str(start_exon) + "-" + str(pos)]
        if len(circ_rna_starts) == 0:
            circ_rna_starts.append(int(start))
            circ_rna_exons[start] = [int(end), chrom + ":" + str(start) + "-" + str(end)]

        circ_rna_seq_pos[key] = []
        sorted_exon_starts = sorted(circ_rna_starts)
        prev_exon_end = 0
        for exon_start in sorted_exon_starts:
            if exon_start > prev_exon_end:
                prev_exon_end = circ_rna_exons[exon_start][0]
                circ_rna_seq_pos[key].append(circ_rna_exons[exon_start][1])

    with open(working_dir + "linear_seq.fasta", "w") as circ_seq_out:
        for key in circ_rna_seq_pos.keys():
            if len(circ_rna_seq_pos[key]) > 0:
                fasta_header = ">" + key
                full_fasta_seq = ""
                strand = key.split(":")[2]
                if strand == "+":
                    for exon_pos in circ_rna_seq_pos[key]:
                        chromosome = exon_pos.split(":")[0]
                        position = exon_pos.split(":")[1]
                        start_pos = int(position.split("-")[0])
                        end_pos = int(position.split("-")[1])
                        fasta_seq = str(genes.fetch(chromosome, start_pos, end_pos))
                        full_fasta_seq += ''.join(fasta_seq)
                elif strand == "-":
                    for exon_pos in circ_rna_seq_pos[key]:
                        chromosome = exon_pos.split(":")[0]
                        position = exon_pos.split(":")[1]
                        start_pos = int(position.split("-")[0])
                        end_pos = int(position.split("-")[1])
                        fasta_seq = str(genes.fetch(chromosome, start_pos, end_pos).reverse.complement)
                        full_fasta_seq = ''.join((fasta_seq, full_fasta_seq))
                circ_seq_out.write("\n" + fasta_header + "\n")
                circ_seq_out.write(full_fasta_seq.upper())
    os.system("sed -i \'1d\' " + working_dir + "linear_seq.fasta")

    with open(working_dir + "linear_seq.fasta", "r") as circ:
        with open(working_dir + "pseudo_circular_seq.fasta", "w") as seq_out:
            line_count = 0
            for line1 in circ:
                header = line1[:-1]
                line2 = next(circ)
                circ_line = line2[:-1]
                seq = circ_line[-25:] + circ_line + circ_line[:25]
                if line_count != 0:
                    seq_out.write("\n" + header + "\n" + seq)
                else:
                    seq_out.write(header + "\n" + seq)
                    line_count += 1

    with open(working_dir + "linear_seq.fasta", "r") as circ:
        with open(working_dir + "multi_cycle_seq.fasta", "w") as seq_out:
            line_count = 0
            for line1 in circ:
                header = line1[:-1]
                line2 = next(circ)
                seq = 4 * line2[:-1]
                if line_count != 0:
                    seq_out.write("\n" + header + "\n" + seq)
                else:
                    seq_out.write(header + "\n" + seq)
                    line_count += 1


# miRNA binding analysis with miranda on the pseudo circular sequence #
def mirna_analysis(working_dir, mirna_run):

    output_dir = working_dir
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


# use fimo on the pseudo circular circRNA sequence #
# filter results for q-val < 0.1 (default) #
def rbp_analysis_circ(working_dir, rbp_db, qval):
    output_dir = working_dir
    rbp_file = rbp_db
    fimo_cmd = "fimo -o " + output_dir + "fimo_circ_out/ " + rbp_file + " " + output_dir + "pseudo_circular_seq.fasta"
    fimo_file = output_dir + "fimo_circ_out/fimo.tsv"
    # use fimo to get rbp binding on circ exon seq first
    os.system(fimo_cmd)

    # filter the fimo results for q-val < 0.1 (default)
    with open(output_dir + "filtered_fimo_circ_res.txt", "w") as fimo_out:
        with open(fimo_file, "r") as fimo_in:
            for line in fimo_in:
                line_content = line.split()
                if len(line_content) == 10 and line_content[0] != "#":
                    if line_content[0] == "motif_id":
                        fimo_out.write(line)
                    else:
                        if float(line_content[8]) <= qval:
                            fimo_out.write(line)

    fimo_res_dict = {}
    with open(output_dir + "filtered_fimo_circ_res.txt", "r") as rbp_in:
        for line in rbp_in:
            line_content = line.split()
            if line_content[1] in fimo_res_dict:
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
def rbp_analysis_bsj(working_dir, rbp_db, qval):
    output_dir = working_dir
    rbp_file = rbp_db
    fimo_cmd = "fimo -o " + output_dir + "fimo_bsj_out/ " + rbp_file + " " + output_dir + "circ_bsj_seq.fasta"
    fimo_file = output_dir + "fimo_bsj_out/fimo.tsv"
    os.system(fimo_cmd)

    # filter the fimo results for q-val < 0.1 (default)
    with open(output_dir + "filtered_fimo_bsj_res.txt", "w") as fimo_out:
        with open(fimo_file, "r") as fimo_in:
            for line in fimo_in:
                line_content = line.split()
                if len(line_content) == 10 and line_content[0] != "#":
                    if line_content[0] == "motif_id":
                        fimo_out.write(line)
                    else:
                        if float(line_content[8]) <= qval:
                            fimo_out.write(line)

    fimo_res_dict = {}
    with open(output_dir + "filtered_fimo_bsj_res.txt", "r") as rbp_in:
        for line in rbp_in:
            line_content = line.split()
            if line_content[1] in fimo_res_dict:
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


def translate(seq):
    protein_code_dict = {'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M', 'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
                         'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K', 'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
                         'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
                         'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
                         'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V', 'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
                         'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E', 'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
                         'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S', 'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
                         'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*', 'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'}
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += protein_code_dict[codon]
        return protein


def find_amino_acid(pep_seq, amino_acid):
    return [i for i, ltr in enumerate(pep_seq) if ltr == amino_acid]


def orf_detection(working_dir, min_aa):
    translation_cycles = ["linear_seq", "pseudo_circular_seq", "multi_cycle_seq"]
    for cycle in translation_cycles:
        out_file = working_dir + str(cycle) + "_complete_orf.pep"
        circ_seq = working_dir + str(cycle) + ".fasta"
        orf_finder(circ_seq, out_file, min_aa)


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


def longest_orf_filtering(working_dir):
    translation_cycles = ["linear_seq", "pseudo_circular_seq", "multi_cycle_seq"]
    cycle_count = 0
    circ_orf_dict = {}
    unique_orf_dict = {}
    seq_type_dict = {}
    for cycle in translation_cycles:
        in_file = working_dir + str(cycle) + "_complete_orf.pep"
        with open(in_file, "r") as orf_in:
            pep_seq = ""
            for line in orf_in:
                if line[0] == ">":
                    if len(pep_seq) > 0 and "partial" not in header_line:
                        # check if orf for this circRNA is unique or redundant (multi cycle squence short orfs)
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
    with open(working_dir + "circ_orfs.tab", "w") as orf_res_out:
        orf_res_out.write("circID\tlinear_orf\tpseudo_circular_orf\tmulti_cycle_orf")
        for key in circ_orf_dict.keys():
            out_line = "\n" + key + "\t" + str(circ_orf_dict[key][0][0]) + ":" + str(
                circ_orf_dict[key][0][1]) + "\t" + str(circ_orf_dict[key][1][0]) + ":" + str(
                circ_orf_dict[key][1][1]) + "\t" + str(circ_orf_dict[key][2][0]) + ":" + str(circ_orf_dict[key][2][1])
            orf_res_out.write(out_line)

    with open(working_dir + "orf_seq.pep", "w") as pep_out:
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


def ires_m6a_prediction(working_dir):
    circ_seq_file = working_dir + "multi_cycle_seq.fasta"
    circ_pep_file = working_dir + "orf_seq.pep"

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

    with open(working_dir + "ires_like_binding_prediction.tsv", "w") as ires_out:
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

    with open(working_dir + "drach_site_prediction.tsv", "w") as drach_out:
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


def unique_peptides_analysis(working_dir, pep_ref, ires_m6a_dict):
    # cite justin murtagh from schulz lab
    file_path = working_dir
    db_file = pep_ref
    # "/home/andre/mouse_heart_data/EC_data/gencode.vM25.pc_translations.fa"
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


# collect all different analyses results and write important end results in the final calcifer output #
def final_output(working_dir):
    # path to all different results
    circ_file = working_dir + "circrna_name_list.tsv"
    orf_file = working_dir + "circ_orfs.tab"
    miranda_file = working_dir + "analysed_miranda_circ_res.txt"
    unique_pep_file = working_dir + "unique_circ_pep.tab"
    rbp_circ_analysis_file = working_dir + "rbp_analysis_circ_res.tab"
    rbp_bsj_analysis_file = working_dir + "rbp_analysis_bsj_res.tab"

    result_output = {}

    # read in general circRNA results
    with open(circ_file, "r") as circ_in:
        next(circ_in)
        for line in circ_in:
            line_content = line.split()
            circ_id = line_content[1]
            parental_gene = line_content[0]
            result_output[circ_id] = [parental_gene]

    # get circRNA sequence length
    circ_seq_len_dict = {}
    with open(circ_seq_file, "r") as seq_in:
        for line in seq_in:
            header = line[1:]
            header = header.strip("\n")
            seq_len = next(seq_in)
            seq_len = len(seq_len.strip("\n"))
            circ_seq_len_dict[header] = str(seq_len)

    # read in ORF results
    orf_dict = {}
    with open(orf_file, "r") as orf_in:
        next(orf_in)
        for line in orf_in:
            line_content = line.split()
            circ_id = line_content[0]
            orf_dict[circ_id] = [line_content[1], line_content[2], line_content[3]]

    miranda_res_dict = {}
    with open(miranda_file, "r") as miranda_in:
        for line in miranda_in:
            line_content = line.split()
            circ_id = line_content[0]
            mirna_bs_density = line_content[3]
            max_mirna_results = line_content[4]
            miranda_res_dict[circ_id] = [mirna_bs_density, max_mirna_results]

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

        if circ_rna in circ_seq_len_dict:
            result_output[circ_rna].append(circ_seq_len_dict[circ_rna])

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
            result_output[circ_rna].append(3 * "0:0")

        if circ_rna in pep_dict:
            result_output[circ_rna].append(pep_dict[circ_rna])
        else:
            result_output[circ_rna].append("non_unique")

    # calcifer out write in working dir
    with open(working_dir + "calcifer_output.tab", "w") as calcifer_out:
        file_header = "circID\tparental_gene\tcirc_seq_len\tmirna_binding_site_density\tmost_mirna" \
                      "\trbp_circ_binding\trbp_bsj_binding\tlinear_seq_orf\tpseudo_circular_seq_orf" \
                      "\tmulti_cycle_seq_orf\tunique_region"
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
        move_mirna_res_cmd = "mv " + working_dir + "all_circs/*mir* " + working_dir + "all_circs/sub_result_dir/mirna_miranda_res_files/"
        os.system(move_mirna_res_cmd)
    
    dir_exists = os.path.exists(working_dir + "all_circs/sub_result_dir/rbp_fimo_res_files/")
    if not dir_exists:
        rbp_res_dir = "mkdir " + working_dir + "all_circs/sub_result_dir/rbp_fimo_res_files"
        os.system(rbp_res_dir)
    else:
        remove_old_res = "rm -r " + working_dir + "all_circs/sub_result_dir/rbp_fimo_res_files/*fimo*"
        os.system(remove_old_res)
    move_rbp_res_cmd_1 = "mv " + working_dir + "all_circs/*fimo* " + working_dir + "all_circs/sub_result_dir/rbp_fimo_res_files/"
    move_rbp_res_cmd_2 = "mv " + working_dir + "all_circs/*rbp* " + working_dir + "all_circs/sub_result_dir/rbp_fimo_res_files/"
    os.system(move_rbp_res_cmd_1)
    os.system(move_rbp_res_cmd_2)
    
    dir_exists = os.path.exists(working_dir + "all_circs/sub_result_dir/orf_res_files/")
    if not dir_exists:
        orf_res_dir = "mkdir " + working_dir + "all_circs/sub_result_dir/orf_res_files"
        os.system(orf_res_dir)
    move_orf_res_cmd_1 = "mv " + working_dir + "all_circs/*seq* " + working_dir + "all_circs/sub_result_dir/orf_res_files/"
    move_orf_res_cmd_2 = "mv " + working_dir + "all_circs/*orf* " + working_dir + "all_circs/sub_result_dir/orf_res_files/"
    move_orf_res_cmd_3 = "mv " + working_dir + "all_circs/*prediction* " + working_dir + "all_circs/sub_result_dir/orf_res_files/"
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


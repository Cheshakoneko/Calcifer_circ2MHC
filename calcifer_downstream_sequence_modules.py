#!/usr/bin/env python

import os
from pyfaidx import Faidx


# module for functions regarding retrieving the circRNA sequence #


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
                # also get exon numbers for circRNA naming
                elif information_type == "exon":
                    if line_content[0] + ":" + str(feature_start) in exon_annotation:
                        if line_content[0] + ":" + str(feature_end) not in exon_annotation[line_content[0] + ":" +
                                                                                           str(feature_start)][0]:
                            exon_annotation[line_content[0] + ":" + str(feature_start)][0].append(line_content[0] + ":"
                                                                                                  + str(feature_end))
                    else:
                        exon_annotation[line_content[0] + ":" + str(feature_start)] = [[line_content[0] + ":" +
                                                                                       str(feature_end)],
                                                                                       line_content[17]]
                    exon_endings[line_content[0] + ":" + str(feature_end)] = [str(feature_start), line_content[17]]

    return cds_annotation, three_utr_annotation, exon_annotation, exon_endings


# get all exons which are included in all found circRNAs #
# create fasta-files with the linear sequence, psuedo_circular sequence #
# (+ 25 bp from end to start and start to end) and multi cycle sequence (4 * linear seq) #
# also get sequence around bsj (250 bp + 25 bp into circRNA) for RBP analysis #
def circ_exon_seq(working_dir, gene_fasta, exon_anno, exon_endings):
    genes = Faidx(gene_fasta)
    circ_file = working_dir + "all_circs/two_unique_filtered.txt"
    circ_rnas = {}
    # keys are start or end position, every start/end can have multiple circRNAs (isoforms) #
    # read in circRNAs with identifier as key #
    with open(circ_file, "r") as circs_in:
        with open(working_dir + "all_circs/circ_bsj_seq.fasta", "w") as bsj_seq:
            for line in circs_in:
                line_content = line.split()
                all_pos = line_content[0].split(":")
                chromosome = all_pos[0]
                strand = line_content[0].split(":")[2]
                positions = all_pos[1].split("-")
                positions[0] = str(int(positions[0]) + 1)
                circ_rnas[line_content[0]] = [positions[0], positions[1], line_content[1:]]
                # get sequence around bsj +/- 250bp
                header = "\n>" + line_content[0]
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
                    fasta_seq_up = str(
                        genes.fetch(chromosome, (int(positions[0]) - 250), (int(positions[0]) + 24)).reverse.complement)
                    full_fasta_seq_up = ""
                    full_fasta_seq_up += ''.join(fasta_seq_up)
                    fasta_seq_down = str(
                        genes.fetch(chromosome, (int(positions[1]) - 24), (int(positions[1]) + 250)).reverse.complement)
                    full_fasta_seq_down = ""
                    full_fasta_seq_down += ''.join(fasta_seq_down)
                bsj_seq.write(header + ":side1\n")
                bsj_seq.write(full_fasta_seq_up)
                bsj_seq.write(header + ":side2\n")
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
                    real_exon_start = exon_endings[ident_pos_end][0]
                    exon_number = exon_endings[ident_pos_end][1]
                    if int(real_exon_start) < pos:
                        start_exon = start
                        end_point = pos
                        circ_rna_starts.append(int(start_exon))
                        # add exon number here (and also if in exon start)
                        circ_rna_exons[start_exon] = [int(end_point), chrom + ":" + str(start_exon) + "-" +
                                                      str(end_point) + ":" + exon_number]
            end_point = -1
            ident_pos = chrom + ":" + str(pos)
            # check for full exons in circRNA
            if ident_pos in exon_anno:
                start_exon = pos
                end_point = 0
                exon_ending_list = exon_anno[ident_pos][0]
                exon_number = exon_anno[ident_pos][1]
                # only take largest version of a given exon (i.e. largest end point inside a circ)
                for i in exon_ending_list:
                    end_pos = int(i.split(":")[1])
                    if end >= end_pos > end_point:
                        end_point = end_pos
                if end_point > pos:
                    circ_rna_starts.append(int(pos))
                    circ_rna_exons[pos] = [int(end_point), chrom + ":" + str(pos) + "-" + str(end_point) + ":"
                                           + exon_number]
        # add partial exon sequence at the end of circRNA
        if end_point == 0:
            circ_rna_starts.append(int(start_exon))
            circ_rna_exons[start_exon] = [int(pos), chrom + ":" + str(start_exon) + "-" + str(pos) + ":" + exon_number]
        if len(circ_rna_starts) == 0:
            circ_rna_starts.append(int(start))
            circ_rna_exons[start] = [int(end), chrom + ":" + str(start) + "-" + str(end) + ":NA"]

        circ_rna_seq_pos[key] = []
        sorted_exon_starts = sorted(circ_rna_starts)
        prev_exon_end = 0
        circ_rna_end = int(key.split(":")[1].split("-")[1])
        for exon_start in sorted_exon_starts:
            putative_end = circ_rna_exons[exon_start][0]
            if exon_start > prev_exon_end:
                prev_exon_end = circ_rna_exons[exon_start][0]
                circ_rna_seq_pos[key].append(circ_rna_exons[exon_start][1])
            # there is the rare case that an exon in front of the circRNA end overlaps with the real last exon
            # e.g. the backspliced exon. In this case we only consider the exon which was also involved in the
            # back-splicing
            elif putative_end == circ_rna_end and prev_exon_end != circ_rna_end:
                circ_rna_seq_pos[key] = circ_rna_seq_pos[key][:-1]
                prev_exon_end = circ_rna_exons[exon_start][0]
                circ_rna_seq_pos[key].append(circ_rna_exons[exon_start][1])

    # get linear circRNA sequence and exon numbers between start and end of each circRNA
    with open(working_dir + "all_circs/circ_naming.txt", "w") as naming:
        with open(working_dir + "all_circs/linear_seq.fasta", "w") as circ_seq_out:
            for key in circ_rna_seq_pos.keys():
                if len(circ_rna_seq_pos[key]) > 0:
                    exon_number_list = []
                    fasta_header = ">" + key
                    full_fasta_seq = ""
                    strand = key.split(":")[2]
                    if strand == "+":
                        for exon_pos in circ_rna_seq_pos[key]:
                            chromosome = exon_pos.split(":")[0]
                            position = exon_pos.split(":")[1]
                            start_pos = int(position.split("-")[0])
                            end_pos = int(position.split("-")[1])
                            if len(exon_pos.split(":")[2]) == 0:
                            	exon_number = "NA"
                            elif "NA" not in exon_pos.split(":")[2]:
                            	exon_number = int(exon_pos.split(":")[2][1:-2])
                            else:
                            	exon_number = "NA"
                            # double check sequence content in comparing the exon number with all already added exons
                            # the sequence is only added if exon number is greater then all others and not already
                            # included in the sequence
                            exon_number_list.append(exon_number)
                            fasta_seq = str(genes.fetch(chromosome, start_pos, end_pos))
                            full_fasta_seq += ''.join(fasta_seq)
                    elif strand == "-":
                        for exon_pos in circ_rna_seq_pos[key]:
                            chromosome = exon_pos.split(":")[0]
                            position = exon_pos.split(":")[1]
                            start_pos = int(position.split("-")[0])
                            end_pos = int(position.split("-")[1])
                            if len(exon_pos.split(":")[2]) == 0:
                            	exon_number = "NA"
                            elif "NA" not in exon_pos.split(":")[2]:
                            	exon_number = int(exon_pos.split(":")[2][1:-2])
                            else:
                            	exon_number = "NA"
                            # double check sequence content in comparing the exon number with all already added exons
                            # the sequence is only added if exon number is smaller then all others and not already
                            # included in the sequence
                            exon_number_list.append(exon_number)
                            fasta_seq = str(genes.fetch(chromosome, start_pos, end_pos).reverse.complement)
                            full_fasta_seq = ''.join((fasta_seq, full_fasta_seq))
                    circ_seq_out.write("\n" + fasta_header + "\n")
                    circ_seq_out.write(full_fasta_seq.upper())
                    final_exon_number_string = "(" + ','.join(str(e) for e in exon_number_list) + ")"
                    naming.write("\n" + key + "\t" + final_exon_number_string)
    os.system("sed -i \'1d\' " + working_dir + "all_circs/linear_seq.fasta")
    os.system("sed -i \'1d\' " + working_dir + "all_circs/circ_naming.txt")

    with open(working_dir + "all_circs/linear_seq.fasta", "r") as circ:
        with open(working_dir + "all_circs/pseudo_circular_seq.fasta", "w") as seq_out:
            line_count = 0
            for line1 in circ:
                header = line1.replace("\n", "")
                line2 = next(circ)
                circ_line = line2.replace("\n", "")
                seq = circ_line[-25:] + circ_line + circ_line[:25]
                if line_count != 0:
                    seq_out.write("\n" + header + "\n" + seq)
                else:
                    seq_out.write(header + "\n" + seq)
                    line_count += 1

    with open(working_dir + "all_circs/linear_seq.fasta", "r") as circ:
        with open(working_dir + "all_circs/multi_cycle_seq.fasta", "w") as seq_out:
            line_count = 0
            for line1 in circ:
                header = line1.replace("\n", "")
                line2 = next(circ)
                seq = 4 * line2.replace("\n", "")
                if line_count != 0:
                    seq_out.write("\n" + header + "\n" + seq)
                else:
                    seq_out.write(header + "\n" + seq)
                    line_count += 1

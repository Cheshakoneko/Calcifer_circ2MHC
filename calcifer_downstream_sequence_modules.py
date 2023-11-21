#!/usr/bin/env python

import os
import random
from pyfaidx import Fasta
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


# unused function to get random 3'utr and cds annotations to test for miRNA binding density #
# compare these two densities against miRNA density on circRNAs
def random_annotation(working_dir, cds_annotation, three_utr_annotation, gene_fasta):
    random.seed(13)
    rnd_cds = random.sample(cds_annotation.values(), k=1000)
    rnd_three_utr = random.sample(three_utr_annotation.values(), k=1000)
    with open(working_dir + "all_circs/rnd_cds_seq.fasta", "w") as cds_seq_out:
        for seq_pos in rnd_cds:
            sequence = os.popen('samtools faidx ' + gene_fasta + ' ' + seq_pos).read().split()[1:]
            out_sequence = ''.join(sequence)
            header = ">CDS_" + seq_pos[1:-1]
            cds_seq_out.write("\n" + header + "\n")
            cds_seq_out.write(str(out_sequence))
    os.system("sed -i \'1d\' " + working_dir + "all_circs/rnd_cds_seq.fasta")
    with open(working_dir + "all_circs/rnd_three_utr_seq.fasta", "w") as three_utr_seq_out:
        for seq_pos in rnd_three_utr:
            sequence = os.popen('samtools faidx ' + gene_fasta + ' ' + seq_pos).read().split()[1:]
            out_sequence = ''.join(sequence)
            header = ">3UTR_" + seq_pos[1:-1]
            three_utr_seq_out.write("\n" + header + "\n")
            three_utr_seq_out.write(str(out_sequence))
    os.system("sed -i \'1d\' " + working_dir + "all_circs/rnd_three_utr_seq.fasta")


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

    with open(working_dir + "all_circs/linear_seq.fasta", "w") as circ_seq_out:
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
    os.system("sed -i \'1d\' " + working_dir + "all_circs/linear_seq.fasta")

    with open(working_dir + "all_circs/linear_seq.fasta", "r") as circ:
        with open(working_dir + "all_circs/pseudo_circular_seq.fasta", "w") as seq_out:
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
                    
    with open(working_dir + "all_circs/linear_seq.fasta", "r") as circ:
        with open(working_dir + "all_circs/multi_cycle_seq.fasta", "w") as seq_out:
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


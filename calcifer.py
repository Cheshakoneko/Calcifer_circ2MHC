#!/usr/bin/env python

# import argument parser #
# enables the usage of calcifer like a command line tool #
import argparse

# import calcifer modules #
import calcifer_general_modules
import calcifer_filtering_modules
import calcifer_circexplorer2_modules
import calcifer_ciri2_modules
import calcifer_downstream_results_modules
import calcifer_downstream_sequence_modules
import calcifer_downstream_countmatrix_modules
import calcifer_downstream_mirna_modules
import calcifer_downstream_orf_modules
import calcifer_downstream_rbp_modules
import calcifer_list_mode_module


# Main Methods for every mode #


# Main method for circexplorer2 #
# runs ce2 to detect circRNAs #
def circexplorer2(args):

    file_path = args.path
    datasets = args.data.split(",")
    star_index = args.star
    genome_file = args.genome
    gene_pred_file = args.gene_pred
    read_type = args.read_type
    cores = args.threads

    for dataset in datasets:
    # change functions based on type of reads #

        if read_type == "se":
            working_dir = calcifer_general_modules.se_file_structure(file_path, dataset)
            unzip_trimmed_data = calcifer_general_modules.se_flexbar(dataset, working_dir)
            calcifer_circexplorer2_modules.se_star_mapping(working_dir, unzip_trimmed_data, star_index, cores)

        elif read_type == "pe":
            working_dir = calcifer_general_modules.pe_file_structure(file_path, dataset)
            unzip_trimmed_data_1, unzip_trimmed_data_2 = calcifer_general_modules.pe_flexbar(dataset, working_dir)
            calcifer_circexplorer2_modules.pe_star_mapping(working_dir, unzip_trimmed_data_1, unzip_trimmed_data_2,
                                                           star_index, cores)

        # function calls for ce2 #
        calcifer_circexplorer2_modules.ce2_parse(working_dir)
        calcifer_circexplorer2_modules.ce2_annotate(working_dir, gene_pred_file, genome_file)
        calcifer_circexplorer2_modules.ce2_initial_filter(working_dir)


# Main method for ciri2 #
# runs ciri2 to detect circRNAs #
def ciri2(args):

    file_path = args.path
    datasets = args.data.split(",")
    bwa_index = args.bwaindex
    genome_path = args.genome
    gtf_file = args.gtf
    ciri_path = args.cpath
    read_type = args.read_type
    cores = args.threads

    for dataset in datasets:
    # change functions based on type of reads #

        if read_type == "se":
            working_dir = calcifer_general_modules.se_file_structure(file_path, dataset)
            unzip_trimmed_data = calcifer_general_modules.se_flexbar(dataset, working_dir)
            calcifer_ciri2_modules.se_bwa_mapping(working_dir, unzip_trimmed_data, bwa_index, cores)

        elif read_type == "pe":
            working_dir = calcifer_general_modules.pe_file_structure(file_path, dataset)
            unzip_trimmed_data_1, unzip_trimmed_data_2 = calcifer_general_modules.pe_flexbar(dataset, working_dir)
            calcifer_ciri2_modules.pe_bwa_mapping(working_dir, unzip_trimmed_data_1, unzip_trimmed_data_2, bwa_index, cores)

        # function call for ciri2 #
        calcifer_ciri2_modules.ciri2_find(ciri_path, working_dir, gtf_file, genome_path)


# Main method for downstream analysis of ce2 and ciri2 results #
# general downstream analysis for all circRNAs detected with both tools #
def downstream(args):

    datasets = args.data.split(",")
    conditions = args.condition.split(",")
    condition_name = args.condition_names.split(",")
    working_dir = args.path
    genome_fasta = args.genome_fasta
    read_type = args.read_type
    gtf_file = args.gtf
    mirna_run = args.mirna
    pep_ref = args.pep
    rbp_db = args.rbp
    min_aa = args.minlen
    qval = args.qval
    ubsjr_filter = args.ubsjr

    # create suitable data structure for the results #
    calcifer_downstream_results_modules.data_structure_filtering(working_dir)
    # using chimeric junctions from star as ground truth and filter for canonical splice sites #
    calcifer_filtering_modules.chimeric_filtering(working_dir, datasets, genome_fasta, gtf_file)
    # merge results based on conditions and all results per default #
    calcifer_filtering_modules.merging_results(working_dir, datasets, conditions, ubsjr_filter)
    # annotate all results remaining after the strict filters #
    all_filtered_circs = working_dir + "all_circs/two_unique_filtered.txt"

    calcifer_downstream_countmatrix_modules.deseq2_analysis(working_dir, datasets, conditions, condition_name,
                                                            read_type, gtf_file)

    cds_anno, three_utr_anno, exon_anno, exon_endings = calcifer_downstream_sequence_modules.mirna_annotation(gtf_file)
    calcifer_downstream_sequence_modules.circ_exon_seq(working_dir, genome_fasta, exon_anno, exon_endings)
    
    calcifer_downstream_mirna_modules.mirna_analysis(working_dir, mirna_run)
    
    calcifer_downstream_orf_modules.orf_detection(working_dir)
    calcifer_downstream_orf_modules.longest_orf_filtering(working_dir)
    ires_m6a_dict = calcifer_downstream_orf_modules.ires_m6a_prediction(working_dir)
    calcifer_downstream_orf_modules.unique_peptides_analysis(working_dir, pep_ref, ires_m6a_dict)
    
    calcifer_downstream_rbp_modules.rbp_analysis_circ(working_dir, rbp_db, qval)
    calcifer_downstream_rbp_modules.rbp_analysis_bsj(working_dir, rbp_db, qval)

    calcifer_downstream_results_modules.final_output(working_dir, conditions, mirna_run)
    calcifer_downstream_results_modules.clean_up(working_dir, mirna_run)


# main method for the list mode with just circRNA positions as input
def list_mode(args):
    working_dir = args.path
    crna_list_file = args.circ_list
    genome_fasta = args.genome_fasta
    gtf_file = args.gtf
    mirna_run = args.mirna
    pep_ref = args.pep
    rbp_db = args.rbp
    min_aa = args.minlen
    qval = args.qval

    calcifer_list_mode_module.parental_gene_name(working_dir, crna_list_file, gtf_file)

    cds_anno, three_utr_anno, exon_anno, exon_endings = calcifer_list_mode_module.mirna_annotation(gtf_file)
    calcifer_list_mode_module.circ_exon_seq(working_dir, genome_fasta, exon_anno, exon_endings)
    
    calcifer_list_mode_module.mirna_analysis(working_dir, mirna_run)
    
    calcifer_list_mode_module.orf_detection(working_dir)
    calcifer_list_mode_module.longest_orf_filtering(working_dir)
    ires_m6a_dict = calcifer_list_mode_module.ires_m6a_prediction(working_dir)
    calcifer_list_mode_module.unique_peptides_analysis(working_dir, pep_ref, ires_m6a_dict)
    
    calcifer_list_mode_module.rbp_analysis_circ(working_dir, rbp_db, qval)
    calcifer_list_mode_module.rbp_analysis_bsj(working_dir, rbp_db, qval)

    calcifer_list_mode_module.final_output(working_dir)

    calcifer_list_mode_module.clean_up(working_dir, mirna_run)
    

# Main method for a full run containing everything from the previous main functions #
# the general mode to do all analyses in one go #
def full_run(args):

    file_path = args.path
    datasets = args.data.split(",")
    conditions = args.condition.split(",")
    condition_name = args.condition_names.split(",")
    star_index = args.star
    genome_file = args.genome
    read_type = args.read_type
    bwa_index = args.bwaindex
    ciri_path = args.cpath
    gene_pred_file = args.gene_pred
    genome_fasta = args.genome_fasta
    gtf_file = args.gtf
    mirna_run = args.mirna
    pep_ref = args.pep
    rbp_db = args.rbp
    cores = args.threads
    min_aa = args.minlen
    qval = args.qval
    ubsjr_filter = args.ubsjr

    for dataset in datasets:
        # bwa and star mapping of the dataset #
        if read_type == "se":
            working_dir = calcifer_general_modules.se_file_structure(file_path, dataset)
            unzip_trimmed_data = calcifer_general_modules.se_flexbar(dataset, working_dir)
            calcifer_circexplorer2_modules.se_star_mapping(working_dir, unzip_trimmed_data, star_index, cores)
            calcifer_ciri2_modules.se_bwa_mapping(working_dir, unzip_trimmed_data, bwa_index, cores)
        elif read_type == "pe":
            working_dir = calcifer_general_modules.pe_file_structure(file_path, dataset)
            unzip_trimmed_data_1, unzip_trimmed_data_2 = calcifer_general_modules.pe_flexbar(dataset, working_dir)
            calcifer_circexplorer2_modules.pe_star_mapping(working_dir, unzip_trimmed_data_1, unzip_trimmed_data_2,
                                                           star_index, cores)
            calcifer_ciri2_modules.pe_bwa_mapping(working_dir, unzip_trimmed_data_1, unzip_trimmed_data_2, bwa_index,
                                                  cores)

        # running ce2 on dataset #
        calcifer_circexplorer2_modules.ce2_parse(working_dir)
        calcifer_circexplorer2_modules.ce2_annotate(working_dir, gene_pred_file, genome_file)
        calcifer_circexplorer2_modules.ce2_initial_filter(working_dir)

        # running ciri2 on dataset #
        calcifer_ciri2_modules.ciri2_find(ciri_path, working_dir, gtf_file, genome_file)

    # running downstream analysis on full setup #
    working_dir = file_path
    # create suitable data structure for the results #
    calcifer_downstream_results_modules.data_structure_filtering(working_dir)
    # using chimeric junctions from star as ground truth and filter for canonical splice sites #
    calcifer_filtering_modules.chimeric_filtering(working_dir, datasets, genome_fasta, gtf_file)
    # merge results based on conditions and all results per default #
    calcifer_filtering_modules.merging_results(working_dir, datasets, conditions, ubsj_filter)
    # annotate all results remaining after the strict filters #
    all_filtered_circs = working_dir + "all_circs/two_unique_filtered.txt"

    calcifer_downstream_countmatrix_modules.deseq2_analysis(working_dir, datasets, conditions, condition_name,
                                                            read_type, gtf_file)

    cds_anno, three_utr_anno, exon_anno, exon_endings = calcifer_downstream_sequence_modules.mirna_annotation(gtf_file)
    calcifer_downstream_sequence_modules.circ_exon_seq(working_dir, genome_fasta, exon_anno, exon_endings)
    
    calcifer_downstream_mirna_modules.mirna_analysis(working_dir, mirna_run)
    
    calcifer_downstream_orf_modules.orf_detection(working_dir)
    calcifer_downstream_orf_modules.longest_orf_filtering(working_dir)
    ires_m6a_dict = calcifer_downstream_orf_modules.ires_m6a_prediction(working_dir)
    calcifer_downstream_orf_modules.unique_peptides_analysis(working_dir, pep_ref, ires_m6a_dict)
    
    calcifer_downstream_rbp_modules.rbp_analysis_circ(working_dir, rbp_db, qval)
    calcifer_downstream_rbp_modules.rbp_analysis_bsj(working_dir, rbp_db, qval)

    calcifer_downstream_results_modules.final_output(working_dir, conditions, mirna_run)
    calcifer_downstream_results_modules.clean_up(working_dir, mirna_run)


# Multiple Parser for every mode #
# Selection of the different mode options of the pipeline #
# Can be used like a commandline tool #

# add args.parse as parser for the main function #
parser = argparse.ArgumentParser()
# add subparsers to the main parser for different modes #
subparsers = parser.add_subparsers()

# specify different modes for the pipeline with needed arguments and helps #
circexplorer2_parser = subparsers.add_parser('circexplorer2',
                                             usage="calcifer.py circexplorer2 -path [path] -data [name] -star [index] "
                                                   "-genome [fasta] -gene_pred [txt] -rt [se/pe] -t [threads]")
circexplorer2_parser.add_argument("-path", action="store", dest="path", help="Input a path to the files", required=True)
circexplorer2_parser.add_argument("-data", action="store", dest="data", help="Input name of datasets", required=True)
circexplorer2_parser.add_argument("-star", action="store", dest="star", help="Path to star-index", required=True)
circexplorer2_parser.add_argument("-genome", action="store", dest="genome", help="Path to genome-fasta", required=True)
circexplorer2_parser.add_argument("-gene_pred", action="store", dest="gene_pred", help="Path to ref gene_pred.txt "
                                                                                       "(generate before run!)",
                                  required=True)
circexplorer2_parser.add_argument("-rt", action="store", dest="read_type",
                                  help="Type of the reads, se for single-end and pe for paired-end accepted",
                                  required=True)
circexplorer2_parser.add_argument("-t", action="store", dest="threads", help="Amount of threads to use for the mapping",
                                  required=True)
circexplorer2_parser.set_defaults(func=circexplorer2)


ciri2_parser = subparsers.add_parser('ciri2', usage="calcifer.py ciri2 -path [path] -data [name] -bwa [index] -genome "
                                     "[fasta] -ref [gff] -cpath [ciri2] -rt [se/pe] -t [threads]")
ciri2_parser.add_argument("-path", action="store", dest="path", help="Input a path to the files", required=True)
ciri2_parser.add_argument("-data", action="store", dest="data", help="Input name of dataset", required=True)
ciri2_parser.add_argument("-bwa", action="store", dest="bwaindex", help="Path to bwa-index", required=True)
ciri2_parser.add_argument("-genome", action="store", dest="genome", help="Path to genome-fasta", required=True)
ciri2_parser.add_argument("-gtf", action="store", dest="gtf", help="Path to ref genome primary annotation",
                          required=True)
ciri2_parser.add_argument("-cpath", action="store", dest="cpath", help="Path to ciri2.pl", required=True)
ciri2_parser.add_argument("-rt", action="store", dest="read_type",
                          help="Type of the reads, se for single-end and pe for paired-end accepted", required=True)
ciri2_parser.add_argument("-t", action="store", dest="threads", help="Amount of threads to use for the mapping",
                          required=True)
ciri2_parser.set_defaults(func=ciri2)


downstream_parser = subparsers.add_parser('downstream', usage="calcifer.py downstream -path [path] -data [data] -con "
                                                              "[list] -con_names [list] -genome_fasta [fasta] -rt "
                                                              "[se/pe] -gtf [path] -mirna [path] -pep [path] "
                                                              "-rbp [path]")
downstream_parser.add_argument("-path", action="store", dest="path", help="Input a path to the files", required=True)
downstream_parser.add_argument("-data", action="store", dest="data", help="Input list of names of the datasets",
                               required=True)
downstream_parser.add_argument("-con", action="store", dest="condition",
                               help="List of amount of datasets per condition, same order as dataset names",
                               required=True)
downstream_parser.add_argument("-con_names", action="store", dest="condition_names",
                               help="List of condition names, same order as dataset names and conditions",
                               required=True)
downstream_parser.add_argument("-genome_fasta", action="store", dest="genome_fasta",
                               help="Path to fasta-file of ref genome", required=True)
downstream_parser.add_argument("-rt", action="store", dest="read_type",
                               help="Type of reads, se for single-end and pe for paired-end accepted", required=True)
downstream_parser.add_argument("-gtf", action="store", dest="gtf",
                               help="Path to ensembl gtf-file",
                               required=True)
downstream_parser.add_argument("-mirna", action="store", dest="mirna",
                               help="Path to miRNA database", required=True)
downstream_parser.add_argument("-pep", action="store", dest="pep",
                               help="Path to fasta-file with all pc-transcripts", required=True)
downstream_parser.add_argument("-rbp", action="store", dest="rbp", help="Path to rbp db file", required=True)
downstream_parser.add_argument("-min", action="store", dest="minlen", help="Minimal amino acids lenghts for circRNA "
                                                                           "peptides, minlen < 4 may lead to error",
                               default=10)
downstream_parser.add_argument("-rbp_cutoff", action="store", dest="qval", help="q value threshold for FIMO results",
                               default=0.1)
downstream_parser.add_argument("-ubsjr_filter", action="store", dest="ubsjr", help="Minimum unique backsplice junction "
                                                                                   "supporting reads for high conf.",
                               default=2)
downstream_parser.set_defaults(func=downstream)


list_parser = subparsers.add_parser('list', usage="calcifer.py list -path [path] -circ_list [data] -genome_fasta "
                                                  "[path] -gtf [path] -mirna [path] -pep [path] -rbp [path]")
list_parser.add_argument("-path", action="store", dest="path", help="Input a path to the files", required=True)
list_parser.add_argument("-circ_list", action="store", dest="circ_list", help="Path to circRNA list", required=True)
list_parser.add_argument("-genome_fasta", action="store", dest="genome_fasta", help="Path to fasta-file of ref genome",
                         required=True)
list_parser.add_argument("-gtf", action="store", dest="gtf", help="Path to ensembl gtf-file", required=True)
list_parser.add_argument("-mirna", action="store", dest="mirna", help="Path to miRNA database", required=True)
list_parser.add_argument("-pep", action="store", dest="pep", help="Path to fasta-file with all pc-transcripts",
                         required=True)
list_parser.add_argument("-rbp", action="store", dest="rbp", help="Path to rbp db file", required=True)
list_parser.add_argument("-min", action="store", dest="minlen", help="Minimal amino acids lenghts for circRNA peptides,"
                                                                     " minlen < 4 may lead to error", default=10)
list_parser.add_argument("-rbp_cutoff", action="store", dest="qval", help="q value threshold for FIMO results",
                         default=0.1)
list_parser.set_defaults(func=list_mode)


full_run_parser = subparsers.add_parser('full_run', usage="calcifer.py full_run -path [path] -data [data] "
                                                          "-star [index] -genome [fasta] -gtf [gtf] -gene_pred [txt] "
                                                          "-rt [se/pe] -con [list] -cpath [ciri2] -bwa [index] "
                                                          "-genome_fasta [fasta] -mirna [path] -pep [path] -rbp [path]"
                                                          " -t [threads]")
full_run_parser.add_argument("-path", action="store", dest="path", help="Input a path to the files", required=True)
full_run_parser.add_argument("-data", action="store", dest="data", help="Input list of names of the datasets",
                             required=True)
full_run_parser.add_argument("-star", action="store", dest="star", help="Path to star-index", required=True)
full_run_parser.add_argument("-genome", action="store", dest="genome", help="Path to genome-fasta", required=True)
full_run_parser.add_argument("-gtf", action="store", dest="gtf", help="Path to ensembl gtf-file", required=True)
full_run_parser.add_argument("-gene_pred", action="store", dest="gene_pred", help="Path to ref gene_pred.txt",
                             required=True)
full_run_parser.add_argument("-rt", action="store", dest="read_type",
                             help="Type of the reads, se for single-end and pe for paired-end accepted", required=True)
full_run_parser.add_argument("-con", action="store", dest="condition", help="List of amount of datasets per condition",
                             required=True)
full_run_parser.add_argument("-con_names", action="store", dest="condition_names",
                             help="List of condition names, same order as dataset names and conditions",
                             required=True)
full_run_parser.add_argument("-cpath", action="store", dest="cpath", help="Path to ciri2.pl", required=True)
full_run_parser.add_argument("-bwa", action="store", dest="bwaindex", help="Path to bwa-index", required=True)
full_run_parser.add_argument("-genome_fasta", action="store", dest="genome_fasta",
                             help="Path to fasta-file of ref genome", required=True)
full_run_parser.add_argument("-mirna", action="store", dest="mirna",
                             help="Path to miRNA database", required=True)
full_run_parser.add_argument("-pep", action="store", dest="pep", help="Path to fasta-file with all pc-transcripts",
                             required=True)
full_run_parser.add_argument("-rbp", action="store", dest="rbp", help="Path to rbp db file", required=True)
full_run_parser.add_argument("-t", action="store", dest="threads", help="Amount of threads to use for the mapping",
                             required=True)
full_run_parser.add_argument("-min", action="store", dest="minlen", help="Minimal amino acids lenghts for circRNA "
                                                                         "peptides, minlen < 4 may lead to error",
                             default=10)
full_run_parser.add_argument("-rbp_cutoff", action="store", dest="qval", help="q value threshold for FIMO results",
                             default=0.1)
full_run_parser.add_argument("-ubsjr_filter", action="store", dest="ubsjr", help="Minimum unique backsplice junction "
                                                                                 "supporting reads for high conf.",
                             default=2)
full_run_parser.set_defaults(func=full_run)


# main function is passed to arg.parser #
# at the user's choice the appropriate main function for each mode will be executed #
if __name__ == '__main__':
    args = parser.parse_args()
    try:
        args.func(args)
    except AttributeError:
        print("Show usage: calcifer.py [mode] -h")

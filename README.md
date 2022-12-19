# Calcifer: A workflow for circRNA detection and analysis
Author: Andre Brezski, Kathi Zarnack    

---

# 1. Introduction
Calcifer is a workflow for highly automated detection and analysis of circRNAs in RNA-Seq datasets. It allows the evaluation of RNA-Seq read data up to a list of characterized circRNA isoforms, as well as the prediction of possible functions.

# 2. Overview

![Calcifer workflow](/assets/images/calcifer_new_workflow.png)

# 3. CircRNA detection
In a first step the RNA-seq data is aligned with STAR [1] and bwa [2] against the reference genome. The resulting bam-files are used as input for CIRCexplorer2 [3] and CIRI2 [4] respectively. These both tools yield an unfiltered list of putative circRNAs, which is further processed.

# 4. CircRNA filtering
The raw circRNA results from both tools are then further filtered for canonical splice sites, length of the circRNA, encompassing junctions and uniquely mapped back-splice junction reads. At the end only circRNAs which have at least 2 uniquely back-splice junction reads are considered for further analysis.

# 5. Downstream analysis
Based on the filtered circRNA list there is a broad downstream analysis. To shed light on putative biogenesis and function the downstream analysis consists out of detection of putative miRNA binding, RBP binding and open reading frames. 

#### 5.1. Linear and circular count data
In general a count matrix for linear and circular mapped reads is created. If needed these can be utilized by the user for differential expression analysis for multiple condition datasets. A rmd-file for a baseline DESeq2 analysis is included in the major workflow folder.

#### 5.2. miRNA binding site detection
MiRNA binding sites are detected by miRanda [5] on the circular exonic sequence for each circRNA. To enable miRNA binding analysis over the back-splice junction sequence, the linear sequence is extended by 25 bp from the opposite end respectively.

#### 5.3. RBP binding site detection
The RBP binding prediction is performed with FIMO [6] on the same sequence (back-splice junction extended linear exon sequence) and additional on the not included sequence around the back-splice junction. CircRNA biogenesis can be enabled by RBP binding close to the back-splice junction. Also a putative function of circRNAs is the direct binding of RBPs.

#### 5.4. ORF prediction
The ORF prediction is performed on the linear- as well as the 2-fold, 3-fold and 4-fold exonic circRNA sequence. These enables the prediction of longer ORFs, which span over the back-splice junction as well as multiple reading frames. For the prediction of the longest ORF TransDecoder [7] is used.

---

# Example output
CircRNA isoform: 
> chr21:29321220-29329693:+	

Parental gene:
> BACH1	

CircRNA type:
> exon	

Unique back-splice junction reads:
> 21

MiRNA binding site density [binding sites/bp]:
> 0.09532	

Most bound single miRNA (name:binding sites count:proportion to all miRNA binding sites):
> hsa-miR-548am-3p:2:0.01143	

Most binding RBP on circRNA sequence (linear sequence with 25 bp added on both sites to simulate circularization):
> RNCMPT00043:1	

Most binding RBP around the back-splice junction (250 bp up- and downstream of the back-splice junction and 25 bp into the circRNA on both sites):
> no_rbp_binding	

Linear sequence ORF:
> 3prime_partial:593	

1 cycle ORF:
> complete:597	

2 cycle ORF:
> complete:597	

3 cycle ORF:
> complete:597	

Uniqueness of putative ORF peptide:
> non_unique	

Circular-to-linear ratio:
> 0.28852
 

---

# Citing Calcifer

> Manuscript in preparation

---

# Usage
Calcifer can be executed from the command line. The dependencies should be installed beforehand in a suitable environment (e.g. miniconda). A proper conda environment can be imported with the calcifer_env.yml file, which will be present in the main folder.

The proper usage can be shown with calcifer --help.

A small test dataset will be included to verify the functionality.

---

# Literature
[1] Dobin, Alexander, et al. "STAR: ultrafast universal RNA-seq aligner." Bioinformatics 29.1 (2013): 15-21.

[2] Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows–Wheeler transform. Bioinformatics, 25(14), 1754–1760.

[3] Zhang, Xiao-Ou, et al. "Diverse alternative back-splicing and alternative splicing landscape of circular RNAs." Genome research 26.9 (2016): 1277-1287.

[4] Gao, Yuan, Jinyang Zhang, and Fangqing Zhao. "Circular RNA identification based on multiple seed matching." Briefings in bioinformatics 19.5 (2018): 803-810.

[5] John, Bino, et al. "Human microRNA targets." PLoS biology 2.11 (2004): e363.

[6] Grant, Charles E., Timothy L. Bailey, and William Stafford Noble. "FIMO: scanning for occurrences of a given motif." Bioinformatics 27.7 (2011): 1017-1018.

[7] Haas, B., and A. J. G. S. Papanicolaou. "TransDecoder (find coding regions within transcripts)." Google Scholar (2016).

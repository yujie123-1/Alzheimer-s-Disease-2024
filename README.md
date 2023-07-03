# Multi-Omics Factor Analysis in AD using RNA-seq Data
Alzheimer’s disease (AD) is a complex neurodegenerative disorder that lacks effective clinical therapeutics and a detailed understanding of its pathogenic mechanisms. Previous research on the transcriptional-level pathogenic mechanisms of AD has elucidated the distinct characteristics of gene expression, RNA splicing, and post-transcriptional modifications in the brains of patients with AD. However, as these studies primarily focused on exploring AD pathogenesis from independent angles of transcriptional regulation, they did not comprehensively explore the integrated effects arising from diverse transcriptional layers. Therefore, we assessed the influence of 3 transcriptional layers—gene expression, alternative splicing, and alternative polyadenylation—on AD by integrating analysis using RNA sequencing data from patients with AD and controls. 
# Bioinformatic Methods
In this project, it mainly includes four aspects, which are:
* Gene expression
* Alternative splicing
* Alternative polyadenylation
* Multi-omics factor analys
## Gene Expression
Changes in gene expression provide insights into disease mechanisms and drug development. In this project, we utilized `DESeq2` for differential gene expression analysis.
## Alternative Splicing
Dysregulated alternative splicing has emerged as a mechanism broadly implicated in the pathogenesis of AD. In this project, we used `LeafCutter` for differential alternative splicing analysis.
## Alternative Polyadenylation
Dysregulation of alternative polyadenylation can result in inefficient gene expression, potentially contributing to the development of AD. In this project, we employed `Dapars` for differential alternative polyadenylation analysis.
## Multi-Omics Factor Analysis (MOFA)
We used `MOFA2`, which is a factor analysis model that provides a general framework for the integration of multi-omic data sets in an unsupervised fashion, to integrate gene expression, alternative splicing, and alternative polyadenylation datasets of patients with AD and control to explore the complex pathogenesis of AD and the interrelated factors that contribute to its onset and progression.
# Contact
If you have any comments, suggestions, questions, bug reports, etc, feel free to contact Yujie Yang [(yj.yang1@siat.ac.cn)]().
# License


![](https://img.shields.io/badge/Python-3.10.2-blue) 
# The principle of different alternative splicing events discovered by LeafCutter.
LeafCutter gathers all mapped reads in a study and detects overlapping introns indicated by split reads. Next, it builds a graph linking all overlapping introns with shared donor or acceptor splice sites. These interconnected sections of the graph create clusters, representing alternative intron excision events. 
LeafCutter can pinpoint different alternative splicing events, including       
(1) Exon skipping.     
(2) 5' or 3' exon extension.     
(3) Complex splicing.     
(4) Exon skipping or alternative start or end.      

![image](https://github.com/yujie123-1/Alzheimer-s-Disease-2024/assets/74124083/9a4a3c9d-24a2-4dd4-9720-03b18196c247)

# The principle of alternative splicing annotation.
the annotation is based on clusters. It involves three steps: 
* Sorting all introns within a cluster according to their starting positions.
* Checking for exons between the minimum and maximum positions in the cluster.
* Generating plot of all exons and intron positions in the cluster.

## The files involved in the annotation process:
* Human exon information file
which should look something like this:
```swift
chr1	11869	12227	+	DDX11L1
chr1	12613	12721	+	DDX11L1
chr1	13221	14409	+	DDX11L1
chr1	12010	12057	+	DDX11L1
chr1	29534	29570	-	WASH7P
chr1	24738	24891	-	WASH7P
chr1	18268	18366	-	WASH7P
chr1	17915	18061	-	WASH7P
chr1	17606	17742	-	WASH7P
```
Each column contains the following information: chromosome number, exon start position, exon end position, plus/minus strand, and gene name.
* Intron clustering file generated by LeafCutter.
which should look something like this:
```swift
chr1:17055:17233:clu_1 21 13 18 20 17 12 11 8 15 25
chr1:17055:17606:clu_1 4 11 12 7 2 0 5 2 4 4
chr1:17368:17606:clu_1 127 132 128 55 93 90 68 43 112 137
chr1:668593:668687:clu_2 3 11 1 3 4 4 8 1 5 16
chr1:668593:672093:clu_2 11 16 23 10 3 20 9 6 23 31
```
Each column corresponds to a different sample (original bam file) and each row to an intron, which is identified as chromosome:intron_start:intron_end:cluster_id.

# Contact
If you have any comments, suggestions, questions, bug reports, etc, feel free to contact Yujie Yang ([yj.yang1@siat.ac.cn]()).

# License
LicenseFinder is released under the [MIT](http://www.opensource.org/licenses/mit-license) License. 

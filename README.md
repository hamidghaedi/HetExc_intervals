# genomic-interval-enrichment-analysis-for-extreme-heterozygote-excess-HetExc-variants

## Intro

Deviations from Hardy-Weinberg Equilibrium (HWE) are not solely rooted in biological or genetic factors; sequencing/mapping errors play a significant role. Previous studies have underscored extreme heterozygote excess (HetExc) for specific variants in unstable genomic regions, suggesting the influence of genomic instability on observed deviations from HWE. Given the importance of variants occurring in coding regions for germline and somatic variant classification, we aimed at genomic intervals within the human exome enriched for HetExc variants. The identified regions could serve as indicators of potential spurious mutations, providing essential guidance for variant classification in genomic research.

## Methods

### Defining genomic intervals with HetExc variants

**Obtaining gencode.v44.annotation and processing** The file was downloaded from [GENCODE](https://www.gencodegenes.org/human/), and further processed to create a subset of exons from protein-coding Ensembl canonical transcripts. The processing details are outlined in the `gtf2bed.ipynb` notebook, in the scripts directory.

**Obtaining vcf files and processing**

We obtained exome and genome VCF files from gnomAD v4 (accessed on 2023-11-16) and processed them on HPC clusters provided by the Digital Research Alliance of Canada (the Alliance). The processing included the following steps: 1. Obtaining VCF files for chromosomes 1-22 and X.

2.  Filtering variants to retain only HetExc variants (indicated by an inbreeding coefficient ≤ -0.3).

3.  Intersecting VCF files with exon coordinates obtained from GENCODE.

4.  Uniting proximal intervals (+/- 300 bp), then filtering united intervals with a variant count of less than two.

The Bash script used for processing exome VCF files can be found in the scripts directory under the name gnomad_proc.sh.

**Annotating the identified intervals**

## Results

The results section presents details regarding our analysis of exome VCF files. If genome VCF files were considered for analysis, we explicitly indicate that. Therefore, the results are primarily based on exome VCF files unless stated otherwise.

### HetExc interval count

We ordered genes based on the number of HetExc intervals identified in them to see which genes have the highest level of HetExc intervals. The table below represents our findings, where it can be seen that *MUC3A* is the top one, harboring a total of eight intervals residing in four exons.

| Gene Symbol      | Transcript ID      | HetExc interval count | Affected exon IDs                                                                                                |
|------------------|--------------------|-----------------------|------------------------------------------------------------------------------------------------------------------|
| *MUC3A*          | ENST00000379458.9  | 8                     | ENSE00003733255.1, ENSE00001760877.2, ENSE00003728907.1, *ENSE00003711593.1*                                     |
| ENSG00000241489  | ENST00000651111.1  | 7                     | ENSE00003850317.1                                                                                                |
| *FGF13*          | ENST00000315930.11 | 7                     | ENSE00003765119.2, ENSE00001124987.8                                                                             |
| *IDS*            | ENST00000340855.11 | 7                     | ENSE00001368349.7                                                                                                |
| *ANKRD36C*       | ENST00000295246.7  | 6                     | ENSE00003650405.1, ENSE00003564919.1, ENSE00002490625.1, ENSE00002460612.1, ENSE00002471049.1, ENSE00002526875.1 |
| *FCGBP*          | ENST00000616721.6  | 6                     | ENSE00003727252.1, ENSE00003736302.1, ENSE00003726191.1, ENSE00003713265.1, ENSE00003717014.1, ENSE00003742840.1 |
| *MBNL3*          | ENST00000370853.8  | 6                     | ENSE00002203352.2                                                                                                |
| *MUC4*           | ENST00000463781.8  | 6                     | ENSE00001854802.1                                                                                                |
| *NBPF1*          | ENST00000430580.6  | 6                     | ENSE00003524745.1, ENSE00003641270.1, ENSE00003610699.1, ENSE00001784848.1, ENSE00003477087.1, ENSE00001765538.1 |
| *AFF2*           | ENST00000370460.7  | 5                     | ENSE00001452771.2                                                                                                |
| *ANKRD36*        | ENST00000420699.9  | 5                     | ENSE00002456047.1, ENSE00002485219.1, ENSE00002463900.1, ENSE00002496065.1, ENSE00003626869.1                    |
| *CNTNAP3B*       | ENST00000377561.7  | 5                     | ENSE00003468443.1, ENSE00003676979.1, ENSE00003669769.1, ENSE00003522083.1, ENSE00001844388.3                    |
| *ARMCX5-GPRASP2* | ENST00000652409.1  | 4                     | ENSE00003842799.1                                                                                                |
| *DCX*            | ENST00000636035.2  | 4                     | ENSE00003900373.1                                                                                                |
| *GPRASP3*        | ENST00000457056.6  | 4                     | ENSE00001458491.1                                                                                                |
| *HRNR*           | ENST00000368801.4  | 4                     | ENSE00001447986.4                                                                                                |
| *KMT2C*          | ENST00000262189.11 | 4                     | ENSE00003680506.1, ENSE00003668340.1, ENSE00003571980.1, ENSE00001783379.1                                       |
| *LANCL3*         | ENST00000378619.4  | 4                     | ENSE00001957816.2                                                                                                |
| *LILRA2*         | ENST00000391738.8  | 4                     | ENSE00003624947.1, ENSE00003459519.1, ENSE00003895381.1, ENSE00003894905.1                                       |
| *MAP2K3*         | ENST00000342679.9  | 4                     | ENSE00003566030.1, ENSE00003656597.1, ENSE00003544515.1, ENSE00003553037.1                                       |

### HetExc interval variant density

A total of 631 HetExc intervals were identified, with a median size of 89 bp (interquartile range: 227), harboring 2 to 113 HetExc variants. There are cases with HetExc interval length \> 1 kb, with a significant number of variants in them. The table below shows cases where the HetExc interval length is over 1 kb, along with the number of variants identified in them.

To further prioritize the intervals, we calculated a weighted variant density score. The approach involves filtering intervals with a length less than 20bp then the weighted variant density ($wvd_i$) for interval $i$ is calculated using the formula:

$$ wvd_i = \log_{10} \left( \frac{v_i}{\sqrt{\frac{l_i}{v_i^2}}} \right) $$

where $v_i$ represents the variant count in interval $i$ and $l_i$ is the length of the interval in base pairs (bp).The calculated weighted variant density prioritizes intervals with larger sizes while simultaneously considering a higher number of variants. This logarithmic transformation ensures that the metric reflects the joint impact of interval size and variant count, providing a comprehensive measure of variant density across genomic intervals.

| HetExc interval          | HetExc var. count | HetExc int. length(bp) | Weighted HetExc var. density | GeneSymb.(Transcript ID)      | Exon ID           |
|--------------------------|-------------------|------------------------|------------------------------|-------------------------------|-------------------|
| chr11:1016608-1018523    | 113               | 1915                   | 2.47                         | *MUC6* (ENST00000421673.7)    | ENSE00001739408.1 |
| chr3:75736810-75738889   | 77                | 2079                   | 2.11                         | *ZNF717* (ENST00000652011.2)  | ENSE00003849189.2 |
| chr11:56700404-56701554  | 59                | 1150                   | 2.01                         | *OR9G1* (ENST00000642097.1)   | ENSE00003812474.1 |
| chr11:56375721-56376331  | 48                | 610                    | 1.97                         | *OR8U1* (ENST00000302270.1)   | ENSE00001142881.1 |
| chr7:100959267-100960873 | 60                | 1606                   | 1.95                         | *MUC3A* (ENST00000379458.9)   | ENSE00003733255.1 |
| chr7:100959267-100960873 | 60                | 1606                   | 1.95                         | *MUC3A* (ENST00000379458.9)   | ENSE00001760877.2 |
| chr17:21415290-21416967  | 56                | 1677                   | 1.88                         | *KCNJ12* (ENST00000583088.6)  | ENSE00002732101.3 |
| chr7:100953518-100954302 | 43                | 784                    | 1.82                         | *MUC3A* (ENST00000379458.9)   | ENSE00003733255.1 |
| chr7:100951864-100952631 | 40                | 767                    | 1.76                         | *MUC3A* (ENST00000379458.9)   | ENSE00003733255.1 |
| chr3:75665550-75666394   | 34                | 844                    | 1.6                          | *FRG2C* (ENST00000308062.8)   | ENSE00003920067.1 |
| chr19:54631530-54632196  | 27                | 666                    | 1.45                         | *LILRB1* (ENST00000324602.12) | ENSE00003671569.1 |
| chr19:54631530-54632196  | 27                | 666                    | 1.45                         | *LILRB1* (ENST00000324602.12) | ENSE00002470335.1 |
| chr7:100956990-100957266 | 21                | 276                    | 1.42                         | *MUC3A* (ENST00000379458.9)   | ENSE00003733255.1 |
| chr7:142773266-142773518 | 19                | 252                    | 1.36                         | *PRSS2* (ENST00000539842.6)   | ENSE00003783126.1 |
| chr10:46549129-46550723  | 28                | 1594                   | 1.29                         | *GPRIN2* (ENST00000374314.6)  | ENSE00003923399.1 |
| chr19:54633043-54633269  | 16                | 226                    | 1.23                         | *LILRB1* (ENST00000324602.12) | ENSE00002482998.1 |
| chr19:54278966-54279050  | 12                | 84                     | 1.2                          | *LILRB2* (ENST00000314446.10) | ENSE00003604374.1 |
| chr19:39885702-39886439  | 20                | 737                    | 1.17                         | *FCGBP* (ENST00000616721.6)   | ENSE00003736302.1 |
| chr19:39885702-39886439  | 20                | 737                    | 1.17                         | *FCGBP* (ENST00000616721.6)   | ENSE00003726191.1 |
| chr19:54574781-54575019  | 15                | 238                    | 1.16                         | *LILRA2* (ENST00000391738.8)  | ENSE00003895381.1 |
| chr13:25097061-25097260  | 13                | 199                    | 1.08                         | *PABPC3* (ENST00000281589.5)  | ENSE00001425120.4 |
| chr7:100956319-100956539 | 13                | 220                    | 1.06                         | *MUC3A* (ENST00000379458.9)   | ENSE00003733255.1 |
| chr17:21703079-21703539  | 15                | 460                    | 1.02                         | *KCNJ18* (ENST00000567955.3)  | ENSE00002608500.3 |

#### Are HetExc intervals depleted for other class of variants?

# genomic-interval-enrichment-analysis-for-extreme-heterozygote-excess-HetExc-variants

## Intro

Deviations from Hardy-Weinberg Equilibrium (HWE) are not solely rooted in biological or genetic factors; sequencing/mapping errors play a significant role. Previous studies have underscored extreme heterozygote excess (HetExc) for specific variants in unstable genomic regions, suggesting the influence of genomic instability on observed deviations from HWE. Given the importance of variants occurring in coding regions for germline and somatic variant classification, we aimed at genomic intervals within the human exome enriched for HetExc variants. The identified regions could serve as indicators of potential spurious mutations, providing essential guidance for variant classification in genomic research.

## Methods

### Defining genomic intervals with HetExc variants

**Obtaining gencode.v44.annotation and processing** The file was downloaded from [GENCODE](https://www.gencodegenes.org/human/), and further processed to create a subset of exons from protein-coding Ensembl canonical transcripts. The processing details are outlined in the `gtf2bed.ipynb` notebook, in the scripts directory.

**Obtaining vcf files and processing**

We obtained exome and genome vcf files from gnomAD v4 (accessed on 2023-11-16) and processed them on HPC clusters provided by the Digital Research Alliance of Canada (the Alliance). The processing included (i) obtaining vcf files for chromosomes 1-22 and X, (ii) filtering variants to retain only HetExc variants (indicated by an inbreeding_coeff \<= -0.3), (iii) Intersecting VCF files with exon coordinates obtained from GENECODE, (iv) uniting proximal intervals (+/- 300 bp), then filtering united intervals with variant count less than two. The bash script used for exome vcf files processing can be found in the scripts directory under `gnomad_proc.sh` name.

**Annotating the identified intervals**

## Results

### HetExc interval count

We ordered genes based on the number of HetExc intervals identified in them to see which genes have the highest level of HetExc intervals. The table below represents our findings, where it can be seen that *MUC3A* is the top one, harboring a total of eight intervals residing in four exons.

| Gene Symbol      | Transcript ID      | HetExc interval count | Affected exon IDs                                                                                                |
|---------------|---------------|---------------|-----------------------------|
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

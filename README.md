### OVERVIEW

STRsearch is an end-to-end pipeline for targeted profiling of short tandem repeats (STRs) in massively parallel sequencing data. It is implemented using Python, supporting both version Python2 and Python3. 



### Installation

To obtain STRsearch, use:

git clone https://github.com/AnJingwd/STRsearch.git



### PREREQUISITE

​      The following linux utilities are needed and the full path of  them on your  local machine  should be provided in conf.py file

1. bwa  (v1.7 or higher) 

2. samtools  (v1.7 or higher) 

3. bamToFastq (v2.17.0 or higher)

4. seqtk  (v1.2 or higher)

5. usearch (v11 or higher)

   

### BASIC USAGE



### Configuration file format:

​      The first step for  STR analysis with STRsearch is to create a configuration file with your custom set of STR loci. One way to do this is by referring to the most up-to-date revised forensic STR sequence guide  and a worksheet can be downloaded from [link](https://strider.online/bundles/strbaseclient/downloads/Forensic_STR_Sequence_Structure_Guide_v5.xlsx). You will need to make a configuration file with the following columns present: 

| Chr   | Start     | End       | Period | Reference allele | Marker  | STR       | STR sequence structure                                       | Stand | 5' Flanking  sequence | 3' Flanking  sequence |
| ----- | --------- | --------- | ------ | ---------------- | ------- | --------- | ------------------------------------------------------------ | ----- | --------------------- | --------------------- |
| chr1  | 7442891   | 7442934   | 4      | 11               | Marker1 | D1GATA113 | [GATA]n                                                      | +     | ACTTGCTTCCTAGAT       | TTCCTATAGCCTCAA       |
| chr21 | 20554291  | 20554417  | 4      | 29               | Marker2 | D21S11    | [TCTA]n  [TCTG]n [TCTA]n ta [TCTA]n tca [TCTA]n tccata [TCTA]n TA [TCTA]n | +     | CCAAGTGAATTGCCT       | TCGTCTATCTATCCA       |
| chrX  | 149710971 | 149711038 | 4      | 15               | Marker3 | DXS7423   | [TGGA]n aggacaga  [TGGA]n                                    | +     | AAATGAATGAGTATG       | TGGGGAGGAAATCTG       |
| chrY  | 15752608  | 15752715  | 3      | 27               | Marker4 | DYS612    | [CCT]n CTT [TCT]n  CCT [TCT]n                                | +     | AGGTTCAGAGGTTTG       | GTCACTTTTCCAAAT       |
| chrY  | 20842518  | 20842573  | 4      | 14               | Marker5 | DYS385a   | [TTTC]n                                                      | -     | TCCTTTCTTTTTCTC       | CCTTCCTTCCTTCCT       |

1. Column 1 : chromosome  (must)
2. Column 2 : start coordinate of the STR  (must)
3. Column 3  : end coordinate of the STR  (must)
4. Column 4  : period of the STR  (must)
5. Column 5 : reference copy number (option)
6. Column 6 : Marker name (option) 
7. Column 7 : STR name (option) 
8. Column 8 :   Reference Sequence repeat region sequence structure summary  (must)
9. Column 9 : stand ("+" means positive stand;"-" means negative stand)  (must)
10. Column 10 : 5' flanking sequence of repeat region  (must)
11. Column 11 : 3' flanking sequence of repeat region  (must)

Note some columns are not used. You can put any value in the non-required columns, just make sure there are at least 11 columns with the required information listed above. Importanly, flanking sequences are necessarily adjacent to STR repeat region.  



 ### Algorithm description: 

​         Briefly, STRsearch employs an iterative algorithm to obtain the longest continuous interval composed by all motifs of STR sequence structure without a priori assumptions on allele size. The actual STR region is determined by comparing the position of repeat patterns with the best matching location of flanking sequences in reads. Ultimately, allele size is calculated not only for repeat patterns, but also indels that are actually in the STR region.   



### OPTIONS

| Option | Value Type | Default | Summary |
| ------ | ---------- | ------- | ------- |
|        |            |         |         |
|        |            |         |         |
|        |            |         |         |
|        |            |         |         |
|        |            |         |         |
|        |            |         |         |
|        |            |         |         |



### REFERENCE

 The STRsearch publication is available here: 



### CONTACT

Developer:  Dong Wang


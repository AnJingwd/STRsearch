#!/usr/bin/env python
# -*- coding: utf-8 -*-

# @Time    : 2019/10/12 15:23
# @Author  : Dong Wang
# @FileName: STR_parse.py
# @Software: STRsearch
# @github    ：https://github.com/AnJingwd/STRsearch


from __future__ import division
import os,re
from . import fastq
import numpy as np


def get_allele(input_file, full_seq):
    with open(input_file, "r") as f1:
        for line in f1:
            line = line.strip()
            if line.split("\t")[0].strip() == full_seq.strip():
                allele = line.split("\t")[1]
                return allele
                break


def str_info_count(input_file,full_seq):
    distance_left,distance_right = [],[]
    STR_origin = open(input_file,"r")
    next(STR_origin)
    for line in STR_origin:
        line = line.strip()
        content_list = line.split("\t")
        if content_list[0].strip() != full_seq.strip():
            pass
        else:
            distance_left.append(int(content_list[3].split(",")[0]))
            distance_right.append(int(content_list[3].split(",")[1]))
    distance_mean_5,distance_mean_3 = str(int(np.argmax(np.bincount(distance_left)))),str(int(np.argmax(np.bincount(distance_right))))
    distance_min_5,distance_max_5 = str(min(distance_left)),str(max(distance_left))
    distance_min_3,distance_max_3 =str(min(distance_right)),str(max(distance_right))
    genotypes_info_stat = [distance_min_5,distance_max_5,distance_min_3,distance_max_3,distance_mean_5,distance_mean_3]
    
    return genotypes_info_stat
        

def all_allele_sorted_dict(input_file):
    allele_reads = {}
    STR_origin = open(input_file,"r")
    next(STR_origin)
    for line in STR_origin:
        line = line.strip()
        content_list = line.split("\t")
        full_seq = content_list[0]
        if full_seq not in allele_reads.keys():
            allele_reads[full_seq]  = 1
        else:
            allele_reads[full_seq] +=1
    support_reads_sorted = sorted(allele_reads.items(),key=lambda item:item[1],reverse=True)
    return support_reads_sorted


def STR_parse(sample,sex,pos,marker_name,input_file,marker,position,pattern,reads_threshold,stutter_ratio):
    support_reads_sorted = all_allele_sorted_dict(input_file)
    if len(support_reads_sorted)==1:
        full_seq_1 = full_seq_2 = support_reads_sorted[0][0]
        allele1 = allele2= get_allele(input_file,full_seq_1)
        allele1_reads = allele2_reads = str(support_reads_sorted[0][1])
    else:
        full_seq_1,full_seq_2 = support_reads_sorted[0][0],support_reads_sorted[1][0]
        allele1= get_allele(input_file,full_seq_1)
        allele2= get_allele(input_file,full_seq_2)
        allele1_reads,allele2_reads =str(support_reads_sorted[0][1]),str(support_reads_sorted[1][1])
    merge = [str(i) + ", "+ str(j) for i, j in zip([allele1,allele1_reads,full_seq_1],[allele2,allele2_reads,full_seq_2])]
    if int(allele1_reads)+int(allele2_reads)< int(reads_threshold):   ####  cov cut off
        allele_adjust_info = "allele reads too low"
    elif int(allele2_reads)/int(allele1_reads)<float(stutter_ratio):    ## stutter cut off
        allele_adjust_info = allele1 + ", " + allele1
    else:
        allele_adjust_info = "NA"
    result_part1 = [sample,marker_name,marker,position,pattern]
    if sex =="male":
        if pos.split("-")[0]=="chrX" or pos.split("-")[0]=="chrY":
            info = "NA" if int(allele1_reads)>reads_threshold else "allele reads too low"
            result_part2 = [allele1,info,allele1_reads,full_seq_1]
        else:
            merge.insert(1,allele_adjust_info)
            result_part2 = merge
    else:
        merge.insert(1,allele_adjust_info)
        result_part2 = merge
    result_list =result_part1+result_part2
    return result_list,support_reads_sorted


def count_lines(file):
    '''
    count number of lines in txt file
    '''
    myfile = open(file,"r")  
    lines = len(myfile.readlines())  
    return int(lines)


def qual_stat(qstr):
    q20 = 0
    q30 = 0
    for q in qstr:
        qual = ord(q) - 33
        if qual >= 30:
            q30 += 1
            q20 += 1
        elif qual >= 20:
            q20 += 1
    return q20, q30


def stat(fq_file):
    reader = fastq.Reader(fq_file)
    total_count = 0
    q20_count = 0
    q30_count = 0
    n=0
    while True:
        n+=1
        read = reader.nextRead()
        if read == None:
            break
        total_count += len(read[3])
        q20, q30 = qual_stat(read[3])
        q20_count += q20
        q30_count += q30
    q20_percents=100 * float(q20_count)/float(total_count+0.00001)
    q30_percents=100 * float(q30_count)/float(total_count+0.00001)
    total_bases,N_reads,Q20,Q30 = str(round(total_count/1000000,2)),str(n-1),str(round(q20_percents,2)),str(round(q30_percents,2))
    fastq_info_stat = [total_bases,N_reads,Q20,Q30]
    return fastq_info_stat




def write_all_alleles(prefix,str_file,support_reads_sorted,result_file):
    with open(result_file, "a") as fx:
        reads_total = sum(list(zip(*support_reads_sorted))[1])
        for (nomenclature,reads) in support_reads_sorted:
            if float(reads/reads_total)<0.01:
                pass
            else:
                allele= get_allele(str_file,nomenclature)
                if allele:
                    all_allele_results ="\t".join( prefix+[nomenclature,allele,str(reads)])
                else:
                    all_allele_results = prefix+"\tNA\tNA\tNA"
                fx.write(all_allele_results+"\n")
        

def write_qc_matrix(prefix,fq_file,str_file,STR_results_list,result_file):
    fastq_info_stat=stat(fq_file)
    
    genotypes_coverage,genotypes_seq, = STR_results_list[7],STR_results_list[8]
    allele_origin = STR_results_list[5].split(", ")
    genotype_allele = allele_origin if len(STR_results_list[6].split(", "))!=2 else STR_results_list[6].split(", ")
    if genotypes_seq =="NA":
        genotypes_info_stat=["NA"]*12
        genotypes_coverage_new,genotype_allele_new = ["NA"]*2,["NA"]*2
    elif len(genotypes_seq.split(", "))==1:
        genotypes_info_stat = str_info_count(str_file,genotypes_seq)+["NA"]*6
        genotypes_coverage_new,genotype_allele_new = [genotypes_coverage,"NA"],genotype_allele+["NA"]
    else:
        genotype_seq_1,genotype_seq_2 = genotypes_seq.split(", ")[0],genotypes_seq.split(", ")[1]
        genotypes_info_stat= str_info_count(str_file,genotype_seq_1)+str_info_count(str_file,genotype_seq_2)
        genotypes_coverage_new,genotype_allele_new = genotypes_coverage.split(", "),genotype_allele
    qc_matrix_list = prefix+fastq_info_stat+genotypes_info_stat+genotypes_coverage_new+genotype_allele_new
    
    with open(result_file, "a") as ft:
        ft.write("\t".join(qc_matrix_list)+"\n")
    

def main(sample,sex,fastq_dir,STRsearch_dir,ref_bed,genotypes,multiple_alleles,qc_matrix,reads_threshold,stutter_ratio):

    ### write header
    with open(genotypes, "w") as f1:
        f1.write("Sample\tMarker\tSTR\tPosition\tSTR sequence sturucture\tAlleles (a1, a2)\tAlleles correction (a1, a2)\tSupporting reads (a1, a2)\tAllele sequences (a1, a2)\n")
    with open(multiple_alleles, "w") as f2:
        f2.write("Sample\tMarker\tSTR\tAllele sequences\tAllele\tSupporting reads\n")
    with open(qc_matrix, "w") as f3:
            f3.write("Sample\tMarker\tSTR\tTotal_bases\tNum_reads\tQ20\tQ30\tDis1_min_5\tDis1_max_5\tDis1_min_3\tDis1_max_3\tDis1_mean_5\tDis1_mean_3\tDis2_min_5\tDis2_max_5\tDis2_min_3\tDis2_max_3\tDis2_mean_5\tDis2_mean_3\tSupp_reads1\tSupp_reads2\tAllele1\tAllele2\n")

    bed = open(ref_bed,"r")
    next(bed)
    for line in bed:
        content_list = line.strip().split("\t")
        marker_name,official_name = content_list[5],content_list[6]
        chr,pos,nomenclature = content_list[0],content_list[0] +": "+content_list[1] +"-"+content_list[2],content_list[7]
        
        prefix = [sample,marker_name,official_name]
        str_file = os.path.join(STRsearch_dir,marker_name+"_results_"+sample+".txt")
        fq_file = os.path.join(fastq_dir,marker_name+"_reads_"+sample+"_merge.fastq")

        if count_lines(str_file)<=1:
            STR_results_list = prefix +[pos,nomenclature]+["NA"]*4
        else:
            STR_results_list,support_reads_sorted = STR_parse(sample,sex,pos,marker_name,str_file,official_name,pos,nomenclature,reads_threshold,stutter_ratio)
            ## write all alleles results
            write_all_alleles(prefix,str_file,support_reads_sorted,multiple_alleles)
        
        with open(genotypes, "a") as f:
            f.write(re.sub("NA","-","\t".join(STR_results_list))+"\n")
        ## write qc_matrix results
        write_qc_matrix(prefix,fq_file,str_file,STR_results_list,qc_matrix)        

    print("STR parse finished !")


if __name__=='__main__':
    sys.exit(main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6],sys.argv[7], sys.argv[8],sys.argv[9], sys.argv[10]))

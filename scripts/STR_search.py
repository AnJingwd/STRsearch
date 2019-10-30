#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# @Time    : 2019/10/12 15:23
# @Author  : Dong Wang
# @FileName: STR_search.py
# @Software: STRsearch
# @github    ：https://github.com/AnJingwd/STRsearch

from __future__ import division
import os,re
import numpy as np
from functools import reduce
from multiprocessing import Pool as Pool


def STR_nomenclature_trans(period,STR_Repeat_Motif):
    '''
    Transfor STR motif to STR unit list ,
    if stand is "-", get Inverse complementary sequence of STR unit list,
    And record which unit needs to be counted for STR allele using a list.

    input: like [TCTA]1TCTG[TCTA]14
    output: like (['TCTA', 'TCTG', 'TCTA'], [1, 0, 1])
    '''
    STR_units = re.findall(r"([a-zA-Z]+)", re.sub("n","",STR_Repeat_Motif))
    repetitive_motifs = [1 if re.findall(r'[[](.*?)[]]',i) else 0 for i in STR_Repeat_Motif.split(" ")]
    STR_numbers = []
    non_period_unit_len = [len(m) for m in re.findall(r"([A-Z]+)", STR_Repeat_Motif) if len(m) < period]
    if sum(non_period_unit_len) == period:
        weight = 1/len(non_period_unit_len)
    elif sum(non_period_unit_len) > period:
        weight = (sum(non_period_unit_len) // period + sum(non_period_unit_len) % period * 0.1) / len(
            non_period_unit_len)
    else:
        weight = 0.1*sum(non_period_unit_len)

    for i in STR_units:
        if re.match(r"[a-z]+", i):
            STR_numbers.append(0)
        else:
            if len(i) < period:
                STR_numbers.append(weight)
            else:
                STR_numbers.append(1)
    STR_units_upper = [j.upper() for j in STR_units]
    return STR_units_upper, repetitive_motifs,STR_numbers

def find_lcseque(seq,STR_unit):
    '''
     seach STR unit in reads a sequence using Longest common subsequence algorithm
     input: 
           seq: a sequence
           STR_unit: the sequence of a STR unit
    output:
            local_max_list：local max score list
            p_list: position list   
    '''
    m  = np.zeros([len(seq) + 1,len(STR_unit) + 1])
    local_max = 0  
    p = 0  
    local_max_list = []
    p_list = []
    for i in range(len(seq)):
        for j in range(len(STR_unit)):
            if seq[i] == STR_unit[j]:
                m[i + 1][j + 1] = m[i][j] + 1
                if m[i + 1][j + 1] == len(STR_unit): 
                    mmax = m[i + 1][j + 1]
                    p = i + 1
                    local_max_list.append(int(mmax))
                    p_list.append(int(p))
    return (local_max_list,p_list)


def match_flank(seq,flank):
    '''
    match flanking sequences using dynamic programming, and output the  location of the least mismatched
    input:
          seq: reads sequnece
          flank: flank sequence of 15 bp
    output: the  location of the least mismatched and the number of mismatched
    '''
    length = len(flank)
    resultMissmatchCount = length
    seqdict = {}
    for index, s in enumerate(seq[:-length]):
        missmatch = 0
        for j, k in zip(seq[index:index + length], flank):  # [(s1[0],s2[0]),(s1[1],s2[1]),...]
            if j != k:
                missmatch += 1
        if missmatch <= resultMissmatchCount:
            seqdict[missmatch] = seq[index:index + length]
            resultMissmatchCount = missmatch
    minkey = min(seqdict.keys())
    result = seqdict[minkey]
    start,end = seq.index(result),seq.index(result)+length
    return start,end,resultMissmatchCount

def merge_intervals(intervals):
    ''''
    Merge interval list
    input: like [[1,3],[2,6],[8,10],[15,18]]
    output: like [[1,6],[8,10],[15,18]]
    '''
    intervals.sort(key=lambda x: x[0])
    merged = []
    for interval in intervals:
        if not merged or merged[-1][-1] < interval[0]:
            merged.append(interval)
        else:
            merged[-1][-1] = max(merged[-1][-1], interval[-1])
    return merged

def get_STR_unit_region(score_list,pos_list,start_point = 0):
    '''
    Compute the intervals of STR unit and union
    input: 
    output:
    '''
    intervals_list = []
    for i in range(0,len(score_list)):
        start = pos_list[i]-score_list[i] + start_point
        end = pos_list[i] + start_point
        intervals_list.append([start,end])
    if score_list == []:
        return []
    elif score_list[0]==2:  ## when a unit is part of other unit ,don’t merge intervals
        return intervals_list
    else:
        intervals_union = merge_intervals(intervals_list)  ## get union intervals of STR units 
        return intervals_union

def list_depth(items):
    '''
    Calculate the depth of nested lists
    input: list
    output: max depth of list
    '''
    max_depth = 1 if isinstance(items, list) else 0
    if max_depth:
        for item in items:
            if isinstance(item, list):
                max_depth = max(max_depth, list_depth(item) + 1)
    else:
        return max_depth
    return max_depth

def get_interval_max(interval_list):
    if list_depth(interval_list) == 1:  ## like [1,4]
        return interval_list
    else:
        if len(interval_list)==1:  # depth=2 and len = 1 like [[-1,-1,-1,-1]]
            return interval_list[0]
        else:   ## depth=2 and len > 1 like [[-1,-1,-1,-1],[-1,-1,3,10]]
            max_length = 0
            max_index = 0
            for i in range(0, len(interval_list)):
                interval_new = [j for j in interval_list[i] if j != -1]
                if interval_new == []:
                    pass
                else:
                    (start, end) = min(interval_new), max(interval_new)
                    length = end - start
                    if length > max_length:
                        max_length = length
                        max_index = i
            return interval_list[max_index]

def get_continue_region(region1,region2):
    ## all region [)
    if region1==region2 == []:
        return [[-1,-1,-1,-1]]
    elif region1== [] or region2 == []:
        return [[-1,-1]+[k[0],k[1]] for k in region2] if region1== [] else [get_interval_max(region1)+[-1, -1]]
    else:
        continue_region = []
        for i in range(0,len(region1)):
                for j in range(0,len(region2)):
                        if region1[i][-1] == region2[j][0]:
                                region = [k for k in region1[i]]+[m for m in region2[j]]
                                continue_region.append(region)
                        elif region1[i][-1] > region2[j][0] and region1[i][-1] < region2[j][1]:
                                region = [k for k in region1[i]]+[region1[i][-1],region2[j][1]]
                                continue_region.append(region)
                        else:
                                pass
        continue_region.append(get_interval_max(region1) +[-1,-1])
        continue_region +=[[-1,-1]*(int(len(region1[0])/2))+[k[0],k[1]] for k in region2]
        return continue_region


def get_STR_seq(max_order_region, STR_units_upper,repetitive_motifs,STR_numbers,seq,period,flank_left_region,flank_right_region):

    mySTR_start,mySTR_end = int(min([i for i in max_order_region if i!=-1])),int(max(max_order_region))
    start_index, end_index = int(max_order_region.index(mySTR_start)/2),int(max_order_region.index(mySTR_end)/2)

    max_order_region_len = [max_order_region[i+1]- max_order_region[i] for i in range(0,len(max_order_region),2)]
    max_order_region_len_sub = max_order_region_len[start_index:end_index+1]
    
    STR_units_sub = STR_units_upper[start_index:end_index+1]
    STR_units_sub_len = [len(j) for j in STR_units_sub]
    
    repetitive_motifs_sub  = repetitive_motifs[start_index:end_index+1]
    STR_numbers_sub = STR_numbers[start_index:end_index+1]
    
    STR_unit_rep_count = np.array(max_order_region_len_sub)/np.array(STR_units_sub_len)*np.array(STR_numbers_sub)
    STR_unit_rep_count_new = [int(STR_unit_rep_count[i]) if repetitive_motifs_sub[i]==1 else STR_unit_rep_count[i] for i in range(len(repetitive_motifs_sub))]
    
    rep_num = sum(STR_unit_rep_count_new)
    STR_to_flank_distance = len(seq[flank_left_region[1]:mySTR_start]) + len(seq[mySTR_end:flank_right_region[0]])
    allele_sum = rep_num+(STR_to_flank_distance)//period+(STR_to_flank_distance)%period*0.1
    allele = int(allele_sum) if (allele_sum-int(allele_sum))<0.1 else allele_sum

    STR_result_pat = ""
    for k in range(0,len(STR_units_sub)):
        if STR_numbers_sub[k] == 1:
            if repetitive_motifs_sub[k]==1:
                unit = "[" + STR_units_sub[k] + "]"+str(STR_unit_rep_count_new[k])
            else:
                unit = STR_units_sub[k]
        elif STR_numbers_sub[k] == 0:
            unit = STR_units_sub[k].lower()
        else:
            unit = STR_units_sub[k]
        STR_result_pat +=unit+" "

    STR_result_pat_full = seq[flank_left_region[1]:mySTR_start].upper()+" " + STR_result_pat +seq[mySTR_end:flank_right_region[0]].upper()
    flank_mismatchs = str(flank_left_region[2])+","+str(flank_right_region[2])
    distances = str(mySTR_start)+","+str(len(seq)-mySTR_end)
    return STR_result_pat_full,allele,max_order_region_len,rep_num,flank_mismatchs,distances

def STR_search(STR_Repeat_Motif,period,seq,flank_seq_left,flank_seq_right):
    flank_left_region,flank_right_region = match_flank(seq,flank_seq_left),match_flank(seq,flank_seq_right)
    STR_units_upper, repetitive_motifs,STR_numbers= STR_nomenclature_trans(period,STR_Repeat_Motif)

    region_list_all = []
    for STR_unit in STR_units_upper:
        (local_max_score_list,p_list_left)  = find_lcseque(seq,STR_unit)
        region_list = get_STR_unit_region(local_max_score_list, p_list_left, 0)
        region_list_all.append(region_list)

    region_list_all_not_empty = [i for i in region_list_all if i !=[]]
    if len(region_list_all_not_empty)==0: # if all units mismatched, the read will be skipped
        return -1
    else:
        max_order_region_list = reduce(get_continue_region,region_list_all)
        max_order_region_list_new = [j for j in max_order_region_list if len(set(j))!=1]
        max_order_region= get_interval_max(max_order_region_list_new)
        (STR_result_pat_full,allele,max_order_region_len,rep_num,flank_mismatchs,distances) = get_STR_seq(max_order_region, STR_units_upper,repetitive_motifs,STR_numbers, seq,period,flank_left_region,flank_right_region)
        if sum(list(map(lambda a,b:a*b,repetitive_motifs,max_order_region_len))) ==0:
            return -1
        else:
            return (STR_result_pat_full,allele,rep_num,flank_mismatchs,distances)

def STR_search_one_locus(args_list):
    marker_name,STR_Repeat_Motif,period,STR_fastq_file,result_file,flank_seq_left,flank_seq_right= args_list
    results = open(result_file,"w")
    results.write("STR_sequence_structure\tAllele\tFlank_mismatchs(5',3')\tDistance_to_reads_ends(5',3')\n")
    N_line = 0
    with open(STR_fastq_file,"r") as f1:
        for line in f1:
            N_line+=1
            if N_line % 4 != 2:
                pass
            else:
                myseq = line.strip()
                res =STR_search(STR_Repeat_Motif,period,myseq,flank_seq_left,flank_seq_right)
                if res == -1:
                    continue
                else:
                    (STR_result_pat_full,allele,rep_num,flank_mismatchs,distances) =res
                    if rep_num<=1:
                        pass
                    else:
                        results.write(str(STR_result_pat_full)+"\t"+ str(allele) +"\t"+flank_mismatchs+ "\t"+distances+"\n")
    results.close()
    print("{} has been decoded completely! ".format(marker_name))


def count_reads(fq_file):
    '''
    count number of reads from fastq file
    '''
    count = 0
    with open(fq_file,"r") as f1:
        for line in f1:
            count+=1
    return int(count/4)



###################  STR genotype for locus listed in bed file one by one  #####################

def main(working_path,sample,sex,fastq_dir,ref_bed,reads_threshold,num_processors):
    ## create STR genotyping path
    STR_results_dir=os.path.join(working_path,"STRsearch")
    if not os.path.exists(STR_results_dir):
        os.makedirs(STR_results_dir)

    info_list = []
    bed = open(ref_bed,"r")
    next(bed)
    N = 0
    for line in bed:
        line = line.strip()
        content_list = line.split("\t")
        chr = content_list[0]
        period = int(content_list[3])
        marker_name = content_list[5]
        STR_Repeat_Motif = content_list[7]
        stand = content_list[-3]
        flank_seq_left,flank_seq_right = content_list[-2],content_list[-1]
        STR_fastq_merge_file = os.path.join(fastq_dir, marker_name +"_reads_"+sample+"_merge.fastq")
        result_file = os.path.join(STR_results_dir,marker_name+"_results_{}.txt".format(sample))
        if count_reads(STR_fastq_merge_file)<= int(reads_threshold): # if reads < 30 in fastq, the STR locus won't be genotyped , but create empty file
            os.system("touch {0}".format(result_file))
            print("{} has been decoded completely!".format(marker_name))
        else:
            if sex =="female" and chr == "chrY":
                os.system("touch {0}".format(result_file))
                print("{} has been decoded completely!".format(marker_name))
            else:
                info_list.append([marker_name,STR_Repeat_Motif,period,STR_fastq_merge_file,result_file,flank_seq_left,flank_seq_right])
        N+=1
    bed.close()

    pool = Pool(num_processors)
    pool.imap(STR_search_one_locus, info_list)
    pool.close()
    pool.join()
    print('Finished searching for total {} STR locus'.format(N))


if __name__=='__main__':
    sys.exit(main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6],sys.argv[7]))



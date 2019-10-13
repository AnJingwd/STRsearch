#!/usr/bin/env python
# -*- coding: utf-8 -*-

# @Time    : 2019/10/12 15:23
# @Author  : Dong Wang
# @FileName: get_STR_fastq.py
# @Software: STRsearch
# @github    ：https://github.com/AnJingwd/STRsearch

from __future__ import division
import os,re,argparse
from pathlib import Path

from multiprocessing import Pool
try:
    import configparser
except:
    from six.moves import configparser

parser = argparse.ArgumentParser(description='Get fastq file of each STR locus from bam file')
parser.add_argument('--work_dir', help='(must) assign working path for STR genotyping',required=True)
parser.add_argument('--sample', help='(must) assign the sample name',required=True)
parser.add_argument('--bam', help='(must) The input bam file',required=True)
parser.add_argument('--ref_bed', help='(must) The configuration file of STR',required=True)
parser.add_argument('--conf', help='(must) The arguments file of linux utilities',required=True)
parser.add_argument('--type', help='(option) The sequencing type (single or paired) (Default value:paired)',default="paired")
parser.add_argument('--merge_pairs', help='(option) if true, paired-end reads are merged (Default value:false)',default="false")
parser.add_argument('--num_processors', help='(option) The number of multiprocess (Default value:4)',type = int,default=4)
args = parser.parse_args()


def conf_parse(conf_file):
    cf = configparser.ConfigParser()
    cf.read(conf_file)
    conf_dict = {}
    for section in cf.sections():
        for (key,value) in cf.items(section):
            conf_dict[key] = cf.get(section, key)
    return conf_dict

config = conf_parse(args.conf)
mysamtools=eval(config["mysamtools"])
mybamToFastq=eval(config["mybamtofastq"])
myseqtk=eval(config["myseqtk"])
myusearch=eval(config["myusearch"])

def get_STR_fq_from_bam(args_list):
    '''
    (1)Get STR bam forward stand and reverse stand bam from  bam file
    (2)then transfor to fastq file using bamToFastq
    (3)get the inverse complementary sequence from reverse fastq file
    (4)combined forward and reverse stand fastq files 
    
    '''
    sample,type,pos,marker_name,bam_file,result_dir,stand,merge_pairs= args_list
    COMMAND_index="{0} index -b {1}"
    if not os.path.exists(bam_file+".bai"):
        os.system(COMMAND_index.format(mysamtools,bam_file))

    STR_bam = os.path.join(result_dir,marker_name+"_reads_"+sample+".bam")
    STR_bam_sort = os.path.join(result_dir,marker_name+"_reads_"+sample+"_sortByname.bam")
    fq_R1 = os.path.join(result_dir,marker_name+"_reads_"+sample+"_R1.fastq")
    fq_R2 = os.path.join(result_dir,marker_name+"_reads_"+sample+"_R2.fastq")
    merge_fq = os.path.join(result_dir,marker_name+"_reads_"+sample+"_merge.fastq")
    merge_fq_forward =os.path.join(result_dir,marker_name+"_reads_"+sample+"_merge_forward.fastq") 
    COMMAND_samtools = "{0} view -b1 {1} {2}>{3}"
    COMMAND_sort = "{0} sort -n {1}>{2}"
    COMMAND_bamToFastq = "{0} -i {1} -fq {2} -fq2 {3}"
    COMMAND_reverse = "{0} seq -r {1}>{2}"
    COMMAND_merge = "{0} -fastq_mergepairs {1} -reverse {2} -fastqout {3}"
	
    os.system(COMMAND_samtools.format(mysamtools,bam_file,pos,STR_bam))
    os.system(COMMAND_sort.format(mysamtools,STR_bam,STR_bam_sort))
    if type=="paired":
        os.system(COMMAND_bamToFastq.format(mybamToFastq,STR_bam_sort,fq_R1,fq_R2))
        if stand == "+":
            fq_R2_reverse = os.path.join(result_dir,marker_name+"_reads_"+sample+"_R2_trans.fastq")
            os.system(COMMAND_reverse.format(myseqtk,fq_R2,fq_R2_reverse))
            if merge_pairs =="false":
                os.system("cat {0} {1}>{2}".format(fq_R1,fq_R2_reverse,merge_fq))
            else:
                os.system(COMMAND_merge.format(myusearch,fq_R1,fq_R2,merge_fq))
        else:
            fq_R1_reverse = os.path.join(result_dir,marker_name+"_reads_"+sample+"_R1_trans.fastq")
            os.system(COMMAND_reverse.format(myseqtk,fq_R1,fq_R1_reverse))
            if merge_pairs =="false":
                os.system("cat {0} {1}>{2}".format(fq_R1_reverse,fq_R2,merge_fq))
            else:
                os.system(COMMAND_merge.format(myusearch,fq_R1,fq_R2,merge_fq_forward))
                os.system(COMMAND_reverse.format(myseqtk,merge_fq_forward,merge_fq))
    else:
        COMMAND_bamToFastq = "{0} -i {1} -fq {2}"
        os.system(COMMAND_bamToFastq.format(mybamToFastq,STR_bam_sort,fq_R1))
        if stand == "+":
            os.system("mv {0} {1}".format(fq_R1,merge_fq))
        else:
            os.system(COMMAND_reverse.format(myseqtk,fq_R1,fq_R1_reverse))
            os.system("mv {0} {1}".format(fq_R1_reverse,merge_fq))
    print("{}: finished getting STR fq from bam!".format(marker_name))

### check bam file ###
my_file = Path(args.bam)
try:
    bam = my_file.resolve()
except FileNotFoundError:
    print("Bam file doesn't exist!!")

### check fastq_dir ###
STR_fastq_dir=os.path.join(args.work_dir,"STRfq")
if not os.path.exists(STR_fastq_dir):
    os.makedirs(STR_fastq_dir)

info_list = []
bed = open(args.ref_bed,"r")
next(bed)
N = 0
for line in bed:
    line = line.strip()
    mylist = line.split("\t")
    pos = mylist[0]+":"+mylist[1]+"-"+mylist[2]
    marker_name,stand = mylist[5],mylist[8]
    info_list.append([args.sample,args.type,pos,marker_name,args.bam,STR_fastq_dir,stand,args.merge_pairs])
    N+=1

pool = Pool(args.num_processors)
pool.imap(get_STR_fq_from_bam, info_list)
pool.close()
pool.join()
print('Finished getting STR fq from bam for total {} markers!'.format(N))







#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# @Time    : 2019/10/12 15:23
# @Author  : Dong Wang
# @FileName: bwa_align.py
# @Software: STRsearch
# @github    ：https://github.com/AnJingwd/STRsearch


import os,sys,argparse
from pathlib import Path

try:
    import configparser
except:
    from six.moves import configparser


parser = argparse.ArgumentParser(description='Align fq with ref genome fasta for getting bam file')
parser.add_argument('--sample', help='(must) The sample name',required=True)
parser.add_argument('--fq1', help='(must) The input R1_fastq.gz file',required=True)
parser.add_argument('--fq2', help='(must) The input R2_fastq.gz file',required=True)
parser.add_argument('--runID', help='(must) provide the runID information for fq',required=True)
parser.add_argument('--lane', help='(must) provide the lane information for fq',required=True)
parser.add_argument('--ref', help='(must) The reference genome fasta file and index file in the same dir',required=True)
parser.add_argument('--output_dir', help='(must) Assign a path for output',required=True)
parser.add_argument('--conf', help='(must) The arguments file of linux utilities',required=True)
parser.add_argument('--num_threads', help='(option) The number of multiple threads (Default value:5)',type = int,default=5)
args = parser.parse_args()

## Assign softwares path
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
mybwa=eval(config["mybwa"])


############ Create output subdirectories  ############



sambam_dir = os.path.join(args.output_dir, "alignments")
if not os.path.exists(sambam_dir):
    os.makedirs(sambam_dir)

def check_file(file):
    my_file = Path(file)
    try:
        fq_file = my_file.resolve()
    except FileNotFoundError:
            print("fastq.gz file doesn't exist!!")    
    
def make_metadata_string(metadata):
    return r'"@RG\tID:%s\tSM:%s\tPL:%s"' % (metadata['ID'], metadata['SM'], metadata['PL'])


def bwaPE(sample,runID, lane, fq1,fq2,ref,output_dir,num_threads):
    """
    Aligns two paired-end fastq files to a reference genome to produce a sam file.
    """
    readgroup_metadata = { 'PL': 'ILLUMINA',
                       'SM': sample,
                       'ID': "%s_%s_Lane%d" % (sample, runID, int(lane)) }
    metadata_str = make_metadata_string(readgroup_metadata)
    sam_file = os.path.join(output_dir,sample+".sam")
    COMMAND_BWA = "{0} mem -t ${1} -M -R {2} {3} {4} {5}>{6}"
    os.system(COMMAND_BWA.format(mybwa,num_threads,metadata_str,ref,fq1,fq2,sam_file))
    
    """
    Convert sam to bam and sort, using Picard.
    """
    bam_file = os.path.join(output_dir,sample+".bam")
    COMMAND_sortTobam = "{0} sort {1}>{2}"
    COMMAND_index = "{0} index {1}"
    os.system(COMMAND_sortTobam.format(mysamtools,sam_file,bam_file))
    os.system(COMMAND_index.format(mysamtools,bam_file))
    os.system("rm {0}".format(sam_file))

check_file(args.fq1)
check_file(args.fq2)
bwaPE(args.sample,args.runID, args.lane, args.fq1,args.fq2,args.ref,sambam_dir,args.num_threads)

    
    







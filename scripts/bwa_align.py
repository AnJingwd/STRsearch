#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# @Time    : 2019/10/12 15:23
# @Author  : Dong Wang
# @FileName: bwa_align.py
# @Software: STRsearch
# @github    ：https://github.com/AnJingwd/STRsearch

import os,sys
from pathlib import Path


try:
    import configparser
except:
    from six.moves import configparser


## Assign softwares path
def conf_parse(conf_file):
        cf = configparser.ConfigParser()
        cf.read(conf_file)
        conf_dict = {}
        for section in cf.sections():
                for (key,value) in cf.items(section):
                        conf_dict[key] = cf.get(section, key)
        return conf_dict


def check_file(file):
    my_file = Path(file)
    try:
        fq_file = my_file.resolve()
    except FileNotFoundError:
            print("fastq.gz file doesn't exist!!")    
    
def make_metadata_string(metadata):
    return r'"@RG\tID:%s\tSM:%s\tPL:%s"' % (metadata['ID'], metadata['SM'], metadata['PL'])


def bwaPE(sample,fq1,fq2,runID, lane,ref,alignments_path,mysamtools,mybwa,num_threads):
    """
    Aligns two paired-end fastq files to a reference genome to produce a sam file.
    """
    readgroup_metadata = { 'PL': 'ILLUMINA',
                       'SM': sample,
                       'ID': "%s_%s_Lane%d" % (sample, runID, int(lane)) }
    metadata_str = make_metadata_string(readgroup_metadata)
    sam_file = os.path.join(alignments_path,sample+".sam")
    COMMAND_BWA = "{0} mem -t ${1} -M -R {2} {3} {4} {5}>{6}"
    os.system(COMMAND_BWA.format(mybwa,num_threads,metadata_str,ref,fq1,fq2,sam_file))
    
    """
    Convert sam to bam and sort, using Picard.
    """
    bam_file = os.path.join(alignments_path,sample+".bam")
    COMMAND_sortTobam = "{0} sort {1}>{2}"
    COMMAND_index = "{0} index {1}"
    os.system(COMMAND_sortTobam.format(mysamtools,sam_file,bam_file))
    os.system(COMMAND_index.format(mysamtools,bam_file))
    os.system("rm {0}".format(sam_file))

def main(sample,fq1,fq2,runID, lane,ref,working_path,num_threads):
    conf_path = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
    config = conf_parse(os.path.join(conf_path,"conf.py"))
    mysamtools =eval(config["mysamtools"])
    mybwa=eval(config["mybwa"])

    alignments_path = os.path.join(working_path, "alignments")
    if not os.path.exists(alignments_path):
        os.makedirs(alignments_path)
    check_file(fq1)
    check_file(fq2)
    bwaPE(sample,fq1,fq2,runID, lane,ref,alignments_path,mysamtools,mybwa,num_threads)


if __name__=='__main__':
    sys.exit(main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6],sys.argv[7], sys.argv[8]))

    
    

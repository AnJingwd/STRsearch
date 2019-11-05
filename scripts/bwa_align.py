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
    file_name = Path(file)
    if not file_name.exists():
        print("{0} doesn't exist!".format(file_name.name))
        sys.exit()
    else:
        if file_name.suffix not in [".fastq",".fq",".fastq.gz",".fq.gz"]:
            print("FASTQ file isn't correct!")
            sys.exit()
    

def bwa_mapping(sample,type,fq1,ref,alignments_path,mysamtools,mybwa,num_threads,fq2=None):
    '''
    Aligns fastq files to a reference genome to produce a sam file.
    '''
    sam_file = os.path.join(alignments_path,sample+".sam")
    if type =="single":
        COMMAND_BWA = "{0} mem -t {1} {2} {3}>{4}"
        os.system(COMMAND_BWA.format(mybwa,num_threads,ref,fq1,sam_file))
    else:
        COMMAND_BWA = "{0} mem -t {1} {2} {3} {4}>{5}"
        os.system(COMMAND_BWA.format(mybwa,num_threads,ref,fq1,fq2,sam_file))
    
    """
    Convert sam to bam and sort, using Picard.
    """
    bam_file = os.path.join(alignments_path,sample+".bam")
    COMMAND_sortTobam = "{0} sort {1}>{2}"
    COMMAND_index = "{0} index {1}"
    os.system(COMMAND_sortTobam.format(mysamtools,sam_file,bam_file))
    os.system(COMMAND_index.format(mysamtools,bam_file))
    os.system("rm {0}".format(sam_file))

def main(sample,type,fq1,ref,working_path,num_threads,fq2):
    conf_path = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
    config = conf_parse(os.path.join(conf_path,"conf.py"))
    mysamtools =eval(config["mysamtools"])
    mybwa=eval(config["mybwa"])

    alignments_path = os.path.join(working_path, "alignments")
    if not os.path.exists(alignments_path):
        os.makedirs(alignments_path)
    check_file(fq1)
    if fq2:
        check_file(fq2)
    bwa_mapping(sample,type,fq1,ref,alignments_path,mysamtools,mybwa,num_threads,fq2)


if __name__=='__main__':
    sys.exit(main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6],sys.argv[7]))

    
    

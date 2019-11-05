#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# @Time    : 2019/10/12 15:23
# @Author  : Dong Wang
# @FileName: pipeline.py
# @Software: STRsearch
# @github    ：https://github.com/AnJingwd/STRsearch

import os,sys,argparse
from scripts.bwa_align import main as bwa_align
from scripts.get_STR_fastq import main as get_STR_fastq
from scripts.STR_search import main as STR_search
from scripts.STR_parse import main as STR_parse


# create the top-level parser
parser = argparse.ArgumentParser(prog='STRsearch')

parser.add_argument('--type',type=str,default="paired",choices=['paired', 'single'], help='(option) The sequencing type (default: paired)')
parser.add_argument('--assemble_pairs', type =bool,default=False,choices=[False, True],help='(option) if True, paired-end reads are assembled (default: False)')
parser.add_argument('--reads_threshold',type = int,default=30, help='(option) The analytical threshold for reads(default: 30)')
parser.add_argument('--stutter_ratio',type=float,default=0.5,help='(option) The stutter ratio(default: 0.5)')
parser.add_argument('--num_threads',type=int, default=4,help='(option) The number of multiple threads (default: 4)')
parser.add_argument('--num_processors',type = int,default=4, help='(option) The number of multiprocess (default: 4)')

subparsers = parser.add_subparsers(title="commands", dest="command",help='start to run STRsearch with BAM-file or fastq file')

# create the parser for the "from_bam" command
parser_a = subparsers.add_parser('from_bam',help='run STRsearch with BAM-file directly')
parser_a.add_argument('--working_path', type=str,help='(must) The working path')
parser_a.add_argument('--sample',type=str,help='(must) The sample name')
parser_a.add_argument('--sex', type=str,choices=['male', 'female'],help='(must) The sample sex')
parser_a.add_argument('--bam',type=str, help='(must) The input BAM-file')
parser_a.add_argument('--ref_bed', type=str,help='(must) The configuration file of STRs')
parser_a.add_argument('--genotypes', help='(must) The output for STR genotypes')
parser_a.add_argument('--multiple_alleles', help='(must) The output for multiple alleles')
parser_a.add_argument('--qc_matrix', help='(must) The output for quality control matrix')


# create the parser for the "from_fastq" command
parser_b = subparsers.add_parser('from_fastq',help='map reads to create BAM-file firstly')
parser_b.add_argument('--working_path', type=str,help='(must) The working path')
parser_b.add_argument('--sample',type=str,help='(must) The sample name')
parser_b.add_argument('--fq1',type=str, help='(must) The input R1_fastq.gz file')
parser_b.add_argument('--fq2',type=str, help='(must) The input R2_fastq.gz file')
parser_b.add_argument('--runID',type=str, help='(must) The runID information')
parser_b.add_argument('--lane',type=str, help='(must) The lane information')
parser_b.add_argument('--ref',type=str, help='(must) The reference genome fasta and index file in the same path')
args = parser.parse_args()


if args.command == 'from_bam':
    print("from_bam")
elif args.command == 'from_fastq':
    print("from_fastq")
else:
    print("111")


if args.command == 'from_bam':
    bam_file = args.bam
elif args.command == 'from_fastq':
    bwa_align(args.sample,args.fq1,args.fq2,args.runID, args.lane,args.ref,args.working_path,args.num_threads)
    bam_file=os.path.join(args.working_path,"alignments",args.sample,".bam")
    sys.exit()
get_STR_fastq(args.sample,bam_file,args.working_path,args.ref_bed,args.type,args.assemble_pairs,args.num_processors)

fastq_dir = os.path.join(args.working_path,"STRfq")
STR_search(args.working_path,args.sample,args.sex,fastq_dir,args.ref_bed,args.reads_threshold,args.num_processors)

STRsearch_dir=os.path.join(args.working_path,"STRsearch")
STR_parse(args.sample,args.sex,fastq_dir,STRsearch_dir,args.ref_bed,args.genotypes,args.multiple_alleles,args.qc_matrix,args.reads_threshold,args.stutter_ratio)
    

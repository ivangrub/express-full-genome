#!/usr/bin/env python

from Bio import SeqIO
import pysam as pys
import numpy as np
import argparse
import os

def build_fasta(bw):
	"""Initialize a dictionary of each reference and the relevant sequence"""
	chr = {}
	for x in SeqIO.parse('%s/%s' %(bw,args.g),'fasta'):
		chr[x.id] = x.seq
	return chr 

def MakeSAMHeader(chrom,gen,path):
	"""Print a SAM header that will be used in the final SAM file that is outputted after the projection of the reads back onto the genome"""
	ref = []
	header = {}
	header['HD'] = {'VN':'1.0'}
	for i in chrom.keys():
		ref.append({'LN':len(chrom[i]), 'SN':i})
	header['SQ'] = ref
	samheader = pys.Samfile('%s/Post-eXpress_Header_%s.sam' % (path,gen),'wh',header = header)
	samheader.close()
	
def ReadLength(file):
	"""Determine the length of read outputted from the sequencing"""
	for x in SeqIO.parse(file,'fastq'):
		return len(x.seq)
	
def print_fasta(pos,seq,header,count):
	"""Print the new FASTA file that will be used to build the bowtie index and be the target sequence in eXpress"""
	
	binnum = pos.split('!')
	chr = binnum[1]
	coor1 = int(binnum[2])-1
	coor2 = int(binnum[3])-1
	s = seq[coor1:coor2].upper()
	if len(s) != coor2-coor1:
		m = len(seq[coor1:coor2])
		en = coor2-coor1
		s2 = s + 'N'*(en-m)
		header.write('@SQ\tSN:%s\tLN:%d\n' % (pos,int(args.b)-1))
		NEWFasta.write('>%s\n%s\n\n' % (pos, s2))
		count += 1
	else:
		header.write('@SQ\tSN:%s\tLN:%d\n' % (pos,int(args.b)-1))
		NEWFasta.write('>%s\n%s\n\n' % (pos, s))
		count += 1
	
	return count		
	
parser = argparse.ArgumentParser(description='Create a genome fasta file that will be used to create the bowtie index and used as an input eXpress.')
parser.add_argument('-g',help = 'Enter the name of the fasta file. The default is mm9.fa',default = 'mm9.fa')
parser.add_argument('-f',help = 'Input the fastq file that will be used by bowtie.',default = None)
parser.add_argument('-l',help = 'Input the length of the fasta read',default = None)
parser.add_argument('-b',help = 'Input the length of the sequence bin',default = 1000)
args = parser.parse_args()

if args.f is None and args.l is None:
	error("You need to input a fastq file or at least the length of the read")
	
# Start of the main function
genome = args.g.strip().split('.')
direct = os.environ['EXPRESS_FILES']
bowtie = os.environ['BOWTIE_INDEXES']

# Read in information about the chromosomes (Both sequence and length parameters)
if args.f is None:
	length = int(args.l)
else:
	length = ReadLength(args.f)

NEWFasta = open('%s/eXpress_%sbp_%d.%s.fa' % (direct,args.b,length,genome[0]),'w')

# Output the reference dictionary into chr and print a new SAM header that will be used to make the final SAM file
chr = build_fasta(bowtie)
MakeSAMHeader(chr,genome[0],direct)

# Bin the sequences into smaller fasta format sequences which will be used to create the bowtie index
k = 0
header = open('%s/Header_%s_%s_%d.sam' % (direct,args.g,args.b,length),'w')
header.write('@HD\tVN:1.0\n')
ref = []

chrindex = open('%s/chrindex_%s_%s_%d.txt' % (direct,args.g,args.b,length),'w')
for key in chr.keys():
	left = np.arange(1,len(chr[key]),int(args.b)-length)
	right = np.arange(int(args.b),len(chr[key])+int(args.b),int(args.b)-length)
	
	print 'On %s' % key
	chrindex.write('%s\t%d\n' % (key,k))
	for j in xrange(len(left)):
		bin = 'bin%d!%s!%d!%d' % (k,key,left[j],right[j])
		k = print_fasta(bin,chr[key],header,k)
		

chrindex.close()		
header.close()
NEWFasta.close()			
											

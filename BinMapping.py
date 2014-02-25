#!/usr/bin/env python

import pysam as pys
import argparse
import os

def gen2bin(infile,convbam,read,chrlist,bin,cc):
	#Find the appropriate bin and update the tid in the BAM file
	chrom = infile.getrname(read.tid)
	ind = int(round(read.pos/float(bin-1)))
	chrst = int(chrlist[chrom])
	
	try:
		name = convbam.getrname(chrst+ind)
		x = name.split('!')
		y = chrst+ind
	except ValueError:
		length = convbam.nreferences
		name = convbam.getrname(length -1)
		x = name.split('!')
		y = length-1
	
	off = read.pos - int(x[2])+1
	k = 1
	
	while chrom != x[1] or off < 0:
		name = convbam.getrname(y-k)
		x = name.split('!')
		off = read.pos - int(x[2]) + 1
		k += 1

	if off > int(int(args.b) - 1) or off < 0:
		print 'ERROR',read

	read.pos = off
	read.tid = convbam.gettid(name)

	if (cc % 1000000 == 0):
		print '%d alignments processed' % cc

	convbam.write(read)		

def offset(conv):
	# Establish the offset 
	n1 = conv.getrname(2)
	x = n1.split('!')
	start = int(x[2])
	n2 = conv.getrname(3)
	x = n2.split('!')
	end = int(x[2])
	off = end - start + 1
	return off

# Input paramaters
parser = argparse.ArgumentParser(description='Put mapped reads into BAM format in chromosomal and express coordinates')
parser.add_argument('-r',help ='Read file name. SAM/BAM format.',default = '')
parser.add_argument('-l',default=50)
parser.add_argument('-g',default='mm9.fa')
parser.add_argument('-o',default=None)
parser.add_argument('-b',default=200)
args = parser.parse_args()

path = os.environ['EXPRESS_FILES']
if args.r is '-':
	infile = pys.Samfile('-','r')
else:
	infile = pys.Samfile(args.r,'rb')

# Load binned header
conv = pys.Samfile('%s/Header_%s_%s_%s.sam' % (path,args.g,args.b,args.l),'r')

# Create template for converted BAM
if args.r is '-' and (args.o is '-' or args.o is None):
	convbam = pys.Samfile('-','wb',template = conv)
elif args.r is '-':
	convbam = pys.Samfile('%s_converted.bam' % args.o,'wb',template = conv)
else:
	convbam = pys.Samfile('%s_converted.bam' % args.o,'wb',template = conv)

# Building conversion headers
chrind = open('%s/chrindex_%s_%s_%s.txt' % (path,args.g,args.b,args.l),'r')
chrindex = {}
for line in chrind:
	s = line.strip().split('\t')
	chrindex[s[0]] = s[1]
chrind.close()

# Getting the offset
binning = offset(conv)

# Read through the input alignments
count = 0
for line in infile:
	if line.is_unmapped:
		continue
	gen2bin(infile,convbam,line,chrindex,binning,count)
	count += 1

convbam.close()
conv.close()

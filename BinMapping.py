#!/usr/bin/env python

import pysam as pys
import argparse
import os

def gen2bin(infile,convbam,read,chrlist,bin,headlist,cc):
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
	# indice = [-1,0,1]
	# for ran in indice:
	# 	try:
	# 		if (convbam.getrname(y-k+1+ran) not in headlist):
	# 			headlist.add(convbam.getrname(y-k+1+ran))
	# 			binheader.write('@SQ\tSN:%s\tLN:%d\n' % (convbam.getrname(y-k+1+ran),int(args.b)-1))
	# 		else:
	# 			continue
	# 	except ValueError:
	# 		continue

	if (cc % 1000000 == 0):
		print '%d reads done' % cc

	convbam.write(read)	
	return headlist	

def TrimBAM(conv,nheader):
	trimmed = pys.Samfile('%s_trimmed.bam' % args.o,'wb',template=nheader)
	for read in conv:
		binid = conv.getrname(read.tid)
		index = trimmed.gettid(binid)
		read.tid = index
		trimmed.write(read)
	trimmed.close()

def offset(conv):
	n1 = conv.getrname(2)
	x = n1.split('!')
	start = int(x[2])
	n2 = conv.getrname(3)
	x = n2.split('!')
	end = int(x[2])
	off = end - start + 1
	return off

parser = argparse.ArgumentParser(description='Put mapped reads into BAM format in chromosomal and express coordinates')
parser.add_argument('-r',help ='Read file name. SAM/BAM format.',default = '')
parser.add_argument('-l',default='50')
parser.add_argument('-g',default='mm9_norandom')
parser.add_argument('-o',default='outputted')
parser.add_argument('-b',default=1000)
args = parser.parse_args()

path = os.environ['EXPRESS_FILES']
if args.r is '-':
	infile = pys.Samfile('-','r')
else:
	infile = pys.Samfile(args.r,'rb')

print 'Load binned header'
conv = pys.Samfile('%s/Header_%s_%s_%s.sam' % (path,args.g,args.b,args.l),'r')
if args.r is '-':
	print 'Create template for converted BAM'
	convbam = pys.Samfile('-','wb',template = infile)
else:
	print 'Create template for converted BAM'
	convbam = pys.Samfile('%s_converted.bam' % args.o,'wb',template = conv)

#bindict = {}
#index = 0
#for ref in conv.header['SQ']:
#	bindict[ref["SN"]] = index
#	index += 1

print 'Building conversion headers'
chrind = open('%s/chrindex_%s_%s_%s.txt' % (path,args.g,args.b,args.l),'r')
chrindex = {}
for line in chrind:
	s = line.strip().split('\t')
	chrindex[s[0]] = s[1]
chrind.close()
print 'Getting the offset'
binning = offset(conv)

print 'Binning the alignments - Might take a long time'
binheader = open('newheaders.txt','w')
binheader.write('@HD\tVN:1.0\n')
headlist = set()
count = 0
for line in infile:
	if line.is_unmapped:
		continue
	# if args.r is '-':
	# 	convbam.write(line)
	headlist = gen2bin(infile,convbam,line,chrindex,binning,headlist,count)
	count += 1

convbam.close()
conv.close()
# conv = pys.Samfile('%s_converted.bam' % args.o,'rb')

# binheader.close()
# print 'Making new header'
# newheader = pys.Samfile('newheaders.txt','r')

# print 'Trimming and updating tid'
# TrimBAM(conv,newheader)
# conv.close()
# newheader.close()
# if args.r is '-':
# 	outbam.close()

# infile.close()




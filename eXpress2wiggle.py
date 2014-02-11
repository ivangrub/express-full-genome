#!/usr/bin/env python

# Convert Tommy Kaplan's Grizzly Peak (ChIP-Seq Peak calling program) from MATLAB into 
# python for integration with eXpress.

import numpy as np
import pysam as pys
import argparse
import glob
import pickle
import os

def GenomeLen():
	"""Open the new SAM header file and create a binned array of zeros for each chromosome"""
	
	global sep
	genome = {}
	path = os.environ['EXPRESS_FILES']
	samfile = pys.Samfile('%s/Post-eXpress_Header_%s.sam' % (path,args.g),'r')
	header = samfile.header['SQ']
	
	# Iterate through the reference sequences in the header
	for i in header:
		genome[i['SN']] = np.zeros(np.floor(i['LN']/int(args.w))+sep)
	
	samfile.close()
	return genome,header
	
def AddtoBin(read,ref):
	global chip, sep#, chipF, chipR
	
	# Check the strand
	if read.is_reverse:
		strand = '-'
	else:
		strand = '+'
	
	# Figure out which reference sequence the read mapped to
	
	b = ref[read.tid]
	bin = b.split('!')
	try:	
		chr = bin[1]
		start = int(read.pos) + int(bin[2])
	except IndexError:
		chr = b
		start = int(read.pos)
		
	st = int(np.floor(start/int(args.w)))
		
	# Get the percent likelihood value of the read mapping to that location
	#if 'hits.prob.bam' in args.r:
	try:
		frac_str = read.opt('XP')
		frac = float(frac_str)
	except KeyError:
		frac = 1.
		
	# Create a coordinate index to update the initial zeros vector array in the chip dictionary
	
	#if strand == '+':
	#	pp = np.arange(st,st+sep-1)
	#else:
	#	pp = np.arange(st-sep+1+int(np.floor(len(read.seq)/int(args.w))),st+int(np.floor(len(read.seq)/int(args.w))))
	
	# Check boundary conditions	
	if strand == '+':
		try:
			#chipF[chr][st:st+sep-1] += frac
			chip[chr][st:st+sep-1] += frac
		except IndexError:
			#chipF[chr][st:len(chipF[chr])-1] += frac
			chip[chr][st:len(chip[chr])-1] += frac
	else:
		try:
			#chipR[chr][st-sep+1+int(np.floor(len(read.seq)/int(args.w))):st+int(np.floor(len(read.seq)/int(args.w)))] += frac
			chip[chr][st-sep+1+int(np.floor(len(read.seq)/int(args.w))):st+int(np.floor(len(read.seq)/int(args.w)))] += frac
		except IndexError:
			#chipR[chr][0:st+int(np.floor(len(read.seq)/int(args.w)))] += frac
			chip[chr][0:st+int(np.floor(len(read.seq)/int(args.w)))] += frac		
	
	endread = AdjustRead(read,chr,start)
	return endread	
		
def NormalizeReads(nreads):
	global chip#, chipF, chipR
	# Normalize to 10 million reads
	ratio = nreads/10000000. #1e7
	for i in chip.keys():
		chip[i] = chip[i]/ratio	
		#chipF[i] = chipF[i]/ratio
		#chipR[i] = chipR[i]/ratio

def AdjustRead(read,chrom,st):
	"""Update or write a new SAM file"""
	global chrlist
	read.tid = chrlist.index(chrom)
	read.pos = st
	
	return read
		
def ReadReads(dir):
	"""Pass in the BAM/SAM file reads and project them back onto the genome"""	
	path = os.environ['EXPRESS_FILES']
	samfile = pys.Samfile('%s/Post-eXpress_Header_%s.sam' % (path,args.g),'r')
	
	if args.bo == 'y':
		outfile = pys.Samfile('%s/Post-eXpress-%s.%s.bam' %(dir,args.o,args.g),'wb',template = samfile)
	
	if args.r is not '-':
		if '.bam' in args.r:
			r = pys.Samfile(args.r,'rb')
		else:
			r = pys.Samfile(args.r,'r')
			
		# Create a tuple of the reference sequences in the file
		refdict = r.references
		
		nreads = 0
		for seq in r:	
			if seq.is_unmapped:
				continue
			read = AddtoBin(seq,refdict)
			if args.bo == 'y':
				outfile.write(read)
			
			# Count the number of mapped reads. Keep track of whether it is a multiply aligned read
			if nreads > 0:
				#print prev_seq.qname, seq.qname, nreads
				if prev_seq.qname != seq.qname:
					nreads += 1
				prev_seq = seq
			else:
				nreads += 1
				prev_seq = seq
			
			# Keep track of the reads
			if nreads % 1e6 == 0:
				print '\r%d tags read' % nreads
		r.close()
	else:
		# Create a tuple of the reference sequences in the file
		refdict = r.references
		nreads = 0
		for seq in sys.stdin:
			if seq.is_unmapped:
				continue
			read = AddtoBin(seq,refdict)
			if args.bo == 'y':
				outfile.write(read)
			
			# Count the number of mapped reads. Keep track of whether it is a multiply aligned read
			if nreads > 0:
				if prev_seq.qname != seq.qname:
					nreads += 1
				prev_seq = seq
			else:
				nreads += 1
				prev_seq = seq
			
			# Keep track of the reads
			if nreads % 1e6 == 0:
				print '\r%d tags read' % nreads
	NormalizeReads(nreads)
	if args.bo == 'y':
		outfile.close()
	samfile.close()
	return nreads 
	
def ParseHeader(hd):
	"""Parse the header list of dictionaries to get the chromosome name and length"""
	chr_list = []
	chr_len = []
	for i in hd:
		chr_list.append(i['SN'])
		chr_len.append(i['LN'])
	return chr_list,chr_len	
		
def PrintWiggle(reads,name,dir):
	"""Print a wiggle file of the projected reads so that it can be visualized on the UCSC genome browser"""
	global chrlist,chrlen,readcount
	
	print 'Printing the %s wiggle file' % name
	f = open('%s/%s.%s.%s.wig' % (dir,name,args.o,args.g),'w')
	bed = open('%s/%s.%s.%s.bedgraph' % (dir,name,args.o,args.g),'w')
	f.write('track type=bedGraph name="%s_%s" description="%s_%s" visibility=full graphType=bar\n' % (args.o,name,args.o,name))
	bed.write('track type=bedGraph name="%s_%s" description="%s_%s" visibility=full graphType=bar\n' % (args.o,name,args.o,name))
	for i in reads.keys():
		length = chrlen[chrlist.index(i)]
		lastpos = np.nan
		lastval = np.nan
		X = np.arange(int(args.w),int(args.w)*(len(reads[i])-17),int(args.w))
		dat = reads[i]
		for zz in xrange(len(X)):
			d = abs(dat[zz] - lastval)
			bed.write('%s %d %e\n' % (i,X[zz],dat[zz]*(readcount/10000000.)))
			if not np.isfinite(d) or d > 1e-1:
				if np.isfinite(lastval):
					if X[zz] > length:
						continue
					f.write('%s %d %d %.1f\n' % (i,lastpos,X[zz],lastval))
				lastval = dat[zz]
				lastpos = X[zz]
	f.close()
	bed.close()

parser = argparse.ArgumentParser(description='Project SAM/BAM format reads back onto the appropriate genome and call significant peaks.')
parser.add_argument('-r',help ='Read file name. SAM/BAM format.',default = '')
parser.add_argument('-g',help = 'The genome abbreviation. Mouse is mm9. Human is hg18. Drosophila is dm3. The default is mm9.',default = 'mm9')
parser.add_argument('-w',help = 'The base pair binning. The default = 50bp.',default = 50)
parser.add_argument('-d',help = 'The DNA fragment length. The default = 250bp.',default = 250)
parser.add_argument('-o',help = 'The prefix of the output files.',default = 'test')
parser.add_argument('-p',help = 'Add this to the command line if you want to call peaks using the Grizzly Peak algorithm.',default = None)
parser.add_argument('-bo',help = 'To output newly converted BAM file, (y/n)',default = 'n')
args = parser.parse_args()

# Establish express output directory and print new files there
dir = args.r.split('/')
dirname = '/'.join(dir[:-1])

# Binning
sep = int(np.ceil(int(args.d)/int(args.w)))

# Allow for the summed input along with specific strands
v = ['chip'] #,'chipF','chipR']
chip,samheader = GenomeLen()
#chipF, ign = GenomeLen()
#chipR, ign = GenomeLen()
chrlist,chrlen = ParseHeader(samheader)

if '%s/%s.%s.pkl'% (dirname,args.o,v[0]) not in glob.glob('%s/*pkl' % dirname):
	readcount = ReadReads(dirname)
	for i in v:
		output = open('%s/%s.%s.pkl'% (dirname,args.o,i),'wb')
		pickle.dump(eval(i),output)
		output.close()
else:
	for i in v:
		pkl_file = open('%s/%s.%s.pkl'% (dirname,args.o,i), 'rb')
		vars()[i] = pickle.load(pkl_file)
		pkl_file.close()

for i in v:
	if '%s/%s.%s.%s.wig' % (dirname,i,args.o,args.g) not in glob.glob('%s/*.wig' % dirname):
		c = eval(i)
		PrintWiggle(c,i,dirname)

os.system('rm *pkl')
if args.p is not None:
	GrizzlyPeak(chip)


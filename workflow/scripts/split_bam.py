#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 14:00:02 2020

@author: m102324
"""

import sys,os
import pysam
from optparse import OptionParser
import shutil
import subprocess


__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "GPL"
__version__="1.0.5"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Development"


def create_headers(bamfile, ex_chr='',  prog_name='spiker', prog_ver = __version__, co=[]):
	"""
	Create BAM headers for human and exogenous BAM files. Note, to distinguish
	the chromosomes of human and exogenous species, the chrom ID of
	exogenous chromosomes cannot start with "chr" (in other words, "chr" is
	reserved only for human chromosomes). For example, prefix "dm6_" was added
	to the chromosome IDs of Drosophila melanogaster (fruit fly).
	Parameters
	----------
	bamfile : AlignmentFile object
		a=pysam.AlignmentFile('33_K27Ac-112.sorted.bam','rb')
	Returns
	-------
	None.
	"""
	bam_header = bamfile.header #b am header from the composite genome.
	hs_header = {} # bam header for human genome alignments
	ex_header = {} # bam header for exogenous genome (such as fly) alignments
	# update 'SQ'
	for key in bam_header.keys():
		if key == 'SQ':
			if 'SQ' not in hs_header:
				hs_header['SQ'] = []
			if 'SQ' not in ex_header:
				ex_header['SQ'] = []
			lst = bam_header['SQ']
			for i in lst:
				if i['SN'].startswith('chr'):
					hs_header['SQ'].append(i)
				elif i['SN'].startswith(ex_chr):
					ex_header['SQ'].append(i)
				else:
					continue
		else:
			if key not in hs_header:
				hs_header[key] = bam_header[key]
	# update 'PG'
	if 'PG' in hs_header:
		hs_header['PG'].append( {'ID':prog_name,'PN':prog_name, 'VN':prog_ver})
	else:
		hs_header['PG'] = [{'ID':prog_name,'PN':prog_name, 'VN':prog_ver}]
	if 'PG' in ex_header:
		ex_header['PG'].append( {'ID':prog_name,'PN':prog_name, 'VN':prog_ver})
	else:
		ex_header['PG'] = [{'ID':prog_name,'PN':prog_name, 'VN':prog_ver}]
	# update 'CO'
	for comment in co:
		if 'CO' in hs_header:
			hs_header['CO'].append(comment)
		else:
			hs_header['CO'] = [comment]
		if 'CO' in ex_header:
			ex_header['CO'].append(comment)
		else:
			ex_header['CO'] = [comment]
	return(hs_header, ex_header)

def lookup_refid(hd,chr_name):
	lst = hd['SQ']
	for i in range(len(lst)):
		if chr_name == lst[i]['SN']:
			return i
	return None

def sort_bam(bamfile, threads = 1):

	sorted_bam = bamfile.replace('.bam','.sorted.bam')
	sorted_bam_bai = sorted_bam + '.bai'

	# find samtools command
	samtools_cmd = shutil.which("samtools")
	if samtools_cmd is None:
		sys.exit("\tCannot find the \"samtools\" command!")


	# sort BAM file
	if os.path.exists(sorted_bam) and os.path.getsize(sorted_bam) > 0:
		print("\t\"%s\" exists and non-empty, skip BAM sorting." % sorted_bam)
		pass
	else:
		print("\tSorting %s ..." % bamfile)
		samtools_sort = "%s sort -@ %d %s > %s" % (samtools_cmd, threads, bamfile, sorted_bam)
		subprocess.call(samtools_sort, shell=True)

	# index BAM
	if os.path.exists(sorted_bam_bai) and os.path.getsize(sorted_bam_bai) > 0:
		print("\t\"%s\" exists and non-empty, skip BAM indexing." % sorted_bam_bai)
		pass
	else:
		print("\Indexing %s ..." % sorted_bam)
		samtools_index = "%s index %s" % (samtools_cmd, sorted_bam)
		subprocess.call(samtools_index, shell=True)

 	# remove not sorted bam file
	if os.path.exists(bamfile):
		os.remove(bamfile)


def divided_bam(bam_file, outfile, q_cut=30, chr_prefix='dm6_', threads = 1):

	unmapped_reads = 0
	qcfail_reads = 0
	duplicate_reads = 0
	secondary_reads = 0
	low_maq = 0
	diff_genome = 0	#read pairs mapped to fly and human simutaneously
	fly_reads = 0
	human_reads = 0
	samfile = pysam.AlignmentFile(bam_file,'rb')
	human_header, ex_header = create_headers(samfile, ex_chr = chr_prefix)
	OUT = open(outfile + '.report.txt','w')
	HU = pysam.AlignmentFile(outfile + '_ref.bam', "wb", header=human_header)
	EX = pysam.AlignmentFile(outfile + '_spike.bam', "wb", header=ex_header)
	#PE = False
	try:
		while(1):
			aligned_read = next(samfile)
			if aligned_read.is_unmapped:
				unmapped_reads += 1
				continue
			elif aligned_read.is_qcfail:
				qcfail_reads += 1
				continue
			elif aligned_read.is_duplicate:
				duplicate_reads += 1
				continue
			elif aligned_read.is_secondary:
				secondary_reads += 1
				continue
			elif aligned_read.mapq < q_cut:
				low_maq += 1
				continue
			else:
				if aligned_read.is_paired:
					#PE = True
					read1_chr = samfile.get_reference_name(aligned_read.reference_id)
					read2_chr = samfile.get_reference_name(aligned_read.next_reference_id)

					if read1_chr.startswith(chr_prefix):
						# both reads maped to fly
						if read2_chr.startswith(chr_prefix):
							if aligned_read.is_read1:
								aligned_read.reference_id = lookup_refid(ex_header, read1_chr)
							else:
								aligned_read.reference_id = lookup_refid(ex_header, read2_chr)
							EX.write(aligned_read)
							fly_reads += 1
						
						# read-1 mapped to fly, read-2 mapped to human
						else:
							if aligned_read.is_read1:
								aligned_read.reference_id = lookup_refid(ex_header, read1_chr)
							else:
								aligned_read.reference_id = lookup_refid(human_header, read2_chr)							
							diff_genome += 1
					else:
						# read-1 mapped to human, read-2 mapped to fly
						if read2_chr.startswith(chr_prefix):
							if aligned_read.is_read1:
								aligned_read.reference_id = lookup_refid(human_header, read1_chr)
							else:
								aligned_read.reference_id = lookup_refid(ex_header, read2_chr)	
							diff_genome += 1
						
						# both reads maped to human
						else:
							if aligned_read.is_read1:
								aligned_read.reference_id = lookup_refid(human_header, read1_chr)
							else:
								aligned_read.reference_id = lookup_refid(human_header, read2_chr)							
							HU.write(aligned_read)
							human_reads += 1
				else:
					# single-end reads
					read_chr = samfile.get_reference_name(aligned_read.reference_id)
					# reads maped to fly genome
					if read_chr.startswith(chr_prefix):
						aligned_read.reference_id = lookup_refid(ex_header, read_chr)
						EX.write(aligned_read)
						fly_reads += 1
					# reads maped to human genome
					else:
						aligned_read.reference_id = lookup_refid(human_header, read_chr)
						HU.write(aligned_read)
						human_reads += 1
	except StopIteration:
		print("Done")

	HU.close()
	EX.close()


	readInfo = open(snakemake.log["readsInfo"], 'w')
	readInfo.write("Ref-spike_commonReads:%s\n" % diff_genome)
	readInfo.write("NSampleReads:%s\n" % human_reads)
	readInfo.write("NSpikeReads:%s\n" % fly_reads)
	readInfo.write("UnmappedReads:%s\n" % unmapped_reads)
	readInfo.write("QCFailReads:%s\n" % qcfail_reads)
	readInfo.write("SecondaryReads:%s\n" % secondary_reads)
	readInfo.write("LowMapQReads:%s\n" % low_maq)
	readInfo.close()


	for i in ((outfile + '_ref.bam'), (outfile + '_spike.bam')):
		sort_bam(bamfile = i, threads = threads)


	print("Sample\tn_unmapped\tn_qcFail\tn_duplicate\tn_secondary\tn_low_mapq\tn_both\tn_sample\tn_exogenous", file = OUT)
	print("%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d" % (os.path.basename(bam_file),unmapped_reads, qcfail_reads, duplicate_reads, secondary_reads, low_maq, diff_genome, human_reads, fly_reads), file=OUT)
	OUT.close()

def main():

	sys.stdout = open(snakemake.log["otherLog"], 'w')
	sys.stderr = sys.stdout

	# Use Snakemake's input, output, params, etc.
	bam_file = snakemake.input.sample_mixed
	out_prefix = snakemake.params.outprefix
	chr_prefix = "EXO_"
	map_qual = snakemake.params.map_qual
	n_thread = snakemake.threads
	
	
	divided_bam(bam_file = bam_file, outfile = out_prefix, q_cut= map_qual, chr_prefix= chr_prefix, threads = n_thread)

if __name__=='__main__':
	main()


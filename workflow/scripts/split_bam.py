#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 2024

@author: Davide Bressan
"""

import sys,os
import shutil
import subprocess

def sort_bam(bamfile, threads = 1):

	sorted_bam = bamfile.replace('.bam','.sorted.bam')
	sorted_bam_bai = sorted_bam + '.bai'
	samtools_cmd = shutil.which("samtools")

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

	print("Sample %s sorted and indexed." % bamfile)
 	# remove not sorted bam file
	if os.path.exists(bamfile):
		os.remove(bamfile)


def remove_uncertain_reads(filtered_bam, new_filt_bam, chr_prefix, threads = 1):
	"""
	This function removes uncertain reads from a paired end bam file (i.e  when two mates map to different organisms).
	Activated only if the bam file is paired end.
	"""
    # Define the awk script.
    #if there is not EXO_ in both mate chromosomes we keep the read
    #also, if there is EXO_ in both mate chromosomes we keep the read
    #if the mate chromosome is '=' we keep the read (same chromosome)
	samtools_cmd = "samtools view -h %s" % filtered_bam
	awk_script = "awk '{count = ($3 ~ /%s/) + ($7 ~ /%s/); if (count == 0 || count == 2 || $7 == \"=\") print $0}'" % (chr_prefix, chr_prefix)
	saveFile_script = "samtools view -Sb - > %s" % new_filt_bam
    # Run the command using subprocess
	try:
		subprocess.run("%s | %s | %s" % (samtools_cmd, awk_script, saveFile_script), shell=True)
	except:
		sys.exit("\tError in removing uncertain and singletons reads!")
		
	#count reads that were removed in the previous step
	try:
		result=subprocess.run("samtools view -c %s" % new_filt_bam, shell=True, capture_output=True, text=True)
		reads_left = int(result.stdout.strip())
	except:
		sys.exit("\tError in counting uncertain reads!")

	#sort and index the filtered bam file and remove the unfiltered one
	try:
		sort_bam(bamfile = new_filt_bam, threads = threads)
	except:
		sys.exit("\tError in sorting the filtered (no uncertain reads) bam file!")

	return reads_left

def filter_and_split_bam(bam_file, outfile, q_cut=30, chr_prefix='EXO_', threads = 1):

	low_mapq = 0
	exo_reads = 0
	endo_reads = 0
	total_reads = 0
	filtered_reads = 0
	uncertain_reads = 0
	is_paired_end = False


	ENDO_bam=outfile + '_ref.bam'
	EXO_bam=outfile + '_spike.bam'

	#first we count reads on the full bam file
	if os.path.exists(bam_file) and os.path.getsize(bam_file) > 0:
		print("\t\"%s\" exists and non-empty, proceed with read counting." % bam_file)
		result=subprocess.run("samtools view -c %s" % bam_file, shell=True, capture_output=True, text=True)
		total_reads = int(result.stdout.strip())
		print("\tTotal reads: %d" % total_reads)
	else:
		sys.exit("\tCannot find (or it is empty) the input BAM file \"%s\"!" % bam_file)
	 
	#filter out low mapq reads( Skip alignments with MAPQ smaller than q_cut) and sort
	if (q_cut > 0):
		print("\tFiltering out reads with MAPQ < %d..." % q_cut)
		filtered_bam = bam_file.replace('.bam','.filt.bam')
		subprocess.run("samtools view -b -q %d %s > %s" % (q_cut, bam_file ,filtered_bam), shell=True)
		#sort and index the filtered bam file and remove the unfiltered one
		try:
			sort_bam(bamfile = filtered_bam, threads = threads)
		except:
			sys.exit("\tError in sorting the filtered (low mapq) bam file!")
		filtered_bam = filtered_bam.replace('.bam','.sorted.bam')
	else:
		print("\tNo MAPQ filtering applied.")
		filtered_bam = bam_file
	
	#count reads on the filtered bam file
	if (q_cut > 0):
		print("\tCounting reads left after quality filtering in bam file...")
		if os.path.exists(filtered_bam) and os.path.getsize(filtered_bam) > 0:
			print("\t\"%s\" exists and non-empty, proceed with read counting." % filtered_bam)
			result=subprocess.run("samtools view -c %s" % filtered_bam, shell=True, capture_output=True, text=True)
			filtered_reads = int(result.stdout.strip())
		else:
			print("\tFiltered bam file \"%s\" not found!" % filtered_bam)

		#calculate the number of low mapq reads
		low_mapq = total_reads - filtered_reads
		print("\tLow MAPQ reads: %d" % low_mapq)
		#if no reads are left after filtering we exit
		if (total_reads == filtered_reads):
			sys.exit("\tNo reads left after filtering out low MAPQ reads! Consider lowering the MAPQ threshold.")
	else:
		low_mapq = 0
	
	#Now we add another FILTERING to remove uncertain reads (reads that map to both endogenous and exogenous genomes)
	#only for paired end samples
	#check if the bam file is paired end (we check first 1000 to reduce time)
	cmd = """{ samtools view -H %s ; samtools view %s | head -n1000; } | samtools view -c -f 1""" % (filtered_bam, filtered_bam)
	try:
		result=subprocess.run(cmd, shell=True, capture_output=True, text=True)
	except:
		sys.exit("\tError in checking if the bam file is paired end!")

	#only if paired end we remove uncertain reads and singleton reads	
	if (int(result.stdout.strip()) != 0):
		print("\tPaired end sample detected, checking and removing uncertain reads...")
		is_paired_end = True
		if (q_cut > 0):
			new_filt_bam = filtered_bam.replace('.sorted.bam','.filtPE.bam')
		else:
			new_filt_bam = filtered_bam.replace('.bam','.filtPE.bam')
		reads_after_pe_filtering = remove_uncertain_reads(filtered_bam, new_filt_bam ,chr_prefix, threads)
		#if no reads are left after filtering we exit
		if (reads_after_pe_filtering == 0):
			sys.exit("\tNo reads left after filtering out uncertain reads! Consider checking the chromosome names.")
		#calculate the number of uncertain reads to save for report table
		if (q_cut > 0):
			uncertain_reads = filtered_reads - reads_after_pe_filtering
		else:
			uncertain_reads = total_reads - reads_after_pe_filtering
		print("\tUncertain reads: %d" % uncertain_reads)
		print("Uncertain reads removed successfully.")
		#at this point, only if mapq filtering was applied, we remove the old mapq filtered bam file and rename the new one
		if q_cut > 0:
			if os.path.exists(filtered_bam):
				os.remove(filtered_bam)
				os.remove(filtered_bam + '.bai')
		filtered_bam = new_filt_bam.replace('.bam','.sorted.bam')
	else:
		print("\tSingle end sample detected, skipping uncertain reads removal.")
	#Now we split the bam file into endogenous and exogenous bam files
	#First we need to extract the chromosome names
	#first endogenous (human or mouse usually)
	try:
		endo_chromNames=subprocess.run("samtools idxstats %s | cut -f 1 | grep '^chr' | sed 's/^/ /'  | tr '\n' ' '" % filtered_bam, shell=True, capture_output=True, text=True)
		exo_chromNames=subprocess.run("samtools idxstats %s | cut -f 1 | grep '^%s' | sed 's/^/ /'  | tr '\n' ' '" % (filtered_bam, chr_prefix), shell=True, capture_output=True, text=True)
	except:
		sys.exit("\tError in extracting chromosome names!")

	print("\tEndogenous chromosomes: %s" % endo_chromNames.stdout.strip())
	print("\tExogenous chromosomes: %s" % exo_chromNames.stdout.strip())
	#now we generate the two bam files
	try:
		subprocess.run("samtools view -b %s %s > %s" % (filtered_bam, endo_chromNames.stdout.strip(), ENDO_bam), shell=True)
	except:
		sys.exit("\tError in creating the endogenous bam file!")
	try:
		subprocess.run("samtools view -b %s %s > %s" % (filtered_bam, exo_chromNames.stdout.strip(), EXO_bam), shell=True)
	except:
		sys.exit("\tError in creating the exogenous bam file!")

	#count reads on the endogenous bam file
	if os.path.exists(ENDO_bam) and os.path.getsize(ENDO_bam) > 0:
		print("\t\"%s\" exists and non-empty, proceed with read counting." % ENDO_bam)
		result=subprocess.run("samtools view -c %s" % ENDO_bam, shell=True, capture_output=True, text=True)
		endo_reads = int(result.stdout.strip())
	else:
		sys.exit("\tCannot find (or it is empty) the endogenous BAM file \"%s\"!" % ENDO_bam)
	
	#count reads on the exogenous bam file
	if os.path.exists(EXO_bam) and os.path.getsize(EXO_bam) > 0:
		print("\t\"%s\" exists and non-empty, proceed with read counting." % EXO_bam)
		result=subprocess.run("samtools view -c %s" % EXO_bam, shell=True, capture_output=True, text=True)
		exo_reads = int(result.stdout.strip())
	else:
		sys.exit("\tCannot find (or it is empty) the exogenous BAM file \"%s\"!" % EXO_bam)

	#if qc filtering was applied we remove the filtered bam file
	if (q_cut > 0 or is_paired_end):
		os.remove(filtered_bam)
		os.remove(filtered_bam + '.bai')

	readInfo = open(snakemake.log["readsInfo"], 'w')
	readInfo.write("NSampleReads:%s\n" % endo_reads)
	readInfo.write("NSpikeReads:%s\n" % exo_reads)
	readInfo.write("LowMapQReads:%s\n" % low_mapq)
	readInfo.write("UncertainReads:%s\n" % uncertain_reads)
	readInfo.close()

	for i in ((outfile + '_ref.bam'), (outfile + '_spike.bam')):
		sort_bam(bamfile = i, threads = threads)

	print("Sample %s succesfully divided into endogenous and exogenous bam files" % bam_file)


def main():

	sys.stdout = open(snakemake.log["otherLog"], 'w')
	sys.stderr = sys.stdout

	# Use Snakemake's input, output, params, etc.
	bam_file = snakemake.input.sample_mixed
	out_prefix = snakemake.params.outprefix
	chr_prefix = "EXO_"
	map_qual = snakemake.params.map_qual
	n_thread = snakemake.threads
	
	# find samtools command
	samtools_cmd = shutil.which("samtools")
	if samtools_cmd is None:
		sys.exit("\tCannot find the \"samtools\" command!")

	#split bam into two separate bams (endo and exo) and sort 
	filter_and_split_bam(bam_file = bam_file, outfile = out_prefix, q_cut= map_qual, chr_prefix= chr_prefix, threads = n_thread)

if __name__=='__main__':
	main()


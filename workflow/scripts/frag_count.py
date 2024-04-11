#!/usr/bin/env python

'''
Calculate fragment counts (for paired-end) or read counts (for single-end) for each peak.

Example:
python3 frag_count.py -b peaks.bed  -i '*.bam' -o output_counts.tsv
or
python3 frag_count.py -b peaks.bed  -i sample1.bam,sample2.bam,sample3.bam  -o output_counts.tsv
'''

#import built-in modules
import os,sys
if sys.version_info[0] != 3:
        print("\nYou are using python" + str(sys.version_info[0]) + '.' + str(sys.version_info[1]) + " This program needs python3!\n", file=sys.stderr)
        sys.exit()

import pandas as pd
import glob
import pysam
from optparse import OptionParser
from time import strftime

#import third-party modules
#from bx.bitset import *
#from bx.bitset_builders import *
#from bx.intervals import *


__author__ = "Liguo Wang"
__copyright__ = "Copyleft"
__credits__ = []
__license__ = "MIT"
__version__="1.0.5"
__maintainer__ = "Liguo Wang"
__email__ = "wang.liguo@mayo.edu"
__status__ = "Production"



def overlap(a_start, a_end, b_start, b_end):
        return max(0, min(a_end, b_end) - max(a_start, b_start))


def main():
        usage="%prog [options]" + '\n' + __doc__ + "\n"
        parser = OptionParser(usage,version="%prog " + __version__)
        parser.add_option("-i","--bam",action="store",type="string",dest="bam_files", help="comma separated BAM file(s). [required]")
        parser.add_option("-o","--out-prefix",action="store",type="string",dest="output", help="Output counts file. [required]")
        parser.add_option("-b","--bed",action="store",type="string",dest="peak_bed", help="Peak file in bed format (at least 3 columns). [required]")
        parser.add_option("-u","--skip-multi-hits",action="store_true",dest="skip_multi", help="How to deal with multiple hit reads. Presence this option renders program to skip multiple hits reads.")
        parser.add_option("-q","--mapq",action="store",type="int",dest="map_qual",default=30, help="The minimum mapping quality (phred scaled). default=%default")
        parser.add_option("--leng",action="store",type="int",dest="overlap_len",default=1, help="The minimum *ChIP fragment (or read) length* required to overlap with a peak. default=%default")
        parser.add_option("--frac",action="store",type="float",dest="overlap_frac",default=0.0, help="The minimum *fraction of ChIP fragment (or read)* required to overlap with a peak. default=%default")

        (options,args)=parser.parse_args()

        if not (options.output and options.bam_files and options.peak_bed):
                parser.print_help()
                sys.exit(0)

        if not os.path.exists(options.peak_bed):
                print(options.peak_bed + " does NOT exists" + '\n', file=sys.stderr)
                sys.exit(0)

        #get all bam files
        bam_files = []
        if ',' in options.bam_files:
                bam_files = options.bam_files.replace(' ','').split(',')
        else:
                bam_files = glob.glob(options.bam_files)
        bam_files = sorted(bam_files)

        for file in bam_files:
                if not os.path.exists(file):
                        print(file + " does NOT exists" + '\n', file=sys.stderr)
                        sys.exit(0)
                if not os.path.exists(file + '.bai'):
                        print(file + '.bai' + " does NOT exists" + '\n', file=sys.stderr)
                        sys.exit(0)
        print ("Number of BAM files: %d" % len(bam_files), file=sys.stderr)
        #get all genomic regions
        print ("reading BED file ...", file=sys.stderr)
        genomic_regions = []
        chromosomes_all = []
        start_all = []
        end_all = []
        for line in open(options.peak_bed,'r'):
                if line.startswith(('#','track','browser')):
                        continue
                fields = line.split()
                if len(fields) < 3:
                        continue
                chrom     = fields[0]
                start  = fields[1]
                end    = fields[2]
                genomic_regions.append(chrom + ':' + start + '-' + end)
                chromosomes_all.append(chrom)
                start_all.append(int(start))
                end_all.append(int(end))

        counts_dict = {}
        for bam_file in bam_files:
                print ("reading BAM file %s" % bam_file, file=sys.stderr)
                sample_name = os.path.basename(bam_file).replace('_ref.sorted.bam','')
                obj = pysam.AlignmentFile(bam_file, 'rb')

                frag_counts = []
                for region in genomic_regions:
                        frag_count = 0
                        #we need to check if chrom name contains also ':' (liek this chr10:14000000-25000000)
                        # Use rsplit to split from the right, limiting to 2 splits
                        parts = region.rsplit(':', 1)

                        # The first part of the split will contain the chromosome and its range
                        chrom = parts[0]
                        # The second and third parts are the start and end, respectively
                        start = int(parts[1].split('-')[0])
                        end =  int(parts[1].split('-')[1])
                        try:
                                alignedReads = obj.fetch(chrom,start,end)
                        except:
                                frag_counts.append(0)
                                continue
                        for aligned_read in alignedReads:
                                if aligned_read.is_qcfail:continue                      #skip low quanlity
                                if aligned_read.is_duplicate:continue           #skip duplicate read
                                if aligned_read.is_secondary:continue           #skip non primary hit
                                if aligned_read.is_supplementary:continue
                                if options.skip_multi:
                                        if aligned_read.mapq < options.map_qual:
                                                continue
                                # pair-end sequencing
                                if aligned_read.is_paired:
                                        if not aligned_read.is_proper_pair:
                                                continue
                                        if aligned_read.is_read2:
                                                continue
                                        frag_start = aligned_read.reference_start
                                        frag_end = aligned_read.next_reference_start + aligned_read.query_length
                                        if frag_start > frag_end:
                                                (frag_start, frag_end) = (frag_end, frag_start)
                                        frag_length = frag_end - frag_start
                                        if frag_length <= 0:
                                                continue
                                        overlap_length = overlap(start, end, frag_start, frag_end)
                                        if  (overlap_length >= options.overlap_len) and (overlap_length/overlap_length >= options.overlap_frac):
                                                frag_count += 1
                                #single-end sequencing
                                else:
                                        read_length = aligned_read.query_length
                                        if read_length <= 0:
                                                continue
                                        read_start = aligned_read.reference_start
                                        read_end = read_start + read_length
                                        overlap_length = overlap(start, end, read_start, read_end)
                                        if  (overlap_length >= options.overlap_len) and (overlap_length/read_length >= options.overlap_frac):
                                                frag_count += 1
                        frag_counts.append(frag_count)
                counts_dict[sample_name] = frag_counts

        counts_df = pd.DataFrame(counts_dict, index=genomic_regions)
        counts_df.insert(0, 'chr', chromosomes_all)
        counts_df.insert(1, 'start', start_all)
        counts_df.insert(2, 'end', end_all)
        print ("writing counts to %s" % options.output, file=sys.stderr)
        counts_df.to_csv(options.output, sep="\t", index_label='region')

if __name__ == '__main__':
        main()

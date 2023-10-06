import pysam
import sys

#------ Functions -------#
def collect_read_names(bamfile):
    """Collects the names of all reads in a BAM file."""
    return {read.query_name for read in bamfile.fetch()}


def write_unique_reads(bamfile1, bamfile2):
    # This function will open the bam files and write to a new ones just those reads that aren't shared by both bam files
    bam1_before = pysam.AlignmentFile(bamfile1,'rb') #ref
    bam2_before = pysam.AlignmentFile(bamfile2,'rb') #spike

    ids_1= collect_read_names(bam1_before)
    ids_2= collect_read_names(bam2_before)

    # Find the common reads between ref and spike
    common_reads = ids_1.intersection(ids_2)
    print("Ref-spike_commonReads:%s" % len(common_reads))

    # Open the new bam files where were going to write the outputs
    bam1_cleaned = pysam.AlignmentFile(snakemake.output['sample_ref'],'wb', template=bam1_before)
    bam2_cleaned = pysam.AlignmentFile(snakemake.output['sample_spike'],'wb', template=bam2_before)

    # Print the bam lines that are not shared by sample and spike to new bam files
    #ref
    counter_ref=0
    for read in bam1_before.fetch():
        if read.query_name not in common_reads:
            bam1_cleaned.write(read)
            counter_ref+=1
    
    print("NSampleReads:%s" % counter_ref)
    bam1_before.close()
    bam1_cleaned.close()

    #spike
    counter_spike=0
    for read in bam2_before.fetch():
        if read.query_name not in common_reads:
            bam2_cleaned.write(read)
            counter_spike+=1
            
    print("NSpikeReads:%s" % counter_spike)
    bam2_before.close()
    bam2_cleaned.close()

    pysam.index(snakemake.output['sample_ref'])
    pysam.index(snakemake.output['sample_spike'])


#------ Main -------#

sys.stdout = open(snakemake.log[0], 'w')
sys.stderr = sys.stdout


bam1, bam2 = snakemake.input[0], snakemake.input[1]
write_unique_reads(bam1, bam2)

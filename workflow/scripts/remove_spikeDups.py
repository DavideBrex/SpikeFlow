import pysam
import sys

def collect_read_names(bamfile):
    """Collects the names of all reads in a BAM file."""
    return {read.query_name for read in bamfile.fetch()}


def write_unique_reads(bamfile1, bamfile2):
    # This function will open the bam files and write to a new ones just those reads that aren't shared by both bam files
    bam1_open = pysam.AlignmentFile(bamfile1,'rb')
    bam2_open = pysam.AlignmentFile(bamfile2,'rb')

    ids_1= collect_read_names(bam1_open)
    ids_2= collect_read_names(bam2_open)

    common_reads = ids_1.intersection(ids_2)

    # Open the new bam files where were going to write the output
    bam1_outfile = pysam.AlignmentFile(bamfile1 + '.temporary', 'wb', template=bam1_open)
    bam2_outfile = pysam.AlignmentFile(bamfile2 + '.temporary', 'wb', template=bam2_open)

    # Print the bam lines that are not shared by sample and spike to new bam files
    counter = 0
    for read in bam1_open.fetch():
        print("test")
        if read.query_name not in common_reads:
            bam1_outfile.write(read)
    
    bam1_open.close()
    bam1_outfile.close()

    for read in bam2_open.fetch():
        if read.query_name not in common_reads:
            bam2_outfile.write(read)
        else:
            counter += 1
            
    bam2_open.close()
    bam2_outfile.close()
    print("Removed %s reads" % (counter))

if __name__ == "__main__":
    bam1, bam2 = sys.argv[1], sys.argv[2]
    write_unique_reads(bam1, bam2)

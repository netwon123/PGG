import pysam

# read the bam file
bamfile = pysam.AlignmentFile("alignment.bam", "rb")

# get the rate
alignment_rate = bamfile.mapped / (bamfile.mapped + bamfile.unmapped)

print("Alignment rate: %.2f%%" % (alignment_rate * 100))



#Consistent with the result producted from "samtools flagstat"

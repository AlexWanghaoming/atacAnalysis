from Bio import SeqIO, bgzf
from gzip import open as gzopen
import random

fq1 = SeqIO.parse(gzopen("/home/wanghm/whm/ATAC/S0821_05A_CHG036758-Lane41-PH1-7d-ACAGTGGT_L001_R1.fastq.gz","rt"), format="fastq")
fq2 = SeqIO.parse(gzopen("/home/wanghm/whm/ATAC/S0821_05A_CHG036758-Lane41-PH1-7d-ACAGTGGT_L001_R2.fastq.gz","rt"), format="fastq")

handle_out_rep1_r1 = bgzf.BgzfWriter("/home/wanghm/whm/ATAC/split/sex_rep1_r1.fastq.gz", "ab")
handle_out_rep1_r2 = bgzf.BgzfWriter("/home/wanghm/whm/ATAC/split/sex_rep1_r2.fastq.gz", "ab")
rep1_count = 0
handle_out_rep2_r1 = bgzf.BgzfWriter("/home/wanghm/whm/ATAC/split/sex_rep2_r1.fastq.gz", "ab")
handle_out_rep2_r2 = bgzf.BgzfWriter("/home/wanghm/whm/ATAC/split/sex_rep2_r2.fastq.gz", "ab")
rep2_count = 0

ll = [[handle_out_rep1_r1, handle_out_rep1_r2, rep1_count], [handle_out_rep2_r1, handle_out_rep2_r2, rep2_count]]

reads_count = 42304520      # reads count in fastq file
for seq in zip(fq1, fq2):

    tmp_repo = random.choice(ll)    # choose random element in list 
    tmp_repo[2] = tmp_repo[2] + 1
    r1 = seq[0]
    r2 = seq[1]

    if rep1_count < reads_count/2 and rep2_count < reads_count/2:
        SeqIO.write(sequences=r1, handle= tmp_repo[0], format="fastq")
        SeqIO.write(sequences=r2, handle= tmp_repo[1], format="fastq")

    if rep1_count >= reads_count/2 and rep2_count < reads_count/2:
        SeqIO.write(sequences=r1, handle= handle_out_rep2_r1, format="fastq")
        SeqIO.write(sequences=r2, handle= handle_out_rep2_r2, format="fastq")

    if rep1_count < reads_count/2 and rep2_count >= reads_count/2:
        SeqIO.write(sequences=r1, handle= handle_out_rep1_r1, format="fastq")
        SeqIO.write(sequences=r2, handle= handle_out_rep1_r2, format="fastq")



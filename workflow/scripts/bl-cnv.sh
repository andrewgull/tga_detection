# check reads length
seqkit watch --fileds ReadLen data/DA61218.fq.gz

# remove reads shorter than AR (3.5 kb)
filtlong --min_length 3500 data/DA61218.fq.gz > data/DA61218_filt.fq

# convert fastq to fasta
seqkit fq2fa data/DA61218_filt.fq > data/DA61218_filt.fasta

# prepare raw reads for alignment using QAlign
FASTA="data/DA61218_filt.fasta"
OUTDIR="reads_quant"
python QAlign/qalign/main.py convert --input_fasta $FASTA --outdir $OUTDIR --qlevel 2 --rc 1 --kmerpath QAlign/qalign/kmermap_models/r9.4_6mer_nucleotide_template_model.txt

# convert ref gbk to fasta
python -c "from Bio import SeqIO; SeqIO.convert('examples/plasmid_DA61218.gbk', 'genbank', 'DA61218_ref.fasta', 'fasta');"

# map quantized reads with minimap
minimap2 -a data/DA61218_ref.fasta reads_quant/DA61218.q2.fasta > results/DA61218_q2.sam
# sort bam to look at coverage
samtools view -b results/DA61218_q2.sam | samtools sort -o results/DA61218_q2.bam -O BAM -@ 2 && samtools index results/DA61218_q2.bam
samtools coverage -m results/DA61218_q2.bam # no reads mapped!!!


# map raw ONT reads
minimap2 -ax map-ont data/DA61218_ref.fasta data/DA61218_filt.fq > results/DA61218.sam
# sort bam to look at coverage
samtools view -b results/DA61218.sam | samtools sort -o results/DA61218.bam -O BAM -@ 2 && samtools index results/DA61218.bam && samtools coverage -m results/DA61218.bam # 100% bases mapped: there is a coverage break at 31k


# map raw (un-quantized) reads with GraphMap
graphmap2 align -r data/DA61218_ref.fasta -d data/DA61218_filt.fq -o results/DA61218_gm.sam --in-fmt fastq --out-fmt sam
# sort bam to look at coverage
samtools view -b results/DA61218_gm.sam | samtools sort -o results/DA61218_gm.bam -O BAM -@ 2 && samtools index results/DA61218_gm.bam && samtools coverage -m results/DA61218_gm.bam # NOTHING


##############

## MAPPING ON PLASMID W 3 BLA
python -c "from Bio import SeqIO; SeqIO.convert('examples/plasmid_DA61218_3_copies_of_amplified_unit.gbk', 'genbank', 'DA61218_ref_3AU.fasta', 'fasta');"

# map filtered reads with minimap2 (the working strategy from above)
# map raw ONT reads
minimap2 -ax map-ont data/DA61218_ref_3AU.fasta data/DA61218_filt.fq > results/DA61218_3AU.sam
# convert
samtools view -b results/DA61218_3AU.sam | samtools sort -o results/DA61218_3AU.bam -O BAM -@ 2 && samtools index results/DA61218_3AU.bam && samtools coverage -m results/DA61218_3AU.bam
# to get a table for coverage histogram
bedtools coverage -a DA61218_3AU.bed -b results/DA61218_3AU.bam -d > results/DA61218_3AU_cov.tab


# check R script for vizualization of the coverage break with both sets of reads

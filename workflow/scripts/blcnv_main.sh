# path to read set: reads/FAU...
READS=$1

# minimal read length: 3500
MINLEN=$2

# output prefix: FAU.../DA...
OUTPUT_PREFIX=$3

FILTFQ="reads/"$OUTPUT_PREFIX"_filt.fq"
FILTFA="reads/"$OUTPUT_PREFIX"_filt.fasta"

# flanking region sequence: reads/flanking_region_1.fasta
FR1=$4
FR2=$5

# threads
THREADS=$6

# remove reads shorter than AR (3.5 kb)
filtlong --min_length $MINLEN $READS > $FILTFQ

# convert fastq to fasta
seqkit fq2fa -j $TREADS $FILTFQ > $FILTFA

# count number of reads in fasta file
NREADS=$(grep -c ">" $FILTFA)
echo "number of records in fasta file: $NREADS"

# BLAST
blastn -query $FR1 -subject $FILTFA -outfmt 6 -num_threads $THREADS -num_alignments $NREADS > "tables/"$OUTPUT_PREFIX"_FR1_blast.tsv"
blastn -query $FR2 -subject $FILTFA -outfmt 6 -num_threads $THREADS -num_alignments $NREADS > "tables/"$OUTPUT_PREFIX"_FR2_blast.tsv"

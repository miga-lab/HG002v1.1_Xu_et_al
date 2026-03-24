## dimelo_density.py 

This script takes a BAM file with modified base tags (from DiMeLo-seq experiments) and computes the methylation density across genomic regions defined in a BED file. It:

Tiles regions into fixed-size windows (default 1000 bp)
For each window, calculates the fraction of modified bases (either 6mA or CpG 5mC) per read
Filters out low-confidence modification calls below a threshold
Outputs a TSV with columns: chrom, [start, end], density, coverage

It supports parallel processing across chromosomes.

calling command:

python script_name.py \
  --bam /path/to/input.bam \
  --bed /path/to/regions.bed \
  --ref /path/to/reference.fasta \
  --mod-tag A \
  --threshold 0.5 \
  --output /path/to/output.tsv \
  --threads 4 \
  --window-size 1000

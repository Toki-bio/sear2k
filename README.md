# sear2k - Short sequence search tool

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A tool for finding short sequences (like SINEs) in genomes using `ssearch36` with fragment-based searching.

Fast version; consider sear2kloop for more thorough but slower algorithm.

## Features

- Efficient genome-wide searching by splitting into overlapping fragments.
  Prevents ssearch36 from missing hits in long sequences.
  For my puprose 2500 bp fragments with 500 bp overlap are good choice, can be adjusted in the code to match for longer sequences.
- Adjustable search thresholds:
   Minimum hit length: 90% of query length (default: 0.9)
   Minimum similarity: 65% identity (default: 65)
- Multiple search rounds with depletion of the original bank after each round for comprehensive detection. 
- **BED** and **FASTA** output formats.
  Bitscore value and query similarity length are reported in bed for each hit.
  Overlapping hits are merged by coordinates.
  FASTA output bank contains added 50 bp flanks on both sides.
- **Reproducibility note**: Due to heuristic algorithms in FASTA, results for high-copy-number elements may vary between runs.

## Dependencies

- [`ssearch36`](https://fasta.bioch.virginia.edu/fasta_www2/fasta_list2.shtml) (from the FASTA package)
- [`bedtools`](https://bedtools.readthedocs.io/)
- [`seqkit`](https://bioinf.shenwei.me/seqkit/)
- [`samtools`](http://www.htslib.org/)


### Install sear2k

git clone https://github.com/yourusername/sear2k.git
cd sear2k
chmod +x sear2k

## Limitations

**Underscores in FASTA names**: Sequence names with underscores are not supported. Preprocess your genome file like so:
awk '{gsub(/_/,"")}1' genome.fna > genome.bnk

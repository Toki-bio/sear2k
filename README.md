# sear2k - short sequence search in genomes

`sear` scans a genome FASTA with short query sequences using `ssearch36`, splitting the genome into overlapping fragments to avoid boundary misses. It merges fragment hits back into genome coordinates, extracts sequence with flanks, and can optionally keep the top-N best hits and align them with MAFFT.

## Highlights

- Fragmented scanning (2500 bp windows, 2000 bp step)
- Strand-aware merging and coordinate reconstruction
- Optional best-hit selection + MAFFT alignment
- Safe handling of long headers and cleanup of intermediates
- Optional slop control for output BED/FASTA

## Requirements

Required:
- `ssearch36` (FASTA package)
- `bedtools`
- `seqkit`
- `samtools`

Optional (best-hit mode):
- `mafft`

## Install

Make the script executable and ensure dependencies are on your `PATH`:

```bash
chmod +x sear
```

## Usage

```bash
sear [options] <query.fasta> <bank.fasta> [Slf=0.8] [Homology=65] [Flank=50] [KeepSplits=0]
```

Options:
- `-b, --best [N]`: keep top N unique hits by bitscore (default 100) and run MAFFT
- `-m, --minus`: enable minus-bank depletion (multi-round per split)
- `-s, --stop [N]`: stop after N total hits (default 1000)
- `-k, --keep`: keep split files for reuse
- `--slop L,R`: override flank sizes for output bed/bnk and best-hit flanks

Notes:
- `--slop` overrides the positional `Flank` value.
- If your query FASTA does not end with a newline, `sear` inserts one before MAFFT input to avoid header leakage.

## Examples

Basic run:

```bash
sear query.fa genome.fa
```

Best-hit mode:

```bash
sear --best 200 query.fa genome.fa
```

Multi-round search with minus bank:

```bash
sear --minus query.fa genome.fa
```

Custom slop for outputs:

```bash
sear --slop 10,20 query.fa genome.fa
```

## Outputs

Let:
- `qname` = query filename without extension
- `taxname` = first 3 characters of the working genome filename

Main outputs:
- `taxname-qname.bed`: merged hit intervals with stats
- `taxname-qname.bnk`: extracted FASTA with flanks (or `--slop` values)

Best-hit mode outputs:
- `taxname-qname.bestN.mafft.fa`: MAFFT alignment of query + top-N hits

Intermediates are cleaned automatically unless `--keep` is used for split reuse.

## License

MIT
